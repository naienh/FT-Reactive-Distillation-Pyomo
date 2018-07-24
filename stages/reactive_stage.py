# 2nd level block Structure: Stage Block
from global_sets.component import m
from utility.model_utility import select_MPCC
from pyomo import environ as pe

# stage construction rules
from physics.kinetics.kinetics_reactive import kinetic_block_rule
from physics.energy.energy_reactive import energy_block_rule
from physics.VLE.VLE_reactive_MPCC_P import VLE_block_rule
from physics.MPCC.MPCC_P import P_NCP_block_rule, P_Reg_block_rule, P_pf_block_rule

def reactive_stage_rule(block,j):
    #-----------------------------------SETS-----------------------------------

    # local sets that will only be used in reactive stage
    block.inlet = pe.Set(initialize=['in'])
    block.outlet = pe.Set(initialize=['out','P'])
    block.stream = block.inlet | block.outlet

    #---------------------------------VARIABLES---------------------------------

    # Tray Inlet/Outlet Variable
    block.x_ = pe.Var(block.inlet,m.COMP_TOTAL,within=pe.NonNegativeReals)
    block.y_ = pe.Var(block.inlet,m.COMP_TOTAL,within=pe.NonNegativeReals)
    block.x = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals)
    block.y = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals)
    block.z = pe.Var(m.COMP_FEED,within=pe.NonNegativeReals)

    block.L = pe.Var(block.stream,within=pe.NonNegativeReals)
    block.V = pe.Var(block.stream,within=pe.NonNegativeReals)
    block.F = pe.Var(within=pe.NonNegativeReals)

    block.H_L_ = pe.Var(block.inlet,within=pe.Reals)
    block.H_V_ = pe.Var(block.inlet,within=pe.Reals)
    block.H_L = pe.Var(within=pe.Reals)
    block.H_V = pe.Var(within=pe.Reals)
    block.H_F = pe.Var(within=pe.Reals)

    # State Variable
    block.T = pe.Var(within=pe.NonNegativeReals,bounds=(100+273.15,350+273.15)) # K
    block.T_F = pe.Var(within=pe.NonNegativeReals) # K
    block.P = pe.Var(within=pe.NonNegativeReals,bounds=(10,30)) # Bar

    block.f_V = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals,initialize=1e-20)
    block.f_L = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals,initialize=1e-20)

    block.cat = pe.Var(within=pe.NonNegativeReals,initialize=3000) # kg
    block.Q_main = pe.Var(within=pe.Reals) # MW
    block.r_total_comp = pe.Var(m.COMP_TOTAL,within=pe.Reals) # kmol/s

    block.PR_L = pe.Var(within=pe.NonNegativeReals,bounds=(0,1))
    block.PR_V = pe.Var(within=pe.NonNegativeReals,bounds=(0,1))

    print('>','Importing Reactive Stage......')
    print('>','Adding the following local variable:')
    print('-'*36)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)

    print('-'*36)
    print('')
    #---------------------------------Equations---------------------------------
    #with HiddenPrints():
    # Kinetics Block
    block.kinetics_block = pe.Block(rule=kinetic_block_rule)

    # Energy Block
    block.energy_block = pe.Block(rule=energy_block_rule)

    # VLE block
    block.VLE_block = pe.Block(rule=VLE_block_rule)

    # Mass Balance
    def mass_balance_main_rule(block,i):
        if i in m.COMP_FEED:
            return block.F*block.z[i] + sum(block.L[s]*block.x_[s,i] + block.V[s]*block.y_[s,i] for s in block.inlet)\
            + block.r_total_comp[i] - sum(block.L[s]*block.x[i] + block.V[s]*block.y[i] for s in block.outlet) == 0
        else:
            return sum(block.L[s]*block.x_[s,i] + block.V[s]*block.y_[s,i] for s in block.inlet)\
            + block.r_total_comp[i] - sum(block.L[s]*block.x[i] + block.V[s]*block.y[i] for s in block.outlet) == 0
    block.mass_balance_main_con = pe.Constraint(m.COMP_TOTAL,rule=mass_balance_main_rule)

    # Equilibrium
    def VL_equil_rule(block,i):
        return block.f_V[i] == block.f_L[i]
    block.VL_equil_con = pe.Constraint(m.COMP_TOTAL,rule=VL_equil_rule)

    # MPCC formation
    block.MPCC_P_pf = pe.Block(rule = P_pf_block_rule)
    block.MPCC_P_NCP = pe.Block(rule = P_NCP_block_rule)
    block.MPCC_P_Reg = pe.Block(rule = P_Reg_block_rule)

    # by default deactivated, can be switched after the block is constructed
    select_MPCC(block,'pf')

    # Summation
    def summation_x_y_rule(block):
        return sum(block.x[i] for i in m.COMP_TOTAL) == sum(block.y[i] for i in m.COMP_TOTAL)
    block.summation_x_y_con = pe.Constraint(rule=summation_x_y_rule)

    def summation_total_mass_rule(block):
        return block.F + sum(block.L[s] + block.V[s] for s in block.inlet) + sum(block.r_total_comp[i] for i in m.COMP_TOTAL)\
                - sum(block.L[s] + block.V[s] for s in block.outlet) == 0
    block.summation_total_mass_con = pe.Constraint(rule=summation_total_mass_rule)

    # Heat Balance
    def heat_balance_main_rule(block):
        return block.F*block.H_F + sum(block.L[s]*block.H_L_[s] + block.V[s]*block.H_V_[s] for s in block.inlet) \
                + block.Q_main - sum(block.L[s]*block.H_L + block.V[s]*block.H_V for s in block.outlet) == 0
    block.heat_balance_main_con = pe.Constraint(rule=heat_balance_main_rule)

    # product / out ratio
    def PR_L_ratio_rule(block):
        return block.L['P'] == block.PR_L * sum(block.L[s] for s in block.outlet)
    block.PR_L_con = pe.Constraint(rule=PR_L_ratio_rule)

    def PR_V_ratio_rule(block):
        return block.V['P'] == block.PR_V * sum(block.V[s] for s in block.outlet)
    block.PR_V_con = pe.Constraint(rule=PR_V_ratio_rule)

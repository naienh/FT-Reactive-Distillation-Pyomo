# 2nd level block Structure: Stage Block
from global_sets.component import m
from pyomo import environ as pe

# stage construction rules
from physics.energy_bounded import energy_block_rule
from physics.VLE_bounded import VLE_block_rule

def non_reactive_stage_rule(block):
    #-----------------------------------SETS-----------------------------------

    # local sets that will only be used in reactive stage
    block.inlet = pe.Set(initialize=['in'])
    block.outlet = pe.Set(initialize=['out','P'])
    block.stream = block.inlet | block.outlet

    #---------------------------------VARIABLES---------------------------------

    block.T_F = pe.Var(within=pe.NonNegativeReals) # K
    block.P = pe.Var(within=pe.NonNegativeReals,bounds=(10,30)) # Bar
    block.Q_main = pe.Var(within=pe.Reals) # MW
    # Tray Inlet/Outlet Variable
    block.x_ = pe.Var(block.inlet,m.COMP_TOTAL,within=pe.NonNegativeReals,bounds=(0,1))
    block.y_ = pe.Var(block.inlet,m.COMP_TOTAL,within=pe.NonNegativeReals,bounds=(0,1))
    block.x = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals,bounds=(0,1))
    block.y = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals,bounds=(0,1))
    block.z = pe.Var(m.COMP_FEED,within=pe.NonNegativeReals,bounds=(0,1))

    block.L = pe.Var(block.stream,within=pe.NonNegativeReals)
    block.V = pe.Var(block.stream,within=pe.NonNegativeReals)
    block.F = pe.Var(within=pe.NonNegativeReals)

    block.H_L_ = pe.Var(block.inlet,within=pe.Reals)
    block.H_V_ = pe.Var(block.inlet,within=pe.Reals)
    block.H_L = pe.Var(within=pe.Reals)
    block.H_V = pe.Var(within=pe.Reals)

    block.T = pe.Var(within=pe.NonNegativeReals,bounds=(200+273.15,300+273.15)) # K
    block.H_F = pe.Var(within=pe.Reals)
    block.f_V = pe.Var(m.COMP_TOTAL,within=pe.PositiveReals,initialize=1e-20)
    block.f_L = pe.Var(m.COMP_TOTAL,within=pe.PositiveReals,initialize=1e-20)

    print('>','Importing Non Reactive Stage......')
    print('>','Adding the following local variable:')
    print('-'*36)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)

    print('-'*36)
    print('')
    #---------------------------------Equations---------------------------------
    #with HiddenPrints():

    # Energy Block
    block.energy_block = pe.Block(rule=energy_block_rule)

    # VLE block
    block.VLE_block = pe.Block(rule=VLE_block_rule)

    # Mass Balance
    def mass_balance_main_rule(block,i):
        if i in m.COMP_FEED:
            return block.F*block.z[i] + sum(block.L[s]*block.x_[s,i] + block.V[s]*block.y_[s,i] for s in block.inlet)\
             - sum(block.L[s]*block.x[i] + block.V[s]*block.y[i] for s in block.outlet) == 0
        else:
            return sum(block.L[s]*block.x_[s,i] + block.V[s]*block.y_[s,i] for s in block.inlet)\
             - sum(block.L[s]*block.x[i] + block.V[s]*block.y[i] for s in block.outlet) == 0
    block.mass_balance_main_con = pe.Constraint(m.COMP_TOTAL,rule=mass_balance_main_rule)

    # Equilibrium
    def VL_equil_rule(block,i):
        return block.f_V[i] == block.f_L[i]
    block.VL_equil_con = pe.Constraint(m.COMP_TOTAL,rule=VL_equil_rule)

    # Summation
    def summation_x_main_rule(block):
        return sum(block.x[i] for i in m.COMP_TOTAL) == 1
    block.summation_x_main_con = pe.Constraint(rule=summation_x_main_rule)

    def summation_y_main_rule(block):
        return sum(block.y[i] for i in m.COMP_TOTAL) == 1
    block.summation_y_main_con = pe.Constraint(rule=summation_y_main_rule)

    # Heat Balance
    def heat_balance_main_rule(block):
        return block.F*block.H_F + sum(block.L[s]*block.H_L_[s] + block.V[s]*block.H_V_[s] for s in block.inlet) \
                + block.Q_main - sum(block.L[s]*block.H_L + block.V[s]*block.H_V for s in block.outlet) == 0
    block.heat_balance_main_con = pe.Constraint(rule=heat_balance_main_rule)

    def total_mass_balance_rule(block):
        return block.F + sum(block.V[s] + block.L[s] for s in block.inlet) \
         == sum(block.V[s] + block.L[s] for s in block.outlet)
    block.total_mass_balance_con = pe.Constraint(rule=total_mass_balance_rule)
    block.total_mass_balance_con.deactivate()

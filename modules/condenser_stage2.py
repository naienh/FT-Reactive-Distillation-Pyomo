# 2nd level block Structure: Stage Block
from modules.global_set import m
from pyomo import environ as pe
from data.utility import HiddenPrints

# stage construction rules
from modules.energy import energy_block_rule
from modules.VLLE import VLLE_block_rule

# the condenser below represents only VLE, water is removed by fraction
def condenser_stage_rule(block):
    #-----------------------------------SETS-----------------------------------

    # local sets that will only be used in reactive stage
    block.inlet = pe.Set(initialize=['in'])
    block.outlet = pe.Set(initialize=['out','P'])
    block.stream = block.inlet | block.outlet
    block.COMP_WATER = pe.Set(initialize=['H2O'])

    #---------------------------------VARIABLES---------------------------------

    block.T = pe.Var(within=pe.NonNegativeReals,bounds=(20+273.15,40+273.15)) # K
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
    block.W = pe.Var(within=pe.NonNegativeReals)
    block.V = pe.Var(block.stream,within=pe.NonNegativeReals)
    block.F = pe.Var(within=pe.NonNegativeReals)

    block.H_L_ = pe.Var(block.inlet,within=pe.Reals)
    block.H_V_ = pe.Var(block.inlet,within=pe.Reals)
    block.H_L = pe.Var(within=pe.Reals)
    block.H_V = pe.Var(within=pe.Reals)

    block.H_F = pe.Var(within=pe.Reals)
    block.f_V = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals)
    block.f_L = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals)

    print('|','Importing Condenser Stage......')
    print('|','Adding the following local variable:')
    print('-'*36)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)

    print('-'*36)
    print('')

    #---------------------------------Equations---------------------------------

    # Energy Block
    block.energy_block = pe.Block(rule=energy_block_rule)

    # VLE block
    block.VLE_block = pe.Block(rule=VLLE_block_rule)

    # Mass Balance
    def mass_balance_main_rule(block,i):
        if i in m.COMP_FEED:
            return block.F*block.z[i] + sum(block.L[s]*block.x_[s,i] + block.V[s]*block.y_[s,i] for s in block.inlet)\
            - sum(block.L[s]*block.x[i] + block.V[s]*block.y[i] for s in block.outlet) == 0
        elif i in block.COMP_WATER:
            return sum(block.L[s]*block.x_[s,i] + block.V[s]*block.y_[s,i] for s in block.inlet)\
            - sum(block.L[s]*block.x[i] + block.V[s]*block.y[i] for s in block.outlet) - block.W == 0
        else:
            return sum(block.L[s]*block.x_[s,i] + block.V[s]*block.y_[s,i] for s in block.inlet)\
            - sum(block.L[s]*block.x[i] + block.V[s]*block.y[i] for s in block.outlet) == 0
    block.mass_balance_main_con = pe.Constraint(m.COMP_TOTAL,rule=mass_balance_main_rule)

    # Equilibrium
    def VL_equil_rule(block,i):
        return block.f_V[i] == block.f_L[i]
    block.VL_equil_con = pe.Constraint(m.COMP_TOTAL-block.COMP_WATER,rule=VL_equil_rule)

    # Water phase
    def L_water_rule(block,i):
        return block.x[i] == pe.exp(-0.66037 - 7.1130*(539.1/block.T) - 0.67885*(1-block.T/539.1)**(1/3) -1.43381*(1-block.T/539.1))
    block.L_water_con = pe.Constraint(block.COMP_WATER,rule=L_water_rule)

    def V_water_rule(block,i):
        return block.y[i]*block.P == pe.exp(5.11564 - 1687.537/(block.T+230.17-273.15))
    block.V_water_con = pe.Constraint(block.COMP_WATER,rule=V_water_rule)

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

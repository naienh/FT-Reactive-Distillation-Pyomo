# 1st level Model Structure: Equation Block
# import sys
# sys.path.append('..')
# this module define the rules for constructing a energy block in the master block
# this is the global component set import, so that all modules uses the same set
from global_sets.component import m

from pyomo import environ as pe

# define MPCC_beta_NCP rule
def P_NCP_block_rule(block):

    #------------------------------LOCAL VARIABLES------------------------------
    block.s_L = pe.Var(within=pe.NonNegativeReals,initialize=0)#,bounds=(0,1))
    block.s_V = pe.Var(within=pe.NonNegativeReals,initialize=0)

    #-----------------------------LOCAL parameters------------------------------
    block.epi = pe.Param(initialize=1e-4,mutable=True)

    print('>','Importing MPCC_P_NCP Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)
    for i in block.component_objects(pe.Param,active=True):
        print('|',i)

    print('-'*50)

    print('>','Adding complementarity constraint, spliting pressure used in VLE')
    print('')

    #------------------------------MPCC equations-------------------------------
    def s_L_complementarity_rule(block):
        return (sum(block.parent_block().L[s] for s in block.parent_block().outlet) + block.s_L) \
        == ((sum(block.parent_block().L[s] for s in block.parent_block().outlet) - block.s_L)**2+block.epi)**0.5
    block.s_L_complementarity_con = pe.Constraint(rule=s_L_complementarity_rule)

    def s_V_complementarity_rule(block):
        return (sum(block.parent_block().V[s] for s in block.parent_block().outlet) + block.s_V) \
        == ((sum(block.parent_block().V[s] for s in block.parent_block().outlet) - block.s_V)**2+block.epi)**0.5
    block.s_V_complementarity_con = pe.Constraint(rule=s_V_complementarity_rule)

    #-----------------------------Global equations------------------------------
    block.parent_block().VLE_block.del_component(block.parent_block().VLE_block.pressure_equal_con)
    def pressure_equal_rule(block):
        return block.parent_block().VLE_block.P_VLE - block.parent_block().P == block.s_L - block.s_V
    block.pressure_equal_con = pe.Constraint(rule=pressure_equal_rule)

# define MPCC_beta_Reg rule
def P_Reg_block_rule(block):

    #------------------------------LOCAL VARIABLES------------------------------
    block.s_L = pe.Var(within=pe.NonNegativeReals,initialize=0)#,bounds=(0,1))
    block.s_V = pe.Var(within=pe.NonNegativeReals,initialize=0)

    #-----------------------------LOCAL parameters------------------------------
    block.epi = pe.Param(initialize=1e-4,mutable=True)

    print('>','Importing MPCC_P_Reg Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)
    for i in block.component_objects(pe.Param,active=True):
        print('|',i)

    print('-'*50)

    print('>','Adding complementarity constraint, spliting pressure used in VLE')
    print('')

    #------------------------------MPCC equations-------------------------------
    def s_L_complementarity_rule(block):
        return sum(block.parent_block().L[s] for s in block.parent_block().outlet) * block.s_L \
        <= block.epi
    block.s_L_complementarity_con = pe.Constraint(rule=s_L_complementarity_rule)

    def s_V_complementarity_rule(block):
        return sum(block.parent_block().V[s] for s in block.parent_block().outlet) * block.s_V <= block.epi
    block.s_V_complementarity_con = pe.Constraint(rule=s_V_complementarity_rule)

    #-----------------------------Global equations------------------------------
    block.parent_block().VLE_block.del_component(block.parent_block().VLE_block.pressure_equal_con)
    def pressure_equal_rule(block):
        return block.parent_block().VLE_block.P_VLE - block.parent_block().P == block.s_L - block.s_V
    block.pressure_equal_con = pe.Constraint(rule=pressure_equal_rule)

# define MPCC_beta_penalty_function rule
def P_pf_block_rule(block):

    #------------------------------LOCAL VARIABLES------------------------------
    block.s_L = pe.Var(within=pe.NonNegativeReals,initialize=0)#,bounds=(0,1))
    block.s_V = pe.Var(within=pe.NonNegativeReals,initialize=0)
    block.pf = pe.Var(within=pe.NonNegativeReals,initialize=0)

    #-----------------------------LOCAL parameters------------------------------
    block.rho = pe.Param(initialize=1,mutable=True)
    block.epi = pe.Param(initialize=1e-4,mutable=True)

    print('>','Importing MPCC_P_pf Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)
    for i in block.component_objects(pe.Param,active=True):
        print('|',i)

    print('-'*50)
    print('>','Spliting pressure used in VLE')
    print('')

    #------------------------------MPCC equations-------------------------------
    def penalty_rule(block):
        return block.pf >= block.rho*(sum(block.parent_block().L[s] for s in block.parent_block().outlet)*block.s_L \
        + sum(block.parent_block().V[s] for s in block.parent_block().outlet)*block.s_V)
    block.penalty_con = pe.Constraint(rule=penalty_rule)

    #-----------------------------Global equations------------------------------
    block.parent_block().VLE_block.del_component(block.parent_block().VLE_block.pressure_equal_con)
    def pressure_equal_rule(block):
        return block.parent_block().VLE_block.P_VLE - block.parent_block().P == block.s_L - block.s_V
    block.pressure_equal_con = pe.Constraint(rule=pressure_equal_rule)

    # def s_L_complementarity_rule(block):
    #     return (sum(block.parent_block().L[s] for s in block.parent_block().outlet) + block.s_L) \
    #     <= ((sum(block.parent_block().L[s] for s in block.parent_block().outlet) - block.s_L)**2+block.epi)**0.5
    # block.s_L_complementarity_con = pe.Constraint(rule=s_L_complementarity_rule)
    #
    # def s_V_complementarity_rule(block):
    #     return (sum(block.parent_block().V[s] for s in block.parent_block().outlet) + block.s_V) \
    #     <= ((sum(block.parent_block().V[s] for s in block.parent_block().outlet) - block.s_V)**2+block.epi)**0.5
    # block.s_V_complementarity_con = pe.Constraint(rule=s_V_complementarity_rule)

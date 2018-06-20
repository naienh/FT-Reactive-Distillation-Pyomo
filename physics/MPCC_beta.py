# 1st level Model Structure: Equation Block
# import sys
# sys.path.append('..')
# this module define the rules for constructing a energy block in the master block
# this is the global component set import, so that all modules uses the same set
from global_sets.component import m

from pyomo import environ as pe

# define MPCC_beta_NCP rule
def beta_NCP_block_rule(block):

    #------------------------------LOCAL VARIABLES------------------------------
    block.beta = pe.Var(within=pe.NonNegativeReals,initialize=1)
    block.s_L = pe.Var(within=pe.NonNegativeReals,initialize=0,bounds=(0,1))
    block.s_V = pe.Var(within=pe.NonNegativeReals,initialize=0)

    #-----------------------------LOCAL parameters------------------------------
    block.epi = pe.Param(initialize=1e-4,mutable=True)

    print('>','Importing MPCC_beta_NCP Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)
    for i in block.component_objects(pe.Param,active=True):
        print('|',i)

    print('-'*50)

    print('>','Adjusting the f_V = f_L bounds to f_V = beta * f_L')
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

    def beta_rule(block):
        return block.beta == 1 - block.s_L + block.s_V
    block.beta_con = pe.Constraint(rule=beta_rule)

    #-----------------------------Global equations------------------------------
    block.parent_block().del_component(block.parent_block().VL_equil_con)
    def VL_equil_rule(model,i):
        return model.f_V[i] == block.beta * model.f_L[i]
    block.parent_block().VL_equil_con = pe.Constraint(m.COMP_TOTAL,rule=VL_equil_rule)

# define MPCC_beta_Reg rule
def beta_Reg_block_rule(block):

    #------------------------------LOCAL VARIABLES------------------------------
    block.beta = pe.Var(within=pe.NonNegativeReals,initialize=1)
    block.s_L = pe.Var(within=pe.NonNegativeReals,initialize=0,bounds=(0,1))
    block.s_V = pe.Var(within=pe.NonNegativeReals,initialize=0)

    #-----------------------------LOCAL parameters------------------------------
    block.epi = pe.Param(initialize=1e-6,mutable=True)

    print('>','Importing MPCC_beta_Reg Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)
    for i in block.component_objects(pe.Param,active=True):
        print('|',i)

    print('-'*50)

    print('>','Adjusting the f_V = f_L bounds to f_V = beta * f_L')
    print('')

    #------------------------------MPCC equations-------------------------------
    def s_L_complementarity_rule(block):
        return sum(block.parent_block().L[s] for s in block.parent_block().outlet) * block.s_L <= block.epi
    block.s_L_complementarity_con = pe.Constraint(rule=s_L_complementarity_rule)

    def s_V_complementarity_rule(block):
        return sum(block.parent_block().V[s] for s in block.parent_block().outlet) * block.s_V <= block.epi
    block.s_V_complementarity_con = pe.Constraint(rule=s_V_complementarity_rule)

    def beta_rule(block):
        return block.beta == 1 - block.s_L + block.s_V
    block.beta_con = pe.Constraint(rule=beta_rule)

    #-----------------------------Global equations------------------------------
    block.parent_block().del_component(block.parent_block().VL_equil_con)
    def VL_equil_rule(model,i):
        return model.f_V[i] == block.beta * model.f_L[i]
    block.parent_block().VL_equil_con = pe.Constraint(m.COMP_TOTAL,rule=VL_equil_rule)

# define MPCC_beta_penalty_function rule
def beta_pf_block_rule(block):

    #------------------------------LOCAL VARIABLES------------------------------
    block.beta = pe.Var(within=pe.NonNegativeReals,initialize=1)
    block.s_L = pe.Var(within=pe.NonNegativeReals,initialize=0,bounds=(0,1))
    block.s_V = pe.Var(within=pe.NonNegativeReals,initialize=0)
    block.pf = pe.Var(within=pe.NonNegativeReals,initialize=0)

    #-----------------------------LOCAL parameters------------------------------
    block.rho = pe.Param(initialize=1,mutable=True)

    print('>','Importing MPCC_beta_Reg Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)
    for i in block.component_objects(pe.Param,active=True):
        print('|',i)

    print('-'*50)
    print('')

    #------------------------------MPCC equations-------------------------------
    def penalty_rule(block):
        return block.pf == block.rho*(sum(block.parent_block().L[s] for s in block.parent_block().outlet)*block.s_L \
        + sum(block.parent_block().V[s] for s in block.parent_block().outlet)*block.s_V)
    block.s_L_complementarity_con = pe.Constraint(rule=s_L_complementarity_rule)

    def beta_rule(block):
        return block.beta == 1 - block.s_L + block.s_V
    block.beta_con = pe.Constraint(rule=beta_rule)

    #-----------------------------Global equations------------------------------
    block.parent_block().del_component(block.parent_block().VL_equil_con)
    def VL_equil_rule(model,i):
        return model.f_V[i] == block.beta * model.f_L[i]
    block.parent_block().VL_equil_con = pe.Constraint(m.COMP_TOTAL,rule=VL_equil_rule)

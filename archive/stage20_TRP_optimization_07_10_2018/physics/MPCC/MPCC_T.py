# 1st level Model Structure: Equation Block
# import sys
# sys.path.append('..')
# this module define the rules for constructing a dew/bubble point block in the master block
# this is the global component set import, so that all modules uses the same set
from global_sets.component import m

# data import and pre-processing
from pyomo import environ as pe
from physics.VLE_reboiler_MPCC_T import VLE_block_rule

# define MPCC_beta_NCP rule
def dew_block_rule(block):

    #------------------------------LOCAL VARIABLES------------------------------

    block.x = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals)
    block.y = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals)
    block.T = pe.Var(within=pe.NonNegativeReals,initialize=240+273.15,bounds=(170,350+273.15)) # K
    block.P = pe.Var(within=pe.NonNegativeReals,bounds=(10,30)) # Bar
    block.f_V = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals,initialize=1e-20)
    block.f_L = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals,initialize=1e-20)
    block.beta = pe.Var(within=pe.NonNegativeReals,initialize=1)
    # block.s_L = pe.Var(within=pe.NonNegativeReals,initialize=0)
    block.s_V = pe.Var(within=pe.NonNegativeReals,initialize=0)

    #-----------------------------LOCAL parameters------------------------------
    block.epi = pe.Param(initialize=1e-3,mutable=True)

    print('>','Importing dew Blocks......')
    print('>','Adding the following local variable:')
    print('-'*36)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)
    for i in block.component_objects(pe.Param,active=True):
        print('|',i)

    print('-'*36)
    print('')

    #---------------------------------Equations---------------------------------
    block.VLE_block = pe.Block(rule=VLE_block_rule)
    block.VLE_block.n_ave.fix(20)

    # Get information from parent block
    def P_equal_rule(block):
        return block.P == block.parent_block().P
    block.P_equal_con = pe.Constraint(rule=P_equal_rule)

    # calculate the pre-separation mixture
    def mixing_rule(block,i):
        return sum(block.parent_block().L[s]*block.parent_block().x[i] + \
            block.parent_block().V[s]*block.parent_block().y[i] for s in block.parent_block().outlet) == \
            block.y[i] * sum(block.parent_block().V[s] + block.parent_block().L[s] for s in block.parent_block().outlet)
        # return block.y[i] == block.parent_block().y[i]
    block.mixing_con = pe.Constraint(m.COMP_TOTAL,rule=mixing_rule)

    # Phase Equilibrium
    def VL_rule(block,i):
        return block.f_V[i] == block.beta * block.f_L[i]
    block.VL_con = pe.Constraint(m.COMP_TOTAL,rule=VL_rule)

    # summation
    def summation_rule(block):
        return sum(block.x[i] for i in m.COMP_TOTAL) == 1
    block.summation_con = pe.Constraint(rule=summation_rule)

    #------------------------------MPCC equations-------------------------------
    def s_L_complementarity_rule(block):
        return (350+273.15 - block.T)*block.s_V <= block.epi
    block.s_L_complementarity_con = pe.Constraint(rule=s_L_complementarity_rule)

    # def s_V_complementarity_rule(block):
    #     return (block.T - 200-273.15)*block.s_L <= block.epi
    # block.s_V_complementarity_con = pe.Constraint(rule=s_V_complementarity_rule)

    def beta_rule(block):
        return block.beta == 1 + block.s_V #- block.s_L
    block.beta_con = pe.Constraint(rule=beta_rule)

    #-----------------------------Global equations------------------------------
    block.parent_block().VLE_block.del_component(block.parent_block().VLE_block.temperature_equal_con)
    def temperature_equal_rule(block):
        return block.parent_block().VLE_block.T_VLE == 0.5*(block.parent_block().T + block.T - ((block.parent_block().T - block.T)**2 + 100*block.epi)**0.5)
    block.temperature_equal_con = pe.Constraint(rule=temperature_equal_rule)

# 1st level Model Structure: Equation Block
# import sys
# sys.path.append('..')
# this module define the rules for constructing a dew/bubble point block in the master block
# this is the global component set import, so that all modules uses the same set
from global_sets.component import m

# data import and pre-processing
from pyomo import environ as pe
from physics.VLE_bounded import VLE_block_rule

# define MPCC_beta_NCP rule
def dew_block_rule(block):

    #------------------------------LOCAL VARIABLES------------------------------

    block.x = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals)
    block.y = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals)
    block.T = pe.Var(within=pe.NonNegativeReals,bounds=(200+273.15,300+273.15)) # K
    block.P = pe.Var(within=pe.NonNegativeReals,bounds=(10,30)) # Bar
    block.f_V = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals,initialize=1e-20)
    block.f_L = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals,initialize=1e-20)

    print('>','Importing dew Blocks......')
    print('>','Adding the following local variable:')
    print('-'*36)

    for i in block.component_objects(pe.Var,active=True):
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
        return sum(block.parent_block().V[s] for s in block.parent_block().outlet)\
         * block.parent_block().y[i] + sum(block.parent_block().L[s] \
         for s in block.parent_block().outlet) * block.parent_block().x[i] == \
         block.y[i] * (sum(block.parent_block().V[s] for s in block.parent_block().outlet)\
         + sum(block.parent_block().L[s] for s in block.parent_block().outlet))
    block.mixing_con = pe.Constraint(m.COMP_TOTAL,rule=mixing_rule)

    # Phase Equilibrium
    def VL_rule(block,i):
        return block.f_V[i] == block.f_L[i]
    block.VL_con = pe.Constraint(m.COMP_TOTAL,rule=VL_rule)

    # summation
    def summation_rule(block):
        return sum(block.x[i] for i in m.COMP_TOTAL) == sum(block.y[i] for i in m.COMP_TOTAL)
    block.summation_con = pe.Constraint(rule=summation_rule)

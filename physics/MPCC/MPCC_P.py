# 1st level Model Structure: Equation Block
# import sys
# sys.path.append('..')
# this module define the rules for constructing a energy block in the master block
# this is the global component set import, so that all modules uses the same set
from global_sets.component import m
from data import kinetic_data as k
from pyomo import environ as pe

# define MPCC_beta_NCP rule
def P_NCP_block_rule(block):

    #------------------------------LOCAL VARIABLES------------------------------
    block.s_L = pe.Var(within=pe.NonNegativeReals,initialize=0)#,bounds=(0,1))
    block.s_V = pe.Var(within=pe.NonNegativeReals,initialize=0)

    #-----------------------------LOCAL parameters------------------------------
    block.epi = pe.Param(initialize=1e-3,mutable=True)

    print('>','Importing MPCC_P_NCP Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)
    for i in block.component_objects(pe.Param,active=True):
        print('|',i)

    print('-'*50)

    print('>','Adding complementarity constraint, spliting pressure used in VLE')

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
    if block.parent_block().VLE_block.find_component('pressure_equal_con'):
        block.parent_block().VLE_block.del_component(block.parent_block().VLE_block.pressure_equal_con)
        print('> Deleted original P_equal constraint')
    else:
        print('> No constraint to delete')
    print('')

    def pressure_equal_rule(block):
        return block.parent_block().VLE_block.P_VLE - block.parent_block().P == block.s_L - block.s_V
    block.pressure_equal_con = pe.Constraint(rule=pressure_equal_rule)

    if 'kinetics_block' in [i.local_name for i in block.parent_block().block_data_objects()]:
        if block.parent_block().kinetics_block.find_component('f_V_MPCC'):
            print('> Already replaced f_V_MPCC')
        else:
            block.parent_block().kinetics_block.del_component(block.parent_block().kinetics_block.r_FT_total_con)
            block.parent_block().kinetics_block.del_component(block.parent_block().kinetics_block.r_WGS_con)

            print('> Deleted kinetics rates constraints')

            block.parent_block().kinetics_block.f_V_MPCC = pe.Var(m.COMP_INORG,within=pe.NonNegativeReals,initialize=1e-20)

            def f_V_MPCC_rule(block,i):
                return block.f_V_MPCC[i] == block.parent_block().P * block.parent_block().y[i]
            block.parent_block().kinetics_block.f_V_MPCC_con = pe.Constraint(m.COMP_INORG,rule=f_V_MPCC_rule)

            def r_FT_total_MPCC_rule(block):
                return block.r_FT_total * (1+k.c_FT*(block.f_V_MPCC['CO'])**0.5+k.d_FT*(block.f_V_MPCC['H2'])**0.5)**2 == \
                (block.parent_block().cat*block.k_FT*block.f_V_MPCC['CO']**0.5*block.f_V_MPCC['H2']**0.5)
            block.parent_block().kinetics_block.r_FT_total_con = pe.Constraint(rule=r_FT_total_MPCC_rule)

            def r_WGS_MPCC_rule(block):
                return block.r_WGS * (block.Ke_WGS*block.f_V_MPCC['H2O'])\
                 == 0.5*block.parent_block().cat*block.k_WGS*block.parent_block().P**0.75*(block.f_V_MPCC['CO']* \
                            (block.Ke_WGS*block.f_V_MPCC['H2O']) - (block.f_V_MPCC['H2']*block.f_V_MPCC['CO2']) )
            block.parent_block().kinetics_block.r_WGS_con = pe.Constraint(rule=r_WGS_MPCC_rule)

            print('> Added f_V_MPCC, rates constraints')

    else:
        print('> Find no kinetics_block')
    print('')

# define MPCC_beta_Reg rule
def P_Reg_block_rule(block):

    #------------------------------LOCAL VARIABLES------------------------------
    block.s_L = pe.Var(within=pe.NonNegativeReals,initialize=0)#,bounds=(0,1))
    block.s_V = pe.Var(within=pe.NonNegativeReals,initialize=0)

    #-----------------------------LOCAL parameters------------------------------
    block.epi = pe.Param(initialize=1e-3,mutable=True)

    print('>','Importing MPCC_P_Reg Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)
    for i in block.component_objects(pe.Param,active=True):
        print('|',i)

    print('-'*50)

    print('>','Adding complementarity constraint, spliting pressure used in VLE')

    #------------------------------MPCC equations-------------------------------
    def s_L_complementarity_rule(block):
        return sum(block.parent_block().L[s] for s in block.parent_block().outlet) * block.s_L \
        <= block.epi
    block.s_L_complementarity_con = pe.Constraint(rule=s_L_complementarity_rule)

    def s_V_complementarity_rule(block):
        return sum(block.parent_block().V[s] for s in block.parent_block().outlet) * block.s_V <= block.epi
    block.s_V_complementarity_con = pe.Constraint(rule=s_V_complementarity_rule)

    #-----------------------------Global equations------------------------------
    if block.parent_block().VLE_block.find_component('pressure_equal_con'):
        block.parent_block().VLE_block.del_component(block.parent_block().VLE_block.pressure_equal_con)
        print('Deleted original P_equal constraint')
    else:
        print('No constraint to delete')
    print('')

    def pressure_equal_rule(block):
        return block.parent_block().VLE_block.P_VLE - block.parent_block().P == block.s_L - block.s_V
    block.pressure_equal_con = pe.Constraint(rule=pressure_equal_rule)

    if 'kinetics_block' in [i.local_name for i in block.parent_block().block_data_objects()]:
        if block.parent_block().kinetics_block.find_component('f_V_MPCC'):
            print('> Already replaced f_V_MPCC')
        else:
            block.parent_block().kinetics_block.del_component(block.parent_block().kinetics_block.r_FT_total_con)
            block.parent_block().kinetics_block.del_component(block.parent_block().kinetics_block.r_WGS_con)

            print('> Deleted kinetics rates constraints')

            block.parent_block().kinetics_block.f_V_MPCC = pe.Var(m.COMP_INORG,within=pe.NonNegativeReals,initialize=1e-20)

            def f_V_MPCC_rule(block,i):
                return block.f_V_MPCC[i] == block.parent_block().P * block.parent_block().y[i]
            block.parent_block().kinetics_block.f_V_MPCC_con = pe.Constraint(m.COMP_INORG,rule=f_V_MPCC_rule)

            def r_FT_total_MPCC_rule(block):
                return block.r_FT_total * (1+k.c_FT*(block.f_V_MPCC['CO'])**0.5+k.d_FT*(block.f_V_MPCC['H2'])**0.5)**2 == \
                (block.parent_block().cat*block.k_FT*block.f_V_MPCC['CO']**0.5*block.f_V_MPCC['H2']**0.5)
            block.parent_block().kinetics_block.r_FT_total_con = pe.Constraint(rule=r_FT_total_MPCC_rule)

            def r_WGS_MPCC_rule(block):
                return block.r_WGS * (block.Ke_WGS*block.f_V_MPCC['H2O'])\
                 == 0.5*block.parent_block().cat*block.k_WGS*block.parent_block().P**0.75*(block.f_V_MPCC['CO']* \
                            (block.Ke_WGS*block.f_V_MPCC['H2O']) - (block.f_V_MPCC['H2']*block.f_V_MPCC['CO2']) )
            block.parent_block().kinetics_block.r_WGS_con = pe.Constraint(rule=r_WGS_MPCC_rule)

            print('> Added f_V_MPCC, rates constraints')

    else:
        print('> Find no kinetics_block')
    print('')

# define MPCC_beta_penalty_function rule
def P_pf_block_rule(block):

    #------------------------------LOCAL VARIABLES------------------------------
    block.s_L = pe.Var(within=pe.NonNegativeReals,initialize=0)#,bounds=(0,1))
    block.s_V = pe.Var(within=pe.NonNegativeReals,initialize=0)
    block.pf = pe.Var(within=pe.NonNegativeReals,initialize=0)
    block.epi = pe.Param(initialize=1e-3,mutable=True)
    #-----------------------------LOCAL parameters------------------------------
    block.rho = pe.Param(initialize=100000,mutable=True)

    print('>','Importing MPCC_P_pf Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)
    for i in block.component_objects(pe.Param,active=True):
        print('|',i)

    print('-'*50)
    print('>','Spliting pressure used in VLE')

    #------------------------------MPCC equations-------------------------------
    def penalty_rule(block):
        return block.pf >= block.rho*(sum(block.parent_block().L[s] for s in block.parent_block().outlet)*block.s_L \
        + sum(block.parent_block().V[s] for s in block.parent_block().outlet)*block.s_V)
    block.penalty_con = pe.Constraint(rule=penalty_rule)

    #-----------------------------Global equations------------------------------
    if block.parent_block().VLE_block.find_component('pressure_equal_con'):
        block.parent_block().VLE_block.del_component(block.parent_block().VLE_block.pressure_equal_con)
        print('Deleted original P_equal constraint')
    else:
        print('No constraint to delete')
    print('')

    def pressure_equal_rule(block):
        return block.parent_block().VLE_block.P_VLE - block.parent_block().P == block.s_L - block.s_V
    block.pressure_equal_con = pe.Constraint(rule=pressure_equal_rule)

    if 'kinetics_block' in [i.local_name for i in block.parent_block().block_data_objects()]:
        if block.parent_block().kinetics_block.find_component('f_V_MPCC'):
            print('> Already replaced f_V_MPCC')
        else:
            block.parent_block().kinetics_block.del_component(block.parent_block().kinetics_block.r_FT_total_con)
            block.parent_block().kinetics_block.del_component(block.parent_block().kinetics_block.r_WGS_con)

            print('> Deleted kinetics rates constraints')

            block.parent_block().kinetics_block.f_V_MPCC = pe.Var(m.COMP_INORG,within=pe.NonNegativeReals,initialize=1e-20)

            def f_V_MPCC_rule(block,i):
                return block.f_V_MPCC[i] == block.parent_block().P * block.parent_block().y[i]
            block.parent_block().kinetics_block.f_V_MPCC_con = pe.Constraint(m.COMP_INORG,rule=f_V_MPCC_rule)

            def r_FT_total_MPCC_rule(block):
                return block.r_FT_total * (1+k.c_FT*(block.f_V_MPCC['CO'])**0.5+k.d_FT*(block.f_V_MPCC['H2'])**0.5)**2 == \
                (block.parent_block().cat*block.k_FT*block.f_V_MPCC['CO']**0.5*block.f_V_MPCC['H2']**0.5)
            block.parent_block().kinetics_block.r_FT_total_con = pe.Constraint(rule=r_FT_total_MPCC_rule)

            def r_WGS_MPCC_rule(block):
                return block.r_WGS * (block.Ke_WGS*block.f_V_MPCC['H2O'])\
                 == 0.5*block.parent_block().cat*block.k_WGS*block.parent_block().P**0.75*(block.f_V_MPCC['CO']* \
                            (block.Ke_WGS*block.f_V_MPCC['H2O']) - (block.f_V_MPCC['H2']*block.f_V_MPCC['CO2']) )
            block.parent_block().kinetics_block.r_WGS_con = pe.Constraint(rule=r_WGS_MPCC_rule)

            print('> Added f_V_MPCC, updated rates constraints')

    else:
        print('> Find no kinetics_block')
    print('')

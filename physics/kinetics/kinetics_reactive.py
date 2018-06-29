# 1st level Model Structure: Equation Block

# this module define the rules for constructing a kinetics block in the master block
# this is the global component set import, so that all modules uses the same set
from global_sets.component import m
from physics.bounds import kinetic_bounds

# data import
from data import kinetic_data as k
from utility.data_utility import cal_op
from pyomo import environ as pe

# import mean
from statistics import mean

# pre-processing
op_ratio = cal_op(k.op_ratio)
h2_consumption = [(2*(i+1)+2)/2*op_ratio[i] + (2*(i+1))/2*(1-op_ratio[i]) for i in range(len(op_ratio))]

# defile knietic block rule
def kinetic_block_rule(block):
    #-----------------------------------SETS-----------------------------------

    # local sets that will only be used in kinetics model
    block.C_NUMBER = pe.RangeSet(1,56)

    #-----------------------------GLOBAL VARIABLES-----------------------------

    # global variables
    # print('\t'*2,'Importing Kinetics Block......')
    # print('\t'*2,'Using the following parent variable:')
    # print('\t'*2,'-'*36)
    # print('\t'*2,block.parent_block().T.name)
    # print('\t'*2,block.parent_block().P.name)
    # print('\t'*2,block.parent_block().cat.name)
    # print('\t'*2,block.parent_block().f_V.name)
    # print('\t'*2,block.parent_block().r_total_comp.name)
    # print('\t'*2,'-'*36)
    # print('')

    #-----------------------------VARIABLES Bounds------------------------------
    def k_FT_bounds(model):
        lower = min(kinetic_bounds['k_FT'])
        lower = lower - abs(lower)*0.1
        upper = max(kinetic_bounds['k_FT'])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def g0_FT_bounds(model):
        lower = min(kinetic_bounds['g0_FT'])
        lower = lower - abs(lower)*0.1
        upper = max(kinetic_bounds['g0_FT'])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def alpha_bounds(model):
        lower = min(kinetic_bounds['alpha'])
        lower = lower - abs(lower)*0.1
        upper = max(kinetic_bounds['alpha'])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def k_WGS_bounds(model):
        lower = min(kinetic_bounds['k_WGS'])
        lower = lower - abs(lower)*0.1
        upper = max(kinetic_bounds['k_WGS'])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def Ke_WGS_bounds(model):
        lower = min(kinetic_bounds['Ke_WGS'])
        lower = lower - abs(lower)*0.1
        upper = max(kinetic_bounds['Ke_WGS'])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    #------------------------------LOCAL VARIABLES------------------------------

    # FT Reaction
    block.k_FT = pe.Var(within=pe.PositiveReals,initialize=5e-4,bounds=k_FT_bounds)
    block.r_FT_total = pe.Var(within=pe.PositiveReals)
    block.g0_FT = pe.Var(within=pe.Reals,bounds=g0_FT_bounds)
    block.alpha = pe.Var(within=pe.PositiveReals,initialize=0.7,bounds=alpha_bounds)
    block.r_FT_cnum = pe.Var(block.C_NUMBER,within=pe.PositiveReals)
    block.r_FT_comp = pe.Var(m.COMP_TOTAL,within=pe.Reals)         # kmol/s

    # WGS Reaction
    block.k_WGS = pe.Var(within=pe.PositiveReals,initialize=2e-4,bounds=k_WGS_bounds)
    block.Ke_WGS = pe.Var(within=pe.PositiveReals,initialize=3,bounds=Ke_WGS_bounds)
    block.r_WGS = pe.Var(within=pe.PositiveReals)
    block.r_WGS_comp = pe.Var(m.COMP_INORG,within=pe.Reals)         # kmol/s

    # initialize these variable: 1/2(ub+lb)
    block.k_FT = mean(kinetic_bounds['k_FT'])
    block.g0_FT = mean(kinetic_bounds['g0_FT'])
    block.alpha = mean(kinetic_bounds['alpha'])
    block.k_WGS = mean(kinetic_bounds['k_WGS'])
    block.Ke_WGS = mean(kinetic_bounds['Ke_WGS'])

    print('>','Importing Kinetics Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)

    print('-'*50)
    print('')


    #---------------------------------Equations---------------------------------

    # FT Reaction
    def r_FT_total_rule(block):
        return block.r_FT_total * (1+k.c_FT*(block.parent_block().f_V['CO'])**0.5+k.d_FT*(block.parent_block().f_V['H2'])**0.5)**2 == \
            (block.parent_block().cat*block.k_FT*block.parent_block().f_V['CO']**0.5*block.parent_block().f_V['H2']**0.5)
    block.r_FT_total_con = pe.Constraint(rule=r_FT_total_rule)

    def k_FT_rule(block):
        return block.k_FT == k.k0_FT * pe.exp(-k.E_FT*1e3/(k.R*block.parent_block().T))
    block.k_KT_con = pe.Constraint(rule=k_FT_rule)

    # ASF Distribution
    def r_FT_cnum_rule1(block,i):
        if i == 1:
            return pe.Constraint.Skip
        else:
            return block.r_FT_cnum[i] == block.alpha*block.r_FT_cnum[i-1]
    block.r_FT_cnum_con1 = pe.Constraint(block.C_NUMBER, rule=r_FT_cnum_rule1)

    def r_FT_cnum_rule2(block):
        return sum(i*block.r_FT_cnum[i] for i in block.C_NUMBER) == block.r_FT_total
    block.r_FT_cnum_con2 = pe.Constraint(rule=r_FT_cnum_rule2)

    def g0_FT_rule(block):
        return block.g0_FT == k.g0_inter_FT + k.g0_slope_FT * block.parent_block().T
    block.g0_FT_con = pe.Constraint(rule=g0_FT_rule)

    def alpha_rule(block):
        return block.alpha**2 == (1-block.alpha) * pe.exp(block.g0_FT*1e3/(k.R*block.parent_block().T))
    block.alpha_con = pe.Constraint(rule=alpha_rule)

    # Apply O/P Ratio
    def r_FT_para_rule(block,i):
        k = m.COMP_PARAFFIN.ord(i)
        return block.r_FT_comp[i] == op_ratio[k-1] * block.r_FT_cnum[k]
    block.r_FT_para_con = pe.Constraint(m.COMP_PARAFFIN, rule=r_FT_para_rule)

    def r_FT_ole_rule(block,i):
        k = m.COMP_OLEFIN.ord(i)+1
        return block.r_FT_comp[i] == (1-op_ratio[k-1]) * block.r_FT_cnum[k]
    block.r_FT_ole_con = pe.Constraint(m.COMP_OLEFIN, rule=r_FT_ole_rule)

    def r_FT_inorg_rule(block,i):
        if i == 'CO': return block.r_FT_comp[i] == -block.r_FT_total
        if i == 'H2O': return block.r_FT_comp[i] == block.r_FT_total
        if i == 'H2': return block.r_FT_comp[i] == -sum(h2_consumption[n-1]*block.r_FT_cnum[n] for n in block.C_NUMBER) - block.r_FT_total
        if i == 'CO2': return block.r_FT_comp[i] == 0
    block.r_FT_inorg_con = pe.Constraint(m.COMP_INORG, rule=r_FT_inorg_rule)

    # WGS Reaction
    def r_WGS_rule(block):
        return block.r_WGS * (block.Ke_WGS*block.parent_block().f_V['H2O']) == 0.5*block.parent_block().cat*block.k_WGS*block.parent_block().P**0.75*(block.parent_block().f_V['CO']* \
                    (block.Ke_WGS*block.parent_block().f_V['H2O']) - (block.parent_block().f_V['H2']*block.parent_block().f_V['CO2']) )
    block.r_WGS_con = pe.Constraint(rule=r_WGS_rule)

    def k_WGS_rule(block):
        return block.k_WGS == k.k0_WGS * pe.exp(-k.E_WGS*1e3/(k.R*block.parent_block().T))
    block.k_WGS_con = pe.Constraint(rule=k_WGS_rule)

    def Ke_WGS_rule(block):
        return block.Ke_WGS == 5*pe.exp((k.s1_WGS/block.parent_block().T + k.s2_WGS + k.s3_WGS*block.parent_block().T + k.s4_WGS*(block.parent_block().T**2))/k.R) * \
                                          block.parent_block().T**(k.s5_WGS/k.R)
    block.Ke_WGS_con = pe.Constraint(rule=Ke_WGS_rule)

    def r_WGS_comp_rule(block,i):
        if i == 'CO' or i == 'H2O': return block.r_WGS_comp[i] == -block.r_WGS
        if i == 'CO2' or i == 'H2': return block.r_WGS_comp[i] == block.r_WGS
    block.r_WGS_comp_con = pe.Constraint(m.COMP_INORG, rule=r_WGS_comp_rule)

    # Combine both reactions
    def r_total_comp_rule(block,i):
        if i in m.COMP_INORG: return block.parent_block().r_total_comp[i] == block.r_FT_comp[i] + block.r_WGS_comp[i]
        return block.parent_block().r_total_comp[i] == block.r_FT_comp[i]
    block.r_total_comp_con = pe.Constraint(m.COMP_TOTAL,rule=r_total_comp_rule)

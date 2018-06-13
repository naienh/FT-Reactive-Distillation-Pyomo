# 1st level Model Structure: Equation Block

# this module define the rules for constructing a VLE block in the master block
# this is the global component set import, so that all modules uses the same set
from modules.global_set import m

# data import and pre-processing
from data import VLE_data as e
from data.utility import cal_cnumber
from pyomo import environ as pe

# defile knietic block rule
def VLE_block_rule(block):
    #-----------------------------------SETS-----------------------------------

    # local sets that will only be used in VLE block
    block.COMP_HENRY = pe.Set(initialize=['H2','CO','CO2','H2O','C1H4','C2H6','C3H8','C2H4','C3H6'])
    block.COMP_NONHENRY = m.COMP_TOTAL - block.COMP_HENRY

    #-----------------------------GLOBAL VARIABLES-----------------------------

    # global variables
    # print('\t'*2,'Importing VLE Block......')
    # print('\t'*2,'Using the following parent variable:')
    # print('\t'*2,'-'*36)
    # print('\t'*2,block.parent_block().T.name)
    # print('\t'*2,block.parent_block().P.name)
    # print('\t'*2,block.parent_block().x.name)
    # print('\t'*2,block.parent_block().y.name)
    # print('\t'*2,block.parent_block().f_L.name)
    # print('\t'*2,block.parent_block().f_V.name)
    # print('\t'*2,'-'*36)
    # print('')

    #------------------------------LOCAL VARIABLES------------------------------

    # teared n_ave, initial guess try to converge to calculated average
    block.n_ave = pe.Var(within=pe.NonNegativeReals,bounds=(7,40))
    block.n_ave_cal = pe.Var(within=pe.NonNegativeReals)

    # fugacity variable
    block.Hen = pe.Var(block.COMP_HENRY,within=pe.NonNegativeReals)  # Bar
    block.Hen0 = pe.Var(block.COMP_HENRY,within=pe.Reals)
    block.gamma = pe.Var(block.COMP_NONHENRY,within=pe.NonNegativeReals)
    block.P_sat = pe.Var(block.COMP_NONHENRY,within=pe.NonNegativeReals)  # Bar
    block.P_sat_dY_inf = pe.Var(within=pe.Reals)
    block.P_sat_dY0 = pe.Var(within=pe.Reals)

    block.Hen_ref = pe.Var(within=pe.NonNegativeReals)
    block.Hen0_ref = pe.Var(within=pe.Reals)
    block.gamma_ref = pe.Var(within=pe.NonNegativeReals)

    print('>','Importing VLE Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)

    print('-'*50)
    print('')

    #---------------------------------Equations---------------------------------

    # henry component
    def f_L_HENRY_rule(block,i):
        return block.parent_block().f_L[i] == block.Hen[i]*block.parent_block().x[i]
    block.f_L_HENRY_con = pe.Constraint(block.COMP_HENRY,rule=f_L_HENRY_rule)

    # henry's constant
    def Hen_rule(block,i):
        return block.Hen[i]*pe.exp(block.n_ave*e.henry.dHen[i]) == pe.exp(block.Hen0[i])
    block.Hen_con = pe.Constraint(block.COMP_HENRY,rule=Hen_rule)

    def Hen0_rule(block,i):
        return block.Hen0[i] == e.henry.A[i] + e.henry.B[i]/block.parent_block().T + e.henry.C[i]*pe.log(block.parent_block().T) + \
                    e.henry.D[i]*(block.parent_block().T**2) + e.henry.E[i]/(block.parent_block().T**2)
    block.Hen0_con = pe.Constraint(block.COMP_HENRY,rule=Hen0_rule)

    def Hen_ref_rule(block):
        return block.Hen_ref*pe.exp(block.n_ave*e.henry.dHen['C6H14']) == pe.exp(block.Hen0_ref)
    block.Hen_ref_con = pe.Constraint(rule=Hen_ref_rule)

    def Hen0_ref_rule(block):
        return block.Hen0_ref == e.henry.A['C6H14'] + e.henry.B['C6H14']/block.parent_block().T + e.henry.C['C6H14']*pe.log(block.parent_block().T) + \
                    e.henry.D['C6H14']*(block.parent_block().T**2) + e.henry.E['C6H14']/(block.parent_block().T**2)
    block.Hen0_ref_con = pe.Constraint(rule=Hen0_ref_rule)

    # non-henry component
    def f_L_NONHENRY_rule(block,i):
        return block.parent_block().f_L[i] == block.gamma[i]*block.P_sat[i]*block.parent_block().x[i]
    block.f_L_NONHENRY_con = pe.Constraint(block.COMP_NONHENRY,rule=f_L_NONHENRY_rule)

    # acticity coefficient gamma
    def gamma_rule(block,i):
        return pe.log(block.gamma[i])*(block.n_ave - cal_cnumber('C6H14')) == \
                pe.log(block.gamma_ref)*(block.n_ave - cal_cnumber(i))
    block.gamma_con = pe.Constraint(block.COMP_NONHENRY,rule=gamma_rule)

    def gamma_ref_rule(block):
        return block.gamma_ref*block.P_sat['C6H14'] == block.Hen_ref
    block.gamma_ref_con = pe.Constraint(rule=gamma_ref_rule)

    def n_ave_rule(block):
        return block.n_ave_cal == sum(block.parent_block().x[i]*cal_cnumber(i) for i in m.COMP_PARAFFIN | \
                m.COMP_OLEFIN)/(1-sum(block.parent_block().x[i] for i in m.COMP_INORG))
    block.n_ave_con = pe.Constraint(rule=n_ave_rule)

    # saturated pressure
    def P_sat_rule(block,i):
        if i in m.COMP_PARAFFIN:
            n_n0 = cal_cnumber(i)-e.nonhenry.n0_paraffin
        elif i in m.COMP_OLEFIN:
            n_n0 = cal_cnumber(i)-e.nonhenry.n0_olefin
        return pe.log(block.P_sat[i]) == e.nonhenry.Y_inf_0 + block.P_sat_dY_inf*(n_n0) \
                    - block.P_sat_dY0*pe.exp(-e.nonhenry.beta*(n_n0)**e.nonhenry.gamma)
    block.P_sat_con = pe.Constraint(block.COMP_NONHENRY,rule=P_sat_rule)

    def P_sat_dY_inf_rule(block):
        return e.nonhenry.dY_inf.A + e.nonhenry.dY_inf.B/block.parent_block().T + e.nonhenry.dY_inf.C*pe.log(block.parent_block().T) + \
                e.nonhenry.dY_inf.D*(block.parent_block().T)**2 + e.nonhenry.dY_inf.E/(block.parent_block().T)**2 == block.P_sat_dY_inf
    block.P_sat_dY_inf_con = pe.Constraint(rule=P_sat_dY_inf_rule)

    def P_sat_dY0_rule(block):
        return e.nonhenry.dY0.A + e.nonhenry.dY0.B/block.parent_block().T + e.nonhenry.dY0.C*pe.log(block.parent_block().T) + \
                e.nonhenry.dY0.D*(block.parent_block().T)**2 + e.nonhenry.dY0.E/(block.parent_block().T)**2 == block.P_sat_dY0
    block.P_sat_dY0_con = pe.Constraint(rule=P_sat_dY0_rule)

    # gas phase assume ideal
    def f_V_rule(block,i):
        return block.parent_block().P*block.parent_block().y[i] == block.parent_block().f_V[i]
    block.f_V_con = pe.Constraint(m.COMP_TOTAL,rule=f_V_rule)

# 1st level Model Structure: Equation Block

# this module define the rules for constructing a VLE block in the master block
# this is the global component set import, so that all modules uses the same set
from global_sets.component import m
from physics.bounds import VLE_bounds

# data import and pre-processing
from data import VLE_data as e
from utility.data_utility import cal_cnumber
from pyomo import environ as pe

# import mean
from statistics import mean

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

    #-----------------------------VARIABLES Bounds------------------------------
    def Hen_bounds(model,i):
        lower = min(VLE_bounds['Hen[{}]'.format(i)])
        lower = lower - abs(lower)*0.1
        upper = max(VLE_bounds['Hen[{}]'.format(i)])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def Hen0_bounds(model,i):
        lower = min(VLE_bounds['Hen0[{}]'.format(i)])
        lower = lower - abs(lower)*0.1
        upper = max(VLE_bounds['Hen0[{}]'.format(i)])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def gamma_bounds(model,i):
        lower = min(VLE_bounds['gamma[{}]'.format(i)])
        lower = lower - abs(lower)*0.1
        upper = max(VLE_bounds['gamma[{}]'.format(i)])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def P_sat_bounds(model,i):
        lower = min(VLE_bounds['P_sat[{}]'.format(i)])
        lower = lower - abs(lower)*0.1
        upper = max(VLE_bounds['P_sat[{}]'.format(i)])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def P_sat_Y_bounds(model,i):
        lower = min(VLE_bounds['P_sat_Y[{}]'.format(i)])
        lower = lower - abs(lower)*0.1
        upper = max(VLE_bounds['P_sat_Y[{}]'.format(i)])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def P_sat_dY_inf_bounds(model):
        lower = min(VLE_bounds['P_sat_dY_inf'])
        lower = lower - abs(lower)*0.1
        upper = max(VLE_bounds['P_sat_dY_inf'])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def P_sat_dY0_bounds(model):
        lower = min(VLE_bounds['P_sat_dY0'])
        lower = lower - abs(lower)*0.1
        upper = max(VLE_bounds['P_sat_dY0'])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def Hen_ref_bounds(model):
        lower = min(VLE_bounds['Hen_ref'])
        lower = lower - abs(lower)*1
        upper = max(VLE_bounds['Hen_ref'])
        upper = upper + abs(upper)*1
        return (lower,upper)

    def Hen0_ref_bounds(model):
        lower = min(VLE_bounds['Hen0_ref'])
        lower = lower - abs(lower)*1
        upper = max(VLE_bounds['Hen0_ref'])
        upper = upper + abs(upper)*1
        return (lower,upper)

    def gamma_ref_bounds(model):
        lower = min(VLE_bounds['gamma_ref'])
        lower = lower - abs(lower)*1
        upper = max(VLE_bounds['gamma_ref'])
        upper = upper + abs(upper)*1
        return (lower,upper)

    def V_L_bounds(model,i):
        lower = min(VLE_bounds['V_L[{}]'.format(i)])
        lower = lower - abs(lower)*0.1
        upper = max(VLE_bounds['V_L[{}]'.format(i)])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def V_L_dY_inf_bounds(model):
        lower = min(VLE_bounds['V_L_dY_inf'])
        lower = lower - abs(lower)*0.1
        upper = max(VLE_bounds['V_L_dY_inf'])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def V_L_dY0_bounds(model):
        lower = min(VLE_bounds['V_L_dY0'])
        lower = lower - abs(lower)*0.5
        upper = max(VLE_bounds['V_L_dY0'])
        upper = upper + abs(upper)*1
        return (lower,upper)

    def poynting_bounds(model,i):
        lower = min(VLE_bounds['poynting[{}]'.format(i)])
        lower = lower - abs(lower)*0.1
        upper = max(VLE_bounds['poynting[{}]'.format(i)])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    #------------------------------LOCAL VARIABLES------------------------------

    # teared n_ave, initial guess try to converge to calculated average
    block.n_ave = pe.Var(within=pe.NonNegativeReals,bounds=(7,58))
    block.n_ave_cal = pe.Var(within=pe.NonNegativeReals)

    # fugacity variable
    block.Hen = pe.Var(block.COMP_HENRY,within=pe.NonNegativeReals,bounds=Hen_bounds)  # Bar
    block.Hen0 = pe.Var(block.COMP_HENRY,within=pe.Reals,initialize=5,bounds=Hen0_bounds)
    block.gamma = pe.Var(block.COMP_NONHENRY,within=pe.PositiveReals,initialize=0.1,bounds=gamma_bounds)
    block.P_sat = pe.Var(block.COMP_NONHENRY,within=pe.PositiveReals,initialize=1e-10,bounds=P_sat_bounds)  # Bar
    block.P_sat_Y = pe.Var(block.COMP_NONHENRY,within=pe.Reals,bounds=P_sat_Y_bounds)
    block.P_sat_dY_inf = pe.Var(within=pe.Reals,bounds=P_sat_dY_inf_bounds)
    block.P_sat_dY0 = pe.Var(within=pe.Reals,bounds=P_sat_dY0_bounds)

    # block.Hen_ref = pe.Var(within=pe.NonNegativeReals,initialize=14,bounds=Hen_ref_bounds)
    # block.Hen0_ref = pe.Var(within=pe.Reals,initialize=3.6,bounds=Hen0_ref_bounds)
    # block.gamma_ref = pe.Var(within=pe.PositiveReals,initialize=0.2,bounds=gamma_ref_bounds)

    block.Hen_ref = pe.Var(within=pe.NonNegativeReals,initialize=14)
    block.Hen0_ref = pe.Var(within=pe.Reals,initialize=3.6)
    block.gamma_ref = pe.Var(within=pe.PositiveReals,initialize=0.2)

    # molar volume variable
    block.V_L = pe.Var(m.COMP_TOTAL,within=pe.Reals,bounds=V_L_bounds) # cm3/mole
    block.V_L_dY_inf = pe.Var(within=pe.Reals,bounds=V_L_dY_inf_bounds)
    block.V_L_dY0 = pe.Var(within=pe.Reals,bounds=V_L_dY0_bounds)

    # poynting facotor variable
    block.poynting = pe.Var(m.COMP_TOTAL,within=pe.Reals,bounds=poynting_bounds)

    # initialize these variable: 1/2(ub+lb)
    for i in block.COMP_HENRY:
        block.Hen[i] = mean(VLE_bounds['Hen[{}]'.format(i)])
        block.Hen0[i] = mean(VLE_bounds['Hen0[{}]'.format(i)])

    # for i in block.COMP_NONHENRY:
        # block.gamma[i] = mean(VLE_bounds['gamma[{}]'.format(i)])
        # block.P_sat[i] = mean(VLE_bounds['P_sat[{}]'.format(i)])
        # block.P_sat_Y[i] = mean(VLE_bounds['P_sat_Y[{}]'.format(i)])
    #
    block.P_sat_dY_inf = mean(VLE_bounds['P_sat_dY_inf'])
    block.P_sat_dY0 = mean(VLE_bounds['P_sat_dY0'])
    block.V_L_dY_inf = mean(VLE_bounds['V_L_dY_inf'])
    block.V_L_dY0 = mean(VLE_bounds['V_L_dY0'])

    # block.Hen_ref = mean(VLE_bounds['Hen_ref'])
    # block.Hen0_ref = mean(VLE_bounds['Hen0_ref'])
    # block.gamma_ref = mean(VLE_bounds['gamma_ref'])

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
        return block.parent_block().f_L[i] == block.Hen[i]*block.parent_block().x[i]*block.poynting[i]
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
        return block.parent_block().f_L[i] == block.gamma[i]*block.P_sat[i]*block.parent_block().x[i]*block.poynting[i]
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
    def P_sat_rule1(block,i):
        if i in m.COMP_PARAFFIN:
            n_n0 = cal_cnumber(i)-e.nonhenry.n0_paraffin
        elif i in m.COMP_OLEFIN:
            n_n0 = cal_cnumber(i)-e.nonhenry.n0_olefin
        return block.P_sat_Y[i] == e.nonhenry.Y_inf_0 + block.P_sat_dY_inf*(n_n0) \
                    - block.P_sat_dY0*pe.exp(-e.nonhenry.beta*(n_n0)**e.nonhenry.gamma)
    block.P_sat_con = pe.Constraint(block.COMP_NONHENRY,rule=P_sat_rule1)

    def P_sat_rule2(block,i):
        return block.P_sat[i] == pe.exp(block.P_sat_Y[i])
    block.P_sat_con2 = pe.Constraint(block.COMP_NONHENRY,rule=P_sat_rule2)

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

    # molar volume Equation
    def V_L_nonHen_rule(block,i):
        if i in m.COMP_PARAFFIN:
            n_n0 = cal_cnumber(i)-e.V_L_nonhenry.n0_paraffin
            n_n0_ = cal_cnumber(i)+e.V_L_nonhenry.n0_paraffin
        elif i in m.COMP_OLEFIN:
            n_n0 = cal_cnumber(i)-e.V_L_nonhenry.n0_olefin
            n_n0_ = cal_cnumber(i)+e.V_L_nonhenry.n0_olefin
        return block.V_L[i] == e.V_L_nonhenry.Y_inf_0 + block.V_L_dY_inf*(n_n0) \
                    - block.V_L_dY0*pe.exp(-e.V_L_nonhenry.beta*(n_n0_)**e.V_L_nonhenry.gamma)
    block.V_L_nonHen_con = pe.Constraint(block.COMP_NONHENRY | {'C3H8','C3H6'},rule=V_L_nonHen_rule)

    def V_L_dY_inf_rule(block):
        return e.V_L_nonhenry.dY_inf.A + e.V_L_nonhenry.dY_inf.B*block.parent_block().T + e.V_L_nonhenry.dY_inf.C*(block.parent_block().T)**2 + \
                e.V_L_nonhenry.dY_inf.D*(block.parent_block().T)**3 == block.V_L_dY_inf
    block.V_L_dY_inf_con = pe.Constraint(rule=V_L_dY_inf_rule)

    def V_L_dY0_rule(block):
        return e.V_L_nonhenry.dY0.A + e.V_L_nonhenry.dY0.B*block.parent_block().T + e.V_L_nonhenry.dY0.C*(block.parent_block().T)**2 + \
                e.V_L_nonhenry.dY0.D*(block.parent_block().T)**3 == block.V_L_dY0
    block.V_L_dY0_con = pe.Constraint(rule=V_L_dY0_rule)

    def V_L_Hen_rule(block,i):
        return block.V_L[i] == e.V_L_henry.A[i] + e.V_L_henry.B[i]*block.parent_block().T + block.n_ave*e.V_L_henry.dV[i]
    block.V_L_Hen_con = pe.Constraint(block.COMP_HENRY - {'C3H8','C3H6','H2O'},rule=V_L_Hen_rule)

    def V_L_H2O_rule(block):
        return block.V_L['H2O'] == e.V_L_water.A + e.V_L_water.B*block.parent_block().T + e.V_L_water.C*(block.parent_block().T)**2 \
                    + e.V_L_water.D*(block.parent_block().T)**3
    block.V_L_H2O_con = pe.Constraint(rule=V_L_H2O_rule)

    # poynting factor equation
    def poynting_rule(block,i):
        return block.poynting[i] == pe.exp(0.1*block.V_L[i]*(block.parent_block().P)/(e.R*block.parent_block().T))
    block.poynting_con = pe.Constraint(m.COMP_TOTAL,rule=poynting_rule)

# 1st level Model Structure: Equation Block
# import sys
# sys.path.append('..')
# this module define the rules for constructing a energy block in the master block
# this is the global component set import, so that all modules uses the same set
from global_sets.component import m
from physics.bounds import energy_bounds3 as energy_bounds

# data import and pre-processing
from data import thermal_data as h
from pyomo import environ as pe

# import mean
from statistics import mean

# defile knietic block rule
def energy_block_rule(block):
    #-----------------------------------SETS-----------------------------------

    # No local SETS

    #-----------------------------GLOBAL VARIABLES-----------------------------

    # global variables
    # print('\t'*2,'Importing Energy Block......')
    # print('\t'*2,'Using the following parent variable:')
    # print('\t'*2,'-'*36)
    # print('\t'*2,block.parent_block().T.name)
    # print('\t'*2,block.parent_block().T_F.name)
    # print('\t'*2,block.parent_block().x.name)
    # print('\t'*2,block.parent_block().y.name)
    # print('\t'*2,block.parent_block().z.name)
    # print('\t'*2,block.parent_block().H_F.name)
    # print('\t'*2,block.parent_block().H_V.name)
    # print('\t'*2,block.parent_block().H_L.name)
    # print('\t'*2,'-'*36)
    # print('')
    #-----------------------------VARIABLES Bounds------------------------------
    def dH_V_bounds(model,i):
        lower = min(energy_bounds['dH_V[{}]'.format(i)])
        lower = lower - abs(lower)*0.1
        upper = max(energy_bounds['dH_V[{}]'.format(i)])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def dH_L_bounds(model,i):
        lower = min(energy_bounds['dH_L[{}]'.format(i)])
        lower = lower - abs(lower)*0.1
        upper = max(energy_bounds['dH_L[{}]'.format(i)])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    def dH_vap_bounds(model,i):
        lower = min(energy_bounds['dH_vap[{}]'.format(i)])
        lower = lower - abs(lower)*0.1
        upper = max(energy_bounds['dH_vap[{}]'.format(i)])
        upper = upper + abs(upper)*0.1
        return (lower,upper)

    #------------------------------LOCAL VARIABLES------------------------------

    # Molar Enthalpy for gas, liquid and vaporization at temperature
    block.dH_F = pe.Var(m.COMP_FEED,within=pe.Reals)  # kJ/mol
    block.dH_V = pe.Var(m.COMP_TOTAL,within=pe.Reals,bounds=dH_V_bounds)
    block.dH_L = pe.Var(m.COMP_TOTAL,within=pe.Reals,bounds=dH_L_bounds)
    block.dH_vap = pe.Var(m.COMP_TOTAL,within=pe.Reals,bounds=dH_vap_bounds)

    # initialize these variable: 1/2(ub+lb)
    for i in m.COMP_TOTAL:
        block.dH_V[i] = mean(energy_bounds['dH_V[{}]'.format(i)])
        block.dH_L[i] = mean(energy_bounds['dH_L[{}]'.format(i)])
        block.dH_vap[i] = mean(energy_bounds['dH_vap[{}]'.format(i)])

    print('>','Importing Energy Blocks......')
    print('>','Adding the following local variable:')
    print('-'*50)

    for i in block.component_objects(pe.Var,active=True):
        print('|',i)

    print('-'*50)
    print('')

    #---------------------------------Equations---------------------------------

    # Vapor
    def dH_V_rule(block,i):
        return block.dH_V[i] == h.Hf0[i] + 1e-3*(h.a0[i]*(block.parent_block().T-h.T0_f) + 1/2*h.a1[i]*(block.parent_block().T**2-h.T0_f**2) + \
                1/3*h.a2[i]*(block.parent_block().T**3-h.T0_f**3) + 1/4*h.a3[i]*(block.parent_block().T**4-h.T0_f**4) + 1/5*h.a4[i]*(block.parent_block().T**5-h.T0_f**5))
    block.dH_V_con = pe.Constraint(m.COMP_TOTAL,rule=dH_V_rule)

    def H_V_rule(block):
        return block.parent_block().H_V == sum(block.parent_block().y[i] * block.dH_V[i] for i in m.COMP_TOTAL)
    block.H_V_con = pe.Constraint(rule=H_V_rule)

    # Liquid
    def dH_vap_rule(block,i):
        tmp = (h.Tc[i]-block.parent_block().T)/(h.Tc[i]-h.Tb[i])
        return block.dH_vap[i] == h.HV[i] * (0.5*tmp + 0.5*(tmp**2+1e-6)**0.5)**0.38
    block.dH_vap_con = pe.Constraint(m.COMP_TOTAL,rule=dH_vap_rule)

    def dH_L_rule(block,i):
        return block.dH_L[i] == block.dH_V[i] - block.dH_vap[i]
    block.dH_L_con = pe.Constraint(m.COMP_TOTAL,rule=dH_L_rule)

    def H_L_rule(block):
        return block.parent_block().H_L == sum(block.parent_block().x[i] * block.dH_L[i] for i in m.COMP_TOTAL)
    block.H_L_con = pe.Constraint(rule=H_L_rule)

    # Feed
    def dH_F_rule(block,i):
        return block.dH_F[i] == h.Hf0[i] + 1e-3*(h.a0[i]*(block.parent_block().T_F-h.T0_f) + 1/2*h.a1[i]*(block.parent_block().T_F**2-h.T0_f**2) + \
                1/3*h.a2[i]*(block.parent_block().T_F**3-h.T0_f**3) + 1/4*h.a3[i]*(block.parent_block().T_F**4-h.T0_f**4) + 1/5*h.a4[i]*(block.parent_block().T_F**5-h.T0_f**5))
    block.dH_F_con = pe.Constraint(m.COMP_FEED,rule=dH_F_rule)

    def H_F_rule(block):
        return block.parent_block().H_F == sum(block.parent_block().z[i] * block.dH_F[i] for i in m.COMP_FEED)
    block.H_F_con = pe.Constraint(rule=H_F_rule)

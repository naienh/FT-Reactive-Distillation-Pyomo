from pyomo import environ as pe
from data.utility import HiddenPrints
from modules.global_set import m

from modules.kinetics import kinetic_block_rule

model = pe.ConcreteModel()

model.P = pe.Var(within=pe.NonNegativeReals,bounds=(10,30)) # Bar
model.cat = pe.Var(within=pe.NonNegativeReals) # kg

model.T = pe.Var(within=pe.NonNegativeReals,bounds=(200+273.15,300+273.15)) # K
model.f_V = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals)

model.r_total_comp = pe.Var(m.COMP_TOTAL,within=pe.Reals) # kmol/s

with HiddenPrints():
    model.kinetics_block = pe.Block(rule=kinetic_block_rule)

model.obj = pe.Objective(expr = 1,sense=pe.maximize)
opt = pe.SolverFactory('ipopt')

model.f_V.fix(0)

def u_cal_R(T,P,cat,h2,co,h2o,co2):
    model.P.fix(P)
    model.cat.fix(cat)
    model.T.fix(T)
    model.f_V['H2'].fix(h2)
    model.f_V['CO'].fix(co)
    model.f_V['H2O'].fix(h2o)
    model.f_V['CO2'].fix(co2)
    opt.solve(model,tee=False)
    return {i:model.r_total_comp[i].value for i in m.COMP_TOTAL}

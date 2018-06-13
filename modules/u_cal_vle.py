from pyomo import environ as pe
# from data.utility import HiddenPrints, check_DOF
import data.utility as u
from modules.global_set import m

from modules.VLE import VLE_block_rule

model = pe.ConcreteModel()

model.P = pe.Var(within=pe.NonNegativeReals,bounds=(10,30)) # Bar
model.x = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals,bounds=(0,1))
model.y = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals,bounds=(0,1))
model.T = pe.Var(within=pe.NonNegativeReals,bounds=(20+273.15,350+273.15)) # K
model.f_V = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals)
model.f_L = pe.Var(m.COMP_TOTAL,within=pe.NonNegativeReals)


with u.HiddenPrints():
    model.VLE_block = pe.Block(rule=VLE_block_rule)

model.obj = pe.Objective(expr = 1,sense=pe.maximize)
opt = pe.SolverFactory('ipopt')

model.x.fix(0.1)
# model.y.fix(1)

model.VLE_block.n_ave_con.deactivate()

def u_cal_vle(T,P,n_ave):
    model.T.fix(T)
    model.P.fix(P)
    model.VLE_block.n_ave.fix(n_ave)
    opt.solve(model,tee=False)
    return {i:model.f_L[i].value/P for i in m.COMP_TOTAL}

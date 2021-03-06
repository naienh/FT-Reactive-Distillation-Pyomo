# a file hosts all global referred parameters/sets/etcself.
from pyomo import environ as pe
from data.thermal_data import Tb

m = pe.ConcreteModel()

m.COMP_OLEFIN = pe.Set(initialize=['C{0}H{1}'.format(i,2*i) for i in range(2,21)],ordered=True)
m.COMP_PARAFFIN = pe.Set(initialize=['C{0}H{1}'.format(i,2*i+2) for i in range(1,57)],ordered=True)
m.COMP_INORG = pe.Set(initialize=['H2','CO','CO2','H2O'],ordered=True)
m.COMP_ORG = m.COMP_OLEFIN | m.COMP_PARAFFIN
m.COMP_TOTAL = m.COMP_INORG | m.COMP_OLEFIN | m.COMP_PARAFFIN

m.COMP_FEED = pe.Set(initialize=['H2','CO','C30H62'],ordered=True)
# m.COMP_FEED = m.COMP_INORG | m.COMP_OLEFIN | m.COMP_PARAFFIN

'''
Sort components based on boiling point
'''
product_boiling_C = [(i,Tb[i]) for i in m.COMP_ORG]
def pick_comp(x):
    return x[1]
m.COMP_SORTED_BP = pe.Set(initialize=[i for i,T in sorted(product_boiling_C,key =pick_comp )],ordered=True)

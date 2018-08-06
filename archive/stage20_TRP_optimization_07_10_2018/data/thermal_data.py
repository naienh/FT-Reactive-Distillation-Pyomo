# lets start with reading the excel file
# Cp = a0 + a1*T + a2*T^2 + a3*T^3 + a4*T^4 (J/mol/K)
# Hf0 (kJ/mol)
# Gf0 (kJ/mol)
# dHv (kJ/mol)

import xlrd
from utility.data_utility import readcol
try:
    workbook = xlrd.open_workbook("data/thermo_data.xlsx")
except:
    try:
        workbook = xlrd.open_workbook("../../data/thermo_data.xlsx")
    except:
        workbook = xlrd.open_workbook("../data/thermo_data.xlsx")
sheet = workbook.sheet_by_index(0)
data_name = [sheet.cell_value(0,i) for i in range(sheet.ncols)]

def readonecol(name):
    return readcol(sheet,name)

Tb = readonecol('Tb')
Tc = readonecol('Tc')
Hf0 = readonecol('Hf0')
HV = readonecol('HV')
a0 = readonecol('a0')
a1 = readonecol('a1')
a2 = readonecol('a2')
a3 = readonecol('a3')
a4 = readonecol('a4')

T0_f = 298.15 # K

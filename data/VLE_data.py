# lets start with reading the excel file
# Henry H (bar)
# ln (H) = H0 - n_ave*dHi
# dHi = A + B/T + C*ln(T)+DT^2 + E/T^2 (T: K)

import xlrd
from data.utility import readcol
try:
    workbook = xlrd.open_workbook("data/VLE_data.xlsx")
except:
    try:
        workbook = xlrd.open_workbook("VLE_data.xlsx")
    except:
        workbook = xlrd.open_workbook("../data/VLE_data.xlsx")


class data_object(object):
    def __init__(self, name):
        self.name = name

sheet_Henry = workbook.sheet_by_index(0)
henry = data_object('henry')
henry.A = readcol(sheet_Henry,'A')
henry.B = readcol(sheet_Henry,'B')
henry.C = readcol(sheet_Henry,'C')
henry.D = readcol(sheet_Henry,'D')
henry.E = readcol(sheet_Henry,'E')
henry.dHen = readcol(sheet_Henry,'dHen')

nonhenry = data_object('nonhenry')
nonhenry.n0_paraffin = 1.126231
nonhenry.n0_olefin = 1.281405
nonhenry.Y_inf_0 = 2.72709
nonhenry.beta = 0.619226
nonhenry.gamma = 0.416321
nonhenry.dY0 = data_object('dY0')
nonhenry.dY0.A = -5.75509
nonhenry.dY0.B = -7.56568
nonhenry.dY0.C = 0.0857734
nonhenry.dY0.D = -1.41964e-5
nonhenry.dY0.E = 2.67209e5
nonhenry.dY_inf = data_object('dY_inf')
nonhenry.dY_inf.A = 15.8059
nonhenry.dY_inf.B = -1496.56
nonhenry.dY_inf.C = -2.17342
nonhenry.dY_inf.D = 7.27763e-7
nonhenry.dY_inf.E = 37876.2

epi=1e-6

R = 8.3144598 # J/K/mol

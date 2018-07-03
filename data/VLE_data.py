# lets start with reading the excel file
# Henry H (bar)
# ln (H) = H0 - n_ave*dHi
# dHi = A + B/T + C*ln(T)+DT^2 + E/T^2 (T: K)

import xlrd
from utility.data_utility import readcol
try:
    workbook = xlrd.open_workbook("data/VLE_data.xlsx")
except:
    try:
        workbook = xlrd.open_workbook("../../data/VLE_data.xlsx")
    except:
        workbook = xlrd.open_workbook("../data/VLE_data.xlsx")


class data_object(object):
    def __init__(self, name):
        self.name = name

sheet_Henry = workbook.sheet_by_index(0)
henry = data_object('henry')
henry.A = readcol(sheet_Henry,'hen_A')
henry.B = readcol(sheet_Henry,'hen_B')
henry.C = readcol(sheet_Henry,'hen_C')
henry.D = readcol(sheet_Henry,'hen_D')
henry.E = readcol(sheet_Henry,'hen_E')
henry.dHen = readcol(sheet_Henry,'dHen')

V_L_henry = data_object('V_L_henry')
V_L_henry.dV = readcol(sheet_Henry,'dV')
V_L_henry.A = readcol(sheet_Henry,'V_A')
V_L_henry.B = readcol(sheet_Henry,'V_B')

V_L_water = data_object('V_L_water')
V_L_water.A = 17.8937
V_L_water.B = 0.00574149
V_L_water.C = -5.19630e-5
V_L_water.D = 1.12264e-7

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

V_L_nonhenry = data_object('V_L_nonhenry')
V_L_nonhenry.n0_paraffin = -1.388524
V_L_nonhenry.n0_olefin = -1.061318
V_L_nonhenry.Y_inf_0 = 0
V_L_nonhenry.beta = 5.519846
V_L_nonhenry.gamma = 0.0570632
V_L_nonhenry.dY0 = data_object('dY0')
V_L_nonhenry.dY0.A = 8592.30
V_L_nonhenry.dY0.B = -85.7292
V_L_nonhenry.dY0.C = 0.280284
V_L_nonhenry.dY0.D = -4.48451e-4
V_L_nonhenry.dY_inf = data_object('dY_inf')
V_L_nonhenry.dY_inf.A = 12.7924
V_L_nonhenry.dY_inf.B = 0.0150627
V_L_nonhenry.dY_inf.C = -1.30794e-5
V_L_nonhenry.dY_inf.D = 1.59611e-8

epi=1e-6

R = 8.3144598 # J/K/mol

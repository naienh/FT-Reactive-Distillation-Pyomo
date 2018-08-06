# this contains the data necessary for modelling FT and WGS kinetics

epi=1e-6

R = 8.3144598 # J/K/mol
k0_FT = 3.8276e4 # kmol/kgcat/s/bar
E_FT = 79.9 # kJ/mol
c_FT = 0.242
d_FT = 0.185

k0_WGS = 8.0695e10 # kmol/kgcat/s/bar
E_WGS = 149.7 # kJ/mol
s1_WGS = 9998.22
s2_WGS = -10.213
s3_WGS = 2.7465e-3
s4_WGS = -4.53e-6
s5_WGS = -0.201

g0_inter_FT = 80.065 # kJ/mol
g0_slope_FT = -0.1471

op_ratio = [[1,0],[2,1.86],[3,4],[4,3.17],[5,1.7],[11,0.67],[21,0],[57,0]]

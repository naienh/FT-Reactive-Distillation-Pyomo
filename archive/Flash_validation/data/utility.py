from modules.global_set import m
from pyomo.core.base.block import generate_cuid_names
import json
import os, sys
import pickle
import numpy as np

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = None

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout
#####----------------------------------------------------------------------#####
def trans_product(dic):
    c1 = [i for i in m.COMP_ORG if cal_cnumber(i) == 1]
    c2 = [i for i in m.COMP_ORG if cal_cnumber(i) == 2]
    c3 = [i for i in m.COMP_ORG if cal_cnumber(i) == 3]
    c4 = [i for i in m.COMP_ORG if cal_cnumber(i) == 4]
    napha = [i for i in m.COMP_ORG if cal_cnumber(i) >= 5 and cal_cnumber(i) <= 7]
    gasoline = [i for i in m.COMP_ORG if cal_cnumber(i) >= 8 and cal_cnumber(i) <= 12]
    diesel = [i for i in m.COMP_ORG if cal_cnumber(i) >= 13 and cal_cnumber(i) <= 18]
    heavy = [i for i in m.COMP_ORG if cal_cnumber(i) >= 19 and cal_cnumber(i) <= 56]

    dataset = {}
    for i in c1:
        tmp = np.array([np.array(dic[i]) for i in c1])
        dataset['c1'] = np.sum(tmp,0)

    for i in c2:
        tmp = np.array([np.array(dic[i]) for i in c2])
        dataset['c2'] = np.sum(tmp,0)

    for i in c3:
        tmp = np.array([np.array(dic[i]) for i in c3])
        dataset['c3'] = np.sum(tmp,0)

    for i in c4:
        tmp = np.array([np.array(dic[i]) for i in c4])
        dataset['c4'] = np.sum(tmp,0)

    for i in napha:
        tmp = np.array([np.array(dic[i]) for i in napha])
        dataset['napha'] = np.sum(tmp,0)

    for i in gasoline:
        tmp = np.array([np.array(dic[i]) for i in gasoline])
        dataset['gasoline'] = np.sum(tmp,0)

    for i in diesel:
        tmp = np.array([np.array(dic[i]) for i in diesel])
        dataset['diesel'] = np.sum(tmp,0)

    for i in heavy:
        tmp = np.array([np.array(dic[i]) for i in heavy])
        dataset['heavy'] = np.sum(tmp,0)

    return dataset

def trans_product_scaled(dic):
    tmp = trans_product(dic)
    scale_factor = 1/sum(tmp[i] for i in tmp.keys())
    for i in tmp.keys():
        tmp[i] = tmp[i]*scale_factor
    return tmp

#####----------------------------------------------------------------------#####

def cal_MW(name):
    if name == 'H2O': return 18
    if name == 'CO2': return 44
    if name == 'CO': return 28
    if name == 'H2': return 2
    start = name.index('C')+1
    end = name.index('H')
    MW = int(name[start:end])*12 + int(name[end+1:])
    return MW

#####----------------------------------------------------------------------#####
# print data from any of the data modules imported
def print_pkg(m):
    key = [v for v in vars(m).keys() if not v.startswith('_')]
    for i in key: print(i)
#####----------------------------------------------------------------------#####
# generate entire list of op ratio based on data, format [#C,ratio]
def cal_op(op_ratio):
    breakpoints = len(op_ratio);i = 1; paraffin_ratio = []; olefin_ratio = [];
    for k in range(breakpoints-1):
        while i < op_ratio[k+1][0]:
            paraffin_ratio.append(1/(op_ratio[k][1]+1)); i += 1;
    return paraffin_ratio

def cal_cnumber(i):
    if i in m.COMP_OLEFIN: return m.COMP_OLEFIN.ord(i)+1
    if i in m.COMP_PARAFFIN: return m.COMP_PARAFFIN.ord(i)

#####----------------------------------------------------------------------#####

class result_analysis_function(object):

    def gas_percentage(self,m,j,n):
        value = sum(m.y[j,i].value for i in m.COMP_TOTAL if cal_cnumber(i) == n)
        return value

    def liquid_percentage(self,m,j,n):
        value = sum(m.x[j,i].value for i in m.COMP_TOTAL if cal_cnumber(i) == n)
        return value

    def print_main(self,pyomo,m):
        for j in m.TRAY:
            print('\nTray: ',j)

            print('-'*100)
            print('Block Variable\n')
            print('T (K):\t\t{:7.1f}\t|\tF (mol/s):\t{:7.4f}\t\t|\tR_FT (kmol/s):\t{:7.4f}'.format(m.T[j].value,m.F[j].value,m.r_FT_total[j].value))
            print('P (bar):\t{:7.1f}\t|\tV (mol/s):\t{:7.4f}\t\t|\tR_WGS (kmol/s):\t{:7.4f}'.format(m.P[j].value,m.V[j].value,m.r_WGS[j].value))
            print('Q (MW):\t\t{:7.1f}\t|\tL (mol/s):\t{:7.4f}\t\t|\t'.format(m.Q_main[j].value,m.L[j].value))
            print('')

            cog = m.y[j,'CO'].value
            h2g = m.y[j,'H2'].value
            co2g = m.y[j,'CO2'].value
            h2og = m.y[j,'H2O'].value

            c1g = self.gas_percentage(m,j,1)
            c2g = self.gas_percentage(m,j,2)
            c3g = self.gas_percentage(m,j,3)
            c4g = self.gas_percentage(m,j,4)
            c57g = sum(self.gas_percentage(m,j,i) for i in range(5,8))
            c812g = sum(self.gas_percentage(m,j,i) for i in range(8,13))
            c1318g = sum(self.gas_percentage(m,j,i) for i in range(13,19))
            c19pg = sum(self.gas_percentage(m,j,i) for i in range(19,57))

            print('-'*100)
            print('Gas Phase\n')
            print('CO:\t\t{:.2%}\t|\tC1:\t\t{:.2%}\t\t|\tC5-C7:\t\t{:.2%}'.format(cog,c1g,c57g))
            print('H2:\t\t{:.2%}\t|\tC2:\t\t{:.2%}\t\t|\tC8-12:\t\t{:.2%}'.format(h2g,c2g,c812g))
            print('CO2:\t\t{:.2%}\t|\tC3:\t\t{:.2%}\t\t|\tC13-18:\t\t{:.2%}'.format(co2g,c3g,c1318g))
            print('H2O:\t\t{:.2%}\t|\tC4:\t\t{:.2%}\t\t|\tC19+:\t\t{:.2%}'.format(h2og,c4g,c19pg))


            col = m.x[j,'CO'].value
            h2l = m.x[j,'H2'].value
            co2l = m.x[j,'CO2'].value
            h2ol = m.x[j,'H2O'].value

            c1l = self.liquid_percentage(m,j,1)
            c2l = self.liquid_percentage(m,j,2)
            c3l = self.liquid_percentage(m,j,3)
            c4l = self.liquid_percentage(m,j,4)
            c57l = sum(self.liquid_percentage(m,j,i) for i in range(5,8))
            c812l = sum(self.liquid_percentage(m,j,i) for i in range(8,13))
            c1318l = sum(self.liquid_percentage(m,j,i) for i in range(13,19))
            c19pl = sum(self.liquid_percentage(m,j,i) for i in range(19,57))

            print('-'*100)
            print('Liquid Phase\n')
            print('CO:\t\t{:.2%}\t|\tC1:\t\t{:.2%}\t\t|\tC5-C7:\t\t{:.2%}'.format(col,c1l,c57l))
            print('H2:\t\t{:.2%}\t|\tC2:\t\t{:.2%}\t\t|\tC8-12:\t\t{:.2%}'.format(h2l,c2l,c812l))
            print('CO2:\t\t{:.2%}\t|\tC3:\t\t{:.2%}\t\t|\tC13-18:\t\t{:.2%}'.format(co2l,c3l,c1318l))
            print('H2O:\t\t{:.2%}\t|\tC4:\t\t{:.2%}\t\t|\tC19+:\t\t{:.2%}'.format(h2ol,c4l,c19pl))


result = result_analysis_function()
#------------------------------------------------------------------------------
def add_dual(pyomo,model):
    ### Declare all pe.Suffixes
    # Ipopt bound multipliers (obtained from solution)
    model.ipopt_zL_out = pyomo.Suffix(direction=pyomo.Suffix.IMPORT)
    model.ipopt_zU_out = pyomo.Suffix(direction=pyomo.Suffix.IMPORT)
    # Ipopt bound multipliers (sent to solver)
    model.ipopt_zL_in = pyomo.Suffix(direction=pyomo.Suffix.EXPORT)
    model.ipopt_zU_in = pyomo.Suffix(direction=pyomo.Suffix.EXPORT)
    # Obtain dual solutions from first solve and send to warm start
    model.dual = pyomo.Suffix(direction=pyomo.Suffix.IMPORT_EXPORT)
    ###
    print('Created the follow pyomo suffixes:')
    print('ipopt_zL_out, ipopt_zU_out, ipopt_zL_in, ipopt_zU_in, dual')

def update_dual(pyomo,model):
    model.ipopt_zL_in.update(model.ipopt_zL_out)
    model.ipopt_zU_in.update(model.ipopt_zU_out)

#------------------------------------------------------------------------------
def save_solution(pyomo,model,name):
    labels = generate_cuid_names(model)
    solution = {}
    for var in model.component_data_objects(pyomo.Var):
        solution[labels[var]] = var.value
    with open('{}.json'.format(name), "w") as f:
        json.dump(solution, f)
    del solution
    print('Valar Morghulis\nVariable Solution Saved!')

def load_solution(pyomo,model,name):
    with open('{}.json'.format(name)) as f:
        solution = json.load(f)
    for cuid, val in solution.items():
        model.find_component(cuid).value = val
    print('Valar Dohaeris\nVariable Solution Loaded!')

def save_all(model,results,name):
    model.solutions.store_to(results)
    with open('{}.pickle'.format(name), 'wb') as f:
        pickle.dump(results,f)
    print('Valar Morghulis\nEntire Solution (include dual) Saved')
    print('-'*100)
    print(results.Solution.Variable.items())
    print('-'*100)
    print(results.Solution.Constraint.items())

def load_all(model,results,name):
    with open('{}.pickle'.format(name), 'rb') as f:
        results = pickle.load(f)
    model.solutions.load_from(results)
    print('Valar Dohaeris\nEntire Solution (include dual) Loaded')
    print('-'*100)
    print(results.Solution.Variable.items())
    print('-'*100)
    print(results.Solution.Constraint.items())

#------------------------------------------------------------------------------

def beautify(pe,model):
    print('Here comes the result:')
    print('-'*100)

    print('T\t\t\t Q\t\t\t V\t\t\t L')
    print(model.condenser.T.value - 273.15,'\t',model.condenser.Q_main.value,'\t',model.condenser.V['out'].value,'\t',model.condenser.L['out'].value)
    for j in model.TRAY:
        print(model.reactive[j].T.value - 273.15,'\t',model.reactive[j].Q_main.value,'\t',model.reactive[j].V['out'].value,'\t',model.reactive[j].L['out'].value)
    print(model.reboiler.T.value - 273.15,'\t',model.reboiler.Q_main.value,'\t',model.reboiler.V['out'].value,'\t',model.reboiler.L['out'].value)
    print('-'*100)

    print('Top Product V/L/W')
    print(model.condenser.V['out'].value)
    print(model.condenser.L['P'].value)
    print(model.reactive[model.TRAY.first()].V['out'].value*0.95*model.reactive[model.TRAY.first()].y['H2O'].value)
    print('-'*100)

    print('Bottom Product L')
    print(model.reboiler.L['out'].value)
    print('-'*100)

    print('Condenser:\tVapor\t\tLiquid\t\tReboiler\tVapor\t\tLiquid')
    for i in m.COMP_TOTAL:
        if model.condenser.y[i].value > 1e-5 or model.condenser.x[i].value > 1e-5:
            print(i,'\t\t{:6.3%}\t\t{:6.3%}\t\t\t\t{:6.3%}\t\t{:6.3%}'\
                  .format(model.condenser.y[i].value,model.condenser.x[i].value,model.reboiler.y[i].value,model.reboiler.x[i].value))
#------------------------------------------------------------------------------
def beautify2(pe,model):
    print('Here comes the result:')
    print('-'*100)

    print('\t\tT\t\t Q\t\t\t V\t\t\t L')
    print('{:20.2f}'.format(model.condenser.T.value - 273.15),'\t','{:20.10f}'.format(model.condenser.Q_main.value),'\t','{:20.10f}'.format(model.condenser.V['out'].value),'\t','{:20.10f}'.format(model.condenser.L['out'].value))
    for j in model.TRAY:
        print('{:20.2f}'.format(model.reactive[j].T.value - 273.15),'\t','{:20.10f}'.format(model.reactive[j].Q_main.value),'\t','{:20.10f}'.format(model.reactive[j].V['out'].value),'\t','{:20.10f}'.format(model.reactive[j].L['out'].value))
    print('-'*100)

    print('Top Product')
    print('V\t',model.condenser.V['out'].value)
    print('L\t',model.condenser.L['P'].value)
    print('W\t',model.condenser.W.value)
    print('-'*100)

    print('Bottom Product L')
    print(model.reactive[model.TRAY.last()].L['out'].value)
    print('-'*100)

    print('Condenser:\tVapor\t\tLiquid\t\tLast Stage\tVapor\t\tLiquid')
    for i in m.COMP_TOTAL:
        if model.condenser.y[i].value > 1e-5 or model.condenser.x[i].value > 1e-5:
            print(i,'\t\t{:6.3%}\t\t{:6.3%}\t\t'.format(model.condenser.y[i].value,model.condenser.x[i].value),'{:8s}'.format(i),'\t{:6.3%}\t\t{:6.3%}'\
                  .format(model.reactive[model.TRAY.last()].y[i].value,model.reactive[model.TRAY.last()].x[i].value))



#####----------------------------------------------------------------------#####
# from excel read one comlumn, for convenience only
def readcol(sheet,name):
    data_name = [sheet.cell_value(0,i) for i in range(sheet.ncols)]
    data = {}
    j = data_name.index(name)
    for i in range(1,sheet.nrows):
        data[sheet.cell_value(i,0)]=sheet.cell_value(i,j)
    return data

# sometimes bounds are None (inf), needs tools to print bounds reliably
def print_bounds(l,body,u):
    if l is None and u is None:
        print('-Inf','<=',body,'<=','+Inf')
    elif l is None:
        print('-Inf','<=',body,'<=',u)
    elif u is None:
        print(l,'<=',body,'<=','+Inf')
    else:
        print(l,'<=',body,'<=',u)

# given a specific constraint, print its name, expression and current value
def print_constraint(pyomo,c):
    for i in c:
        print('*'*60)
        print('Constraint: ',c[i])
        print(c[i].body)
        print('-'*60)
        print_bounds(c[i].lower,pyomo.value(c[i].body),c[i].upper)

def check_violate_constraint(pyomo,m):
    counter = 0;
    for v in m.component_data_objects(pyomo.Constraint, active=True):
        resL = 0; resU = 0;
        if v.lower is not None: resL = max(0,v.lower-pyomo.value(v.body))
        if v.upper is not None: resU = max(0,pyomo.value(v.body)-v.upper)
        if resL > 1e-6 or resU > 1e-6:
            print('*'*60)
            print('Constraint: ',v)
            print(v.body)
            print('-'*60)
            print_bounds(v.lower,pyomo.value(v.body),v.upper)
            counter += 1
    if counter == 0:
        print('Congratz, all constraint vialation <= 1e-6')

def print_variable(pyomo,v):
    for i in v:
        print('-'*60)
        print('Variable: ',v[i])
        print_bounds(v[i].lb,pyomo.value(v[i]),v[i].ub)

def check_forced_variable(pyomo,m):
    counter = 0;
    for v in m.component_data_objects(pyomo.Var, active=True):
        resL = 1; resU = 1;
        if v.fixed == True: continue
        if v.lb is not None: resL = abs(pyomo.value(v)-v.lb)
        if v.ub is not None: resU = abs(v.ub-pyomo.value(v))
        if resL < 1e-8 or resU < 1e-8:
            print('-'*50)
            print(v)
            print_bounds(v.lb,pyomo.value(v),v.ub)
            counter += 1
    if counter == 0:
        print('Congratz, no variable is at bounds')

def check_DOF(pyomo,m):
    constraint_number = len([i for i in m.component_data_objects(pyomo.Constraint, active=True)])
    variable_number = len([i for i in m.component_data_objects(pyomo.Var, active=True)])
    fixed_variable_number = sum([i.fixed for i in m.component_data_objects(pyomo.Var, active=True)])
    print('Active Constraints:\t',constraint_number)
    print('Active Variables:\t',variable_number)
    print('Fixed Variables:\t',fixed_variable_number)
    print('DOF:\t\t\t',variable_number-fixed_variable_number-constraint_number)

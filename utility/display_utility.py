from global_sets.component import m
import os, sys
import copy
import numpy as np
from utility.data_utility import cal_MW, cal_cnumber

'''-----------------------------------------------------------------------------
This can be used to slient output produced by constructing blocks
Usage:
with HiddenPrints():
    print("This will not be printed")
-----------------------------------------------------------------------------'''
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = None

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout

'''-----------------------------------------------------------------------------
This is used to group component (C10H22) data into product (gasoline) data
1. trans_product_mole: mole, mole%
2. trans_product_mass: mass, mass%
Usage:
1. Prepare the dictionary into the following shape:
    dic = {'C10H20':[...],'C5H12':[...],...}
2. Use the function:
    reaction_data = trans_product_mole(dic)
3. Retuen dictionary:
    reaction_data['unscaled'] = {'gasoline':[...]}
    reaction_data['scaled'] = {'gasoline':[...]}
-----------------------------------------------------------------------------'''
def trans_product_mole(dic):
    product = {}
    product['c1'] = [i for i in m.COMP_ORG if cal_cnumber(i) == 1]
    product['c2'] = [i for i in m.COMP_ORG if cal_cnumber(i) == 2]
    product['c3'] = [i for i in m.COMP_ORG if cal_cnumber(i) == 3]
    product['c4'] = [i for i in m.COMP_ORG if cal_cnumber(i) == 4]
    product['napha'] = [i for i in m.COMP_ORG if cal_cnumber(i) >= 5 and cal_cnumber(i) <= 7]
    product['gasoline'] = [i for i in m.COMP_ORG if cal_cnumber(i) >= 8 and cal_cnumber(i) <= 12]
    product['diesel'] = [i for i in m.COMP_ORG if cal_cnumber(i) >= 13 and cal_cnumber(i) <= 18]
    product['heavy'] = [i for i in m.COMP_ORG if cal_cnumber(i) >= 19 and cal_cnumber(i) <= 56]

    # compute mole flow rate
    dataset = {}
    for c in product.keys():
        for i in product[c]:
            tmp = np.array([np.array(dic[i]) for i in product[c]])
            dataset[c] = np.sum(tmp,0)

    # compute mole % scaled
    dataset_scaled = copy.deepcopy(dataset)
    scale_factor = 1/sum(dataset_scaled[i] for i in dataset_scaled.keys())
    for i in dataset_scaled.keys():
        dataset_scaled[i] = dataset_scaled[i]*scale_factor

    return {'unscaled':dataset,'scaled':dataset_scaled}

def trans_product_mass(dic):
    mass_data = {}
    for i in m.COMP_ORG:
        mass_data[i] = np.array(dic[i])*cal_MW(i)
    return trans_product_mole(mass_data)

'''-----------------------------------------------------------------------------
This can be used to pretty print a reactive distillation solution
Usage:
    1. beautify: with Reboiler
    2. beautify2: without Reboiler
This section updates regularly to reflect the latest need to print solution
-----------------------------------------------------------------------------'''
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

    print('Top')
    print('V\t',model.condenser.V['out'].value)
    print('L\t',model.condenser.L['P'].value)
    print('W\t',model.condenser.W.value)
    print('-'*100)

    print('Bottom')
    print('L\t',model.reactive[model.TRAY.last()].L['out'].value)
    print('-'*100)

    print('Condenser:\tVapor\t\tLiquid\t\tLast Stage\tVapor\t\tLiquid')
    for i in m.COMP_TOTAL:
        if model.condenser.y[i].value > 1e-5 or model.condenser.x[i].value > 1e-5:
            print(i,'\t\t{:6.3%}\t\t{:6.3%}\t\t'.format(model.condenser.y[i].value,model.condenser.x[i].value),'{:8s}'.format(i),'\t{:6.3%}\t\t{:6.3%}'\
                  .format(model.reactive[model.TRAY.last()].y[i].value,model.reactive[model.TRAY.last()].x[i].value))

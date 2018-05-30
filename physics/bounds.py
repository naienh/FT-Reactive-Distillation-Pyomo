import pickle

'''-----------------------------------------------------------------------------
This is used to collect data from solved example to set up variable bounds
A lot of local variables are temperature dependent, so can be safely bounded
-----------------------------------------------------------------------------'''
with open('../saved_solutions/reactive_flash_200C.pickle', 'rb') as f:
    C200 = pickle.load(f)
with open('../saved_solutions/reactive_flash_300C.pickle', 'rb') as f:
    C300 = pickle.load(f)

def collect_bounds1(name):
    return {i.replace(name,''):[C200.solution.Variable[i]['Value'],C300.solution.Variable[i]['Value']]\
     for i in C200.solution.Variable.keys() if i.startswith(name)}

kinetic_bounds = collect_bounds1('kinetics_block.')
energy_bounds = collect_bounds1('energy_block.')
VLE_bounds = collect_bounds1('VLE_block.')

with open('../saved_solutions/condenser_20C.pickle', 'rb') as f:
    C20 = pickle.load(f)
with open('../saved_solutions/condenser_40C.pickle', 'rb') as f:
    C40 = pickle.load(f)

def collect_bounds2(name):
    return {i.replace(name,''):[C20.solution.Variable[i]['Value'],C40.solution.Variable[i]['Value']]\
     for i in C20.solution.Variable.keys() if i.startswith(name)}

energy_bounds2 = collect_bounds2('energy_block.')
VLLE_bounds = collect_bounds2('VLE_block.')
water_x = [C20.solution.Variable['x[H2O]']['Value'],C40.solution.Variable['x[H2O]']['Value']]
water_yp = [C20.solution.Variable['y[H2O]']['Value']*21,\
            C40.solution.Variable['y[H2O]']['Value']*21]

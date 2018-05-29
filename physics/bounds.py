import pickle

with open('../saved_solutions/reactive_flash_200C.pickle', 'rb') as f:
    C200 = pickle.load(f)
with open('../saved_solutions/reactive_flash_300C.pickle', 'rb') as f:
    C300 = pickle.load(f)

def collect_bounds(name):
    return {i.replace(name,''):[C200.solution.Variable[i]['Value'],C300.solution.Variable[i]['Value']]\
     for i in C200.solution.Variable.keys() if i.startswith(name)}

kinetic_bounds = collect_bounds('kinetics_block.')
energy_bounds = collect_bounds('energy_block.')
VLE_bounds = collect_bounds('VLE_block.')

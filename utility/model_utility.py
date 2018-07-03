from global_sets.component import m
import pickle

'''-----------------------------------------------------------------------------
This is used to create dual variables for a model
Usage:
1. Input the local pyomo and model name
2. This function will create ipopt dual variables for you
-----------------------------------------------------------------------------'''
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

'''-----------------------------------------------------------------------------
This is used to delete dual variables for a model
Usage:
1. To add dual variables for a newly added block, you have to delete old dual
variables first, this function does just that
-----------------------------------------------------------------------------'''
def delete_dual(pyomo,model):
    model.del_component(model.ipopt_zL_out)
    model.del_component(model.ipopt_zU_out)
    model.del_component(model.ipopt_zL_in)
    model.del_component(model.ipopt_zU_in)
    model.del_component(model.dual)

'''-----------------------------------------------------------------------------
This is used to update dual variables for a model
Usage:
1. Input the local pyomo and model name
2. This function will update ipopt dual variables for you. Note that after each
    solve, pyomo will put value to .ipopt_zL_out variable, you have to update
    this value to .ipopt_zL_in variable, in order to warm start
-----------------------------------------------------------------------------'''
def update_dual(pyomo,model):
    model.ipopt_zL_in.update(model.ipopt_zL_out)
    model.ipopt_zU_in.update(model.ipopt_zU_out)

'''-----------------------------------------------------------------------------
This is used to save solutions of a solve, there are two types of a solve, one
is to only save the values of variables, the other saves the entire thing
Usage:
1. Input the local pyomo and model name
2. This function will update ipopt dual variables for you. Note that after each
    solve, pyomo will put value to .ipopt_zL_out variable, you have to update
    this value to .ipopt_zL_in variable, in order to warm start
-----------------------------------------------------------------------------'''
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

def load_result(model,name):
    with open('{}.pickle'.format(name), 'rb') as f:
        results = pickle.load(f)
    return results

def load_all(model,results,name):
    model.solutions.load_from(results)
    print('Valar Dohaeris\nEntire Solution (include dual) Loaded')
    print('-'*100)
    print(results.Solution.Variable.items())
    print('-'*100)
    print(results.Solution.Constraint.items())
    return results

'''-----------------------------------------------------------------------------
This is used to debug pyomo models by printing relevant constraints and
variable bounds
-----------------------------------------------------------------------------'''
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

'''-----------------------------------------------------------------------------
This is used to get the current DOF of the model, will be checking for
equality and inequality, as well as fixed variables
-----------------------------------------------------------------------------'''
def check_DOF(pyomo,m):
    equality_number = len([i for i in m.component_data_objects(pyomo.Constraint, active=True)\
                if i.upper is not None and i.lower is not None])
    inequality_number = len([i for i in m.component_data_objects(pyomo.Constraint, active=True)\
                if i.upper is None or i.lower is None])
    variable_number = len([i for i in m.component_data_objects(pyomo.Var, active=True)])
    fixed_variable_number = sum([i.fixed for i in m.component_data_objects(pyomo.Var, active=True)])
    print('Active Equality Constraints:\t',equality_number)
    print('Active Inequality Constraints:\t',inequality_number)
    print('Active Variables:\t\t',variable_number)
    print('Fixed Variables:\t\t',fixed_variable_number)
    print('DOF:\t\t\t\t',variable_number-fixed_variable_number-equality_number)

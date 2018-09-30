#!/usr/bin/env python3
'''
How to Train Your Dragon: V6
New changes, remove stripping section (set it to fixed value in this case)
To be used with multi-start, randomized system operating parameters
Sequentially initialize FT reactive distillation model automatically
Add DDF module and then optimize during one step
'''
# system imports
import sys
import os
import datetime
import shutil
import atexit
sys.path.append(os.path.abspath('..'))
sys.path.append(os.path.abspath('../..'))

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pickle
import dill
from copy import deepcopy
from datetime import timedelta, datetime

# pyomo imports
from pyomo import environ as pe
from global_sets.component import m

from stages.reactive_stage import reactive_stage_rule
from stages.condenser_stage import condenser_stage_rule
from stages.reboiler_stage import reboiler_stage_rule

from utility.display_utility import beautify, beautify_reactive, HiddenLogs, HiddenPrints, \
                                    plot_distribution, plot_product_distribution, check_product_spec
from utility.model_utility import add_dual, update_dual, delete_dual, check_DOF, check_iteration, tray_translator
from utility.model_utility import which_MPCC, select_MPCC, augmented_objective, add_solver, disable_restoration
from utility.time_utility import create_filename_time, log_now, log_end

class data_object(object):
    def __init__(self, name):
        self.name = name

plt.ion()

'''
Constructing the model and logfile
'''
model = pe.ConcreteModel(name='reactive_distillation')
mpcc_type = 'pf'
# digest input arguments
if len(sys.argv) == 1:
    logname = create_filename_time()+'_'+mpcc_type
    log_text_dir = './log/text/mul_onestep_'+logname+'.dat'
    log_figure_dir = './log/figure/mul_onestep_'+logname+'.pdf'
    output_dir = './tmp/'+logname+'.output'
    log_master_dir = './log/master/master_log.txt'
elif len(sys.argv) == 4:
    rand_file_name = sys.argv[1]
    logname = 'Preset_Case:_{}'.format(sys.argv[2])+'_'+mpcc_type
    log_text_dir = './log/text/mul_onestep_'+logname+'.dat'
    log_figure_dir = './log/figure/mul_onestep_'+logname+'.pdf'
    output_dir = './tmp/'+logname+'.output'
    log_master_dir = './log/master/'+sys.argv[3]
else:
    exit('Argument Number DO NOT MATCH')

'''
write its own option files (because it's gonna alter it during solve), clean up after exit
'''

option_dir = logname+'.opt'
shutil.copy('ipopt.opt',option_dir)

def cleanup():
    if option_dir:
        if os.path.exists(option_dir):
            os.remove(option_dir)
atexit.register(cleanup)

os.makedirs('./log/text',exist_ok=True)
os.makedirs('./log/figure',exist_ok=True)
os.makedirs('./log/model',exist_ok=True)
os.makedirs('./log/master',exist_ok=True)
os.makedirs('./tmp',exist_ok=True)

'''
model input, randomized, or input from already generated random files
'''

if len(sys.argv) == 1:
    # The following parameters requires manual input
    tray_number = 20
    non_reactive_flag = [1,2,3,4,5,6,7]

    # The following parameters is generated manually
    # reset random seeds
    np.random.seed()

    # distillate ratio = 0.05 - 0.15 total condenser liquid
    rand = np.random.rand()
    rr_ratio = 0.05 + 0.1*(rand)

    # generate three product flag, randint produces [low,high)
    side_draw_flag = {}
    # intermediate 2 - 4
    rand_int = np.random.randint(low=1,high=3)
    rand = np.random.rand()
    side_draw_flag.update({rand_int:0.01+0.02*rand})
    # gasoline 5 - 9
    rand_int = np.random.randint(low=3,high=10)
    rand = np.random.rand()
    side_draw_flag.update({rand_int:0.1+0.2*rand})
    # diesel 11 - 18
    rand_int = np.random.randint(low=10,high=19)
    rand = np.random.rand()
    side_draw_flag.update({rand_int:0.2+0.3*rand})

    # generating temperature profile
    profile_tmp = sorted(np.random.rand(tray_number - len(non_reactive_flag)))
    temperature_flag = {}; n = 0
    for j in range(1,tray_number+1):
        if j not in non_reactive_flag:
            temperature_flag[j] = 200+100*profile_tmp[n]
            n += 1

    # generating catalyst profile, total = 30000 kg
    profile_tmp = {j:np.random.rand() for j in range(1,tray_number+1) if j not in non_reactive_flag}
    total_tmp = sum(profile_tmp[j] for j in profile_tmp)
    catalyst_flag = {j:100 + (30000-100*len(profile_tmp))*profile_tmp[j]/total_tmp for j in range(1,tray_number+1) if j not in non_reactive_flag}
    catalyst_flag.update({j:0 for j in non_reactive_flag})

    # generating feed profile, total = 10 kmol/s
    profile_tmp = {j:np.random.rand() for j in range(1,tray_number+1) if j not in non_reactive_flag}
    total_tmp = sum(profile_tmp[j] for j in profile_tmp)
    feed_flag = {j:0.01 + (10-0.01*len(profile_tmp))*profile_tmp[j]/total_tmp for j in range(1,tray_number+1) if j not in non_reactive_flag}
    feed_flag.update({j:0 for j in non_reactive_flag})

elif len(sys.argv) == 4:
    with open(rand_file_name,'rb') as f:
        collection = pickle.load(f)
    random_state = collection[int(sys.argv[2])]

    tray_number = random_state.tray_number
    non_reactive_flag = random_state.non_reactive_flag
    rr_ratio = random_state.rr_ratio
    side_draw_flag = random_state.side_draw_flag
    temperature_flag = random_state.temperature_flag
    catalyst_flag = random_state.catalyst_flag
    feed_flag = random_state.feed_flag

with PdfPages(log_figure_dir,keep_empty=False) as pdf, open(log_master_dir,'a') as master_log:

    master_log.write('\n'+'-'*108+'\n')
    master_log.write('ATTENTION: NEW RUN BEGINNING\n')
    master_log.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    master_log.write('\nlog name: {}\n'.format(logname))

    '''
    construct the reactive stages
    '''
    model.TRAY = pe.RangeSet(1,tray_number)
    model.TRAY_nonreactive = pe.Set(initialize=non_reactive_flag)
    model.TRAY_reactive = model.TRAY - model.TRAY_nonreactive
    model.TRAY_total = pe.Set(initialize=['condenser']+[str(i) for i in model.TRAY]+['reboiler'],ordered=True)

    with HiddenPrints():
        model.reactive = pe.Block(model.TRAY,rule=reactive_stage_rule)
        # change reactive flash MPCC type based on model input
        for j in model.reactive:
            select_MPCC(model.reactive[j],mpcc_type)
        add_dual(pe,model)

    '''
    setting stream variables
    '''
    # in/out variable
    for j in model.reactive:
        model.reactive[j].x_.fix(0)
        model.reactive[j].y_.fix(0)
        model.reactive[j].L['in'].fix(0)
        model.reactive[j].V['in'].fix(0)
        model.reactive[j].L['R'].fix(0)
        model.reactive[j].V['R'].fix(0)
        model.reactive[j].H_L_.fix(0)
        model.reactive[j].H_V_.fix(0)

    # operating parameters
    for j in model.reactive:
        model.reactive[j].cat.fix(3000)
        model.reactive[j].P.fix(20)
        model.reactive[j].VLE_block.n_ave.fix(20)

        model.reactive[j].F.fix(1)
        model.reactive[j].T_F.fix(200+273.15)
        model.reactive[j].z['CO'].fix(1/(1+2)-0/2)
        model.reactive[j].z['H2'].fix(2/(1+2)-0/2)
        model.reactive[j].z['C30H62'].fix(0)

        model.reactive[j].PR_L.fix(1)
        model.reactive[j].PR_V.fix(1)

        # model.reactive[j].Q_main.fix(0)
        model.reactive[j].T.setub(220+273.15)
        model.reactive[j].T.setlb(200+273.15)

    model.obj = augmented_objective(pe,model,expr = sum(model.reactive[j].T for j in model.reactive),sense=pe.maximize)

    '''
    solver options
    '''
    disable_restoration(mode = 'enable', option_dir = option_dir)
    opt = add_solver(pe, max_iter = 2000, warm_start = False, output = output_dir, option_dir = option_dir)

    progress = '> First Solve, disconnected reactive stages'

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    with HiddenLogs(log_text_dir,'w'):
        print('\n{}'.format(progress))
        print('-'*108)
        beautify_reactive(pe,model)
        log_now()
        check_iteration(filename=output_dir)
        if results.solver.termination_condition.key != 'optimal':
            print('NOT Optimum: ' + results.solver.termination_condition.key)
            master_log.write('Failed: {}\n'.format(progress))
            exit()


    '''
    connect reactive stages
    '''

    def V_between_rule(model,j):
        if j == model.TRAY.last(): return pe.Constraint.Skip
        return model.reactive[j].V['in'] == model.reactive[j+1].V['out']
    model.V_between_con = pe.Constraint(model.TRAY,rule=V_between_rule)

    def Vy_between_rule(model,j,i):
        if j == model.TRAY.last(): return pe.Constraint.Skip
        return model.reactive[j].y_['in',i] == model.reactive[j+1].y[i]
    model.Vy_between_con = pe.Constraint(model.TRAY,m.COMP_TOTAL,rule=Vy_between_rule)

    def Vh_between_rule(model,j):
        if j == model.TRAY.last(): return pe.Constraint.Skip
        return model.reactive[j].H_V_['in'] == model.reactive[j+1].H_V
    model.Vh_between_con = pe.Constraint(model.TRAY,rule=Vh_between_rule)

    def L_between_rule(model,j):
        if j == model.TRAY.last(): return pe.Constraint.Skip
        return model.reactive[j+1].L['in'] == model.reactive[j].L['out']
    model.L_between_con = pe.Constraint(model.TRAY,rule=L_between_rule)

    def Lx_between_rule(model,j,i):
        if j == model.TRAY.last(): return pe.Constraint.Skip
        return model.reactive[j+1].x_['in',i] == model.reactive[j].x[i]
    model.Ly_between_con = pe.Constraint(model.TRAY,m.COMP_TOTAL,rule=Lx_between_rule)

    def Lh_between_rule(model,j):
        if j == model.TRAY.last(): return pe.Constraint.Skip
        return model.reactive[j+1].H_L_['in'] == model.reactive[j].H_L
    model.Lh_between_con = pe.Constraint(model.TRAY,rule=Lh_between_rule)

    for j in model.reactive:
        if j != model.TRAY.first():
            for i in m.COMP_TOTAL:
                model.reactive[j].x_['in',i].unfix()
            model.reactive[j].H_L_['in'].unfix()
            model.reactive[j].L['in'].unfix()
        if j != model.TRAY.last():
            for i in m.COMP_TOTAL:
                model.reactive[j].y_['in',i].unfix()
            model.reactive[j].V['in'].unfix()
            model.reactive[j].H_V_['in'].unfix()

    '''
    warm start options
    '''
    opt = add_solver(pe, max_iter = 1000, warm_start = True, output = output_dir, option_dir = option_dir, tol=1e-8)

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    not_optimum_counter = 0

    if results.solver.termination_condition.key != 'optimal':
        not_optimum_counter += 1
    else:
        not_optimum_counter = 0

    if not_optimum_counter == 1:
        with HiddenLogs(log_text_dir):
            print('NOT Optimum: ' + results.solver.termination_condition.key)
        exit()

    '''
    The model above has liquid and vapor leave stage as products
    Now let liquid and vapor flow up/down stages
    '''

    PR_range = [0.8,0.6,0.4,0.2,0.1,0]
    with HiddenLogs(log_text_dir):
        print('\n> Connect stages and solve with 0 inter-stage flow')
        log_now()
        check_iteration(output_dir)

    for r in PR_range:
        for j in model.reactive:
            model.reactive[j].PR_L.fix(r)
            model.reactive[j].PR_V.fix(r)

        progress = '> Working on PR ratio = {:.2f}'.format(r)
        results = opt.solve(model,tee=False)
        update_dual(pe,model)

        if results.solver.termination_condition.key != 'optimal':
            not_optimum_counter += 1
        else:
            not_optimum_counter = 0

        with HiddenLogs(log_text_dir):
            print('\n{}'.format(progress))
            print('-'*108)
            beautify_reactive(pe,model)
            log_now()
            check_iteration(output_dir)
            if not_optimum_counter == 1:
                print('NOT Optimum: ' + results.solver.termination_condition.key)
                master_log.write('Failed: {}\n'.format(progress))
                exit()

    '''
    Deactivate reactive part of the model to initialize condenser and reboiler
    '''

    for i in model.block_data_objects():
        if i.name != model.name:
            i.deactivate()
    for i in model.component_objects(pe.Constraint, active=True):
        i.deactivate()

    with HiddenPrints():
        model.condenser = pe.Block(rule=condenser_stage_rule)
        delete_dual(pe,model)
        add_dual(pe,model)

    '''
    setting condenser stream variables to match reactive stages vapor output
    '''

    # in/out variables
    model.condenser.x_.fix(0)
    for i in m.COMP_TOTAL:
        model.condenser.y_['in',i].fix(model.reactive[model.TRAY.first()].y[i].value)
    model.condenser.V['in'].fix(model.reactive[model.TRAY.first()].V['out'].value)
    model.condenser.L['in'].fix(0)
    model.condenser.V['out'].fix(0)
    model.condenser.H_L_.fix(0)
    model.condenser.H_V_.fix(model.reactive[model.TRAY.first()].H_V.value)

    # operating parameters
    model.condenser.P.fix(19)
    model.condenser.T_F.fix(200+273.15)
    model.condenser.F.fix(0)
    model.condenser.z.fix(0)
    model.condenser.VLE_block.n_ave.fix(4)
    model.condenser.PR_L.fix(1)

    model.condenser.T.setub(30+273.15)

    model.del_component(model.obj)
    model.obj = augmented_objective(pe,model, expr = model.condenser.T, sense = pe.maximize)

    '''
    setting solver optopms
    '''
    opt = add_solver(pe, max_iter = 2000, warm_start = False, output = output_dir, option_dir = option_dir)

    progress = '> Initialize condenser'

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    if results.solver.termination_condition.key != 'optimal':
        not_optimum_counter += 1
    else:
        not_optimum_counter = 0

    with HiddenLogs(log_text_dir):
        print('\n{}'.format(progress))
        log_now()
        check_iteration(output_dir)
        if not_optimum_counter == 1:
            print('NOT Optimum: ' + results.solver.termination_condition.key)
            master_log.write('Failed: {}\n'.format(progress))
            exit()

    model.condenser.deactivate()

    with HiddenPrints():
        model.reboiler = pe.Block(rule=reboiler_stage_rule)
        select_MPCC(model.reboiler,mpcc_type)
        delete_dual(pe,model)
        add_dual(pe,model)

    '''
    setting reboiler stream variables to match reactive stages vapor output
    '''

    # in/out variables
    model.reboiler.y_.fix(0)
    for i in m.COMP_TOTAL:
        model.reboiler.x_['in',i].fix(model.reactive[model.TRAY.last()].x[i].value)
    model.reboiler.L['in'].fix(model.reactive[model.TRAY.last()].L['out'].value)
    model.reboiler.V['in'].fix(0)
    model.reboiler.L['out'].fix(0)
    model.reboiler.V['P'].fix(0)
    model.reboiler.H_L_.fix(model.reactive[model.TRAY.last()].H_L.value)
    model.reboiler.H_V_.fix(0)

    # operating parameters
    model.reboiler.P.fix(20)
    model.reboiler.T_F.fix(200+273.15)
    model.reboiler.F.fix(0)
    model.reboiler.z['CO'].fix(1/(1+2)-0/2)
    model.reboiler.z['H2'].fix(2/(1+2)-0/2)
    model.reboiler.z['C30H62'].fix(0)
    model.reboiler.VLE_block.n_ave.fix(20)

    model.reboiler.T.setub(model.reactive[model.TRAY.last()].T.value)

    model.del_component(model.obj)
    model.obj = augmented_objective(pe,model, expr = model.reboiler.T, sense = pe.maximize)

    progress = '> Initialize reboiler'

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    if results.solver.termination_condition.key != 'optimal':
        not_optimum_counter += 1
    else:
        not_optimum_counter = 0

    with HiddenLogs(log_text_dir):
        print('\n{}'.format(progress))
        log_now()
        check_iteration(output_dir)
        if not_optimum_counter == 1:
            print('NOT Optimum: ' + results.solver.termination_condition.key)
            master_log.write('Failed: {}\n'.format(progress))
            exit()

    '''
    Now all the variables are initialized, reactivate everything, time to connect condenser and reboiler
    '''

    for i in model.block_data_objects():
        if i.name != 'reactive_distillation':
            i.activate()
    for i in model.component_objects(pe.Constraint):
        i.activate()
    with HiddenPrints():
        for j in model.reactive:
            select_MPCC(model.reactive[j],mpcc_type)
        select_MPCC(model.reboiler,mpcc_type)


    def V_condenser_rule(model):
        return model.reactive[model.TRAY.first()].V['out'] == model.condenser.V['in']
    model.V_condenser_con = pe.Constraint(rule=V_condenser_rule)

    def Vy_condenser_rule(model,i):
        return model.reactive[model.TRAY.first()].y[i] == model.condenser.y_['in',i]
    model.Vy_condenser_con = pe.Constraint(m.COMP_TOTAL,rule=Vy_condenser_rule)

    def Vh_condenser_rule(model):
        return model.reactive[model.TRAY.first()].H_V == model.condenser.H_V_['in']
    model.Vh_condenser_con = pe.Constraint(rule=Vh_condenser_rule)

    def L_condenser_rule(model):
        return model.reactive[model.TRAY.first()].L['in'] == model.condenser.L['out']
    model.L_condenser_con = pe.Constraint(rule=L_condenser_rule)

    def Lx_condenser_rule(model,i):
        return model.reactive[model.TRAY.first()].x_['in',i] == model.condenser.x[i]
    model.Lx_condenser_con = pe.Constraint(m.COMP_TOTAL,rule=Lx_condenser_rule)

    def Lh_condenser_rule(model):
        return model.reactive[model.TRAY.first()].H_L_['in'] == model.condenser.H_L
    model.Lh_condenser_con = pe.Constraint(rule=Lh_condenser_rule)

    def V_reboiler_rule(model):
        return model.reactive[model.TRAY.last()].V['in'] == model.reboiler.V['out']
    model.V_reboiler_con = pe.Constraint(rule=V_reboiler_rule)

    def Vy_reboiler_rule(model,i):
        return model.reactive[model.TRAY.last()].y_['in',i] == model.reboiler.y[i]
    model.Vy_reboiler_con = pe.Constraint(m.COMP_TOTAL,rule=Vy_reboiler_rule)

    def Vh_reboiler_rule(model):
        return model.reactive[model.TRAY.last()].H_V_['in'] == model.reboiler.H_V
    model.Vh_reboiler_con = pe.Constraint(rule=Vh_reboiler_rule)

    # note the epi here
    def L_reboiler_rule(model):
        return model.reactive[model.TRAY.last()].L['out'] + 1e-4 == model.reboiler.L['in']
    model.L_reboiler_con = pe.Constraint(rule=L_reboiler_rule)

    def Lx_reboiler_rule(model,i):
        return model.reactive[model.TRAY.last()].x[i] == model.reboiler.x_['in',i]
    model.Lx_reboiler_con = pe.Constraint(m.COMP_TOTAL,rule=Lx_reboiler_rule)

    def Lh_reboiler_rule(model):
        return model.reactive[model.TRAY.last()].H_L == model.reboiler.H_L_['in']
    model.Lh_reboiler_con = pe.Constraint(rule=Lh_reboiler_rule)

    model.del_component(model.obj)
    model.obj = augmented_objective(pe,model,expr = sum(model.reactive[j].T for j in model.TRAY) ,sense=pe.maximize)

    with HiddenPrints():
        delete_dual(pe,model)
        add_dual(pe,model)

    '''
    Make sure all stream variables are set correctly
    '''

    # in/out variables
    model.condenser.x_.fix(0)
    for i in m.COMP_TOTAL:
        model.condenser.y_['in',i].unfix()
    model.condenser.V['in'].unfix()
    model.condenser.L['in'].fix(0)
    model.condenser.V['out'].fix(0)
    model.condenser.H_L_.fix(0)
    model.condenser.H_V_.unfix()

    # operating parameters
    model.condenser.P.fix(19)
    model.condenser.T_F.fix(200+273.15)
    model.condenser.F.fix(0)
    model.condenser.z.fix(0)
    model.condenser.VLE_block.n_ave.fix(4)
    model.condenser.PR_L.fix(1)

    model.condenser.T.fix(30+273.15)

    # in/out variables
    model.reboiler.y_.fix(0)
    for i in m.COMP_TOTAL:
        model.reboiler.x_['in',i].unfix()
    model.reboiler.L['in'].unfix()
    model.reboiler.V['in'].fix(0)
    model.reboiler.L['out'].fix(0)
    model.reboiler.V['P'].fix(0)
    model.reboiler.H_L_.unfix()
    model.reboiler.H_V_.fix(0)

    # operating parameters
    model.reboiler.P.fix(20)
    model.reboiler.T_F.fix(200+273.15)
    model.reboiler.F.fix(0)
    model.reboiler.z['CO'].fix(1/(1+2)-0/2)
    model.reboiler.z['H2'].fix(2/(1+2)-0/2)
    model.reboiler.z['C30H62'].fix(0)

    model.reboiler.VLE_block.n_ave.fix(20)

    model.reboiler.T.fix(model.reactive[model.TRAY.last()].T.value)

    # unlock reflux and reboiler vapor
    for j in model.reactive:
        for i in m.COMP_TOTAL:
            model.reactive[j].x_['in',i].unfix()
        model.reactive[j].H_L_['in'].unfix()
        model.reactive[j].L['in'].unfix()
        for i in m.COMP_TOTAL:
            model.reactive[j].y_['in',i].unfix()
        model.reactive[j].V['in'].unfix()
        model.reactive[j].H_V_['in'].unfix()

    for j in model.reactive:
            model.reactive[j].cat.fix(3000)
            model.reactive[j].P.fix(20)
            model.reactive[j].VLE_block.n_ave.fix(20)

            model.reactive[j].F.fix(1)
            model.reactive[j].T_F.fix(200+273.15)
            model.reactive[j].z['CO'].fix(1/(1+2)-0/2)
            model.reactive[j].z['H2'].fix(2/(1+2)-0/2)
            model.reactive[j].z['C30H62'].fix(0)

            model.reactive[j].PR_L.fix(0)
            model.reactive[j].PR_V.fix(0)

            # model.reactive[j].Q_main.fix(0)
            model.reactive[j].T.setub(220+273.15)
            model.reactive[j].T.setlb(200+273.15)


    with HiddenLogs(log_text_dir):
        print('\n> Reactive distillation column assembled, ready to solve, DOF: \n')
        print('-'*108)
        check_DOF(pe,model)

    opt = add_solver(pe, max_iter = 1000, warm_start = True, output = output_dir, option_dir = option_dir, tol=1e-4)

    progress = '> Initialization with no reflux and reboil'

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    if results.solver.termination_condition.key != 'optimal':
        not_optimum_counter += 1
    else:
        not_optimum_counter = 0

    with HiddenLogs(log_text_dir):
        print('\n{}'.format(progress))
        beautify(pe,model)
        log_now()
        check_iteration(output_dir)
        if not_optimum_counter == 1:
            print('NOT Optimum: ' + results.solver.termination_condition.key)
            master_log.write('Failed: {}\n'.format(progress))
            exit()

    # plot_distribution(model,pdf,progress)

    '''
    Introduce reflux, in a gentle way.
    linspace for rr_ratio is not linear in terms of reflux flow, has to be modified
    '''

    R_range = np.arange(0,1/rr_ratio-1,1)
    R_range = np.delete(R_range,0)
    PR_range = 1/(1+R_range)

    initial_reflux_profile = [0.7,0.5]
    PR_range = np.hstack([[r for r in initial_reflux_profile if r > max(PR_range)],PR_range,rr_ratio])


    for r in PR_range:
        model.condenser.PR_L.fix(r)

        progress = '> Working on Reflux, PR ratio = {:.2f}'.format(r)

        results = opt.solve(model,tee=False)
        update_dual(pe,model)

        if results.solver.termination_condition.key != 'optimal':
            not_optimum_counter += 1
        else:
            not_optimum_counter = 0

        with HiddenLogs(log_text_dir):
            print('\n{}'.format(progress))
            print('-'*108)
            beautify(pe,model)
            log_now()
            check_iteration(output_dir)
            if not_optimum_counter == 1:
                print('NOT Optimum: ' + results.solver.termination_condition.key)
                master_log.write('Failed: {}\n'.format(progress))
                exit()

    # plot_distribution(model,pdf,progress)

    '''
    Start the reboiler
    '''

    T_range = np.linspace(model.reboiler.T.value,350+273.15,2)
    for t in T_range[1:]:
        model.reboiler.T.fix(t)

        progress = '> Working on reboiler temperature = {:.2f}'.format(t)

        results = opt.solve(model,tee=False)
        update_dual(pe,model)

        if results.solver.termination_condition.key != 'optimal':
            not_optimum_counter += 1
        else:
            not_optimum_counter = 0

        with HiddenLogs(log_text_dir):
            print('\n{}'.format(progress))
            print('-'*108)
            beautify(pe,model)
            log_now()
            check_iteration(output_dir)
            if not_optimum_counter == 1:
                print('NOT Optimum: ' + results.solver.termination_condition.key)
                master_log.write('Failed: {}\n'.format(progress))
                exit()

    # # plot_distribution(model,pdf,'Reboiler')


    '''
    Extract side-draws
    '''

    for j in side_draw_flag.keys():

        r = side_draw_flag[j]
        model.reactive[j].PR_L.fix(r)

        progress = '> Working on side draw of {:.1%} on stage {}'.format(r,j)

        results = opt.solve(model,tee=False)
        update_dual(pe,model)

        if results.solver.termination_condition.key != 'optimal':
            not_optimum_counter += 1
        else:
            not_optimum_counter = 0

        with HiddenLogs(log_text_dir):
            print('\n{}'.format(progress))
            print('-'*108)
            beautify(pe,model)
            log_now()
            check_iteration(output_dir)
            if not_optimum_counter == 1:
                print('NOT Optimum: ' + results.solver.termination_condition.key)
                master_log.write('Failed: {}\n'.format(progress))
                exit()

    # plot_distribution(model,pdf,'Side draw')

    '''
    Following design specification, remove non-reactive stages' catalyst and feed, feed cat first then Q
    '''
    for j in model.TRAY_nonreactive:

        model.reactive[j].cat.fix(catalyst_flag[j])
        model.reactive[j].F.fix(feed_flag[j])
        # During this step, it is very likely for the temperature to not conform to the given profile
        # Therefore, add an additional layer of protection with temperature lower bounds
        # model.reactive[j].T.setlb(model.reactive[j].T.ub - 1)

        progress = '> Working to adjust catalyst and feed {}:'.format(j)

        results = opt.solve(model,tee=False)
        update_dual(pe,model)

        if results.solver.termination_condition.key != 'optimal':
            not_optimum_counter += 1
        else:
            not_optimum_counter = 0

        with HiddenLogs(log_text_dir):
            print('\n{}'.format(progress))
            print('-'*108)
            beautify(pe,model)
            log_now()
            check_iteration(output_dir)
            if not_optimum_counter == 1:
                print('NOT Optimum: ' + results.solver.termination_condition.key)
                master_log.write('Failed: {}\n'.format(progress))
                exit()

    for alpha in np.arange(0,1,0.04)+0.04:
        for j in model.TRAY_reactive:

            model.reactive[j].cat.fix(catalyst_flag[j]*alpha + 3000*(1-alpha))
            model.reactive[j].F.fix(feed_flag[j]*alpha + 1*(1-alpha))
        # During this step, it is very likely for the temperature to not conform to the given profile
        # Therefore, add an additional layer of protection with temperature lower bounds
        # model.reactive[j].T.setlb(model.reactive[j].T.ub - 1)

        progress = '> Working to adjust catalyst and feed, alpha = {:.2f}:'.format(alpha)

        results = opt.solve(model,tee=False)
        update_dual(pe,model)

        if results.solver.termination_condition.key != 'optimal':
            not_optimum_counter += 1
        else:
            not_optimum_counter = 0

        with HiddenLogs(log_text_dir):
            print('\n{}'.format(progress))
            print('-'*108)
            beautify(pe,model)
            log_now()
            check_iteration(output_dir)
            if not_optimum_counter == 1:
                print('NOT Optimum: ' + results.solver.termination_condition.key)
                master_log.write('Failed: {}\n'.format(progress))
                exit()

    tmp_expr = sum(model.reactive[j].T for j in model.TRAY)
    for j in model.TRAY_nonreactive:

        tmp_expr -= model.reactive[j].T

        model.reactive[j].T.unfix()
        model.reactive[j].cat.fix(0)
        model.reactive[j].F.fix(0)
        model.reactive[j].T.setub(350+273.15)
        model.reactive[j].T.setlb(100+273.15)

        model.del_component(model.obj)
        model.obj = augmented_objective(pe,model,expr = tmp_expr ,sense=pe.maximize)

        # iterativly change Q
        Q_steps = int(abs(model.reactive[j].Q_main.value)/sum(model.reactive[j].L[s].value+\
                    model.reactive[j].V[s].value for s in model.reactive[j].outlet)/2 + 1)
        Q_range = np.linspace(model.reactive[j].Q_main.value,0,Q_steps+1)
        Q_range = np.delete(Q_range,0)

        for q in Q_range:
            model.reactive[j].Q_main.fix(q)

            progress = '> Working on stage {:}, changing Q to {:.2f}:'.format(j,q)

            results = opt.solve(model,tee=False)
            update_dual(pe,model)

            if results.solver.termination_condition.key != 'optimal':
                not_optimum_counter += 1
            else:
                not_optimum_counter = 0

            with HiddenLogs(log_text_dir):
                print('\n{}'.format(progress))
                print('-'*108)
                beautify(pe,model)
                log_now()
                check_iteration(output_dir)
                if not_optimum_counter == 1:
                    print('NOT Optimum: ' + results.solver.termination_condition.key)
                    master_log.write('Failed: {}\n'.format(progress))
                    exit()

    # plot_distribution(model,pdf,'Non-reactive stage + catalyst / feed adjustment')

    '''
    Move temperatures to design point (gently)
    '''

    for j in sorted(temperature_flag,reverse=True):

        t_set = temperature_flag[j] + 273.15
        trange = np.arange(model.reactive[j].T.value,t_set,10)

        if trange.size > 0:
            trange = np.delete(trange,0)
            trange = np.hstack([trange,t_set])

            for t in trange:
                model.reactive[j].T.setub(t)
                # model.reactive[j].T.setlb(model.reactive[j].T.ub - 10)

                progress = '> Working on adjusting stage {} temperature to {:.2f}C '.format(j,t-273.15)

                results = opt.solve(model,tee=False)
                update_dual(pe,model)

                if results.solver.termination_condition.key != 'optimal':
                    not_optimum_counter += 1
                else:
                    not_optimum_counter = 0

                with HiddenLogs(log_text_dir):
                    print('\n{}'.format(progress))
                    print('-'*108)
                    beautify(pe,model)
                    log_now()
                    check_iteration(output_dir)
                    if not_optimum_counter == 1:
                        print('NOT Optimum: ' + results.solver.termination_condition.key)
                        master_log.write('Failed: {}\n'.format(progress))
                        exit()

        trange = np.arange(model.reactive[j].T.value,t_set,-10)

        if trange.size > 0:
            trange = np.delete(trange,0)
            trange = np.hstack([trange,t_set])

            for t in trange:
                model.reactive[j].T.setub(t)
                # model.reactive[j].T.setlb(model.reactive[j].T.ub - 10)

                progress = '> Working on adjusting stage {} temperature to {:.2f}C '.format(j,t-273.15)

                results = opt.solve(model,tee=False)
                update_dual(pe,model)

                if results.solver.termination_condition.key != 'optimal':
                    not_optimum_counter += 1
                else:
                    not_optimum_counter = 0

                with HiddenLogs(log_text_dir):
                    print('\n{}'.format(progress))
                    print('-'*108)
                    beautify(pe,model)
                    log_now()
                    check_iteration(output_dir)
                    if not_optimum_counter == 1:
                        print('NOT Optimum: ' + results.solver.termination_condition.key)
                        master_log.write('Failed: {}\n'.format(progress))
                        exit()

    # plot_distribution(model,pdf,'Finalized Stage Temperatures')

    with HiddenLogs(log_text_dir):
        if results.solver.termination_condition.key != 'optimal':
            print('Initialization NOT Optimum: ')
            master_log.write('Failed: {}\n'.format(progress))
            exit()
        else:
            print('Initialization Complete\nPlease check the logs for details')
            master_log.write('Success: > Initialization\n')

    '''
    Now move on to optimization
    The first thing to do is to add on new components necessary for tray optimization
    The reason is to maintain the above simulator's functionality for exact modelling, DDF is not exact
    '''

    '''
    Redefine MPCC type

    '''
    for j in model.reactive:
        select_MPCC(model.reactive[j],'pf')
    select_MPCC(model.reboiler,'pf')

    model.sigma = pe.Param(initialize=0.5,mutable=True)
    model.epi = pe.Param(initialize=1e-5,mutable=True)
    model.scale_epi = pe.Param(initialize=10,mutable=True)

    model.P_tray = pe.Var(model.TRAY_total,m.PRODUCT,within=pe.NonNegativeReals,initialize=0)
    model.N_tray = pe.Var(m.PRODUCT,within=pe.NonNegativeReals) # extended range from condenser (0) to reboiler (N+1)
    model.P_total = pe.Var(m.PRODUCT,within=pe.NonNegativeReals,initialize=0)
    model.P_total_dry = pe.Var(m.PRODUCT,within=pe.NonNegativeReals,initialize=0)
    model.x_P = pe.Var(m.COMP_TOTAL,m.PRODUCT,within=pe.NonNegativeReals,bounds=(0,1))
    model.x_P_dry = pe.Var(m.COMP_ORG,m.PRODUCT,within=pe.NonNegativeReals,bounds=(0,1))

    for i,j in model.P_tray:
        if j != 'naphtha' and j != 'heavy':
            model.P_tray[i,j].setlb(model.epi)

    for j in model.P_total:
        if j != 'naphtha' and j != 'heavy':
            model.P_total[j].setlb(model.epi*len(model.TRAY_total))
        else:
            continue
            model.P_total[j].setlb(model.epi)

    # sum of liquid draw for all products for each stage
    def stage_sum_product_rule(model,j):
        return tray_translator(model,j).L['P'] == sum(model.P_tray[j,p] for p in m.PRODUCT)
    model.stage_sum_product_con = pe.Constraint(model.TRAY_total,rule=stage_sum_product_rule)

    # liquid product mass balance
    def product_sum_stage_rule(model,j,p):
        if p == 'naphtha' or p == 'heavy':
            return pe.Constraint.Skip
        return pe.log(model.scale_epi + (model.P_tray[j,p] - model.epi) * sum(pe.exp(-(model.TRAY_total.ord(j_)-1-model.N_tray[p])**2/model.sigma) for j_ in model.TRAY_total))\
                == pe.log(model.scale_epi + (model.P_total[p] - model.epi*len(model.TRAY_total)) * pe.exp(-(model.TRAY_total.ord(j)-1-model.N_tray[p])**2/model.sigma))
    model.product_sum_stage_con = pe.Constraint(model.TRAY_total,m.PRODUCT,rule=product_sum_stage_rule)

    # condenser and reboiler
    model.product_sum_stage_con2 = pe.ConstraintList()
    model.product_sum_stage_con2.add(expr = model.P_tray['condenser','naphtha'] == model.P_total['naphtha'])
    model.product_sum_stage_con2.add(expr = model.P_tray['reboiler','heavy'] == model.P_total['heavy'])

    # liquid product component mass balance
    def mass_balance_product_rule(model,i,p):
        return sum(model.P_tray[j,p]*tray_translator(model,j).x[i] for j in model.TRAY_total) == (model.P_total[p])*model.x_P[i,p]
    model.mass_balance_product_con = pe.Constraint(m.COMP_TOTAL,m.PRODUCT,rule=mass_balance_product_rule)

    # dry liquid product component
    def product_sum_dry_rule(model,p):
        return model.P_total_dry[p] == model.P_total[p] * (1 - sum(model.x_P[i,p] for i in m.COMP_INORG))
    model.product_sum_dry_con = pe.Constraint(m.PRODUCT,rule=product_sum_dry_rule)

    # dry liquid product component mass balance
    def mass_balance_dry_rule(model,i,p):
        return model.x_P_dry[i,p] * (1 - sum(model.x_P[i,p] for i in m.COMP_INORG)) == model.x_P[i,p]
    model.mass_balance_dry_con = pe.Constraint(m.COMP_ORG,m.PRODUCT,rule=mass_balance_dry_rule)

    # initialize with values
    for j in model.reactive:
        model.reactive[j].PR_L.unfix();

    for j in model.TRAY_total:
        model.P_tray[j,'naphtha'].fix(0)
        model.P_tray[j,'heavy'].fix(0)

    model.P_tray['condenser','naphtha'].unfix(); # model.P_tray['condenser','naphtha'].setlb(model.epi)
    model.P_tray['condenser','naphtha'].set_value(model.condenser.L['P'].value)
    model.P_tray[str(sorted(side_draw_flag)[0]),'intermediate'].set_value(model.reactive[sorted(side_draw_flag)[0]].L['P'].value)
    model.P_tray[str(sorted(side_draw_flag)[1]),'gasoline'].set_value(model.reactive[sorted(side_draw_flag)[1]].L['P'].value)
    model.P_tray[str(sorted(side_draw_flag)[2]),'diesel'].set_value(model.reactive[sorted(side_draw_flag)[2]].L['P'].value)
    model.P_tray['reboiler','heavy'].unfix(); # model.P_tray['reboiler','heavy'].setlb(model.epi)
    model.P_tray['reboiler','heavy'].set_value(model.reboiler.L['P'].value)

    model.N_tray['naphtha'].fix(0)
    model.N_tray['intermediate'].fix(sorted(side_draw_flag)[0])
    model.N_tray['gasoline'].fix(sorted(side_draw_flag)[1])
    model.N_tray['diesel'].fix(sorted(side_draw_flag)[2])
    model.N_tray['heavy'].fix(21)

    model.P_total['naphtha'].set_value(model.P_tray['condenser','naphtha'].value)
    model.P_total['intermediate'].fix(model.P_tray[str(sorted(side_draw_flag)[0]),'intermediate'].value)
    model.P_total['gasoline'].fix(model.P_tray[str(sorted(side_draw_flag)[1]),'gasoline'].value)
    model.P_total['diesel'].fix(model.P_tray[str(sorted(side_draw_flag)[2]),'diesel'].value)
    model.P_total['heavy'].set_value(model.P_tray['reboiler','heavy'].value)

    for i in m.COMP_TOTAL:
        model.x_P[i,'naphtha'].set_value(model.condenser.x[i].value)
        model.x_P[i,'intermediate'].set_value(model.reactive[sorted(side_draw_flag)[0]].x[i].value)
        model.x_P[i,'gasoline'].set_value(model.reactive[sorted(side_draw_flag)[1]].x[i].value)
        model.x_P[i,'diesel'].set_value(model.reactive[sorted(side_draw_flag)[2]].x[i].value)
        model.x_P[i,'heavy'].set_value(model.reboiler.x[i].value)

    model.P_total_dry['naphtha'].set_value(model.P_tray['condenser','naphtha'].value * (1 - sum(model.condenser.x[i].value for i in m.COMP_INORG)))
    model.P_total_dry['intermediate'].set_value(model.P_tray[str(sorted(side_draw_flag)[0]),'intermediate'].value\
                                                * (1 - sum(model.reactive[sorted(side_draw_flag)[0]].x[i].value for i in m.COMP_INORG)))
    model.P_total_dry['gasoline'].set_value(model.P_tray[str(sorted(side_draw_flag)[1]),'gasoline'].value\
                                            * (1 - sum(model.reactive[sorted(side_draw_flag)[1]].x[i].value for i in m.COMP_INORG)))
    model.P_total_dry['diesel'].set_value(model.P_tray[str(sorted(side_draw_flag)[2]),'diesel'].value\
                                          * (1 - sum(model.reactive[sorted(side_draw_flag)[2]].x[i].value for i in m.COMP_INORG)))
    model.P_total_dry['heavy'].set_value(model.P_tray['reboiler','heavy'].value * (1 - sum(model.reboiler.x[i].value for i in m.COMP_INORG)))

    for i in m.COMP_ORG:
        model.x_P_dry[i,'naphtha'].set_value(model.condenser.x[i].value / (1 - sum(model.condenser.x[i].value for i in m.COMP_INORG)))
        model.x_P_dry[i,'intermediate'].set_value(model.reactive[sorted(side_draw_flag)[0]].x[i].value / (1 - sum(model.reactive[sorted(side_draw_flag)[0]].x[i].value for i in m.COMP_INORG)))
        model.x_P_dry[i,'gasoline'].set_value(model.reactive[sorted(side_draw_flag)[1]].x[i].value / (1 - sum(model.reactive[sorted(side_draw_flag)[1]].x[i].value for i in m.COMP_INORG)))
        model.x_P_dry[i,'diesel'].set_value(model.reactive[sorted(side_draw_flag)[2]].x[i].value / (1 - sum(model.reactive[sorted(side_draw_flag)[2]].x[i].value for i in m.COMP_INORG)))
        model.x_P_dry[i,'heavy'].set_value(model.reboiler.x[i].value / (1 - sum(model.reboiler.x[i].value for i in m.COMP_INORG)))

    model.del_component(model.obj)
    model.obj = augmented_objective(pe,model,expr = sum(model.reactive[j].T for j in model.TRAY_reactive), sense = pe.maximize)

    for j in model.reactive:
        model.reactive[j].MPCC_P_pf.rho = 100000
    model.reboiler.MPCC_P_pf.rho = 100000

    disable_restoration(mode = 'disable', option_dir = option_dir)

    opt = add_solver(pe, max_iter = 1000, warm_start = True, output = output_dir, option_dir = option_dir)

    progress = '> Added DDF formulation - Product'

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    with HiddenLogs(log_text_dir):
        print('\n{}'.format(progress))
        print('-'*108)
        beautify(pe,model)
        check_product_spec(model)
        log_now()
        check_iteration(output_dir)
        if results.solver.termination_condition.key != 'optimal':
            print('NOT Optimum: ' + results.solver.termination_condition.key)
            master_log.write('Failed: {}\n'.format(progress))
            exit()
        else:
            print('DDF Complete\nPlease check the logs for details')
            master_log.write('Success: {}\n'.format(progress))

    '''
    Optimization
    '''

    for j in model.TRAY_reactive:
        model.reactive[j].T.setlb(200+273.15)
        model.reactive[j].T.setub(300+273.15)

    model.condenser.PR_L.unfix()
    model.condenser.PR_L.setlb(0.05)
    model.condenser.PR_L.setub(0.5)

    model.P_total['intermediate'].unfix()
    model.P_total['gasoline'].unfix()
    model.P_total['diesel'].unfix()

    model.quality_spec = pe.Param(m.PRODUCT,initialize={\
                        'naphtha':0.75,'gasoline':0.75,'diesel':0.75,'heavy':0.75},mutable=True)

    def product_spec_rule(model,p):
        if p == 'intermediate':
            return pe.Constraint.Skip
        return sum(model.x_P_dry[i,p] for i in m.PRODUCT_cnumber[p]) >= model.quality_spec[p]
    model.product_spec_con = pe.Constraint(m.PRODUCT,rule=product_spec_rule)

    model.total_feed = pe.Var(within=pe.NonNegativeReals,bounds=(1,30),initialize=10)
    model.total_feed_con = pe.ConstraintList()
    model.total_feed_con.add(expr = sum(model.reactive[j].F for j in model.TRAY_reactive) == model.total_feed);

    model.total_feed.fix(10)
    for j in model.TRAY_reactive:
        model.reactive[j].F.unfix()
        model.reactive[j].F.setlb(1e-3)
        model.reactive[j].F.setub(10)

    for j in model.TRAY_reactive:
        model.reactive[j].cat.unfix()
        model.reactive[j].cat.setlb(10)
        model.reactive[j].cat.setub(30000)

    model.total_cat_con = pe.ConstraintList()
    model.total_cat_con.add(expr = sum(model.reactive[j].cat for j in model.reactive) == 10*3000);

    model.N_tray['gasoline'].unfix();
    model.N_tray['gasoline'].setlb(3)
    model.N_tray['gasoline'].setub(9)

    model.N_tray['diesel'].unfix();
    model.N_tray['diesel'].setlb(10)
    model.N_tray['diesel'].setub(18)

    with HiddenPrints():
        delete_dual(pe,model)
        add_dual(pe,model)

    '''save a version that could be reused
    '''
    # master_model = deepcopy(model)

    model.del_component(model.obj)
    model.obj = augmented_objective(pe,model,expr = \
                                    43*model.P_total['naphtha']+\
                                    20*model.P_total['intermediate']+\
                                    90*model.P_total['gasoline']+\
                                    128*model.P_total['diesel']+\
                                    100*model.P_total['heavy'], sense = pe.maximize)

    progress = '> One-step Optimization - Revenue'

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    with HiddenLogs(log_text_dir):
        print('\n{}'.format(progress))
        print('-'*108)
        print('obj',model.obj())
        beautify(pe,model)
        check_product_spec(model)
        log_now()
        check_iteration(output_dir)
        if results.solver.termination_condition.key != 'optimal':
            print('NOT Optimum: ' + results.solver.termination_condition.key)
            master_log.write('Failed: {}\n'.format(progress))
            exit()
        else:
            print('Optimization Complete\nPlease check the logs for details')
            master_log.write('Success: {}\n'.format(progress))

    # plot_distribution(model,pdf,'Optimized Temperature, reflux, product flow and tray, feed, catalyst')
    # plot_product_distribution(model,pdf)

    '''profit opt
    '''
    # model = deepcopy(master_model)

    '''
    First thing, again, is to transfer reflux to DDF formulation
    '''

    model.sigma_reflux = pe.Param(initialize=0.5,mutable=True)
    model.N_reflux_tray = pe.Var(within=pe.NonNegativeReals,bounds=(1,5),initialize=1) # bounded by upper section non-reactive trays

    model.del_component(model.L_condenser_con)
    model.del_component(model.Lx_condenser_con)
    model.del_component(model.Lh_condenser_con)

    def L_condenser_rule(model,j):
        return pe.log(model.scale_epi + (model.reactive[j].L['R'] - model.epi)*sum(pe.exp(-(model.TRAY_total.ord(str(j_))-1-model.N_reflux_tray)**2\
                /model.sigma_reflux) for j_ in model.TRAY)) == pe.log(model.scale_epi + (model.condenser.L['out'] - model.epi*len(model.TRAY)) * \
                pe.exp(-(model.TRAY_total.ord(str(j))-1-model.N_reflux_tray)**2/model.sigma_reflux))
    model.L_condenser_con = pe.Constraint(model.TRAY,rule=L_condenser_rule)

    def Lx_condenser_rule(model,j,i):
        return model.reactive[j].x_['R',i] == model.condenser.x[i]
    model.Lx_condenser_con = pe.Constraint(model.TRAY,m.COMP_TOTAL,rule=Lx_condenser_rule)

    def Lh_condenser_rule(model,j):
        return model.reactive[j].H_L_['R'] == model.condenser.H_L
    model.Lh_condenser_con = pe.Constraint(model.TRAY,rule=Lh_condenser_rule)

    ''' Initialize
    '''

    model.reactive[model.TRAY.first()].L['in'].fix(0)
    for i in m.COMP_TOTAL:
        model.reactive[model.TRAY.first()].x_['in',i].fix(0)
    model.reactive[model.TRAY.first()].H_L_['in'].fix(0)

    for j in model.reactive:
        model.reactive[j].L['R'].unfix(); model.reactive[j].L['R'].setlb(model.epi.value); model.reactive[j].L['R'].set_value(1.1*model.epi.value)
        for i in m.COMP_TOTAL:
            model.reactive[j].x_['R',i].unfix(); model.reactive[j].x_['R',i].set_value(model.condenser.x[i].value)
        model.reactive[j].H_L_['R'].unfix(); model.reactive[j].H_L_['R'].set_value(model.condenser.H_L.value)

    model.N_reflux_tray.fix(1)
    model.reactive[model.TRAY.first()].L['R'].set_value(model.condenser.L['out'].value)

    progress = '> Added DDF formulation - Reflux'

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    with HiddenLogs(log_text_dir):
        print('\n{}'.format(progress))
        print('-'*108)
        beautify(pe,model)
        check_product_spec(model)
        log_now()
        check_iteration(output_dir)
        if results.solver.termination_condition.key != 'optimal':
            print('NOT Optimum: ' + results.solver.termination_condition.key)
            master_log.write('Failed: {}\n'.format(progress))
            exit()
        else:
            print('Optimization Complete\nPlease check the logs for details')
            master_log.write('Success: {}\n'.format(progress))


    ''' optimize
    '''

    model.N_reflux_tray.unfix();
    model.total_feed.unfix();

    model.del_component(model.obj)
    model.obj = augmented_objective(pe,model,expr = \
                                    43*model.P_total['naphtha']+\
                                    20*model.P_total['intermediate']+\
                                    90*model.P_total['gasoline']+\
                                    128*model.P_total['diesel']+\
                                    100*model.P_total['heavy']+\
                                    1.3*model.condenser.V['P']-\
                                    2.24*model.total_feed+\
                                    0.5*(model.N_reflux_tray-1), sense = pe.maximize)

    progress = '> One-step Optimization - Profit 3-1'

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    with HiddenLogs(log_text_dir):
        print('\n{}'.format(progress))
        print('-'*108)
        print('obj',model.obj())
        beautify(pe,model)
        check_product_spec(model)
        print('Reflux Tray Location:',model.N_reflux_tray.value)
        print('Total Feed:',model.total_feed.value)
        log_now()
        check_iteration(output_dir)

        if results.solver.termination_condition.key != 'optimal':
            print('NOT Optimum: ' + results.solver.termination_condition.key)
            master_log.write('Failed: {}\n'.format(progress))
        else:
            print('Optimization Complete\nPlease check the logs for details')
            master_log.write('Success: {}\n'.format(progress))

    '''3-2
    '''

    model.del_component(model.obj)
    model.obj = augmented_objective(pe,model,expr = \
                                    43*model.P_total['naphtha']+\
                                    20*model.P_total['intermediate']+\
                                    90*model.P_total['gasoline']+\
                                    128*model.P_total['diesel']+\
                                    100*model.P_total['heavy']+\
                                    1.3*model.condenser.V['P']-\
                                    1.5*model.total_feed+\
                                    0.005*(model.N_reflux_tray-1), sense = pe.maximize)

    progress = '> One-step Optimization - Profit 3-2'

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    with HiddenLogs(log_text_dir):
        print('\n{}'.format(progress))
        print('-'*108)
        print('obj',model.obj())
        beautify(pe,model)
        check_product_spec(model)
        print('Reflux Tray Location:',model.N_reflux_tray.value)
        print('Total Feed:',model.total_feed.value)
        log_now()
        check_iteration(output_dir)

        if results.solver.termination_condition.key != 'optimal':
            print('NOT Optimum: ' + results.solver.termination_condition.key)
            master_log.write('Failed: {}\n'.format(progress))
        else:
            print('Optimization Complete\nPlease check the logs for details')
            master_log.write('Success: {}\n'.format(progress))

    '''3-3
    '''

    model.del_component(model.obj)
    model.obj = augmented_objective(pe,model,expr = \
                                    43*model.P_total['naphtha']+\
                                    20*model.P_total['intermediate']+\
                                    90*model.P_total['gasoline']+\
                                    128*model.P_total['diesel']+\
                                    100*model.P_total['heavy']+\
                                    1.3*model.condenser.V['P']-\
                                    3.5*model.total_feed+\
                                    0.005*(model.N_reflux_tray-1), sense = pe.maximize)

    progress = '> One-step Optimization - Profit 3-2'

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    with HiddenLogs(log_text_dir):
        print('\n{}'.format(progress))
        print('-'*108)
        print('obj',model.obj())
        beautify(pe,model)
        check_product_spec(model)
        print('Reflux Tray Location:',model.N_reflux_tray.value)
        print('Total Feed:',model.total_feed.value)
        log_now()
        check_iteration(output_dir)

        if results.solver.termination_condition.key != 'optimal':
            print('NOT Optimum: ' + results.solver.termination_condition.key)
            master_log.write('Failed: {}\n'.format(progress))
        else:
            print('Optimization Complete\nPlease check the logs for details')
            master_log.write('Success: {}\n'.format(progress))

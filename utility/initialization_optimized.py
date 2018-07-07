# %matplotlib inline
# %load ../../utility/initialization.py
'''
How to Train Your Dragon: V2
Sequentially initialize FT reactive distillation model automatically
'''
# system imports
import sys
import os
import datetime
sys.path.append(os.path.abspath('..'))
sys.path.append(os.path.abspath('../..'))

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pickle
import dill
from copy import deepcopy

# pyomo imports
from pyomo import environ as pe
from global_sets.component import m

from stages.reactive_stage import reactive_stage_rule
from stages.condenser_stage import condenser_stage_rule
from stages.reboiler_stage import reboiler_stage_rule

from utility.display_utility import trans_product_mole, trans_product_mass, beautify, \
                                    beautify_reactive, HiddenLogs, HiddenPrints, plot_distribution
from utility.model_utility import add_dual, update_dual, delete_dual, check_DOF, check_violate_constraint
from utility.data_utility import cal_cnumber
from utility.time_utility import create_filename_time, log_now, log_end

'''
Constructing the model and logfile
'''

model = pe.ConcreteModel(name='reactive_distillation')
logname = create_filename_time()
log_text_dir = './log/text/'+logname+'.dat'
log_figure_dir = './log/figure/'+logname+'.pdf'

'''
model input
'''
tray_number = 10
non_reactive_flag = [1,2,3,5,10]
# non_reactive_flag = []
side_draw_flag = {4:0.2,7:0.7}
# default temperature is 220C
temperature_flag = {6:240,7:240,8:245,9:250}
rr_ratio = 0.05

with PdfPages(log_figure_dir,keep_empty=True) as pdf:
    '''
    construct the reactive stages
    '''
    model.TRAY = pe.RangeSet(1,tray_number)
    model.TRAY_nonreactive = pe.Set(initialize=non_reactive_flag)
    model.TRAY_reactive = model.TRAY - model.TRAY_nonreactive

    with HiddenPrints():
        model.reactive = pe.Block(model.TRAY,rule=reactive_stage_rule)
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

    model.obj = pe.Objective(expr = sum(model.reactive[j].T - model.reactive[j].MPCC.pf \
                                        for j in model.reactive),sense=pe.maximize)

    '''
    solver options
    '''
    try:
        opt = pe.SolverFactory('ipopt')

        opt.options['print_user_options'] = 'yes'
        opt.options['linear_solver'] = 'ma86'

        opt.options['linear_system_scaling '] = 'mc19'
        opt.options['linear_scaling_on_demand '] = 'no'

        opt.options['max_iter'] = 7000

        results = opt.solve(model,tee=False)
        update_dual(pe,model)
    except:
        opt = None
        opt = pe.SolverFactory('ipopt')

        opt.options['print_user_options'] = 'yes'
        opt.options['linear_scaling_on_demand '] = 'no'

        opt.options['max_iter'] = 7000

        results = opt.solve(model,tee=False)
        update_dual(pe,model)


    with HiddenLogs(log_text_dir,'w'):
        print('\n> First Solve, disconnected reactive stages')
        print('-'*108)
        beautify_reactive(pe,model)
        log_now()

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
            model.reactive[j].x_.unfix()
            model.reactive[j].H_L_.unfix()
            model.reactive[j].L['in'].unfix()
        if j != model.TRAY.last():
            model.reactive[j].y_.unfix()
            model.reactive[j].V['in'].unfix()
            model.reactive[j].H_V_.unfix()

    '''
    warm start options
    '''

    opt.options['warm_start_init_point'] = 'yes'
    opt.options['warm_start_bound_push'] = 1e-20
    opt.options['warm_start_mult_bound_push'] = 1e-20
    opt.options['mu_init'] = 1e-6

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    '''
    The model above has liquid and vapor leave stage as products
    Now let liquid and vapor flow up/down stages
    '''

    PR_range = np.linspace(0.8,0,5)
    with HiddenLogs(log_text_dir):
        print('\n> Connect stages and introduce inter-stage flow')
        log_now()
    for r in PR_range:
        for j in model.reactive:
            model.reactive[j].PR_L.fix(r)
            model.reactive[j].PR_V.fix(r)

        results = opt.solve(model,tee=False)
        update_dual(pe,model)
        with HiddenLogs(log_text_dir):
            print('\n>','Working on PR ratio = {:.2f}'.format(r))
            print('-'*108)
            beautify_reactive(pe,model)
            log_now()

    '''
    Deactivate reactive part of the model to initialize condenser and reboiler
    '''

    for i in model.block_data_objects():
        if i.name != 'reactive_distillation':
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
    model.condenser.V['P'].fix(0)
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
    model.obj = pe.Objective(expr = model.condenser.T, sense = pe.maximize)

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    model.condenser.deactivate()

    with HiddenPrints():
        model.reboiler = pe.Block(rule=reboiler_stage_rule)
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
    model.reboiler.L['P'].fix(0)
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
    model.obj = pe.Objective(expr = model.reboiler.T - model.reboiler.MPCC.pf, sense = pe.maximize)

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    '''
    Now all the variables are initialized, reactivate everything, time to connect condenser and reboiler
    '''

    for i in model.block_data_objects():
        if i.name != 'reactive_distillation':
            i.activate()
    for i in model.component_objects(pe.Constraint):
        i.activate()

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

    def L_reboiler_rule(model):
        return model.reactive[model.TRAY.last()].L['out'] + 1e-6 == model.reboiler.L['in']
    model.L_reboiler_con = pe.Constraint(rule=L_reboiler_rule)

    def Lx_reboiler_rule(model,i):
        return model.reactive[model.TRAY.last()].x[i] == model.reboiler.x_['in',i]
    model.Lx_reboiler_con = pe.Constraint(m.COMP_TOTAL,rule=Lx_reboiler_rule)

    def Lh_reboiler_rule(model):
        return model.reactive[model.TRAY.last()].H_L == model.reboiler.H_L_['in']
    model.Lh_reboiler_con = pe.Constraint(rule=Lh_reboiler_rule)

    model.del_component(model.obj)
    model.obj = pe.Objective(expr = sum(model.reactive[j].T - model.reactive[j].MPCC.pf for j in model.TRAY)\
                             - model.reboiler.MPCC.pf ,sense=pe.maximize)

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
    model.condenser.V['P'].fix(0)
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
    model.reboiler.L['P'].fix(0)
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
        model.reactive[j].x_.unfix()
        model.reactive[j].H_L_.unfix()
        model.reactive[j].L['in'].unfix()
        model.reactive[j].y_.unfix()
        model.reactive[j].V['in'].unfix()
        model.reactive[j].H_V_.unfix()

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

    results = opt.solve(model,tee=False)
    update_dual(pe,model)

    with HiddenLogs(log_text_dir):
        print('\n> Initialization with no reflux and reboil complete')
        beautify(pe,model)
        log_now()

    plot_distribution(model,pdf,'Initialized with no reflux and reboil')

    '''
    Introduce reflux, in a gentle way.
    linspace for rr_ratio is not linear in terms of reflux flow, has to be modified
    '''

    R_range = np.arange(0,1/rr_ratio-1,3)
    R_range = np.delete(R_range,0)
    PR_range = 1/(1+R_range)

    initial_reflux_profile = [0.7,0.5]
    PR_range = np.hstack([[r for r in initial_reflux_profile if r > max(PR_range)],PR_range])


    for r in PR_range:
        model.condenser.PR_L.fix(r)

        results = opt.solve(model,tee=False)
        update_dual(pe,model)

        with HiddenLogs(log_text_dir):
            print('\n>','Working on Reflux, PR ratio = {:.2f}'.format(r))
            print('-'*108)
            beautify(pe,model)
            log_now()

        plot_distribution(model,pdf,'Refulx PR ratio = {:.2f}'.format(r))

    '''
    Start the reboiler
    '''

    T_range = np.linspace(model.reboiler.T.value,350+273.15,2)
    for t in T_range[1:]:
        model.reboiler.T.fix(t)

        results = opt.solve(model,tee=False)
        update_dual(pe,model)

        with HiddenLogs(log_text_dir):
            print('\n>','Working on reboiler temperature = {:.2f}'.format(t))
            print('-'*108)
            beautify(pe,model)
            log_now()

    plot_distribution(model,pdf,'Reboiler')


    '''
    Following design specification, remove non-reactive stages' catalyst and feed
    '''
    tmp_expr = sum(model.reactive[j].T - model.reactive[j].MPCC.pf for j in model.TRAY) - model.reboiler.MPCC.pf
    for j in model.TRAY_nonreactive:

        tmp_expr -= model.reactive[j].T

        model.reactive[j].T.unfix()
        model.reactive[j].Q_main.fix(0)
        model.reactive[j].cat.fix(0)
        model.reactive[j].F.fix(0)
        model.reactive[j].T.setub(350+273.15)
        model.reactive[j].T.setlb(100+273.15)

        model.del_component(model.obj)
        model.obj = pe.Objective(expr = tmp_expr ,sense=pe.maximize)

        results = opt.solve(model,tee=False)
        update_dual(pe,model)
        with HiddenLogs(log_text_dir):
            print('\n>','Working to remove catalyst and feed from stage {}, unfixing temperature, changing to adiabatic:'.format(j))
            print('-'*108)
            beautify_reactive(pe,model)
            log_now()

        plot_distribution(model,pdf,'Non-reactive stage {}'.format(j))


    '''
    Extract side-draws
    '''

    for j in side_draw_flag.keys():

        r = side_draw_flag[j]
        model.reactive[j].PR_L.fix(r)

        results = opt.solve(model,tee=False)
        update_dual(pe,model)

        with HiddenLogs(log_text_dir):
            print('\n>','Working on side draw of {:.1%} on stage {}'.format(r,j))
            print('-'*108)
            beautify(pe,model)
            log_now()

        plot_distribution(model,pdf,'Side draw of {:.1%} on stage {}'.format(r,j))

    '''
    Move temperatures to design point (gently)
    '''

    for j in temperature_flag.keys():

        t_set = temperature_flag[j] + 273.15
        trange = np.arange(model.reactive[j].T.value,t_set,20)
        trange = np.delete(trange,0)
        trange = np.hstack([trange,t_set])

        for t in trange:
            model.reactive[j].T.setub(t)

            results = opt.solve(model,tee=False)
            update_dual(pe,model)

            with HiddenLogs(log_text_dir):
                print('\n>','Working on adjusting stage {} temperature to {:.2f}C '.format(j,t-273.15))
                print('-'*108)
                beautify(pe,model)
                log_now()

    plot_distribution(model,pdf,'Finalized Stage Temperatures')

    with HiddenLogs(log_text_dir):
        log_end()

    with open('./log/model/{}.pickle'.format(logname),'wb') as f:
        dill.dump(model,f)

    with HiddenLogs(log_text_dir):
        print('Initialization Complete\nPlease check the logs for details')

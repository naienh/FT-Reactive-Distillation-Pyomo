from global_sets.component import m
import os, sys
import matplotlib
if "DISPLAY" not in os.environ:
    matplotlib.use("Agg")
import copy
import numpy as np
from utility.data_utility import cal_MW, cal_cnumber
from matplotlib import pyplot as plt
from utility.model_utility import tray_translator
from matplotlib.backends.backend_pdf import PdfPages

'''-----------------------------------------------------------------------------
This can be used to duplicate the output to an external
Usage:
with HiddenLogs('filename','mode'):
    print("This will be displayed as usual + saved in log")
-----------------------------------------------------------------------------'''

class Logger(object):
    def __init__(self,_original_stdout, logfile):
        self.terminal = _original_stdout
        self.log = logfile

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.log.flush()

class HiddenLogs(object):
    def __init__(self,textlogdir,mode='a'):
        self.filename = textlogdir
        self.mode = mode

    def __enter__(self):
        self.logfile = open(self.filename, self.mode)
        self._original_stdout = sys.stdout
        sys.stdout = Logger(self._original_stdout, self.logfile)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.logfile.close()
        sys.stdout = self._original_stdout

class HiddenPrints(object):
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
1. Prepare the dictionary into the following shapyomo:
    dic = {'C10H20':[...],'C5H12':[...],...}
2. Use the function:
    reaction_data = trans_product_mole(dic)
3. Retuen dictionary:
    reaction_data['unscaled'] = {'gasoline':[...]}
    reaction_data['scaled'] = {'gasoline':[...]}
-----------------------------------------------------------------------------'''
def trans_product_mole(dic):
    # compute mole flow rate
    dataset = {}
    for c in m.PRODUCT_cnumber:
        for i in m.PRODUCT_cnumber[c]:
            tmp = np.array([np.array(dic[i]) for i in m.PRODUCT_cnumber[c]])
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

def check_product_spec(model):
    product_data = {}
    print('The following is considered as dry')
    for p in m.PRODUCT:
        product_data[p] = trans_product_mole({i:model.x_P_dry[i,p].value for i in m.COMP_ORG})
        print('{:<14.14}: {:<14.7}{:<10.10}: {:14.7}{:<10.10}: {:<14.7}{:<10.10}: {:>7}'\
        .format(p,'{:.2f}'.format(product_data[p]['unscaled'][p]),'Wet flow',\
        '{:.4f}'.format(model.P_total[p].value),'Dry flow','{:.4f}'.format(model.P_total_dry[p].value),\
        'N_tray','{:.1f}'.format(model.N_tray[p].value)))

'''-----------------------------------------------------------------------------
This can be used to pretty print a reactive distillation solution
Usage:
    1. beautify: with Reboiler
    2. beautify2: without Reboiler
This section updates regularly to reflect the latest need to print solution
-----------------------------------------------------------------------------'''

def beautify(pyomo,model):
    print('Here comes the result:')
    print('Total Conversion: {:.2%}'.format(cal_total_conversion(model)))
    print('-'*108)
    beautify_condenser(pyomo,model)
    print('')
    beautify_reactive(pyomo,model)
    print('')
    beautify_reboiler(pyomo,model)
    print('-'*108)

#------------------------------------------------------------------------------
def beautify2(pyomo,model):
    print('Here comes the result:')
    print('-'*108)
    beautify_condenser(pyomo,model)
    print('')
    beautify_reactive(pyomo,model)
    print('-'*108)

#------------------------------------------------------------------------------
def beautify_reactive(pyomo,model):
    placeholder = '{:<13.13}'+'{:<5.5}  '*2+' '*5+'{:<5.5}  '*4+' '*5+'{:<6.6}  '*4+' '*5+'{:<6.6}'
    convert_tmp = cal_conversion(model)
    print(placeholder.format('stages','T','Q','r_FT','Conv%','F','cat','V','Re','L','P','P_VLE'))
    for j in model.reactive:
        if model.reactive[j].cat.value >= 1e-6:
            stagename = 'React[{}]'.format(j)
        else:
            stagename = 'NON--[{}]'.format(j)
        temp_num =model.reactive[j].T.value - 273.15, model.reactive[j].Q_main.value,model.reactive[j].kinetics_block.r_FT_total.value,\
        convert_tmp[j-1],model.reactive[j].F.value,model.reactive[j].cat.value,model.reactive[j].V['out'].value,model.reactive[j].L['R'].value,\
        model.reactive[j].L['out'].value,model.reactive[j].L['P'].value,model.reactive[j].VLE_block.P_VLE.value
        temp_string = ['{:.5f}'.format(i) for i in temp_num]
        print(placeholder.format(stagename,*temp_string))
#------------------------------------------------------------------------------
def beautify_condenser(pyomo,model):
    placeholder = '{:<13.13}'+'{:<5.5}  '*2+' '*5+'{:<5.5}  '*4+' '*5+'{:<6.6}  '*4+' '*5+'{:<6.6}'
    print(placeholder.format('stages','T','Q','','','','','V','','L','P','W'))
    temp_num = model.condenser.T.value - 273.15, model.condenser.Q_main.value,model.condenser.V['P'].value,\
    model.condenser.L['out'].value,model.condenser.L['P'].value,model.condenser.W.value
    temp_string = ['{:.5f}'.format(i) for i in temp_num]
    temp_string = ['condenser',temp_string[0],temp_string[1],'','','','',temp_string[2],'',temp_string[3],\
                   temp_string[4],temp_string[5]]
    print(placeholder.format(*temp_string))
#------------------------------------------------------------------------------
def beautify_reboiler(pyomo,model):
    placeholder = '{:<13.13}'+'{:<5.5}  '*2+' '*5+'{:<5.5}  '*4+' '*5+'{:<6.6}  '*4+' '*5+'{:<6.6}'
    print(placeholder.format('stages','T','Q','','','','','V','','L','P','P_VLE'))
    temp_num = model.reboiler.T.value - 273.15, model.reboiler.Q_main.value,model.reboiler.V['out'].value,\
    model.reboiler.L['P'].value,model.reboiler.VLE_block.P_VLE.value
    temp_string = ['{:.5f}'.format(i) for i in temp_num]
    temp_string = ['reboiler',temp_string[0],temp_string[1],'','','','',temp_string[2],'','',\
                   temp_string[3],temp_string[4]]
    print(placeholder.format(*temp_string))


'''-----------------------------------------------------------------------------
This is used to group component (C10H22) data into carbon number data
1. trans_cnumber: mole, mole%
Usage:
1. Prepare the dictionary into the following shapyomo:
    dic = {'C10H20':[...],'C5H12':[...],...}
2. Use the function:
    phase_data = trans_cnumber(dic)
3. Return dictionary:
    phase_data= {'original index':[c1-c56...]}
-----------------------------------------------------------------------------'''
cnumber_range = range(1,57)

def trans_cnumber(dic):
    molefraction = {}
    for i in cnumber_range:
        molefraction[i] = []
    for i in m.COMP_ORG:
        molefraction[cal_cnumber(i)].append(np.array(dic[i]))
    for i in cnumber_range:
        molefraction[i] = np.sum(molefraction[i],0)
    tmp = []
    for i in cnumber_range:
        tmp.append(molefraction[i])
    return tmp

def cal_conversion(model):
    convertion_data = []
    for j in model.reactive:

        inletfeed = sum(model.reactive[j].V[s].value*(model.reactive[j].y_[s,'CO'].value + model.reactive[j].y_[s,'H2'].value) + \
                        model.reactive[j].L[s].value*(model.reactive[j].x_[s,'CO'].value + model.reactive[j].x_[s,'H2'].value) for s in model.reactive[j].inlet) + \
                    model.reactive[j].F.value
        outlet = sum(model.reactive[j].V[s].value*(model.reactive[j].y['CO'].value + model.reactive[j].y['H2'].value) + \
                    model.reactive[j].L[s].value*(model.reactive[j].x['CO'].value + model.reactive[j].x['H2'].value) for s in model.reactive[j].outlet)

        # to ensure other modules display some helpful information when the solution is incorrect
        try:
            tmp = round((inletfeed - outlet)/inletfeed,3)
        except:
            tmp = 0

        if tmp <= 0:
            convertion_data.append(0)
        else:
            convertion_data.append(tmp)
    return convertion_data

def cal_total_conversion(model):
    total_inlet = sum(model.reactive[j].F.value for j in model.reactive) + model.reboiler.F.value
    total_outlet = sum(tray_translator(model,j).V['P'].value*(tray_translator(model,j).y['CO'].value + tray_translator(model,j).y['H2'].value) + \
                    tray_translator(model,j).L['P'].value*(tray_translator(model,j).x['CO'].value + tray_translator(model,j).x['H2'].value) for j in model.TRAY_total)
    # to ensure other modules display properly when the solution is incorrect
    try:
        total_conversion = (total_inlet - total_outlet)/total_inlet
    except:
        total_conversion = 0

    return total_conversion

def plot_distribution(model,open_log_pdf = None,title = None):

    cnumber_range = range(1,57)

    g_data = []; d_data = []; l_out_data = {}; l_P_data = {}; b_data = []
    cd_x_data = []; rf_x_data = {}; rb_x_data = []

    g_data = trans_cnumber({i:model.condenser.V['P'].value*model.condenser.y[i].value for i in m.COMP_TOTAL})
    d_data = trans_cnumber({i:model.condenser.L['P'].value*model.condenser.x[i].value for i in m.COMP_TOTAL})
    for j in model.reactive:
        l_out_data[j] = trans_cnumber({i:model.reactive[j].L['out'].value*model.reactive[j].x[i].value for i in m.COMP_TOTAL})
        l_P_data[j] = trans_cnumber({i:model.reactive[j].L['P'].value*model.reactive[j].x[i].value for i in m.COMP_TOTAL})
    b_data = trans_cnumber({i:model.reboiler.L['P'].value*model.reboiler.x[i].value for i in m.COMP_TOTAL})

    cd_x_data = trans_cnumber({i:model.condenser.x[i].value for i in m.COMP_TOTAL})
    for j in model.reactive:
        rf_x_data[j] = trans_cnumber({i:model.reactive[j].x[i].value for i in m.COMP_TOTAL})
    rb_x_data = trans_cnumber({i:model.reboiler.x[i].value for i in m.COMP_TOTAL})

    # start the plot
    fig = plt.figure(figsize=(16,14))
    gs = plt.GridSpec(4, 3)
    ax1 = plt.subplot(gs[0,:-1])
    ax2 = plt.subplot(gs[1,:-1])
    ax3 = plt.subplot(gs[2,:-1])
    ax4 = plt.subplot(gs[:2,2])
    ax5 = plt.subplot(gs[3,:-1])
    ax6 = plt.subplot(gs[2:,2])
    '''
    ax1, product distribution
    '''

    ax1.plot(cnumber_range,g_data,'C6o-',markeredgecolor='w')
    ax1.plot(cnumber_range,d_data,'C2o-',markeredgecolor='w')

    if model.find_component('epi'):
        for j in model.reactive:
            if model.reactive[j].L['P'].value < 50*model.epi.value:
                ax1.plot(cnumber_range,l_P_data[j],'C7o-',alpha = 0.5,marker = r'${:}$'.format(j),markersize=10)
            else:
                ax1.plot(cnumber_range,l_P_data[j],'C0o-',marker = r'${:}$'.format(j),markersize=14)
    else:
        for j in model.reactive:
            ax1.plot(cnumber_range,l_P_data[j],'C0o-',marker = r'${:}$'.format(j),markersize=14)

    ax1.plot(cnumber_range,b_data,'C1o-',markeredgecolor='w')

    ax1.set_xlim(0, 30)
    ax1.set_yscale("log")
    ax1.set_ylim(1e-6, 1)
    # ax1.legend(['Vapor','Distillate',*['P{:}'.format(i) for i in model.reactive],'Bottom'],fontsize=14,loc=1)
    ax1.set_title('Product Distribution (Mole log)',fontsize=14)

    ax1.set_ylabel('Molar Flow (kmol/s)', color='K',fontsize=14)
    '''
    ax2, tray x composition
    '''
    ax2.plot(cnumber_range,cd_x_data,'C2o-',markeredgecolor='w')

    if model.find_component('epi'):
        for j in model.reactive:
            if model.reactive[j].L['P'].value < 50*model.epi.value:
                ax2.plot(cnumber_range,rf_x_data[j],'C7o-',alpha = 0.5,marker = r'${:}$'.format(j),markersize=10)
            else:
                ax2.plot(cnumber_range,rf_x_data[j],'C0o-',marker = r'${:}$'.format(j),markersize=14)
    else:
        for j in model.reactive:
            if model.reactive[j].L['P'].value < 1e-6:
                ax2.plot(cnumber_range,rf_x_data[j],'C7o-',alpha = 0.5,marker = r'${:}$'.format(j),markersize=10)
            else:
                ax2.plot(cnumber_range,rf_x_data[j],'C0o-',marker = r'${:}$'.format(j),markersize=14)

    ax2.plot(cnumber_range,rb_x_data,'C1o-',markeredgecolor='w')

    ax2.set_xlim(0, 30)
    ax2.set_ylim(0, 0.4)
    ax2.set_title('Liquid Composition (Mole)',fontsize=14)

    ax2.set_ylabel('Molar Fraction', color='K',fontsize=14)

    '''
    ax3, inter stage flows
    '''
    tray_num = len(model.TRAY)
    tray_pos = np.arange(tray_num)

    vapor_flow = [model.reactive[j].V['out'].value for j in model.reactive]
    liquid_flow = [model.reactive[j].L['out'].value for j in model.reactive]
    product_flow = [model.reactive[j].L['P'].value for j in model.reactive]
    ymax = max([model.reactive[j].L['P'].value + model.reactive[j].L['out'].value for j in model.reactive])
    ax3.set_ylim(0,ymax*1.2)

    line1,*trash = ax3.bar(tray_pos,liquid_flow, width = min(5/tray_num,0.1), color='C0')
    line2,*trash = ax3.bar(tray_pos,product_flow, bottom = liquid_flow, width = min(5/tray_num,0.1), color='C1')

    ax3_ = ax3.twinx()
    line3, = ax3_.plot(tray_pos,vapor_flow, 'C2o-',markersize=12,markeredgecolor='w')
    ymax = max(vapor_flow)
    ax3_.set_ylim(0,ymax*1.1)

    ax3.set_xticks(tray_pos)
    ax3.set_xticklabels(['T {:}'.format(j) for j in model.TRAY])
    ax3.set_title('Inter Stage Flow (kmol/s)',fontsize=14)

    ax3.set_ylabel('Liquid Flow (kmol/s)', color='K',fontsize=14)
    ax3_.set_ylabel('Vapor Flow (kmol/s)', color='K',fontsize=14)

    ax3.legend([line1,line2,line3],['Liquid','Product','Vapor'],loc=1,fontsize=12,fancybox=True,framealpha=0.2)

    '''
    ax4, reaction and conversion, also MPCC pf warning is embedded as well
    '''

    tray_pos_reboiler = np.arange(tray_num+1)
    tray_conv = cal_conversion(model)

    ax4_ = ax4.twinx()

    line1, = ax4_.plot([model.reactive[j].kinetics_block.r_FT_total.value for j in model.reactive] + [0],\
            tray_pos_reboiler,'C1o-',markersize=12,markeredgecolor='w')
    line2, *trash = ax4_.barh(tray_pos_reboiler,tray_conv + [0], height = min(20/tray_num,0.3), color='C0')
    if model.find_component('reactive[1].MPCC_P_pf').active and \
        model.find_component('reboiler.MPCC_P_pf').active:
        line3, = ax4.plot([1-model.reactive[j].MPCC_P_pf.pf.value for j in model.reactive] + [1-model.reboiler.MPCC_P_pf.pf.value],\
                tray_pos_reboiler,'C3o-',markersize=8,markeredgecolor='w')

    ax4_.set_xlim(-0.05,1.05)
    ax4_.set_xticklabels(['{:.0%}'.format(x) for x in ax4.get_xticks()])
    ax4_.set_yticks(tray_pos_reboiler)
    ax4_.set_yticklabels(['T {:}'.format(j) for j in model.TRAY] + ['Reboiler'])
    ax4_.invert_yaxis()
    ax4_.invert_xaxis()
    ax4.invert_yaxis()

    ax4.set_yticks([])

    ax4_.set_title('Total Conversion = {:.1%}'.format(cal_total_conversion(model)),fontsize=14)

    if model.find_component('reactive[1].MPCC_P_pf').active and \
        model.find_component('reboiler.MPCC_P_pf').active:
        ax4_.legend([line1,line2,line3],['r_FT','Conversion','Penalty'],loc=9,fontsize=12,fancybox=True,framealpha=0.2)
    else:
        ax4_.legend([line1,line2],['r_FT','Conversion'],loc=9,fontsize=12,fancybox=True,framealpha=0.2)

    '''
    ax5, temperature and Q
    '''
    temp_profile = [model.reactive[j].T.value-273.15 for j in model.reactive]
    Q_profile = [model.reactive[j].Q_main.value for j in model.reactive]

    line1,*trash = ax5.bar(tray_pos,temp_profile, width = min(5/tray_num,0.1), color='C0')

    ax5_ = ax5.twinx()
    line2, = ax5_.plot(tray_pos,Q_profile, 'C2o-',markersize=12,markeredgecolor='w')

    ymax = max(temp_profile+[250])
    ymin = min(temp_profile+[200])

    ax5.set_ylim([ymin*0.9,ymax*1.3])

    ymax = max(Q_profile+[5])
    ymin = min(Q_profile+[-5])
    ax5_.set_ylim([ymin*1.2,ymax*1.2])

    ax5.set_xticks(tray_pos)
    ax5.set_xticklabels(['T {:}'.format(j) for j in model.TRAY])
    ax5.set_title('Stage Temperature and Heat',fontsize=14)

    ax5.set_ylabel('Temperature (C)', color='K',fontsize=14)
    ax5_.set_ylabel('Q_in (MW)', color='K',fontsize=14)

    ax5_.legend([line1,line2],['Temperature','Q'],loc=1,fontsize=12,fancybox=True,framealpha=0.2)

    '''
    ax6, reaction and conversion, also MPCC pf warning is embedded as well
    '''

    ax6_ = ax6.twinx()
    ax6_y = ax6.twiny()

    line1, = ax6_y.plot([model.reactive[j].F.value for j in model.reactive] + [model.reboiler.F.value],\
            tray_pos_reboiler,'C1o-',markersize=12,markeredgecolor='w')
    line2, *trash = ax6_.barh(tray_pos_reboiler,[model.reactive[j].cat.value for j in model.reactive] + [0], height = min(20/tray_num,0.3), color='C0')

    ax6_.set_xticklabels(['{:.0f}'.format(x) for x in ax6.get_xticks()])
    ax6_.set_yticks(tray_pos_reboiler)
    ax6_.set_yticklabels(['T {:}'.format(j) for j in model.reactive] + ['Reboiler'])
    ax6_.invert_yaxis()
    ax6_.invert_xaxis()
    ax6_y.invert_yaxis()
    ax6_y.invert_xaxis()

    ax6.set_yticks([])

    ax6_.legend([line1,line2],['Feed','Catalyst'],loc=9,fontsize=12,fancybox=True,framealpha=0.2)

    '''
    Finalize
    '''

    ax1.grid()
    ax2.grid()
    ax3.grid(axis='y')
    ax4.grid(axis='x')
    ax5.grid(axis='y')
    ax6.grid(axis='x')

    if title:
        fig.suptitle('> '+title+' <', fontsize=18)

    plt.subplots_adjust(
        left  = 0.1,
        right = 0.9,
        bottom = 0.1,
        top = 0.9,
        wspace = 0.3,
        hspace = 0.3)

    if open_log_pdf:
        open_log_pdf.savefig()
    if "DISPLAY" in os.environ:
        plt.show()
    plt.close()

def plot_product_distribution(model,open_log_pdf = None):
    tray_num = len(model.TRAY_total)
    tray_pos = np.arange(tray_num)

    fig, ax = plt.subplots(figsize=(16,4))

    for i in m.PRODUCT:
        product_flow = [model.P_tray[j,i].value for j in model.TRAY_total]
        ax.bar(tray_pos,product_flow,alpha=0.7)

    ax.legend([i for i in m.PRODUCT],fontsize=12)
    ax.set_title('DDF Product Distribution',fontsize=14)
    ax.set_ylabel('Product Flow kmol/s',fontsize=14)

    ax.set_xticks(tray_pos)
    ax.set_xticklabels(['{:}'.format(j) for j in model.TRAY_total])

    ax.grid()

    if open_log_pdf:
        open_log_pdf.savefig()

    if "DISPLAY" in os.environ:
        plt.show()
    plt.close(fig)

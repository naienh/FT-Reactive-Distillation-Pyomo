'''-----------------------------------------------------------------------------
This is used to return what is inside the data module
-----------------------------------------------------------------------------'''
def print_pkg(m):
    key = [v for v in vars(m).keys() if not v.startswith('_')]
    for i in key: print(i)

'''-----------------------------------------------------------------------------
This is used to return the molecular weight of a certain components
-----------------------------------------------------------------------------'''
def cal_MW(name):
    if name == 'H2O': return 18
    if name == 'CO2': return 44
    if name == 'CO': return 28
    if name == 'H2': return 2
    start = name.index('C')+1
    end = name.index('H')
    MW = int(name[start:end])*12 + int(name[end+1:])
    return MW

'''-----------------------------------------------------------------------------
This is used to return the carbon number of a certain components
-----------------------------------------------------------------------------'''
def cal_cnumber(name):
    start = name.index('C')+1
    end = name.index('H')
    return int(name[start:end])

'''-----------------------------------------------------------------------------
This is used to return the olefin and paraffin ratio based on list of points
Usage:
1. Prepare the list of o/p ratio points in the following format:
op_ratio = [[1,0],[2,1.86],[3,4],[4,3.17],[5,1.7],[11,0.67],[21,0],[57,0]]
                                            c5-c11: 1.7, c11-c21: 0.67 ...
2. Return parafin ratios
-----------------------------------------------------------------------------'''
def cal_op(op_ratio):
    breakpoints = len(op_ratio);i = 1; paraffin_ratio = []; olefin_ratio = [];
    for k in range(breakpoints-1):
        while i < op_ratio[k+1][0]:
            paraffin_ratio.append(1/(op_ratio[k][1]+1)); i += 1;
    return paraffin_ratio

'''-----------------------------------------------------------------------------
This is used to return entire column begins with 'name' in a sheet
Usage:
1. workbook = xlrd.open_workbook("VLE_data.xlsx")
2. sheet_Henry = workbook.sheet_by_index(0)
3. henry.A = readcol(sheet_Henry,'A')
-----------------------------------------------------------------------------'''
def readcol(sheet,name):
    data_name = [sheet.cell_value(0,i) for i in range(sheet.ncols)]
    data = {}
    j = data_name.index(name)
    for i in range(1,sheet.nrows):
        data[sheet.cell_value(i,0)]=sheet.cell_value(i,j)
    return data

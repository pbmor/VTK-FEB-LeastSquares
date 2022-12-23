import numpy as np
import os
import csv
import glob
import vtk
import matplotlib.pyplot as plt
from math import cos, sin 
from vtk.util.numpy_support import vtk_to_numpy
from Read_XPLTfuncs import GetFEB, GetMeshInfo, GetData
from ResidualFunction_Time import RunLS, OrderList, GetPressure
from RemeshAll import Remesh_All
from itertools import zip_longest

CountStart = 0
Count = 81

ParamNames = []
Params = []
Params_Avg = []
for i in range(CountStart,Count):
    FileDir = 'ParSweep/Case_'+str(i+1) + '/'
    with open(FileDir+'Parameters.csv') as csv_file:
        XLData = csv.reader(csv_file, delimiter=',')
        next(XLData, None)
        for row in XLData:
            if row[2] != '': 
                if '_' not in row[2]:
                    try:
                        exec(row[2])
                    except NameError:
                        ParamNames.append(row[2])
                        exec(row[2]+'=[]')
                        exec('Params.append([])')
                        exec('Params_Avg.append([])')
                        
                    exec(row[2]+'.append('+row[3]+')')
                    print('Parameter:',row[2],'=', row[3])

Means = [] 
Stds = []    
for j, name in enumerate(ParamNames):
    exec(name+'_mean=np.mean('+name+')')
    exec(name+'_std=np.std('+name+')')
    exec('Params['+str(j)+']='+name)
    exec('Params_Avg['+str(j)+']=np.divide(Params['+str(j)+'],'+name+'_mean)')
    exec('Means.append('+name+'_mean)')
    exec('Stds.append('+name+'_std)')
    print('Mean of '+name+'=',Means[j])
    print('Standard Deviation of '+name+'=',Stds[j])
nCases = len(ParamNames)
CasesRange = list(range(0,nCases))       
plt.figure(1)
plt.boxplot(Params_Avg) 
plt.xticks(np.add(CasesRange,1), ParamNames, rotation='vertical')
plt.ylabel('Average parameter estimation')
plt.show()                     
B_Min = [ 0,    0,  0,  0]
B_Max = [100, 1000, 100, 50] 
for j, name in enumerate(ParamNames):
    plt.figure(j+2)
    exec('plt.plot('+name+')') 
    plt.xlabel(name)
    plt.ylim((B_Min[j],B_Max[j]*1.1))
    plt.ylabel(name+', Parameter Estimation')
    plt.show() 
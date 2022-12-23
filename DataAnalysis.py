import numpy as np
import os
import csv
import glob
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
from ResidualFunction_Time import OrderList
    
#Choose parameter estimation choices
PressureChoice = False              # Choose to vary pressure magnitude
ModelParChoice = True               # Choose to vary model parameters
RunLSChoice    = True               # Choose to run least Squares (or default/initial guess)
FibChoice      = False              # Choose to vary fiber direction, as an angle from the circumferential direction
ProfileChoice  = ['SetWindkessel']  # Choose profile shapes, options are: 'Triangle','Step','SmoothStep','Virtual', 'Fourier','Fitted'm Windkessel'
ResChoice      = ['CellPlane']      # Choose type of residual calculation method: 'P2P', 'CentreLine', 'CellPlane'
ModelChoice    = ['HGO']            # Choose model from 'MR','tiMR','Ogden' and 'Fung',  'HGO'
DataDirChoice  = 'All'              # Choose thewhich directories are to be included: 'Specfic', 'SpecificGroup', 'All_TAVs', 'All_BAVs','All'


# Code to run through all data sets
List_of_Subdirectories = sorted(glob.glob('./Strains/medial_meshes/*'))
ND = len(List_of_Subdirectories)
CommonOfDir = os.path.commonprefix(List_of_Subdirectories)

DataDirs = []
if DataDirChoice == 'Specific':
    DataDirs = ['tav02']
elif DataDirChoice == 'SpecificGroup':
    DataDirs = ['tav02','tav12','tav14','tav16']
elif DataDirChoice == 'All_TAVs':
    for d in List_of_Subdirectories:
        DataDir = d.replace(CommonOfDir,'')
        if DataDir[0] =='t':
            DataDirs.append(DataDir)
elif DataDirChoice == 'All_BAVs':
    for d in List_of_Subdirectories:
        DataDir = d.replace(CommonOfDir,'')
        if DataDir[0] =='b':
            DataDirs.append(DataDir)
elif DataDirChoice == 'All':
    for d in List_of_Subdirectories:
        DataDir = d.replace(CommonOfDir,'')
        DataDirs.append(DataDir)

ParamNames = []
Residuals = []
Residuals_Pts = []
nFs = []
for i, DataDir in enumerate(DataDirs):
    for PC in ProfileChoice:
        for MC in ModelChoice:
            for RC in ResChoice:
                Residuals.append([])
                print('Running Data Set: ',DataDir)
        
                flist = glob.glob('./Remeshed/'+DataDir+'/*')
                nF = len(flist)
                nFs.append(nF)    
                # Start counting parameters, based on model choices
                nC = 0
                Conditions = ''
                if PressureChoice:
                    Conditions += 'PMagT_'
                else:
                    Conditions += 'PMagF_'
                    
                if ModelParChoice:
                    Conditions += 'ModelT_'
                    if MC == 'MR':
                        nC +=4
                    elif MC == 'tiMR':
                        nC +=8
                    elif MC == 'Ogden':
                        nC += 14 
                    elif MC == 'Fung':
                        nC +=12
                    elif MC == 'HGO':
                        nC+=6
                else:
                    print('Model parameters not optimised')
                    Conditions += 'ModelF_'
                    if MC == 'MR':
                        nC +=4
                    elif MC == 'tiMR':
                        nC +=8
                    elif MC == 'Ogden':
                        nC += 14 
                    elif MC == 'Fung':
                        nC +=12
                    elif MC == 'HGO':
                        
                        nC+=6
                        
                if PC == 'Triangle':
                    nC +=1
                elif PC == 'Step':
                    nC +=2
                elif PC == 'SmoothStep':
                    nC +=4
                elif PC == 'Fourier':
                    nC +=4
                elif PC == 'Fitted':
                    nC +=nF-2
                elif ProfileChoice == 'Windkessel':
                    nC +=4
                '''
                elif PC == 'SetTriangle':
                    nP = 1
                    ParamNames  = ['Pressure Peak']
                    ParamValues = [0.5]
                elif PC == 'SetStep':
                    nP = 2
                    ParamNames  = ['Pressure_Increase','Pressure_Decrease']
                    ParamValues = [0.2,0.5]
                elif PC == 'SetSmoothStep':
                    nP = 4
                    ParamNames  = ['Pressure_Increase','Pressure_Decrease','Increase_Angle','Decrease_Angle']
                    ParamValues = [0.2,0.5,50,50]
                elif PC == 'SetFourier':
                    nP = 4
                    ParamNames  = ['Fourier_a1','Fourier_a2','Fourier_a3','Fourier_a4']
                    ParamValues = [3.0,1.0,0.0,0.05]
                elif PC == 'SetFitted':
                    nP = nF-2
                    ParamNames  = ['Pressure_P_'+ str(i) for i in range(nF)]
                    ParamValues = np.zeros(nF)
                elif PC == 'SetWindkessel':
                    nP = 4 
                    ParamNames  = ['qi_mag','Rp','Rd','Cp']
                    ParamValues = [11,1.4,14,0.004]
                
                # Define Timesteps
                TimeData = np.linspace(0,1,nF)
                # Get Pressure Profile
                PressureData = GetPressure(ParamValues,0,nF,CF,PC,TimeData)
                    
                # Get Interpolated Data
                TimeInterp = np.linspace(0,1,101)
                PressureInterp = np.interp(TimeInterp,TimeData,PressureData)
                PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
                '''    
                if FibChoice and MC == 'tiMR':
                    Conditions += 'FibT_'
                else:
                    Conditions += 'FibF_'
                 
                if RunLSChoice:
                    RunDir = 'RunLS'
                else:
                    RunDir = 'RunDefault'
        
                Conditions += RC
                
                FileDir = './NewFiles/'+ DataDir+'/'+MC+'/'+RunDir + '/' + PC+ '/' + Conditions+'/'
                
                flist = glob.glob(FileDir+'/*.vtk')
                nF = len(flist)
                
                #Get case info
                with open('echoframetime.csv') as csv_file:
                    XLData = csv.reader(csv_file, delimiter=',')
                    for row in XLData:
                        if DataDir == row[0]:
                            DataInfo = row
                            
                #Get reference frame number Choose frame before valve opens      
                refN = int(DataInfo[4])-1 
                
                #Define Frame Time Length
                FT = float(DataInfo[1])
                #Define Open Frame and Close Frame
                OF = DataInfo[4]
                CF = DataInfo[5] 
                
                #Get name of first and reference file
                common = os.path.commonprefix(flist)
                for Fname in list(flist):
                    Residuals_Pts.append([])
                    X = Fname.replace(common,'')
                    X = X.replace('.vtk','')
                    X = np.fromstring(X, dtype=int, sep=' ')
                    X=X[0]
                    if X==refN:
                        ref=Fname
                        
                FListOrdered, FId, refN = OrderList(flist, nF, ref)
            
                for X, Fname in enumerate(FListOrdered):                            
                    reader = vtk.vtkPolyDataReader()
                    reader.SetFileName(Fname)
                    reader.ReadAllScalarsOn()
                    reader.ReadAllVectorsOn()
                    reader.Update()
                    
                    polydata = reader.GetOutput()  
                    # Pts = vtk_to_numpy(polydata.GetPoints().GetData())
                    nPts = polydata.GetNumberOfPoints()
                    Res = vtk_to_numpy(polydata.GetPointData().GetArray('Residual'))
                    STJ = vtk_to_numpy(polydata.GetPointData().GetArray('STJ'))
                    VAJ = vtk_to_numpy(polydata.GetPointData().GetArray('VAJ'))
                    
                    BoundIds = np.add(STJ,VAJ)
                    
                    Res_reduced = []
                    for j, check in enumerate(BoundIds):
                        if round(check,1) !=1:
                            Res_reduced.append(Res[j])
                    
                    Residuals[i].append(sum(Res)/nPts)
                    if X!=refN:
                        if np. sum(Residuals_Pts)==0:
                            Residuals_Pts.append(Res_reduced)
                        else:
                            Residuals_Pts[i] += Res_reduced
                # Residuals_Pts[i] /= nFs[i]
                
                with open(FileDir+'Parameters.csv') as csv_file:
                    XLData = csv.reader(csv_file, delimiter=',')
                    next(XLData, None)
                    for row in XLData:
                        if row[2] != '':
                            try:
                                exec(row[2])
                            except NameError:
                                ParamNames.append(row[2])
                                exec(row[2]+'=[]')
                            exec(row[2]+'.append('+row[3]+')')
                            print('Parameter:',row[2],' = ', row[3])
                            
nCases = len(DataDirs)
Means = [] 
Stds = []    
CasesRange = list(range(0,nCases))
for i,name in enumerate(ParamNames):
    exec(name+'_mean = np.mean('+name+')')
    exec(name+'_std = np.std('+name+')')
    exec('Means.append('+name+'_mean)')
    exec('Stds.append('+name+'_std)')
    
    print('Mean of '+name+' is: ',Means[i])
    print('Standard Deviation of '+name+' is: ',Stds[i])
    exec('X = '+name)
    exec('X_mean = '+name+'_mean')
    exec('X_std = '+name+'_std')
    
    plt.figure(i)
    plt.plot(X,'ok')
    plt.plot(X_mean*np.ones(nCases),'k')
    plt.plot((X_mean+X_std)*np.ones(nCases),'--k',alpha = 0.4)
    plt.plot((X_mean-X_std)*np.ones(nCases),'--k',alpha = 0.4)
    plt.xticks(CasesRange, DataDirs, rotation='vertical')
    plt.ylabel(name)
    # plt.boxplot(exec(name))    
    
    
nBAV = 0
nTAV = 0      
if DataDirChoice == 'All':
    for i,DataDir in enumerate(DataDirs):
        if DataDir[0] == 'b':
            nBAV += 1
        if DataDir[0] == 't':
            nTAV += 1
    
    for i,name in enumerate(ParamNames):
        
        exec(name+'_BAVmean = np.mean('+name+'[0:nBAV])')
        exec(name+'_BAVstd = np.std('+name+'[0:nBAV])')
        exec(name+'_TAVmean = np.mean('+name+'[nBAV:])')
        exec(name+'_TAVstd = np.std('+name+'[nBAV:])')
        exec('X = '+name)
        exec('X_BAVmean = '+name+'_BAVmean')
        exec('X_BAVstd = '+name+'_BAVstd')
        exec('X_TAVmean = '+name+'_TAVmean')
        exec('X_TAVstd = '+name+'_TAVstd')
        plt.figure(i)
        plt.plot(CasesRange[0:nBAV],X_BAVmean*np.ones(nBAV),'b')
        plt.fill_between(CasesRange[0:nBAV],(X_BAVmean-X_BAVstd)*np.ones(nBAV),(X_BAVmean+X_BAVstd)*np.ones(nBAV),color='b',alpha=0.2)
        plt.plot(CasesRange[nBAV:],X_TAVmean*np.ones(nTAV),'r')
        plt.fill_between(CasesRange[nBAV:],(X_TAVmean-X_TAVstd)*np.ones(nTAV),(X_TAVmean+X_TAVstd)*np.ones(nTAV),color='r',alpha=0.2)
    

AvgRes = [sum(Residuals[i])/nFs[i] for i in range(nCases)]

plt.figure(i+1)
plt.plot(CasesRange[0:nBAV],AvgRes[0:nBAV],'ob')
plt.plot(CasesRange[nBAV:],AvgRes[nBAV:],'or')
plt.xticks(CasesRange, DataDirs, rotation='vertical')
plt.ylabel('Average Residual (per point per frame)')
plt.figure(i+2)
plt.plot(nFs,'ok')
plt.xticks(CasesRange, DataDirs, rotation='vertical')
plt.ylabel('Number of Frames')
    
Ind_Sort = (np.argsort(nFs)[i] for i in range(nCases))
AvgRes_Sort = []
DataDirs_Sort = []
for i in Ind_Sort:
    AvgRes_Sort.append(AvgRes[i])
    DataDirs_Sort.append(DataDirs[i])  
    
plt.figure(i+3)
plt.plot(AvgRes_Sort,'ok')
plt.xticks(CasesRange, DataDirs_Sort, rotation='vertical')
plt.ylabel('Average Residual (per point per frame)')

CasesRange = np.add(CasesRange,1)
plt.figure(i+4)
plt.boxplot(Residuals_Pts) 
plt.xticks(CasesRange, DataDirs, rotation='vertical')
plt.ylabel('Average Residual (per frame)')

Res_Pts_Avg = []
for j in range(nCases):
    Res_Pts_Avg.append(np.divide(Residuals_Pts[j],nFs[j]))
    
plt.figure(i+5)
plt.boxplot(Res_Pts_Avg) 
plt.xticks(CasesRange, DataDirs, rotation='vertical')
plt.ylabel('Average Residual (per frame)')
plt.show()                     
            
             
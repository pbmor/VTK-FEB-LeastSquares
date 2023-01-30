import numpy as np
import os
import csv
import glob
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
from ResidualFunction_Time import OrderList
from Read_XPLTfuncs import GetFEB
    
#Choose parameter estimation choices
PressureChoice = False              # Choose to vary pressure magnitude
RunLSChoice    = True               # Choose to run least Squares (or default/initial guess)
FibChoice      = False              # Choose to vary fiber direction, as an angle from the circumferential direction
ProfileChoice  = ['SetWindkessel']  # Choose profile shapes, options are: 'Triangle','Step','SmoothStep','Virtual', 'Fourier','Fitted'm Windkessel'
ResChoice      = ['P2P']      # Choose type of residual calculation method: 'P2P', 'CentreLine', 'CellPlane'
ModelChoice    = ['HGO']            # Choose model from 'MR','tiMR','Ogden' and 'Fung',  'HGO'
DataDirChoice  = 'Specific'              # Choose thewhich directories are to be included: 'Specfic', 'SpecificGroup', 'All_TAVs', 'All_BAVs','All'


BonusDetails = ''

# Code to run through all data sets
List_of_Subdirectories = sorted(glob.glob('./Strains/medial_meshes/*'))
ND = len(List_of_Subdirectories)
CommonOfDir = os.path.commonprefix(List_of_Subdirectories)

DataDirs = []
if DataDirChoice == 'Specific':
    DataDirs = ['bav01']
elif DataDirChoice == 'SpecificGroup':
    DataDirs = ['bav01','tav02']
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
elif DataDirChoice == 'AllExcept':
    Except = ['tav23','tav20','tav26']
    for d in List_of_Subdirectories:
        DataDir = d.replace(CommonOfDir,'')
        if DataDir not in Except:
            DataDirs.append(DataDir)

ParamNames = []
Residuals = []
Residuals_Pts = []
ModelTime = []
ModelPressure = []
PMag=[]
nFs = []
FTs = []
for i, DataDir in enumerate(DataDirs):
    for PC in ProfileChoice:
        for MC in ModelChoice:
            for RC in ResChoice:
                Residuals.append([])
                Residuals_Pts.append([])
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
                    
                if MC[0:3]!='Set':
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
                        nC+=4
                else:
                    print('Model parameters not optimised')
                    Conditions += 'ModelF_'
                    if MC[3:] == 'MR':
                        nC +=4
                    elif MC[3:] == 'tiMR':
                        nC +=8
                    elif MC[3:] == 'Ogden':
                        nC += 14 
                    elif MC[3:] == 'Fung':
                        nC +=12
                    elif MC[3:] == 'HGO':
                        nC+=4
                        
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
                FTs.append(FT)
                #Define Open Frame and Close Frame
                OF = DataInfo[4]
                CF = DataInfo[5] 
                
                #Get name of first and reference file
                common = os.path.commonprefix(flist)
                for Fname in list(flist):
                    # Residuals_Pts[i].append([])
                    X = Fname.replace(common,'')
                    X = X.replace('.vtk','')
                    X = np.fromstring(X, dtype=int, sep=' ')
                    X=X[0]
                    if X==refN:
                        ref=Fname
                
                                 
                reader = vtk.vtkPolyDataReader()
                reader.SetFileName(Fname)
                reader.ReadAllScalarsOn()
                reader.ReadAllVectorsOn()
                reader.Update()
                
                polydata = reader.GetOutput()  
                
                nPts = polydata.GetNumberOfPoints()
                nCls = polydata.GetNumberOfCells()
                
                infile = './FEB_Files/'+DataDir+'.log'
                
                important = []
                keep_phrases = ['N O R M A L']
                
                with open(infile) as f:
                    f = f.readlines()
                
                for line in f:
                    for phrase in keep_phrases:
                        if phrase in line:
                            important.append(line)
                            break
                if important == []:
                    print('Error in Febio Run...')
                    if MC in ['tiMR','HGO']:
                        nDoms = nCls
                    else:
                        nDoms=1
                    XPLTname = './FEB_Files/'+DataDir+'.xplt'
                    _,_,nStates, _ = GetFEB(XPLTname,nDoms,False)
                    if nF != nStates:
                        print('Note: there is a discrepancy in the number of frames.')
                        print('The number of simulated frames:', nStates)
                        print('The number of expected frames:', nF)
                    
                FListOrdered, FId, refN = OrderList(flist, nF, ref)
            
                for X, Fname in enumerate(FListOrdered):     
                    print(X)                       
                    reader = vtk.vtkPolyDataReader()
                    reader.SetFileName(Fname)
                    reader.ReadAllScalarsOn()
                    reader.ReadAllVectorsOn()
                    reader.Update()
                    
                    polydata = reader.GetOutput()  
                    Res  = vtk_to_numpy(polydata.GetPointData().GetArray('Residual'))
                    STJ  = vtk_to_numpy(polydata.GetPointData().GetArray('STJ'))
                    VAJ  = vtk_to_numpy(polydata.GetPointData().GetArray('VAJ'))
                    
                    BoundIds = np.add(STJ,VAJ)
                    
                    Res_reduced = []
                    for j, check in enumerate(BoundIds):
                        if round(check,1) !=1:
                            Res_reduced.append(Res[j])
                    
                    Residuals[i].append(sum(Res)/nPts)
                    if X!=refN:
                        if np. sum(Residuals_Pts[i])==0:
                            Residuals_Pts[i] = Res_reduced
                        else:
                            Residuals_Pts[i] += Res_reduced
                
                ModelTime.append([])
                ModelPressure.append([])
                with open(FileDir+'Parameters.csv') as csv_file:
                    XLData = csv.reader(csv_file, delimiter=',')
                    next(XLData, None)
                    for row in XLData:
                        if row[2] != '':
                            if 'Starting' not in row[2]:
                                try:
                                    exec(row[2])
                                except NameError:
                                    
                                    ParamNames.append(row[2])
                                    exec(row[2]+'=[]')
                                exec(row[2]+'.append('+row[3]+')')
                                print('Parameter:',row[2],' = ', row[3])
                        if row[4] == 'Pressure_Magnitude':
                            exec('PMag.append('+row[5]+')')
                        if row[6]!='':
                            exec('ModelTime[i].append('+row[6]+')')
                            exec('ModelPressure[i].append('+row[7]+')')
  
                    
CaseName = Conditions+BonusDetails

figname = 'Figures/'+CaseName +'/'  
os.system('mkdir '+figname)

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
    plt.xticks(CasesRange, DataDirs, rotation='vertical',fontsize=13)
    plt.yticks(fontsize=15)
    plt.ylabel(name,fontsize=20)
    # plt.boxplot(exec(name))  
    
    
    
PMag_mean = np.mean(PMag)
PMag_std  = np.std(PMag)
print('PMag_mean = ',PMag_mean)
print('PMag_std = ',PMag_std)    
    
nBAV = 0
nTAV = 0      
if DataDirChoice == 'AllExcept':
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
    

for i,name in enumerate(ParamNames):
    plt.figure(i,figsize=(8, 6), dpi=300)
    plt.savefig(name+'.png')
    plt.tight_layout()
    os.system('mv '+name+'.png '+figname+name+'.png')
    
AvgRes = [sum(Residuals[i])/nFs[i] for i in range(nCases)]

plt.figure(i+1)
plt.plot(CasesRange[0:nBAV],AvgRes[0:nBAV],'ob')
plt.plot(CasesRange[nBAV:],AvgRes[nBAV:],'or')
plt.xticks(CasesRange, DataDirs, rotation='vertical',fontsize=13)
plt.yticks(fontsize=15)
plt.ylabel('Average Residual (per point per frame)',fontsize=13)
plt.savefig('AverageResidual.png')
os.system('mv AverageResidual.png '+figname+'AverageResidual.png')

plt.figure(i+2)
plt.plot(nFs,'ok')
plt.xticks(CasesRange, DataDirs, rotation='vertical',fontsize=13)
plt.yticks(fontsize=15)
plt.ylabel('Number of Frames',fontsize=20)
    
# Ind_Sort = (np.argsort(nFs)[i] for i in range(nCases))
# AvgRes_Sort = []
# DataDirs_Sort = []
# for j in Ind_Sort:
#     AvgRes_Sort.append(AvgRes[j])
#     DataDirs_Sort.append(DataDirs[j])  
    
plt.figure(i+3)
plt.plot(AvgRes,'ok')
plt.xticks(CasesRange, DataDirs, rotation='vertical',fontsize=13)
plt.yticks(fontsize=15)
plt.ylabel('Average Residual',fontsize=20)


Res_Pts_Avg = []
for j in range(nCases):
    Residuals_Pts[i]
    Res_Pts_Avg.append(np.divide(Residuals_Pts[j],nFs[j]))

CasesRange = np.add(CasesRange,1)
plt.figure(i+4)
plt.boxplot(Residuals_Pts) 
plt.xticks(CasesRange, DataDirs, rotation='vertical',fontsize=13)
plt.yticks(fontsize=15)
plt.ylabel('Residual',fontsize=20)

    
plt.figure(i+5)
plt.boxplot(Res_Pts_Avg) 
plt.xticks(CasesRange, DataDirs, rotation='vertical',fontsize=13)
plt.yticks(fontsize=15)
plt.ylabel('Residual (per frame)',fontsize=20) 
plt.savefig('ResidualBoxplots.png')     
            
for j in range(nCases):
    if DataDirs[j][0] == 'b':
        plt.figure(i+6)
        plt.plot(np.multiply(ModelTime[j],np.multiply(nFs[j],FTs[j])),ModelPressure[j],label=DataDirs[j])
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('Time (ms)',fontsize=20)
        plt.ylabel('Pressure',fontsize=20)
        plt.legend()
        
        plt.figure(i+7)
        # plt.figure(figsize=(4,3))
        plt.plot(np.multiply(ModelTime[j],np.multiply(nFs[j],FTs[j])),np.multiply(ModelPressure[j],PMag[j]),label=DataDirs[j])
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('Time (ms)',fontsize=20)
        plt.ylabel('Pressure',fontsize=20)
        plt.legend()
        
    else:
        plt.figure(i+8)
        plt.plot(np.multiply(ModelTime[j],np.multiply(nFs[j],FTs[j])),ModelPressure[j],label=DataDirs[j])
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('Time (ms)',fontsize=20)
        plt.ylabel('Pressure',fontsize=20)
        plt.legend()
        
        plt.figure(i+9)
        # plt.figure(figsize=(4,3))
        plt.plot(np.multiply(ModelTime[j],np.multiply(nFs[j],FTs[j])),np.multiply(ModelPressure[j],PMag[j]),label=DataDirs[j])
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('Time (ms)',fontsize=20)
        plt.ylabel('Pressure',fontsize=20)
        plt.legend()

plt.figure(i+6)
plt.savefig('BAVPressue_excMagnitude.png')
os.system('mv BAVPressue_excMagnitude.png '+figname+'BAVPressue_excMagnitude.png')

plt.figure(i+7)
plt.savefig('BAVPressue_incMagnitude.png')
os.system('mv BAVPressue_incMagnitude.png '+figname+'BAVPressue_incMagnitude.png')

plt.figure(i+8)
plt.savefig('TAVPressue_excMagnitude.png')
os.system('mv TAVPressue_excMagnitude.png '+figname+'TAVPressue_excMagnitude.png')

plt.figure(i+9)
plt.savefig('TAVPressue_incMagnitude.png')
os.system('mv TAVPressue_incMagnitude.png '+figname+'TAVPressue_incMagnitude.png')

for j in range(nCases):
    if DataDirs[j][0] == 'b':
        plt.figure(i+10)
        plt.plot(np.linspace(0,1,nFs[j]+1),np.concatenate([Residuals[j],[Residuals[j][0]]],axis=0),label=DataDirs[j])
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('Time (ms)',fontsize=20)
        plt.ylabel('Residual',fontsize=20)
        plt.legend()
    if DataDirs[j][0] == 't':
        plt.figure(i+11)
        plt.plot(np.linspace(0,1,nFs[j]+1),np.concatenate([Residuals[j],[Residuals[j][0]]],axis=0),label=DataDirs[j])
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('Time (ms)',fontsize=20)
        plt.ylabel('Residual',fontsize=20)
        plt.legend()

plt.figure(i+10)
plt.savefig('BAVResiduals.png')
os.system('mv BAVResiduals.png '+figname+'BAVResiduals.png')

plt.figure(i+11)
plt.savefig('TAVResiduals.png')
os.system('mv TAVResiduals.png '+figname+'TAVResiduals.png')
plt.show()             
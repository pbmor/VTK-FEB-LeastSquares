import numpy as np
import os
import csv
import glob
from Read_XPLTfuncs import GetFEB, GetMeshInfo
import matplotlib.pyplot as plt
from ResidualFunction_Time import RunLS, SaveFiles, OrderList, GetPressure
from RemeshAll import Remesh_All
from itertools import zip_longest

#Choose parameter estimation choices
PressureChoice = True               # Choose to vary pressure magnitude
RunLSChoice    = True               # Choose to run least Squares (or default/initial guess)
FibChoice      = False              # Choose to vary fiber direction, as an angle from the circumferential direction
ProfileChoice  = ['Windkessel']  # Choose profile shapes, options are: 'Triangle','Step','SmoothStep','Virtual', 'Fourier','Fitted'm Windkessel'
ResChoice      = ['P2P']            # Choose type of residual calculation method: 'P2P', 'CentreLine', 'CellPlane'
ModelChoice    = ['SetHGO']            # Choose model from 'MR','tiMR','Ogden' and 'Fung',  'HGO'
DataDirChoice  = 'AllExcept'         # Choose which directories are to be included: 'Specfic', 'SpecificGroup', 'All_TAVs', 'All_BAVs','All','AllExcept'

# Code to run through all data sets
List_of_Subdirectories = sorted(glob.glob('./Strains/medial_meshes/*'))
ND = len(List_of_Subdirectories)
CommonOfDir = os.path.commonprefix(List_of_Subdirectories)

DataDirs = []
if DataDirChoice == 'Specific':
    DataDirs = ['tav02']
elif DataDirChoice == 'SpecificGroup':
    # DataDirs = ['bav01','tav02','bav02','tav04','bav07','tav06','bav08','tav11','bav12','tav12','bav13','tav14','bav23','tav16','tav20','tav23','tav26']
    DataDirs = ['bav01','tav04']#,'bav02','tav06']
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
    Except = ['tav02']
    for d in List_of_Subdirectories:
        DataDir = d.replace(CommonOfDir,'')
        if DataDir not in Except:
            DataDirs.append(DataDir)


for d in List_of_Subdirectories:
    DataDir = d.replace(CommonOfDir,'')
    FList = glob.glob('./Remeshed/'+DataDir+'/*')
    if FList ==[]:
        OriginalFiles = sorted(glob.glob('./Strains/medial_meshes/'+DataDir+'/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/*'))
        for OF in OriginalFiles:
            if OF[-4:]==".vtp"or OF[-4:]==".vtk":
                Remesh_All(DataDir,OF)
PressureMags   = []
PressureParams = []
CostsPressure  = []
CostsFinal     = []
for DId, DataDir in enumerate(DataDirs):
    print('Running Data Set: ',DataDir)
    
    FList = glob.glob('./Remeshed/'+DataDir+'/*')
    
    nF = len(FList)
    
    #Get case info
    with open('echoframetime.csv') as csv_file:
        XLData = csv.reader(csv_file, delimiter=',')
        for row in XLData:
            if DataDir == row[0]:
                DataInfo = row
                
    #Get reference frame number Choose frame before valve opens      
    refN = int(DataInfo[4])-1 
    CF = int(DataInfo[5])-refN

    #Get name of first and reference file
    common = os.path.commonprefix(FList)
    for Fname in list(FList):
        X = Fname.replace(common,'')
        X = X.replace('.vtk','')
        X = np.fromstring(X, dtype=int, sep=' ')
        X=X[0]
        if X==refN:
            ref=Fname
            
    #Order list to put reference frame first
    FListOrdered, FId, refN = OrderList(FList, nF, ref)
    FListOrdered.append(FListOrdered[0])
    FId = np.concatenate([FId,[FId[0]]],axis=0)
    nF+=1
        
    for PC in ProfileChoice:
        for MC in ModelChoice:
            for RC in ResChoice:
                # Initialise empty array for Setting pressure parameters
                SetP=[[],[]]
                # Initialise initial parameter arrays
                InitParams = [[],[],[],[]]
                
                # Choose initial estimate of pressure magnitude
                if PressureChoice:
                    InitParams[0] = [5]
                else:
                    SetP[0]=[5]
                 
                # Choose initial model parameters
                if ModelChoice[0:3] != 'Set':
                    if MC == 'MR':
                        InitParams[1] = [1,200,250]
                    if MC == 'tiMR':
                        InitParams[1] = [10,10,10,10,10,10,10,10]
                    elif MC == 'Ogden':
                        InitParams[1] = [1,10,10,10,10,10,10,10,10,10,10,10,10,10]
                    elif MC == 'Fung':
                        InitParams[1] = [1,1,2,1.5,0.5,0.1,1,1.2,2,1.5,0.5,0.1]
                    elif MC == 'HGO':     
                        InitParams[1] = [ 7.64,  996.6, 524.6, 49.98]
                else:
                    if MC == 'SetMR':
                        ModelPar = [1,10,10,10]
                    if MC == 'SettiMR':
                        ModelPar = [10,10,10,10,10,10,10,10]
                    elif MC == 'SetOgden':
                        ModelPar = [1,10,10,10,10,10,10,10,10,10,10,10,10,10]
                    elif MC == 'SetFung':
                        ModelPar = [1,1,2,1.5,0.5,0.1,1,1.2,2,1.5,0.5,0.1]
                    elif MC == 'SetHGO':     
                        ModelPar =  [ 7.64,  996.6, 524.6, 49.98]
                  
                # Define initial fiber angle
                # Note if the chosen model is not 'tiMR' this parameter array will be made to be empty
                if FibChoice and ModelChoice =='tiMR':
                    InitParams[2] = [0]
            
                #Choose initial paramters for the pressure profile
                if PC == 'Triangle':
                    InitParams[3] = [0.5]
                elif PC == 'Step':
                    InitParams[3] = [0.2, 0.5]
                elif PC == 'SmoothStep':
                    InitParams[3] = [0.2, 0.5,50,50]
                elif PC == 'Virtual':
                    InitParams[3] = []
                elif PC == 'Fourier':
                    InitParams[3] = [3.0, 1.0, 0.0, 0.05]
                elif PC == 'Fitted':
                    InitParams[3] = np.zeros(nF-2)
                elif PC == 'Windkessel':
                    InitParams[3] = [11, 1.4,14, 0.004]
                elif PC == 'SetTriangle':
                    SetP[1] = [0.5]
                elif PC == 'SetStep':
                    SetP[1] = [0.2, 0.5]
                elif PC == 'SetSmoothStep':
                    SetP[1] = [0.2, 0.5,50,50]
                elif PC == 'SetVirtual':
                    SetP[1] = []
                elif PC == 'SetFourier':
                    SetP[1] = [3.0, 1.0, 0.0, 0.05]
                elif PC == 'SetFitted':
                    SetP[1] = np.zeros(nF-2)
                elif PC == 'SetWindkessel':
                    SetP[1] = [11, 1.4,14, 0.004]
                    
                C = np.concatenate(InitParams, axis=0 )
                
                # If C is empty and RunLS is True, change to False
                if C.size == 0 and RunLSChoice:
                    print('Note: no parameters are being fitted, thus RunLSChoice is updated to be False')
                    RunLSChoice = False
                    
                #Run least squares script
                Out, Pts, Disp_Wall, Norm, Circ, Long, CellIds, nCls, STJ_Id, VAJ_Id, Cost = RunLS(DataDir,FListOrdered,FId,ref,CF,C,PressureChoice,PC,RunLSChoice,RC,MC,FibChoice,SetP)
                PressureMags.append(Out[0])
                PressureParams.append(Out[1:])
                CostsPressure.append(Cost)
                # Save new files
                nStates,Residual = SaveFiles(DataDir,FList,ref,Pts,Disp_Wall,Norm,Circ,Long,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,PC,RunLSChoice,RC,MC,FibChoice,Out,SetP)
    
                
                # Initialise Parameter Array to be saved as csv
                Params = [['Condition'],['Choice'],['Model Parameter Name'],['Parameter Value'],['Pressure Parameter Name'],['Parameter Value'],['Model Time'],['Model Pressure'],['Interpolated Time'],['Interpolated Pressure']]
                
                if nF != nStates:
                    print('The number of frames in original dataset: ', nF)
                    print('The number of frames in simulated dataset: ', nStates)
                else:
                    print('The number of frames: ', nF)
                    
                # Start counting parameters, based on model choices
                nC = 0
                Conditions = ''
                if PressureChoice:
                    Params[0].append('Pressure_Magnitude')
                    Params[1].append('True')
                    nP = 1
                    ParamNames = ['Pressure_Magnitude']
                    ParamValues = [Out[0]]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    for i in range(nP):
                        Params[4].append(ParamNames[i]+'_Starting')
                        Params[5].append(InitParams[0][i])
                    nC +=nP
                    Conditions += 'PMagT_'
                else:
                    Params[0].append('Pressure_Magnitude')
                    Params[1].append('False')
                    Params[4].append('Pressure_Magnitude')
                    Params[5].append(SetP[0][0])
                    print('Pressure Magnitude is Default')
                    Conditions += 'PMagF_'
                    
                if MC[0:3] != 'Set':
                    Conditions += 'ModelT_'
                    Params[0].append('Model_Parameters')
                    Params[1].append('True')
                    if MC == 'MR':
                        Params[0].append('Model')
                        Params[1].append('MR')
                        nP = 3
                        ParamNames  = ['density','C1','C2']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[2].append(ParamNames[i])
                            Params[3].append(ParamValues[i])
                        for i in range(nP):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(InitParams[1][i])
                        nC += nP
                    elif MC == 'tiMR':
                        Params[0].append('Model')
                        Params[1].append('tiMR')
                        nP = 8
                        ParamNames  = ['density','C1','C2','C3','C4','C5','k','lam_max']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[2].append(ParamNames[i])
                            Params[3].append(ParamValues[i])
                        for i in range(nP):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(InitParams[1][i])
                        nC +=nP
                    elif MC == 'Ogden':
                        Params[0].append('Model')
                        Params[1].append('Ogden')
                        nP = 14
                        ParamNames  = ['density','k','c1','c2','c3','c4','c5','c6','m1','m2','m3','m4','m5','m6']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[2].append(ParamNames[i])
                            Params[3].append(ParamValues[i])
                        for i in range(nP):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(InitParams[1][i])
                        nC += nP
                    elif MC == 'Fung':
                        Params[0].append('Model')
                        Params[1].append('Fung')
                        nP = 12
                        ParamNames  = ['density','E1','E2','E3','G12','G23','G31','v12','v23','v31','c','k']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[2].append(ParamNames[i])
                            Params[3].append(ParamValues[i])
                        for i in range(nP):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(InitParams[1][i])
                        nC +=nP
                    elif MC == 'HGO':
                        Params[0].append('Model')
                        Params[1].append('HGO')
                        nP = 4
                        ParamNames  = ['c','k1','k2','gamma']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[2].append(ParamNames[i])
                            Params[3].append(ParamValues[i])
                        for i in range(nP):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(InitParams[1][i])
                        nC+=nP
                else:
                    print('Model parameters not optimised')
                    Conditions += 'ModelF_'
                    Params[0].append('Model_Parameters')
                    Params[1].append('True')
                    if MC == 'SetMR':
                        Params[0].append('Model')
                        Params[1].append('MR')
                        for i in range(3):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(ModelPar[i])
                    elif MC == 'SettiMR':
                        Params[0].append('Model')
                        Params[1].append('tiMR')
                        for i in range(8):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(ModelPar[i])
                    elif MC == 'SetOgden':
                        Params[0].append('Model')
                        Params[1].append('Ogden')
                        for i in range(14):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(ModelPar[i])
                    elif MC == 'SetFung':
                        Params[0].append('Model')
                        Params[1].append('Fung')
                        for i in range(12):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(ModelPar[i])
                    elif MC == 'SetHGO':
                        Params[0].append('Model')
                        Params[1].append('HGO')
                        for i in range(4):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(ModelPar[i])
                
                if RunLSChoice:
                    Params[0].append('Run_Least_Squares')
                    Params[1].append('True')
                else:
                    Params[0].append('Run_Least_Squares')
                    Params[1].append('False')
                
                if PC == 'Triangle':
                    nP = 1
                    ParamNames  = ['Pressure Peak']
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nP
                elif PC == 'Step':
                    nP = 2
                    ParamNames  = ['Pressure_Increase','Pressure_Decrease']
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nP
                elif PC == 'SmoothStep':
                    nP = 4
                    ParamNames  = ['Pressure_Increase','Pressure_Decrease','Increase_Angle','Decrease_Angle']
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nP
                elif PC == 'Fourier':
                    nP = 4
                    ParamNames  = ['Fourier_a1','Fourier_a2','Fourier_a3','Fourier_a4']
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nP
                elif PC == 'Fitted':
                    nP = nF-2
                    ParamNames  = ['Pressure_P_'+ str(i) for i in range(nF)]
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nF-2
                elif PC == 'Windkessel':
                    nP = 4
                    ParamNames  = ['qi_mag','Rp','Rd','Cp']
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nP
                elif PC == 'SetTriangle':
                    nP = 1
                    ParamNames  = ['Pressure Peak']
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                elif PC == 'SetStep':
                    nP = 2
                    ParamNames  = ['Pressure_Increase','Pressure_Decrease']
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                elif PC == 'SetSmoothStep':
                    nP = 4
                    ParamNames  = ['Pressure_Increase','Pressure_Decrease','Increase_Angle','Decrease_Angle']
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                elif PC == 'SetFourier':
                    nP = 4
                    ParamNames  = ['Fourier_a1','Fourier_a2','Fourier_a3','Fourier_a4']
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                elif PC == 'SetFitted':
                    nP = nF-2
                    ParamNames  = ['Pressure_P_'+ str(i) for i in range(nF)]
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                elif PC == 'SetWindkessel':
                    nP = 4 
                    ParamNames  = ['qi_mag','Rp','Rd','Cp']
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                
                # Define Timesteps
                TimeData = np.linspace(0,1,nF)
                # Get Pressure Profile
                PressureData = GetPressure(ParamValues,0,nF,CF,PC,TimeData,SetP)
                for i in range(len(TimeData)):
                    Params[6].append(TimeData[i])
                    Params[7].append(PressureData[i])
                    
                # Get Interpolated Data
                TimeInterp = np.linspace(0,1,101)
                PressureInterp = np.interp(TimeInterp,TimeData,PressureData)
                PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
                
                for i in range(len(TimeInterp)):
                    Params[8].append(TimeInterp[i])
                    Params[9].append(PressureInterp[i])
                
                # Save Pressure Profile Choice
                Params[0].append('Pressure_Profile')
                Params[1].append(PC)
                    
                if FibChoice and MC == 'tiMR':
                    Params[0].append('Include_Fiber')
                    Params[1].append('True')
                    nP = 1
                    ParamNames  = ['Fiber_Angle']
                    ParamValues = Out[nC:nC+nP+1]
                    
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC += nP
                    Conditions += 'FibT_'
                else:
                    Params[0].append('Include_Fiber')
                    Params[1].append('False')
                    Conditions += 'FibF_'
                 
                if RunLSChoice:
                    RunDir = 'RunLS'
                else:
                    RunDir = 'RunDefault'
    
                Conditions += RC
                 
                # Save Pressure Profile Choice
                Params[0].append('Residual_Choice')
                Params[1].append(RC)
                
                Params[0].append('Final_Cost')
                Params[1].append(Cost)

                with open('./NewFiles/'+ DataDir+'/'+MC+'/'+RunDir + '/' + PC+ '/' + Conditions+'/Parameters.csv','w') as result_file:
                    wr = csv.writer(result_file, dialect='excel')
                    for values in zip_longest(*Params):
                        wr.writerow(values)
                
                XPLTfilename = './FEB_Files/' + DataDir + '.xplt'

                if ModelChoice in ['MR','SetMR']:
                    nDoms = 1
                else:
                    nDoms = nCls                  
                feb,file_size,nStates, mesh = GetFEB(XPLTfilename,nDoms,False)
                nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
                
                Residual_Mag = np.zeros(nF)
                Residual_VectorMag = np.zeros((nF,nNodes))
                for i in range(nF):
                    for j in range(nNodes):
                        Residual_VectorMag[i,j] = np.sqrt(np.sum(np.power(Residual[i,j],2)))
                    Residual_Mag[i] = np.sum(Residual_VectorMag[i])
                
                plt.figure(1)
                plt.plot(np.linspace(0,1,nF),Residual_Mag/nNodes,label = ProfileChoice)
                plt.xlabel('Time')
                plt.ylabel('Relative Residual (per node)')
                    
    #Choose parameter estimation choices
    PressureChoice = False               # Choose to vary pressure magnitude
    RunLSChoice    = True              # Choose to run least Squares (or default/initial guess)
    FibChoice      = False              # Choose to vary fiber direction, as an angle from the circumferential direction
    ProfileChoice  = ['SetWindkessel']  # Choose profile shapes, options are: 'Triangle','Step','SmoothStep','Virtual', 'Fourier','Fitted'm Windkessel'
    ResChoice      = ['P2P']            # Choose type of residual calculation method: 'P2P', 'CentreLine', 'CellPlane'
    ModelChoice    = ['HGO']            # Choose model from 'MR','tiMR','Ogden' and 'Fung',  'HGO'

                
# for DId, DataDir in enumerate(DataDirs):
    print('Running Data Set: ',DataDir)

    FList = glob.glob('./Remeshed/'+DataDir+'/*')
    
    nF = len(FList)
    
    #Get case info
    with open('echoframetime.csv') as csv_file:
        XLData = csv.reader(csv_file, delimiter=',')
        for row in XLData:
            if DataDir == row[0]:
                DataInfo = row
                
    #Get reference frame number Choose frame before valve opens      
    refN = int(DataInfo[4])-1 
    CF = int(DataInfo[5])-refN

    #Get name of first and reference file
    common = os.path.commonprefix(FList)
    for Fname in list(FList):
        X = Fname.replace(common,'')
        X = X.replace('.vtk','')
        X = np.fromstring(X, dtype=int, sep=' ')
        X=X[0]
        if X==refN:
            ref=Fname
            
    #Order list to put reference frame first
    FListOrdered, FId, refN = OrderList(FList, nF, ref)
    FListOrdered.append(FListOrdered[0])
    FId = np.concatenate([FId,[FId[0]]],axis=0)
    nF+=1
        
    for PC in ProfileChoice:
        for MC in ModelChoice:
            for RC in ResChoice:
                SetP = [[0],[1]]
                # Initialise initial parameter arrays
                InitParams = [[],[],[],[]]
                
                # Choose initial estimate of pressure magnitude
                if PressureChoice:
                    InitParams[0] = PressureMags[DId]
                else:
                    SetP[0] = PressureMags[DId]
                 
                # Choose initial model parameters
                if ModelChoice[0:3] != 'Set':
                    if MC == 'MR':
                        InitParams[1] = [1,200,250]
                    if MC == 'tiMR':
                        InitParams[1] = [10,10,10,10,10,10,10,10]
                    elif MC == 'Ogden':
                        InitParams[1] = [1,10,10,10,10,10,10,10,10,10,10,10,10,10]
                    elif MC == 'Fung':
                        InitParams[1] = [1,1,2,1.5,0.5,0.1,1,1.2,2,1.5,0.5,0.1]
                    elif MC == 'HGO':     
                        InitParams[1] = [ 7.64,  996.6, 524.6, 49.98]
                else:   
                    if ModelChoice == 'SetMR':
                        ModelPar = [1,10,10,10]
                    if ModelChoice == 'SettiMR':
                        ModelPar = [10,10,10,10,10,10,10,10]
                    elif ModelChoice == 'SetOgden':
                        ModelPar = [1,10,10,10,10,10,10,10,10,10,10,10,10,10]
                    elif ModelChoice == 'SetFung':
                        ModelPar = [1,1,2,1.5,0.5,0.1,1,1.2,2,1.5,0.5,0.1]
                    elif ModelChoice == 'SetHGO':     
                        ModelPar =  [ 7.64,  996.6, 524.6, 49.98]
                  
                # Define initial fiber angle
                # Note if the chosen model is not 'tiMR' this parameter array will be made to be empty
                if FibChoice and ModelChoice =='tiMR':
                    InitParams[2] = [0]
                    
                #Choose initial paramters for the pressure profile
                if PC == 'Triangle':
                    InitParams[3] = [0.5]
                elif PC == 'Step':
                    InitParams[3] = [0.2, 0.5]
                elif PC == 'SmoothStep':
                    InitParams[3] = [0.2, 0.5,50,50]
                elif PC == 'Virtual':
                    InitParams[3] = []
                elif PC == 'Fourier':
                    InitParams[3] = [3.0, 1.0, 0.0, 0.05]
                elif PC == 'Fitted':
                    InitParams[3] = np.zeros(nF-2)
                elif PC == 'Windkessel':
                    InitParams[3] = [11, 1.4,14, 0.004]
                elif PC == 'SetTriangle':
                    SetP[1] = PressureParams[DId]
                elif PC == 'SetStep':
                    SetP[1] = PressureParams[DId]
                elif PC == 'SetSmoothStep':
                    SetP[1] = PressureParams[DId]
                elif PC == 'SetVirtual':
                    SetP[1] = []
                elif PC == 'SetFourier':
                    SetP[1] = PressureParams[DId]
                elif PC == 'SetFitted':
                    SetP[1] = np.zeros(nF-2)
                elif PC == 'SetWindkessel':
                    SetP[1] = PressureParams[DId]
                    
                C = np.concatenate(InitParams, axis=0 )
                
                # If C is empty and RunLS is True, change to False
                if C.size == 0 and RunLSChoice:
                    print('Note: no parameters are being fitted, thus RunLSChoice is updated to be False')
                    RunLSChoice = False
                    
                #Run least squares script
                Out, Pts, Disp_Wall, Norm, Circ, Long, CellIds, nCls, STJ_Id, VAJ_Id, Cost = RunLS(DataDir,FListOrdered,FId,ref,CF,C,PressureChoice,PC,RunLSChoice,RC,MC,FibChoice,SetP)
                
                CostsFinal.append(Cost)
                
                # Save new files
                nStates,Residual = SaveFiles(DataDir,FList,ref,Pts,Disp_Wall,Norm,Circ,Long,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,PC,RunLSChoice,RC,MC,FibChoice,Out,SetP)
    
                
                # Initialise Parameter Array to be saved as csv
                Params = [['Condition'],['Choice'],['Model Parameter Name'],['Parameter Value'],['Pressure Parameter Name'],['Parameter Value'],['Model Time'],['Model Pressure'],['Interpolated Time'],['Interpolated Pressure']]
                
                if nF != nStates:
                    print('The number of frames in original dataset: ', nF)
                    print('The number of frames in simulated dataset: ', nStates)
                else:
                    print('The number of frames: ', nF)
                    
                # Start counting parameters, based on model choices
                nC = 0
                Conditions = ''
                if PressureChoice:
                    Params[0].append('Pressure_Magnitude')
                    Params[1].append('True')
                    nP = 1
                    ParamNames = ['Pressure_Magnitude']
                    ParamValues = [Out[0]]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    for i in range(nP):
                        Params[4].append(ParamNames[i]+'_Starting')
                        Params[5].append(InitParams[0][i])
                    nC +=nP
                    Conditions += 'PMagT_'
                else:
                    Params[0].append('Pressure_Magnitude')
                    Params[1].append('False')
                    Params[4].append('Pressure_Magnitude')
                    Params[5].append(SetP[0][0])
                    print('Pressure Magnitude is Default')
                    Conditions += 'PMagF_'
                    
                if MC[0:3] != 'Set':
                    Conditions += 'ModelT_'
                    Params[0].append('Model_Parameters')
                    Params[1].append('True')
                    if MC == 'MR':
                        Params[0].append('Model')
                        Params[1].append('MR')
                        nP = 3
                        ParamNames  = ['density','C1','C2']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[2].append(ParamNames[i])
                            Params[3].append(ParamValues[i])
                        for i in range(nP):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(InitParams[1][i])
                        nC += nP
                    elif MC == 'tiMR':
                        Params[0].append('Model')
                        Params[1].append('tiMR')
                        nP = 8
                        ParamNames  = ['density','C1','C2','C3','C4','C5','k','lam_max']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[2].append(ParamNames[i])
                            Params[3].append(ParamValues[i])
                        for i in range(nP):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(InitParams[1][i])
                        nC +=nP
                    elif MC == 'Ogden':
                        Params[0].append('Model')
                        Params[1].append('Ogden')
                        nP = 14
                        ParamNames  = ['density','k','c1','c2','c3','c4','c5','c6','m1','m2','m3','m4','m5','m6']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[2].append(ParamNames[i])
                            Params[3].append(ParamValues[i])
                        for i in range(nP):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(InitParams[1][i])
                        nC += nP
                    elif MC == 'Fung':
                        Params[0].append('Model')
                        Params[1].append('Fung')
                        nP = 12
                        ParamNames  = ['density','E1','E2','E3','G12','G23','G31','v12','v23','v31','c','k']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[2].append(ParamNames[i])
                            Params[3].append(ParamValues[i])
                        for i in range(nP):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(InitParams[1][i])
                        nC +=nP
                    elif MC == 'HGO':
                        Params[0].append('Model')
                        Params[1].append('HGO')
                        nP = 4
                        ParamNames  = ['c','k1','k2','gamma']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[2].append(ParamNames[i])
                            Params[3].append(ParamValues[i])
                        for i in range(nP):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(InitParams[1][i])
                        nC+=nP
                else:
                    print('Model parameters not optimised')
                    Conditions += 'ModelF_'
                    Params[0].append('Model_Parameters')
                    Params[1].append('True')
                    if MC == 'SetMR':
                        Params[0].append('Model')
                        Params[1].append('MR')
                        for i in range(3):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(ModelPar[i])
                    elif MC == 'SettiMR':
                        Params[0].append('Model')
                        Params[1].append('tiMR')
                        for i in range(8):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(ModelPar[i])
                    elif MC == 'SetOgden':
                        Params[0].append('Model')
                        Params[1].append('Ogden')
                        for i in range(14):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(ModelPar[i])
                    elif MC == 'SetFung':
                        Params[0].append('Model')
                        Params[1].append('Fung')
                        for i in range(12):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(ModelPar[i])
                    elif MC == 'SetHGO':
                        Params[0].append('Model')
                        Params[1].append('HGO')
                        for i in range(4):
                            Params[2].append(ParamNames[i]+'_Starting')
                            Params[3].append(ModelPar[i])
                
                if RunLSChoice:
                    Params[0].append('Run_Least_Squares')
                    Params[1].append('True')
                else:
                    Params[0].append('Run_Least_Squares')
                    Params[1].append('False')
                
                if PC == 'Triangle':
                    nP = 1
                    ParamNames  = ['Pressure Peak']
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nP
                elif PC == 'Step':
                    nP = 2
                    ParamNames  = ['Pressure_Increase','Pressure_Decrease']
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nP
                elif PC == 'SmoothStep':
                    nP = 4
                    ParamNames  = ['Pressure_Increase','Pressure_Decrease','Increase_Angle','Decrease_Angle']
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nP
                elif PC == 'Fourier':
                    nP = 4
                    ParamNames  = ['Fourier_a1','Fourier_a2','Fourier_a3','Fourier_a4']
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nP
                elif PC == 'Fitted':
                    nP = nF-2
                    ParamNames  = ['Pressure_P_'+ str(i) for i in range(nF)]
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nF-2
                elif PC == 'Windkessel':
                    nP = 4
                    ParamNames  = ['qi_mag','Rp','Rd','Cp']
                    ParamValues = Out[nC:nC+nP+1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC +=nP
                elif PC == 'SetTriangle':
                    nP = 1
                    ParamNames  = ['Pressure Peak']
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                elif PC == 'SetStep':
                    nP = 2
                    ParamNames  = ['Pressure_Increase','Pressure_Decrease']
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                elif PC == 'SetSmoothStep':
                    nP = 4
                    ParamNames  = ['Pressure_Increase','Pressure_Decrease','Increase_Angle','Decrease_Angle']
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                elif PC == 'SetFourier':
                    nP = 4
                    ParamNames  = ['Fourier_a1','Fourier_a2','Fourier_a3','Fourier_a4']
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                elif PC == 'SetFitted':
                    nP = nF-2
                    ParamNames  = ['Pressure_P_'+ str(i) for i in range(nF)]
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                elif PC == 'SetWindkessel':
                    nP = 4 
                    ParamNames  = ['qi_mag','Rp','Rd','Cp']
                    ParamValues = SetP[1]
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                
                # Define Timesteps
                TimeData = np.linspace(0,1,nF)
                # Get Pressure Profile
                PressureData = GetPressure(ParamValues,0,nF,CF,PC,TimeData,SetP)
                for i in range(len(TimeData)):
                    Params[6].append(TimeData[i])
                    Params[7].append(PressureData[i])
                    
                # Get Interpolated Data
                TimeInterp = np.linspace(0,1,101)
                PressureInterp = np.interp(TimeInterp,TimeData,PressureData)
                PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
                
                for i in range(len(TimeInterp)):
                    Params[8].append(TimeInterp[i])
                    Params[9].append(PressureInterp[i])
                
                # Save Pressure Profile Choice
                Params[0].append('Pressure_Profile')
                Params[1].append(PC)
                    
                if FibChoice and MC == 'tiMR':
                    Params[0].append('Include_Fiber')
                    Params[1].append('True')
                    nP = 1
                    ParamNames  = ['Fiber_Angle']
                    ParamValues = Out[nC:nC+nP+1]
                    
                    for i in range(nP):
                        print(ParamNames[i],': ',ParamValues[i])
                        Params[4].append(ParamNames[i])
                        Params[5].append(ParamValues[i])
                    nC += nP
                    Conditions += 'FibT_'
                else:
                    Params[0].append('Include_Fiber')
                    Params[1].append('False')
                    Conditions += 'FibF_'
                 
                if RunLSChoice:
                    RunDir = 'RunLS'
                else:
                    RunDir = 'RunDefault'
    
                Conditions += RC
                 
                # Save Pressure Profile Choice
                Params[0].append('Residual_Choice')
                Params[1].append(RC)
                
                # Save Final Cost
                Params[0].append('Cost')
                Params[1].append(Cost)
                
                with open('./NewFiles/'+ DataDir+'/'+MC+'/'+RunDir + '/' + PC+ '/' + Conditions+'/Parameters.csv','w') as result_file:
                    wr = csv.writer(result_file, dialect='excel')
                    for values in zip_longest(*Params):
                        wr.writerow(values)
                
                XPLTfilename = './FEB_Files/' + DataDir + '.xplt'

                if ModelChoice in ['MR','SetMR']:
                    nDoms = 1
                else:
                    nDoms = nCls
                
                feb,file_size,nStates, mesh = GetFEB(XPLTfilename,nDoms,False)
                nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
                
                Residual_Mag = np.zeros(nF)
                Residual_VectorMag = np.zeros((nF,nNodes))
                for i in range(nF):
                    for j in range(nNodes):
                        Residual_VectorMag[i,j] = np.sqrt(np.sum(np.power(Residual[i,j],2)))
                    Residual_Mag[i] = np.sum(Residual_VectorMag[i])
                
                # Plot Final Residuals
                plt.figure(1)
                plt.plot(np.linspace(0,1,nF),Residual_Mag/nNodes,label = ProfileChoice)
                plt.xlabel('Time')
                plt.ylabel('Relative Residual (per node)')
                    
plt.show()
        
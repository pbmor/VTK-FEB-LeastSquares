import numpy as np
import os
import csv
import glob
# from xml.etree import ElementTree as et
from Read_XPLTfuncs import GetFEB, GetMeshInfo
import matplotlib.pyplot as plt
from ResidualFunction_Time import RunLS, SaveFiles, OrderList, GetPressure
from RemeshAll import Remesh_All
from itertools import zip_longest


    
#Choose parameter estimation choices
PressureChoice = False              # Choose to vary pressure magnitude
ModelParChoice = True               # Choose to vary model parameters
RunLSChoice    = True               # Choose to run least Squares (or default/initial guess)
FibChoice      = False              # Choose to vary fiber direction, as an angle from the circumferential direction
ProfileChoice  = ['SetWindkessel']  # Choose profile shapes, options are: 'Triangle','Step','SmoothStep','Virtual', 'Fourier','Fitted'm Windkessel'
ResChoice      = ['CellPlane']      # Choose type of residual calculation method: 'P2P', 'CentreLine', 'CellPlane'
ModelChoice    = ['MR']             # Choose model from 'MR','tiMR','Ogden' and 'Fung',  'HGO'
DataDirChoice  = 'Specific'         # Choose thewhich directories are to be included: 'Specfic', 'SpecificGroup', 'All_TAVs', 'All_BAVs','All'



# Code to run through all data sets
List_of_Subdirectories = sorted(glob.glob('./Strains/medial_meshes/*'))
ND = len(List_of_Subdirectories)
CommonOfDir = os.path.commonprefix(List_of_Subdirectories)

DataDirs = []
if DataDirChoice == 'Specific':
    DataDirs = ['tav02']
elif DataDirChoice == 'SpecificGroup':
    DataDirs = ['tav02','tav04', 'tav06']
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


for d in List_of_Subdirectories:
    DataDir = d.replace(CommonOfDir,'')
    flist = glob.glob('./Remeshed/'+DataDir+'/*')
    if flist ==[]:
        OriginalFiles = sorted(glob.glob('./Strains/medial_meshes/'+DataDir+'/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/*'))
        for OF in OriginalFiles:
            if OF[-4:]==".vtp"or OF[-4:]==".vtk":
                Remesh_All(DataDir,OF)
                
Params = [['Condition'],['Choice'],['Model Parameter Name'],['Parameter Value'],['Pressure Parameter Name'],['Parameter Value'],['Model Time'],['Model Pressure'],['Interpolated Time'],['Interpolated Pressure']]
for DataDir in DataDirs:
        print('Running Data Set: ',DataDir)
        
        # Initialise Parameter Array to be saved as csv
        Params = [['Condition'],['Choice'],['Model Parameter Name'],['Parameter Value'],['Pressure Parameter Name'],['Parameter Value'],['Model Time'],['Model Pressure'],['Interpolated Time'],['Interpolated Pressure']]

        flist = glob.glob('./Remeshed/'+DataDir+'/*')
        
        nF = len(flist)
        
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
        common = os.path.commonprefix(flist)
        for Fname in list(flist):
            X = Fname.replace(common,'')
            X = X.replace('.vtk','')
            X = np.fromstring(X, dtype=int, sep=' ')
            X=X[0]
            if X==refN:
                ref=Fname
                
        #Order list to put reference frame first
        FListOrdered, FId, refN = OrderList(flist, nF, ref)
            
        for PC in ProfileChoice:
            for MC in ModelChoice:
                for RC in ResChoice:
                    
                    # Initialise initial parameter arrays
                    InitParams = [[],[],[],[]]
                    
                    # Choose initial estimate of pressure magnitude
                    if PressureChoice:
                        InitParams[0] = [1]
                     
                    # Choose initial model parameters
                    if ModelParChoice:
                        if MC == 'MR':
                            InitParams[1] = [1,10,10,10]
                        if MC == 'tiMR':
                            InitParams[0] = [10,10,10,10,10,10,10,10]
                        elif MC == 'Ogden':
                            InitParams[0] = [1,10,10,10,10,10,10,10,10,10,10,10,10,10]
                        elif MC == 'Fung':
                            InitParams[0] = [1,1,2,1.5,0.5,0.1,1,1.2,2,1.5,0.5,0.1]
                        elif MC == 'HGO':     
                            InitParams[0] = [ 5,  415,  5, 28, 0.3, 10]
                      
                    # Define initial fiber angle
                    # Note if the chosen model is not 'tiMR' this parameter array will be made to be empty
                    if FibChoice and ModelChoice =='tiMR':
                        InitParams[2] = [0]
                
                    #Choose initial paramters for the pressure profile
                    if PC == 'Triangle':
                        InitParams[3] = [0.5]
                    if PC == 'Step':
                        InitParams[3] = [0.2, 0.5]
                    if PC == 'SmoothStep':
                        InitParams[3] = [0.2, 0.5,50,50]
                    if PC == 'Virtual':
                        InitParams[3] = []
                    if PC == 'Fourier':
                        InitParams[3] = [3.0, 1.0, 0.0, 0.05]
                    if PC == 'Fitted':
                        InitParams[3] = np.zeros(nF-2)
                    if PC == 'Windkessel':
                        InitParams[3] = [11,1.4,14,0.004]
                    
                    C = np.concatenate(InitParams, axis=0 )
                    
                    # If C is empty and RunLS is True, change to False
                    if C.size == 0 and RunLSChoice:
                        print('Note: no parameters are being fitted, thus RunLSChoice is updated to be False')
                        RunLSChoice == False
                    if ModelParChoice and not RunLSChoice:
                        print('Note: Parameters are not being fitted, thus ModelParChoice is updated to be False')
                        ModelParChoice = False
                        
                    #Run least squares script
                    Out, Pts, Disp_Wall, Norm, Circ, Long, CellIds, nCls, STJ_Id, VAJ_Id = RunLS(DataDir,FListOrdered,FId,ref,CF,C,PressureChoice,ModelParChoice,PC,RunLSChoice,RC,MC,FibChoice)
                    
                    # Save new files
                    nStates,Residual = SaveFiles(DataDir,ref,Pts,Disp_Wall,Norm,Circ,Long,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,ModelParChoice,PC,RunLSChoice,RC,MC,FibChoice,Out)
        
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
                        ParamNames = 'Pressure_Magnitude'
                        ParamValues = Out[0]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                        nC +=1
                        Conditions += 'PMagT_'
                    else:
                        Params[0].append('Pressure_Magnitude')
                        Params[1].append('False')
                        Params[4].append('Pressure_Magnitude')
                        Params[5].append('1.0')
                        print('Pressure Magnitude is Default')
                        Conditions += 'PMagF_'
                        
                    if ModelParChoice:
                        Conditions += 'ModelT_'
                        Params[0].append('Model_Parameters')
                        Params[1].append('True')
                        if MC == 'MR':
                            Params[0].append('Model')
                            Params[1].append('MR')
                            nP = 4
                            ParamNames  = ['density','C1','C2','k']
                            ParamValues = Out[nC:nC+nP+1]
                            for i in range(nP):
                                print(ParamNames[i],': ',ParamValues[i])
                                Params[2].append(ParamNames[i])
                                Params[3].append(ParamValues[i])
                            nC +=4
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
                            nC +=8
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
                            nC += 14 
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
                            nC +=12
                        elif MC == 'HGO':
                            Params[0].append('Model')
                            Params[1].append('HGO')
                            nP = 6
                            ParamNames  = ['c','k1','k2','gamma','kappa','k']
                            ParamValues = Out[nC:nC+nP+1]
                            for i in range(nP):
                                print(ParamNames[i],': ',ParamValues[i])
                                Params[2].append(ParamNames[i])
                                Params[3].append(ParamValues[i])
                            nC+=6
                    else:
                        print('Model parameters not optimised')
                        Conditions += 'ModelF_'
                        Params[0].append('Model_Parameters')
                        Params[1].append('True')
                        if MC == 'MR':
                            Params[0].append('Model')
                            Params[1].append('MR')
                            nC +=4
                        elif MC == 'tiMR':
                            Params[0].append('Model')
                            Params[1].append('tiMR')
                            nC +=8
                        elif MC == 'Ogden':
                            Params[0].append('Model')
                            Params[1].append('Ogden')
                            nC += 14 
                        elif MC == 'Fung':
                            Params[0].append('Model')
                            Params[1].append('Fung')
                            nC +=12
                        elif MC == 'HGO':
                            Params[0].append('Model')
                            Params[1].append('HGO')
                            nC+=6
                    
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
                        nC +=1
                    elif PC == 'Step':
                        nP = 2
                        ParamNames  = ['Pressure_Increase','Pressure_Decrease']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                        nC +=2
                    elif PC == 'SmoothStep':
                        nP = 4
                        ParamNames  = ['Pressure_Increase','Pressure_Decrease','Increase_Angle','Decrease_Angle']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                        nC +=4
                    elif PC == 'Fourier':
                        nP = 4
                        ParamNames  = ['Fourier_a1','Fourier_a2','Fourier_a3','Fourier_a4']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                        nC +=4
                    elif PC == 'Fitted':
                        nP = nF-2
                        ParamNames  = ['Pressure_P_'+ str(i) for i in range(nF)]
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                        nC +=nF-2
                    elif ProfileChoice == 'Windkessel':
                        nP = 4
                        ParamNames  = ['qi_mag','Rp','Rd','Cp']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                        nC +=4
                    elif PC == 'SetTriangle':
                        nP = 1
                        ParamNames  = ['Pressure Peak']
                        ParamValues = [0.5]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                    elif PC == 'SetStep':
                        nP = 2
                        ParamNames  = ['Pressure_Increase','Pressure_Decrease']
                        ParamValues = [0.2,0.5]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                    elif PC == 'SetSmoothStep':
                        nP = 4
                        ParamNames  = ['Pressure_Increase','Pressure_Decrease','Increase_Angle','Decrease_Angle']
                        ParamValues = [0.2,0.5,50,50]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                    elif PC == 'SetFourier':
                        nP = 4
                        ParamNames  = ['Fourier_a1','Fourier_a2','Fourier_a3','Fourier_a4']
                        ParamValues = [3.0,1.0,0.0,0.05]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                    elif PC == 'SetFitted':
                        nP = nF-2
                        ParamNames  = ['Pressure_P_'+ str(i) for i in range(nF)]
                        ParamValues = np.zeros(nF)
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                    elif PC == 'SetWindkessel':
                        nP = 4 
                        ParamNames  = ['qi_mag','Rp','Rd','Cp']
                        ParamValues = [11,1.4,14,0.004]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[4].append(ParamNames[i])
                            Params[5].append(ParamValues[i])
                    
                    # Define Timesteps
                    TimeData = np.linspace(0,1,nF)
                    # Get Pressure Profile
                    PressureData = GetPressure(ParamValues,0,nF,CF,PC,TimeData)
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
                        nC += 1
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
    
                    with open('./NewFiles/'+ DataDir+'/'+MC+'/'+RunDir + '/' + PC+ '/' + Conditions+'/Parameters.csv','w') as result_file:
                        wr = csv.writer(result_file, dialect='excel')
                        for values in zip_longest(*Params):
                            wr.writerow(values)
                
                    
                    XPLTfilename = './FEB_Files/' + DataDir + '.xplt'

                    if ModelChoice == 'tiMR':
                        nDoms = nCls
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
plt.show()
        
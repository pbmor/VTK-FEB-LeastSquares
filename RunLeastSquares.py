import numpy as np
import os
# import vtk
import csv
import glob
from xml.etree import ElementTree as et
from Read_XPLTfuncs import GetFEB, GetMeshInfo
# from RemeshAll import Remesh_All
# from vtk2feb_Disp import VTK2Feb_Func
import matplotlib.pyplot as plt
from ResidualFunction_Time import RunLS, SaveFiles, OrderList
from RemeshAll import Remesh_All


# Code to run through all data sets
List_of_Subdirectories = sorted(glob.glob('../RunVTK/Strains/medial_meshes/*'))

# List_of_Subdirectories = '../RunVTK/Strains/medial_meshes/bav02'
ND = len(List_of_Subdirectories)

CommonOfDir = os.path.commonprefix(List_of_Subdirectories)
for d in List_of_Subdirectories:
    DataDir = d.replace(CommonOfDir,'')
    flist = glob.glob('./Remeshed/'+DataDir+'/*')
    if flist ==[]:
        OriginalFiles = sorted(glob.glob('../RunVTK/Strains/medial_meshes/'+DataDir+'/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/*'))
        for OF in OriginalFiles:
            print(OF)
            if OF[-4:]==".vtp"or OF[-4:]==".vtk":
                Remesh_All(DataDir,OF)
                
Params = []
for d in List_of_Subdirectories:
    DataDir = d.replace(CommonOfDir,'')
    # if DataDir[0] == 't':
    if DataDir in ['tav02']:
    # if DataDir in ['tav06','tav11','tav12','tav16']:
    # if DataDir in ['tav02','tav14','tav20','tav23','tav26']:
        print('Running Data Set: ',DataDir)
        
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
        
        #Choose if data needs remeshed
        PressureChoice = False           # Choose to vary pressure magnitude
        ModelParChoice = True            # Choose to vary modelparameters
        RunLSChoice    = True            # Choose to run least Squares (or default/initial guess)
        ProfileChoice  = ['Bio']         #['Triangle','Step','SmoothStep','Bio','Fourier','Fitted'] # Choose profile shapes, options are: 'Triangle','Step','SmoothStep','Bio', 'Fourier','Fitted'
        ResChoice      = ['CellPlane']   # Choose type of residual calculation method: 'P2P', 'CentreLine', 'CellPlane'
        ModelChoice    = ['MR']          #Choose model from 'MR','tiMR','Ogden' and 'Fung'
        
        
        for PC in ProfileChoice:
            for MC in ModelChoice:
                for RC in ResChoice:
                    #Run least squares script
                    Out, Pts, Disp_Wall, Norm, CellIds, nCls, STJ_Id, VAJ_Id, FId, nF = RunLS(DataDir,d,FListOrdered,FId,ref,CF,PressureChoice,ModelParChoice,PC,RunLSChoice,RC,MC)
                    
                    RemeshedFiles = glob.glob('./Remeshed/'+DataDir+'/*')
                    
                    for Fname in list(RemeshedFiles):
                        X = Fname.replace(common,'')
                        X = X.replace('.vtk','')
                        X = np.fromstring(X, dtype=int, sep=' ')
                        X=X[0]
                        if X==refN:
                            RemeshedFile=Fname
                    
                    Params.append([DataDir,Out])
                    nStates,Residual = SaveFiles(DataDir,RemeshedFile,Pts,Disp_Wall,Norm,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,ModelParChoice,PC,RunLSChoice,RC,MC,Out)
                    if nF != nStates:
                        print('The number of frames in original dataset: ', nF)
                        print('The number of frames in simulated dataset: ', nStates)
                    else:
                        print('The number of frames: ', nF)
                        
                    nNodes = len(Residual[0])
                    print('The number of nodes: ',nNodes)
                    nC = 0
                    if PressureChoice:
                        Pressure_Mag = Out[0]
                        print('Pressure Magnitude: ',Pressure_Mag)
                        nC +=1
                    else:
                        print('Pressure Magnitude: ',-0.005332)
                        
                    if ModelParChoice:
                        if MC == 'MR':
                            Model_C1 = Out[0+nC]
                            Model_C2 = Out[1+nC]
                            Model_k  = Out[2+nC]
                            print('Model C1: ',Model_C1)
                            print('Model C2: ',Model_C2)
                            print('Model k:  ',Model_k)
                        elif MC == 'Ogden':
                            Model_k = Out[0+nC]
                            Model_c1   = Out[1+nC]
                            Model_c2   = Out[2+nC]
                            Model_c3   = Out[3+nC]
                            Model_c4   = Out[4+nC]
                            Model_c5   = Out[5+nC]
                            Model_c6   = Out[6+nC]
                            Model_m1   = Out[7+nC]
                            Model_m2   = Out[8+nC]
                            Model_m3   = Out[9+nC]
                            Model_m4   = Out[10+nC]
                            Model_m5   = Out[11+nC]
                            Model_m6   = Out[12+nC]
                            print('Model k: ',Model_k)
                            print('Model C1: ',Model_c1)
                            print('Model C2: ',Model_c2)
                            print('Model C3: ',Model_c3)
                            print('Model C4: ',Model_c4)
                            print('Model C5: ',Model_c5)
                            print('Model C6: ',Model_c6)
                            print('Model m1: ',Model_m1)
                            print('Model m2: ',Model_m2)
                            print('Model m3: ',Model_m3)
                            print('Model m4: ',Model_m4)
                            print('Model m5: ',Model_m5)
                            print('Model m6: ',Model_m6)
                        nC+=3
                    else:
                        print('Model C1: ',1)
                        print('Model C2: ',0)
                        print('Model k: ',10)
               
                    if MC == 'MR':
                        line = ''
                    elif MC == 'tiMR':
                        line = ''
                    elif MC == 'Ogden':
                        line = '--'    
                    elif MC == 'Fung':
                        line = ':'     
                    
                    if PC == 'Triangle':
                        col = 'k'
                    elif PC == 'Step':
                        col = 'r'
                    elif PC == 'SmoothStep':
                        col = 'b'
                    elif PC == 'Bio':
                        col = 'g'
                    elif PC == 'Fourier':
                        col = 'm'
                    elif PC == 'Windkess':
                        col = 'c'
                    elif PC == 'SetWindkess':
                        col = 'c'
                        
                    style = col +line
                    
                    XPLTfilename = './FEB_Files/' + DataDir + '.xplt'
                    feb,file_size,nStates, mesh = GetFEB(XPLTfilename)
                    nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
            
                    TimeData = (np.genfromtxt('TimeProfile.csv', delimiter=','))
                    PressureData = np.genfromtxt('PressureProfile.csv', delimiter=',')
                    TimeInterp = np.linspace(0,1,1001)
                    PressureInterp = np.interp(TimeInterp,TimeData,PressureData)
                    PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
                    
                    Time = np.zeros(nF)
                    Pressure = np.zeros(nF)
                    tree = et.parse(DataDir + '.feb')
                    Pressure_Mag = float(tree.find('Loads').find('surface_load').find('pressure').text)
                    for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
                        for i in range(len(Pt.text)):
                            if Pt.text[i] == ',':
                                Time[idx] = float(Pt.text[0:i])
                                Pressure[idx] = float(Pt.text[i+1:])
                         
                    Residual_Mag = np.zeros(nF)
                    Residual_VectorMag = np.zeros((nF,nNodes))
                    for i in range(nF):
                        for j in range(nNodes):
                            Residual_VectorMag[i,j] = np.sqrt(np.sum(np.power(Residual[i,j],2)))
                        Residual_Mag[i] = np.sum(Residual_VectorMag[i])
                    
                        
                    plt.figure(1)
                    plt.plot(np.linspace(0,1,nF),Residual_Mag/nNodes,label = PC)
                    plt.xlabel('Time')
                    plt.ylabel('Relative Residual')
                    
                    if PC == 'Fitted':
                        plt.figure(2)
                        plt.plot(Time,np.multiply(Pressure,-Pressure_Mag),label = PC)
                        plt.xlabel('Time')
                        plt.ylabel('Pressure')
                
plt.show()
        
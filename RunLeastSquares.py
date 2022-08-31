import numpy as np
import os
import csv
import glob
from xml.etree import ElementTree as et
from Read_XPLTfuncs import GetFEB, GetMeshInfo
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
            if OF[-4:]==".vtp"or OF[-4:]==".vtk":
                Remesh_All(DataDir,OF)
                
Params = [['Name'],['Values']]
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
        ModelParChoice = True            # Choose to vary model parameters
        RunLSChoice    = True            # Choose to run least Squares (or default/initial guess)
        FibChoice      = True            # Choose to vary fiber direction, as an angle from the circumferential direction
        ProfileChoice  = ['Windkess']    # Choose profile shapes, options are: 'Triangle','Step','SmoothStep','Bio', 'Fourier','Fitted'
        ResChoice      = ['P2P']         # Choose type of residual calculation method: 'P2P', 'CentreLine', 'CellPlane'
        ModelChoice    = ['MR','tiMR']        # Choose model from 'MR','tiMR','Ogden' and 'Fung',  'HGO'
        
        
        for PC in ProfileChoice:
            for MC in ModelChoice:
                for RC in ResChoice:
                    #Run least squares script
                    Out, Pts, Disp_Wall, Norm, Circ, Long, CellIds, nCls, STJ_Id, VAJ_Id, FId, nF = RunLS(DataDir,d,FListOrdered,FId,ref,CF,PressureChoice,ModelParChoice,PC,RunLSChoice,RC,MC,FibChoice)
                    
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
                        nP = 1
                        ParamNames = 'Pressure Magnitude'
                        ParamValues = Out[0]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[0].append(ParamNames[i])
                            Params[1].append(ParamValues[i])
                        nC +=1
                        Conditions += 'PMagT_'
                    else:
                        print('Pressure Magnitude is Default')
                        Conditions += 'PMagF_'
                        
                    if ModelParChoice:
                        Conditions += 'ModelT_'
                        if MC == 'MR':
                            nP = 4
                            ParamNames  = ['density','C1','C2','k']
                            ParamValues = Out[nC:nC+nP+1]
                            for i in range(nP):
                                print(ParamNames[i],': ',ParamValues[i])
                                Params[0].append(ParamNames[i])
                                Params[1].append(ParamValues[i])
                            nC +=4
                            Conditions += 'MR_'
                        elif MC == 'tiMR':
                            nP = 8
                            ParamNames  = ['density','C1','C2','C3','C4','C5','k','lam_max']
                            ParamValues = Out[nC:nC+nP+1]
                            for i in range(nP):
                                print(ParamNames[i],': ',ParamValues[i])
                                Params[0].append(ParamNames[i])
                                Params[1].append(ParamValues[i])
                            nC +=8
                            Conditions += 'tiMR_'
                        elif MC == 'Ogden':
                            nP = 14
                            ParamNames  = ['density','k','c1','c2','c3','c4','c5','c6','m1','m2','m3','m4','m5','m6']
                            ParamValues = Out[nC:nC+nP+1]
                            for i in range(nP):
                                print(ParamNames[i],': ',ParamValues[i])
                                Params[0].append(ParamNames[i])
                                Params[1].append(ParamValues[i])
                            nC += 14 
                            Conditions += 'Ogden_'
                        elif MC == 'Fung':
                            nP = 12
                            ParamNames  = ['density','E1','E2','E3','G12','G23','G31','v12','v23','v31','c','k']
                            ParamValues = Out[nC:nC+nP+1]
                            for i in range(nP):
                                print(ParamNames[i],': ',ParamValues[i])
                                Params[0].append(ParamNames[i])
                                Params[1].append(ParamValues[i])
                            nC +=12
                            Conditions += 'Fung_'
                        elif MC == 'HGO':
                            nP = 6
                            ParamNames  = ['c','k1','k2','gamma','kappa','k']
                            ParamValues = Out[nC:nC+nP+1]
                            for i in range(nP):
                                print(ParamNames[i],': ',ParamValues[i])
                                Params[0].append(ParamNames[i])
                                Params[1].append(ParamValues[i])
                            nC+=6
                            Conditions += 'HGO_'
                    else:
                        print('Model parameters not optimised')
                        Conditions += 'ModelF_'
                    
                    if RunLSChoice:
                        Conditions += 'RunLST_'
                    else:
                        Conditions += 'RunLSF_'
                    
                    if PC == 'Triangle':
                        nP = 1
                        ParamNames  = ['Pressure Peak']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[0].append(ParamNames[i])
                            Params[1].append(ParamValues[i])
                        nC +=1
                    elif PC == 'Step':
                        nP = 2
                        ParamNames  = ['Pressure Increase','Pressure Decrease']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[0].append(ParamNames[i])
                            Params[1].append(ParamValues[i])
                        nC +=2
                    elif PC == 'SmoothStep':
                        nP = 4
                        ParamNames  = ['Pressure Increase','Pressure Decrease','Increase Angle','Decrease Angle']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[0].append(ParamNames[i])
                            Params[1].append(ParamValues[i])
                        nC +=4
                    elif PC == 'Fourier':
                        nP = 4
                        ParamNames  = ['Fourier a1','Fourier a2','Fourier a3','Fourier a4']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[0].append(ParamNames[i])
                            Params[1].append(ParamValues[i])
                        nC +=4
                    elif PC == 'Fitted':
                        nP = nF
                        ParamNames  = ['Pressure P_'+ str(i) for i in range(nF)]
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[0].append(ParamNames[i])
                            Params[1].append(ParamValues[i])
                        nC +=nF
                    elif PC == 'Windkess':
                        nP = 4
                        ParamNames  = ['qi_mag','Rp','Rd','Cp']
                        ParamValues = Out[nC:nC+nP+1]
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[0].append(ParamNames[i])
                            Params[1].append(ParamValues[i])
                        nC +=4
                        
                        
                    if FibChoice and MC == 'tiMR':
                        nP = 1
                        ParamNames  = ['Fiber Angle']
                        ParamValues = Out[nC:nC+nP+1]
                        
                        for i in range(nP):
                            print(ParamNames[i],': ',ParamValues[i])
                            Params[0].append(ParamNames[i+nC])
                            Params[1].append(ParamValues[i+nC])
                        nC += 1
                        Conditions += 'FibT_'
                    else:
                        Conditions += 'FibF_'
                        
                    Conditions += PC +'_'+ RC
                        
                    with open('./NewFiles/'+ DataDir+'/'+Conditions+'/Parameters.csv','w') as result_file:
                        wr = csv.writer(result_file, dialect='excel')
                        for i in range(len(Params[0])):
                            wr.writerow([Params[0][i],Params[1][i]])
                        
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
                        
                    # style = col +line
                    style = col
                    
                    XPLTfilename = './FEB_Files/' + DataDir + '.xplt'
                    feb,file_size,nStates, mesh = GetFEB(XPLTfilename,False)
                    nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
            
                    TimeData = (np.genfromtxt('TimeProfile.csv', delimiter=','))
                    PressureData = np.genfromtxt('PressureProfile.csv', delimiter=',')
                    TimeInterp = np.linspace(0,1,1001)
                    PressureInterp = np.interp(TimeInterp,TimeData,PressureData)
                    PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
                    
                    Time = np.zeros(nF)
                    Pressure = np.zeros(nF)
                    tree = et.parse('./FEB_Files/' + DataDir + '.feb')
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
        
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

def SaveFiles(DataDir,RemeshedFile,Pts,Disp_Wall,Norm,Circ_Cls,Long_Cls,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,ProfileChoice,RunLSChoice,ResChoice,ModelChoice,FibChoice,Out,Count):
    '''
    Read xplt file and save to remeshed vtk file and save as a new vtk file.
    
    Inputs:
    DataDir - Name of directory, e.g. 'tav02'
    RemeshedFile - Name of remeshed files
    Pts - Array of points
    Disp_Wall - Array of wall displacements for each point of mesh
    Norm - Array of Normal of points
    CellIds - Array of point Ids in each cell
    nCls - Number of cells
    STJ_Id - Point Ids for the STJ
    VAJ_Id - Point Ids for the VAJ
    FId - Array of the ordered frame Ids, starting the reference frame
    nF - Number of Frames
    CF - Id of the valve closing frame 
    PressureChoice - Choice to include pressure magnitude in least squares optimisation
    ProfileChoice - Choice of pressure profile shape function
    RunLSChoice - Choice to run least quares optimisation
    ResChoice -  Choice of residual calculation method 
    ModelChoice - Choice of constitutive model
    Out - Output of Least squares optimisation
    
    '''
    
    #Name of existing .xplt file that has been created 
    xpltName = './FEB_Files/' + DataDir + '.xplt'
    
    if ModelChoice == 'tiMR':
        nDoms = nCls
    else:
        nDoms = nCls
    
    #Get Febio data tree, feb, and number of states, note the number of frames from Febio may not be the same of the original number of frames
    feb, _,nStates, _ = GetFEB(xpltName,nDoms,False)
    
    #Get number of points , nNodes, and number of elements, nElems, number of variables from febio calc, nVar, the times of each frame, StateTimes, get names of variables and types, VarNames and VarTypes 
    nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
    
    #Collect the data
    displacement = GetData(feb,'displacement',nStates)
    stress = GetData(feb,'stress',nStates)
    
    #Gather stresses and displacements Note this is dependent on the Var Type
    Stress_X, Stress_Y, Stress_Z, Stress_XY, Stress_YZ, Stress_XZ = np.zeros((nStates,nElems)),np.zeros((nStates,nElems)),np.zeros((nStates,nElems)),np.zeros((nStates,nElems)),np.zeros((nStates,nElems)),np.zeros((nStates,nElems))
    Disp_X,Disp_Y,Disp_Z = np.zeros((nStates,nNodes)),np.zeros((nStates,nNodes)),np.zeros((nStates,nNodes))
    Disp_FEB  = np.zeros((nStates,nNodes,3))      
    
    #Restructure displacement arrays and Stress arrays
    for i in range(nStates): 
        for j in range(nNodes):
            Disp_X[i,j]  = displacement[i][j*3+1]
            Disp_Y[i,j]  = displacement[i][j*3+2]
            Disp_Z[i,j]  = displacement[i][j*3+3]
        Disp_FEB[i,:,0]  = Disp_X[i,:]  
        Disp_FEB[i,:,1]  = Disp_Y[i,:]  
        Disp_FEB[i,:,2]  = Disp_Z[i,:]
        
    for i in range(nStates): 
        for j in range(nElems):
            Stress_X[i,j]  = stress[i][j*6]
            Stress_Y[i,j]  = stress[i][j*6+1]
            Stress_Z[i,j]  = stress[i][j*6+2]
            Stress_XY[i,j] = stress[i][j*6+3]
            Stress_YZ[i,j] = stress[i][j*6+4]
            Stress_XZ[i,j] = stress[i][j*6+5]
            
    if nF !=nStates:
        #Warning to to highlight difference in original frames and febio frames
        print('Warning: Number of Frames do not match.')
        
        #Define Timesteps
        TimeInterp = np.linspace(0,1,nF)
        
        #Define empty arrays to be filled with interpolate data arrays
        Disp_X_Interp = np.zeros((nF,nNodes))
        Disp_Y_Interp = np.zeros((nF,nNodes))
        Disp_Z_Interp = np.zeros((nF,nNodes))
        Disp_FEB_Interp = np.zeros(Disp_Wall.shape)
        Stress_X_Interp  = np.zeros((nF,nElems))
        Stress_Y_Interp  = np.zeros((nF,nElems))
        Stress_Z_Interp  = np.zeros((nF,nElems))
        Stress_XY_Interp = np.zeros((nF,nElems))
        Stress_YZ_Interp = np.zeros((nF,nElems))
        Stress_XZ_Interp = np.zeros((nF,nElems))
        for i in range(nNodes):
            Disp_X_Interp[:,i] = np.interp(TimeInterp,StateTimes,Disp_X[:,i])
            Disp_Y_Interp[:,i] = np.interp(TimeInterp,StateTimes,Disp_Y[:,i])
            Disp_Z_Interp[:,i] = np.interp(TimeInterp,StateTimes,Disp_Z[:,i])
            for j in range(3):
                Disp_FEB_Interp[:,i,j] = np.interp(TimeInterp,StateTimes,Disp_FEB[:,i,j])
        for i in range(nElems):
            Stress_X_Interp[:,i]  = np.interp(TimeInterp,StateTimes,Stress_X[:,i])
            Stress_Y_Interp[:,i]  = np.interp(TimeInterp,StateTimes,Stress_Y[:,i])
            Stress_Z_Interp[:,i]  = np.interp(TimeInterp,StateTimes,Stress_Z[:,i])
            Stress_XY_Interp[:,i] = np.interp(TimeInterp,StateTimes,Stress_XY[:,i])
            Stress_YZ_Interp[:,i] = np.interp(TimeInterp,StateTimes,Stress_YZ[:,i])
            Stress_XZ_Interp[:,i] = np.interp(TimeInterp,StateTimes,Stress_XZ[:,i])
    else:
        #When the number of frames are equivalent, interpolation is unnecessary
        Disp_X_Interp  = Disp_X
        Disp_Y_Interp  = Disp_Y
        Disp_Z_Interp  = Disp_Z
        Disp_FEB_Interp  = Disp_FEB
        Stress_X_Interp  = Stress_X
        Stress_Y_Interp  = Stress_Y
        Stress_Z_Interp  = Stress_Z
        Stress_XY_Interp = Stress_XY
        Stress_YZ_Interp = Stress_YZ
        Stress_XZ_Interp = Stress_XZ
    
    #Define array of residuals
    if ResChoice == 'P2P':
        #P2P is 'point to point' and uses explicit differences in point positions
        Residual = np.subtract(Disp_Wall,Disp_FEB_Interp)
    elif ResChoice == 'CentreLine':
        nSTJ = len(STJ_Id)
        STJ_Cent = 0
        for i in STJ_Id:
            STJ_Cent += Pts[i]
        STJ_Cent /= nSTJ
        
        nVAJ = len(VAJ_Id)
        VAJ_Cent = 0
        for i in STJ_Id:
            VAJ_Cent += Pts[i]
        VAJ_Cent /= nVAJ
        Residual = np.zeros((nStates,nNodes))
        for i in range(nF):
            for j in range(nNodes):
                VTK_Pt = np.add(Pts[j],Disp_Wall[i,j])
                FEB_Pt = np.add(Pts[j],Disp_FEB_Interp[i,j])
                VTK2L = (np.linalg.norm(np.cross(np.subtract(VTK_Pt,VAJ_Cent),np.subtract(VTK_Pt,STJ_Cent))))/np.linalg.norm(np.subtract(STJ_Cent,VAJ_Cent))
                FEB2L = (np.linalg.norm(np.cross(np.subtract(FEB_Pt,VAJ_Cent),np.subtract(FEB_Pt,STJ_Cent))))/np.linalg.norm(np.subtract(STJ_Cent,VAJ_Cent))
            
                Residual[i,j] = np.absolute(FEB2L - VTK2L)
    elif ResChoice == 'CellPlane':
        # Define pts of FEB simulation
        VTK_Pts = np.zeros((nF,nNodes,3))
        for i in range(nF):
            for j in range(nNodes):
                VTK_Pts[i,j] = np.add(Pts[j],Disp_Wall[i,j])
                
        FEB_Pts = np.zeros((nF,nNodes,3))
        for i in range(nF):
            for j in range(nNodes):
                FEB_Pts[i,j] = np.add(Pts[j],Disp_FEB_Interp[i,j])
                
        Residual = np.zeros((nF,nNodes))
        CellPtSum    = [[] for i in range(nNodes)]
        CellDistId   = [[] for i in range(nNodes)]
        CellTotalSum = [[] for i in range(nNodes)]
        for i in range(nF):
            for j in range(nNodes): 
                if i ==0:
                    for k in range(nCls):
                        if j in CellIds[k]:
                            CellPtSum[j].append([])
                            CellDistId[j].append(CellIds[k].astype(int))
                            
                    for k in range(len(CellDistId[j])):
                        for l in np.array(CellDistId[j][k]).astype(int):
                            CellPtSum[j][k].append(np.linalg.norm(np.subtract(VTK_Pts[i,j],FEB_Pts[i,l])))
                        CellTotalSum[j].append(np.sum(CellPtSum[j][k]))
                    
                MinCl = np.argmin(CellTotalSum[j])
                
                Pt1, Pt2, Pt3 = FEB_Pts[i,CellDistId[j][MinCl][np.argpartition(CellPtSum[j][MinCl],3)[0]]], FEB_Pts[i,CellDistId[j][MinCl][np.argpartition(CellPtSum[j][MinCl],3)[1]]], FEB_Pts[i,CellDistId[j][MinCl][np.argpartition(CellPtSum[j][MinCl],3)[2]]]
                
                planeNorm = np.cross(np.subtract(Pt1,Pt2),np.subtract(Pt1,Pt3))
                planePt = Pt1
                
                PtNorm = Norm[i,j]
                Pt     = VTK_Pts[i,j]
                
                ndotu = planeNorm.dot(PtNorm)
                if abs(ndotu) < 1e-6:
                    Residual[i,j] =  0.0
                    raise RuntimeError("no intersection or line is within plane")
                else:
                    w = Pt - planePt
                    si = -planeNorm.dot(w) / ndotu
                    Intersect = np.add(np.add(w,np.multiply(si,PtNorm)),planePt)
                    
                    Residual[i,j] =  np.linalg.norm(np.subtract(Pt,Intersect))
     
    if FibChoice and ModelChoice == 'tiMR':
        #Get Fibers with fitted angle              
        nC = 0
        if PressureChoice:
            nC += 1
            
        if ModelChoice == 'MR':
            nC += 3
        if ModelChoice == 'tiMR':
            nC += 8
        elif ModelChoice == 'Ogden':
            nC +=14
        elif ModelChoice == 'Fung':
            nC +=11
        elif ModelChoice == 'HGO':
            nC +=6
        
        theta = Out[0+nC]
        FEB_Fib = np.zeros((nCls,3))
        for i in range(nCls):
            C_tick = np.cross(np.cross(Circ_Cls[i],Long_Cls[i]),Circ_Cls[i])
            FEB_Fib[i] = np.dot(cos(theta),Circ_Cls[i]) + np.dot(sin(theta),C_tick)        
     
    reader = vtk.vtkPolyDataReader()
    
    print('Writing New Files...')
    
    FList = glob.glob('./Remeshed/'+DataDir+'/*')
    
    #Order list to put reference frame first
    FListOrdered, FId, refN = OrderList(FList, nF, RemeshedFile )
        
    #Save a new vtk file for each state
    for j, fid in enumerate((FId.astype(int)-1)):
        
        # Read the source file.
        reader.SetFileName(FListOrdered[j])
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        polydata = reader.GetOutput()
        
        # Add Cell Data to Cells
        CellData = [Stress_X_Interp[j], Stress_Y_Interp[j], Stress_Z_Interp[j], Stress_XY_Interp[j], Stress_YZ_Interp[j], Stress_XZ_Interp[j]]
        CellNames = ['Stress_X','Stress_Y','Stress_Z','Stress_XY','Stress_YZ','Stress_XZ']
        
        for i in range(len(CellNames)) :
            arrayCell = vtk.util.numpy_support.numpy_to_vtk(CellData[i], deep=True)
            arrayCell.SetName(CellNames[i])
            dataCells = polydata.GetCellData()
            dataCells.AddArray(arrayCell)
            dataCells.Modified()
            
        # Convert Cell Data to Point Data
        c2p = vtk.vtkCellDataToPointData()
        c2p.AddInputData(polydata)
        c2p.Update()
        c2p.GetOutput()
        NumOfArr = c2p.GetPolyDataOutput().GetPointData().GetNumberOfArrays()
        
        for i in range(NumOfArr):
            if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_X':
                Stress_X_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
            if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_Y':
                Stress_Y_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
            if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_Z':
                Stress_Z_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
            if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_XY':
                Stress_XY_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
            if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_YZ':
                Stress_YZ_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
            if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Stress_XZ':
                Stress_XZ_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
        
        # Add Point Data
        PointData = [Stress_X_Pt,Stress_Y_Pt,Stress_Z_Pt,Stress_XY_Pt,Stress_YZ_Pt,Stress_XZ_Pt,Disp_X_Interp[j],Disp_Y_Interp[j],Disp_Z_Interp[j],Disp_FEB_Interp[j],Residual[j]]
        PointNames = ['Stress_X_Pt','Stress_Y_Pt','Stress_Z_Pt','Stress_XY_Pt','Stress_YZ_Pt','Stress_XZ_Pt','Disp_X','Disp_Y','Disp_Z','Disp_FEB','Residual']
        
        for i in range(len(PointNames)) :
            arrayPoint = vtk.util.numpy_support.numpy_to_vtk(PointData[i], deep=True)
            arrayPoint.SetName(PointNames[i])
            dataPoints = polydata.GetPointData()
            dataPoints.AddArray(arrayPoint)
            dataPoints.Modified()
            
        if FibChoice and ModelChoice == 'tiMR':
            # Add Vector Data on Points
            VectorData = [FEB_Fib]
            VectorNames = ['Fibers_FEB']
            for i in range(len(VectorNames)) :
                arrayVector = vtk.util.numpy_support.numpy_to_vtk(VectorData[i], deep=True)
                arrayVector.SetName(VectorNames[i])
                dataVectors = polydata.GetCellData()
                dataVectors.AddArray(arrayVector)
                dataVectors.Modified()
            
        Conditions = ''
        if PressureChoice:
            Conditions += 'PMagT_'
        else:
            Conditions += 'PMagF_' 
            
        if ModelChoice[0:3] != 'Set':
            Conditions += 'ModelT_'
        else:
            Conditions += 'ModelF_'
        
        if FibChoice and ModelChoice == 'tiMR': 
            Conditions += 'FibT_'
        else:
            Conditions += 'FibF_'
        
        Conditions += ResChoice
        
        if RunLSChoice:
            RunDir = 'RunLS'
        else:
            RunDir = 'RunDefault'
        
        fname = './ParSweep/Case_'+str(Count)+'/'+DataDir + '_' + str(fid+1)+'.vtk'
        
        directory = os.path.dirname(fname)
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        writer = vtk.vtkDataSetWriter()
        writer.SetFileName(fname)
        writer.SetInputData(polydata)
        writer.Write()
        
        feb_name = './FEB_Files/'+ DataDir+'.feb'
        
        os.system('cp '+ feb_name + ' ./ParSweep/Case_'+str(Count)+'/')
        
    return nStates, Residual
  
    
    
if __name__=='__main__':  
    #Choose parameter estimation choices
    NoVar          = True               # Choose to do parameter variation
    PressureChoice = False              # Choose to vary pressure magnitude
    RunLSChoice    = True               # Choose to run least Squares (or default/initial guess)
    FibChoice      = False              # Choose to vary fiber direction, as an angle from the circumferential direction
    ProfileChoice  = ['SetWindkessel']  # Choose profile shapes, options are: 'Triangle','Step','SmoothStep','Virtual', 'Fourier','Fitted'm Windkessel'
    ResChoice      = ['CellPlane']      # Choose type of residual calculation method: 'P2P', 'CentreLine', 'CellPlane'
    ModelChoice    = ['HGO']            # Choose model from 'MR','tiMR','Ogden' and 'Fung',  'HGO'
    DataDirChoice  = 'Specific'              # Choose thewhich directories are to be included: 'Specfic', 'SpecificGroup', 'All_TAVs', 'All_BAVs','All'
    
    
    # Code to run through all data sets
    List_of_Subdirectories = sorted(glob.glob('./Strains/medial_meshes/*'))
    ND = len(List_of_Subdirectories)
    CommonOfDir = os.path.commonprefix(List_of_Subdirectories)
    
    DataDirs = []
    if DataDirChoice == 'Specific':
        DataDirs = ['tav02']
    elif DataDirChoice == 'SpecificGroup':
        DataDirs = ['tav04', 'tav06','tav11','tav20','tav23','tav26']
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
        FList = glob.glob('./Remeshed/'+DataDir+'/*')
        if FList ==[]:
            OriginalFiles = sorted(glob.glob('./Strains/medial_meshes/'+DataDir+'/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/*'))
            for OF in OriginalFiles:
                if OF[-4:]==".vtp"or OF[-4:]==".vtk":
                    Remesh_All(DataDir,OF)
                    
    for DataDir in DataDirs:
        print('Running Data Set: ',DataDir)
    
        # FList = glob.glob('./Remeshed/'+DataDir+'/*')
        FList = glob.glob('./ParSweep/Case_1/*.vtk')
        
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
            
        for PC in ProfileChoice:
            for MC in ModelChoice:
                for RC in ResChoice:
                    
                    Outs = []
                    Residuals = []
                    Outs.append([])
                    Residuals.append([])
                    with open('./NewFiles/tav02/HGO/RunLS/SetWindkessel/PMagF_ModelT_FibF_CellPlane/Parameters.csv') as csv_file:
                        XLData = csv.reader(csv_file, delimiter=',')
                        next(XLData, None)
                        for row in XLData:
                            if row[2] != '':
                                Outs[0].append(row[3])
                    FList_orig = glob.glob('./NewFiles/tav02/HGO/RunLS/SetWindkessel/PMagF_ModelT_FibF_CellPlane/*.vtk')
                    
                    #Get name of first and reference file
                    common = os.path.commonprefix(FList_orig)
                    for Fname in list(FList_orig):
                        X = Fname.replace(common,'')
                        X = X.replace('.vtk','')
                        X = np.fromstring(X, dtype=int, sep=' ')
                        X=X[0]
                        if X==refN:
                            ref_orig=Fname
                            
                    FListOrdered_orig, _, _ = OrderList(FList_orig, nF, ref_orig)
                
                    for X, Fname in enumerate(FListOrdered_orig):                            
                        reader = vtk.vtkPolyDataReader()
                        reader.SetFileName(Fname)
                        reader.ReadAllScalarsOn()
                        reader.ReadAllVectorsOn()
                        reader.Update()
                        
                        polydata = reader.GetOutput()  
                        # Pts = vtk_to_numpy(polydata.GetPoints().GetData())
                        nPts = polydata.GetNumberOfPoints()
                        Res = vtk_to_numpy(polydata.GetPointData().GetArray('Residual'))
                        
                        Residuals[0].append(Res)
                    Count = 0  
                    Par1 = [0,0.5,1.0]
                    Par2 = [0,0.5,1.0]
                    Par3 = [0,0.5,1.0]
                    Par4 = [0,0.5,1.0]
                    if NoVar:
                        Par1,Par2,Par3,Par4 = [1],[1],[1],[1]
                    for P1 in Par1:
                        for P2 in Par2:
                            for P3 in Par3:
                                for P4 in Par4:
                            
                                    Count+=1
                                    print('Run Count:',Count)
                                    HGOPars = [100,25,90,25]#[100*P1, 1000*P2, 100*P3, 50*P4]
                                            
                                    print('Run Case:',HGOPars)
                                    # Initialise initial parameter arrays
                                    InitParams = [[],[],[],[]]
                                    
                                    # Choose initial estimate of pressure magnitude
                                    if PressureChoice:
                                        InitParams[0] = [0.126]
                                     
                                    # Choose initial model parameters
                                    if ModelChoice[0:3] != 'Set':
                                        if MC == 'MR':
                                            InitParams[1] = [1,10,10,10]
                                        if MC == 'tiMR':
                                            InitParams[0] = [10,10,10,10,10,10,10,10]
                                        elif MC == 'Ogden':
                                            InitParams[0] = [1,10,10,10,10,10,10,10,10,10,10,10,10,10]
                                        elif MC == 'Fung':
                                            InitParams[0] = [1,1,2,1.5,0.5,0.1,1,1.2,2,1.5,0.5,0.1]
                                        elif MC == 'HGO':     
                                            InitParams[0] = HGOPars
                                      
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
                                        RunLSChoice = False
                                        
                                    #Run least squares script
                                    Out, Pts, Disp_Wall, Norm, Circ, Long, CellIds, nCls, STJ_Id, VAJ_Id = RunLS(DataDir,FListOrdered,FId,ref,CF,C,PressureChoice,PC,RunLSChoice,RC,MC,FibChoice)
                                    
                                    # Save new files
                                    nStates,Residual = SaveFiles(DataDir,FList,ref,Pts,Disp_Wall,Norm,Circ,Long,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,PC,RunLSChoice,RC,MC,FibChoice,Out,Count)
                        
                                    Outs.append(Out)
                                    Residuals.append(Residual)
                                    
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
                                        nC +=1
                                        Conditions += 'PMagT_'
                                    else:
                                        Params[0].append('Pressure_Magnitude')
                                        Params[1].append('False')
                                        Params[4].append('Pressure_Magnitude')
                                        Params[5].append('1.0')
                                        print('Pressure Magnitude is Default')
                                        Conditions += 'PMagF_'
                                        
                                    if ModelChoice[0:3] != 'Set':
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
                                            nP = 4
                                            ParamNames  = ['c','k1','k2','gamma','kappa','k']
                                            ParamValues = Out[nC:nC+nP+1]
                                            for i in range(nP):
                                                print(ParamNames[i],': ',ParamValues[i])
                                                Params[2].append(ParamNames[i])
                                                Params[3].append(ParamValues[i])
                                            for i in range(nP):
                                                Params[2].append(ParamNames[i]+'_Starting')
                                                Params[3].append(HGOPars[i])
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
                    
                                    with open('./ParSweep/'+ 'Case_'+str(Count)+'/Parameters.csv','w') as result_file:
                                        wr = csv.writer(result_file, dialect='excel')
                                        for values in zip_longest(*Params):
                                            wr.writerow(values)
                                
                                    # fname = './NewFiles/'+ DataDir+'/'+ModelChoice+'/'+RunDir + '/' + ProfileChoice+ '/' + Conditions+'/*'
                                    # os.system('cp '+ fname + ' ./ParSweep/'+ 'Case_'+str(Count)+'/')
                                    
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
            
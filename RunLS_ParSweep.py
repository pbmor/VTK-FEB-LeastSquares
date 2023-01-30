import numpy as np
import os
import csv
import glob
import vtk
import matplotlib.pyplot as plt
from math import cos, sin 
from xml.etree import ElementTree as et
from vtk.util.numpy_support import vtk_to_numpy
from Read_XPLTfuncs import GetFEB, GetMeshInfo, GetData
from ResidualFunction_Time import RunLS, OrderList, GetPressure
from RemeshAll import Remesh_All
from itertools import zip_longest

def SaveFiles(DataDir,FList,RemeshedFile,Pts,Disp_Wall,Norm,Circ_Cls,Long_Cls,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,ProfileChoice,RunLSChoice,ResChoice,ModelChoice,FibChoice,Out,Count,SetP):
    '''
    Read xplt file and save to remeshed vtk file and save as a new vtk file.
    
    Inputs:
    DataDir - Name of directory, e.g. 'tav02'
    FList -- list of filenames
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
    # Use xml tree to update .feb file
    tree = et.parse('./FEB_Files/' + DataDir + '.feb')
    # Count the number of choices
    nC = 0
    # Define pressure magnitude
    if PressureChoice:
        if ModelChoice in ['MR','SetMR']:
            tree.find('Loads').find('surface_load').find('pressure').text = str(C[0])
        else:
            tree.find('Loads').find('surface_load').find('pressure').text = str(-C[0])
        nC += 1
    
    # Replace material properties with parameter estimates
    if ModelChoice[0:3] != 'Set':
        if ModelChoice == 'MR':
            tree.find('Material').find('material').find('density').text = str(C[0+nC])
            tree.find('Material').find('material').find('c1').text      = str(C[1+nC])
            tree.find('Material').find('material').find('c2').text      = str(C[2+nC])
            # tree.find('Material').find('material').find('k').text       = str(C[3+nC])
            nC += 3
        if ModelChoice == 'tiMR':
            TreeBranches = tree.find('Material').findall('material')
            for i in range(nCls):
                TreeBranches[i].find('density').text = str(C[0+nC])
                TreeBranches[i].find('c1').text = str(C[1+nC])
                TreeBranches[i].find('c2').text = str(C[2+nC])
                TreeBranches[i].find('c3').text = str(C[3+nC])
                TreeBranches[i].find('c4').text = str(C[4+nC])
                TreeBranches[i].find('c5').text = str(C[5+nC])
                TreeBranches[i].find('lam_max').text = str(C[6+nC])
                TreeBranches[i].find('k').text = str(C[7+nC])
            nC += 8
        elif ModelChoice == 'Ogden':
            TreeBranches = tree.find('Material').findall('material')
            for i in range(nCls):
                TreeBranches[i].find('density').text  = str(C[0+nC])
                TreeBranches[i].find('k').text  = str(C[1+nC])
                TreeBranches[i].find('c1').text = str(C[2+nC])
                TreeBranches[i].find('c2').text = str(C[3+nC])
                TreeBranches[i].find('c3').text = str(C[4+nC])
                TreeBranches[i].find('c4').text = str(C[5+nC])
                TreeBranches[i].find('c5').text = str(C[6+nC])
                TreeBranches[i].find('c6').text = str(C[7+nC])
                TreeBranches[i].find('m1').text = str(C[8+nC])
                TreeBranches[i].find('m2').text = str(C[9+nC])
                TreeBranches[i].find('m3').text = str(C[10+nC])
                TreeBranches[i].find('m4').text = str(C[11+nC])
                TreeBranches[i].find('m5').text = str(C[12+nC])
                TreeBranches[i].find('m6').text = str(C[13+nC])
            nC +=14
        elif ModelChoice == 'Fung':
            TreeBranches = tree.find('Material').findall('material')
            for i in range(nCls):
                TreeBranches[i].find('density').text = str(C[0+nC])
                TreeBranches[i].find('E1').text      = str(C[1+nC])
                TreeBranches[i].find('E2').text      = str(C[2+nC])
                TreeBranches[i].find('E3').text      = str(C[3+nC])
                TreeBranches[i].find('G12').text     = str(C[4+nC])
                TreeBranches[i].find('G23').text     = str(C[5+nC])
                TreeBranches[i].find('G31').text     = str(C[6+nC])
                TreeBranches[i].find('v12').text     = str(C[7+nC])
                TreeBranches[i].find('v23').text     = str(C[8+nC])
                TreeBranches[i].find('v31').text     = str(C[9+nC])
                TreeBranches[i].find('c').text       = str(C[10+nC])
                TreeBranches[i].find('k').text       = str(C[11+nC])
            nC +=12
        elif ModelChoice == 'HGO':
            TreeBranches = tree.find('Material').findall('material')
            for i in range(nCls):
                TreeBranches[i].find('c').text     = str(C[0+nC])
                TreeBranches[i].find('k1').text    = str(C[1+nC])
                TreeBranches[i].find('k2').text    = str(C[2+nC])
                TreeBranches[i].find('gamma').text = str(C[3+nC])
            nC +=4
     
        
        if FibChoice and ModelChoice == 'tiMR':
            theta = C[0+nC]
            for i in range(nCls):
                C_tick = np.cross(np.cross(Circ[i],Long[i]),Circ[i])
                New = np.dot(cos(theta),Circ[i]) + np.dot(sin(theta),C_tick)
                TreeBranches[i].find('fiber').text = str(New[0])+str(',')+str(New[1])+str(',')+str(New[2])
            nC += 1
    
    #Define Timesteps
    TimeInterp = np.linspace(0,1,nF)        
    #Get Pressure Profile
    if ProfileChoice[0:3] != 'Set':
        PressureInterp = GetPressure(C,nC,nF,CF,ProfileChoice,TimeInterp,SetP)
        for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
                Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx])
               
    # Rewrite .feb file with paramter updates
    tree.write('./FEB_Files/' + DataDir + '.feb',xml_declaration=True,encoding="ISO-8859-1")
    
    #Run updated .feb file to create new .xplt file
    os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i '+ './FEB_Files/' + DataDir+'.feb >/dev/null 2>&1')
    
    #Name of existing .xplt file that has been created 
    xpltName = './FEB_Files/' + DataDir + '.xplt'
    
    if ModelChoice in ['MR','SetMR']:
        nDoms = 1
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
        Residual = np.zeros((nF,nNodes))
        for i in range(nF):
            for j in range(nNodes):
                Residual[i,j] = np.linalg.norm(np.subtract(Disp_Wall[i,j],Disp_FEB_Interp[i,j]))
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
    #Order list to put reference frame first
    FListOrdered, FId, refN = OrderList(FList, nF-1, RemeshedFile )
    #Save a new vtk file for each state
    for j, fid in enumerate((FId.astype(int))):
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
        
        fname = './ParSweep/'+DataDir+'/Case_'+str(Count)+'/'+DataDir + '_' + str(fid)+'.vtk'
        print(fname)
        
        directory = os.path.dirname(fname)
        if not os.path.exists(directory):
            os.makedirs(directory)
            
        writer = vtk.vtkDataSetWriter()
        writer.SetFileName(fname)
        writer.SetInputData(polydata)
        writer.Write()
        
        feb_name  = './FEB_Files/'+ DataDir+'.feb'
        xplt_name = './FEB_Files/'+ DataDir+'.xplt'
        log_name  = './FEB_Files/'+ DataDir+'.log'
        
        os.system('cp '+ feb_name + ' ./ParSweep/'+DataDir+'/Case_'+str(Count)+'/')
        os.system('cp '+ xplt_name + ' ./ParSweep/'+DataDir+'/Case_'+str(Count)+'/')
        os.system('cp '+ log_name + ' ./ParSweep/'+DataDir+'/Case_'+str(Count)+'/')
        
    return nStates, Residual

    
if __name__=='__main__':  
    #Choose parameter estimation choices
    NoVar          = True               # Choose to do parameter variation
    PressureChoice = False              # Choose to vary pressure magnitude
    RunLSChoice    = True               # Choose to run least Squares (or default/initial guess)
    FibChoice      = False              # Choose to vary fiber direction, as an angle from the circumferential direction
    ProfileChoice  = ['SetWindkessel']  # Choose profile shapes, options are: 'Triangle', 'Step', 'SmoothStep', 'Virtual', 'Fourier','Fitted', 'Windkessel'
    ResChoice      = ['CellPlane']      # Choose type of residual calculation method: 'P2P', 'CentreLine', 'CellPlane'
    ModelChoice    = ['HGO']            # Choose model from 'MR','tiMR','Ogden' and 'Fung',  'HGO'
    DataDirChoice  = 'Specific'         # Choose thewhich directories are to be included: 'Specfic', 'SpecificGroup', 'All_TAVs', 'All_BAVs','All'
    
    
    # Code to run through all data sets
    List_of_Subdirectories = sorted(glob.glob('./Strains/medial_meshes/*'))
    ND = len(List_of_Subdirectories)
    CommonOfDir = os.path.commonprefix(List_of_Subdirectories)
    
    DataDirs = []
    if DataDirChoice == 'Specific':
        DataDirs = ['bav07']
    elif DataDirChoice == 'SpecificGroup':
        DataDirs = ['tav02', 'bav07']
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
                    
    for DataDir in DataDirs:
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
                    
                    Outs = []
                    Residuals = []
                    Outs.append([])
                    Residuals.append([])
                    Count = 0  
                    Par1 = [1/3,2/3,1]
                    Par2 = np.linspace(0,1,3)
                    Par3 = [0.5]
                    Par4 = np.linspace(0,1,3)
                    if NoVar:
                        Par1,Par2,Par3,Par4 = [1],[1],[1],[1]
                    for P1 in Par1:
                        for P2 in Par2:
                            for P3 in Par3:
                                for P4 in Par4:
                                    # Initialise empty array for Setting pressure parameters
                                    SetP=[[],[]]
                            
                                    Count+=1
                                    print('Run Count:',Count)
                                    HGOPars = [10*P1,1000*P2,1000*P3,90*P4] #[ 5,  415,  5, 28]#[100*P1, 1000*P2, 100*P3, 50*P4]
                                            
                                    print('Run Case:',HGOPars)
                                    # Initialise initial parameter arrays
                                    InitParams = [[],[],[],[]]
                                    
                                    # Choose initial estimate of pressure magnitude
                                    if PressureChoice:
                                        InitParams[0] = [5]
                                    else:
                                        SetP[0] = [5]
                                     
                                    # Choose initial model parameters
                                    if ModelChoice[0:3] != 'Set':
                                        if MC == 'MR':
                                            InitParams[1] = [1,10,10,10]
                                        if MC == 'tiMR':
                                            InitParams[1] = [10,10,10,10,10,10,10,10]
                                        elif MC == 'Ogden':
                                            InitParams[1] = [1,10,10,10,10,10,10,10,10,10,10,10,10,10]
                                        elif MC == 'Fung':
                                            InitParams[1] = [1,1,2,1.5,0.5,0.1,1,1.2,2,1.5,0.5,0.1]
                                        elif MC == 'HGO':     
                                            InitParams[1] = HGOPars
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
                                    
                                    # Save new files
                                    nStates,Residual = SaveFiles(DataDir,FList,ref,Pts,Disp_Wall,Norm,Circ,Long,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,PC,RunLSChoice,RC,MC,FibChoice,Out,Count,SetP)
                        
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
                                        Params[5].append(SetP[0][0])
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
                                            nP = 3
                                            # ParamNames  = ['c','k1','k2','gamma','kappa','k']
                                            ParamNames  = ['c','k1','gamma']
                                            ParamValues = Out[nC:nC+nP+1]
                                            for i in range(nP):
                                                print(ParamNames[i],': ',ParamValues[i])
                                                Params[2].append(ParamNames[i])
                                                Params[3].append(ParamValues[i])
                                            for i in range(nP):
                                                Params[2].append(ParamNames[i]+'_Starting')
                                                Params[3].append(HGOPars[i])
                                            nC+=3
                                    else:
                                        print('Model parameters not optimised')
                                        Conditions += 'ModelF_'
                                        Params[0].append('Model_Parameters')
                                        Params[1].append('False')
                                        if MC == 'SetMR':
                                            Params[0].append('Model')
                                            Params[1].append('MR')
                                            nP = 3
                                            ParamNames  = ['density','C1','C2']
                                            for i in range(3):
                                                Params[2].append(ParamNames[i]+'_Starting')
                                                Params[3].append(ModelPar[i])
                                        elif MC == 'SettiMR':
                                            Params[0].append('Model')
                                            Params[1].append('tiMR')
                                            nP = 8
                                            ParamNames  = ['density','C1','C2','C3','C4','C5','k','lam_max']
                                            for i in range(nP):
                                                Params[2].append(ParamNames[i]+'_Starting')
                                                Params[3].append(ModelPar[i])
                                        elif MC == 'SetOgden':
                                            Params[0].append('Model')
                                            Params[1].append('Ogden')
                                            nP = 14
                                            ParamNames  = ['density','k','c1','c2','c3','c4','c5','c6','m1','m2','m3','m4','m5','m6']
                                            for i in range(nP):
                                                Params[2].append(ParamNames[i]+'_Starting')
                                                Params[3].append(ModelPar[i])
                                        elif MC == 'SetFung':
                                            Params[0].append('Model')
                                            Params[1].append('Fung')
                                            nP = 12
                                            ParamNames  = ['density','E1','E2','E3','G12','G23','G31','v12','v23','v31','c','k']
                                            for i in range(nP):
                                                Params[2].append(ParamNames[i]+'_Starting')
                                                Params[3].append(ModelPar[i])
                                        elif MC == 'SetHGO':
                                            Params[0].append('Model')
                                            Params[1].append('HGO')
                                            nP = 4
                                            ParamNames  = ['c','k1','k2','gamma']
                                            for i in range(nP):
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
                
                                    # Save Final Cost
                                    Params[0].append('Cost')
                                    Params[1].append(Cost)
                    
                                    with open('./ParSweep/'+ DataDir+'/Case_'+str(Count)+'/Parameters.csv','w') as result_file:
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


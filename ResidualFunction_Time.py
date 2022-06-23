import numpy as np
import os
import vtk
import csv
import glob
from xml.etree import ElementTree as et
from scipy.optimize import least_squares
from vtk.util.numpy_support import vtk_to_numpy
from Read_XPLTfuncs import GetFEB, GetMeshInfo, GetData
from vtk2feb_Disp import VTK2Feb_Func
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from math import sin, pi
from math import exp as e

def SaveFiles(DataDir,RemeshedFile,Pts,Disp_Wall,Norm,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,ModelParChoice,ProfileChoice,RunLSChoice,ResChoice,ModelChoice,Out):
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
    ModelParChoice - Choice to include model parameters in least squares optimisation
    ProfileChoice - Choice of pressure profile shape function
    RunLSChoice - Choice to run least quares optimisation
    ResChoice -  Choice of residual calculation method 
    ModelChoice - Choice of constitutive model
    Out - Output of Least squares optimisation
    
    '''
    
    #Name of existing .xplt file that has been created 
    xpltName = './FEB_Files/' + DataDir + '.xplt'
    
    #Get Febio data tree, feb, and number of states, note the number of frames from Febio may not be the same of the original number of frames
    feb, _,nStates, _ = GetFEB(xpltName)
    
    #Get number of points , nNodes, and number of elements, nElems, number of variables from febio calc, nVar, the times of each frame, StateTimes, get names of variables and types, VarNames and VarTypes 
    nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
    
    #Collect the data
    displacement = GetData(feb,'displacement',nStates,nVar)
    stress = GetData(feb,'stress',nStates,nVar)
    
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
            Stress_X[i,j]  = stress[i][j*6+1]
            Stress_Y[i,j]  = stress[i][j*6+2]
            Stress_Z[i,j]  = stress[i][j*6+3]
            Stress_XY[i,j] = stress[i][j*6+4]
            Stress_YZ[i,j] = stress[i][j*6+5]
            Stress_XZ[i,j] = stress[i][j*6+6]
            
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
        for i in range(nStates):
            for j in range(nNodes):
                VTK_Pt = np.add(Pts[j],Disp_Wall[i,j])
                FEB_Pt = np.add(Pts[j],Disp_FEB[i,j])
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
                    # raise RuntimeError("no intersection or line is within plane")
                else:
                    w = Pt - planePt
                    si = -planeNorm.dot(w) / ndotu
                    Intersect = np.add(np.add(w,np.multiply(si,PtNorm)),planePt)
                    
                    Residual[i,j] =  np.linalg.norm(np.subtract(Pt,Intersect))
                # if Residual[i,j] >1000:
                #     ax = plt.axes(projection='3d')
                #     # ax.scatter3D(VTK_Pts[i,:,0],VTK_Pts[i,:,1],VTK_Pts[i,:,2],'xk')
                #     # if i ==5 and j==48:
                #     for l in range(len(CellDistId[j])):
                #         ax.scatter3D(VTK_Pts[i,CellDistId[j][l].astype(int),0],VTK_Pts[i,CellDistId[j][l].astype(int),1],VTK_Pts[i,CellDistId[j][l].astype(int),2], c='r', marker='x')
                #         ax.scatter3D(FEB_Pts[i,CellDistId[j][l].astype(int),0],FEB_Pts[i,CellDistId[j][l].astype(int),1],FEB_Pts[i,CellDistId[j][l].astype(int),2], c='g', marker='x')
                #     ax.scatter3D(VTK_Pts[i,j,0],VTK_Pts[i,j,1],VTK_Pts[i,j,2], c='m', marker='x')
                #     ax.scatter3D(Pt1[0],Pt1[1],Pt1[2], c='b', marker='x')
                #     ax.scatter3D(Pt2[0],Pt2[1],Pt2[2], c='b', marker='x')
                #     ax.scatter3D(Pt3[0],Pt3[1],Pt3[2], c='b', marker='x')
                #     ax.scatter3D(Intersect[0],Intersect[1],Intersect[2], c='b', marker='o')
                #     plt.show()
                    
    reader = vtk.vtkPolyDataReader()
    
    print('Writing New Files...')
    #Save a new vtk file for each state
    for j, fid in enumerate((FId.astype(int)-1)):
        
        # Read the source file.
        reader.SetFileName(RemeshedFile)
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
            
        Conditions = ''
        if PressureChoice:
            Conditions += 'PMagT_'
        else:
            Conditions += 'PMagF_' 
            
        if ModelParChoice:
            Conditions += 'ModelT_'
        else:
            Conditions += 'ModelF_'
            
        if ModelChoice == 'MR':
            Conditions += 'MR_'
        elif ModelChoice == 'tiMR':
            Conditions += 'tiMR_'
        elif ModelChoice == 'Ogden':
            Conditions += 'Ogd_'
        elif ModelChoice == 'Fung':
            Conditions += 'Fung_'
            
        if RunLSChoice:
            Conditions += 'RunLST_'
        else:
            Conditions += 'RunLSF_'
        
        
        Conditions += ProfileChoice +'_'+ ResChoice
            
        
        fname = './NewFiles/'+ DataDir+'/'+Conditions+'/'+xpltName[0:5] + '_' + str(fid+1)+'.vtk'
        directory = os.path.dirname(fname)
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        writer = vtk.vtkDataSetWriter()
        writer.SetFileName(fname)
        writer.SetInputData(polydata)
        writer.Write()
        
        np.save('./NewFiles/'+ DataDir+'/'+Conditions+'/'+'Params_'+Conditions,Out)
        
    return nStates, Residual
    
def OrderList(flist,N,ref):
    '''
    Reorders the list so that the first file is the reference file, and returns new list of filenames and their IDs 
    Keyword arguments:
    flist -- list of filenames
    N -- number of files
    ref -- reference file name
    For example: for a list [file1.vtk, file2.vtk, file3.vtk, file4.vtk, file7.vtk, file8.vtk] and ref = file3.vtk
    it will return [file3.vtk file4.vtk file7.vtk file8.vtk file1.vtk file2.vtk], [3 4 7 8 1 2], and 3
    '''
    # Order filenames so that reference frame goes first
    Fno = np.zeros(N)
    FId = np.zeros(N)

    FListOrdered = [None]*N
    common = os.path.commonprefix(flist)
    for i, Fname in enumerate(flist):
        X = Fname.replace(common,'')
        X = X.replace('.vtk','')
        X = np.fromstring(X, dtype=int, sep=' ')
        #Get list of frame labels
        Fno[i] = X
        # Get label of reference frame
        if Fname==ref:
            refN = X
    # Sort fname labels
    Fno.sort()
    
    #Find Id of reference frame
    for i,X in enumerate(Fno):
        if X==refN:
            RefId = i

    # Sort reference area to be first in list
    FId[0:N-RefId] = Fno[RefId:N]
    FId[N-RefId:N]   = Fno[0:RefId]

    # Get list of file names in new order
    for i,F in enumerate(FId):
        for Fname in flist:
            X = Fname.replace(common,'')
            X = X.replace('.vtk','')
            X = np.fromstring(X, dtype=int, sep=' ')
            if X ==F:
                FListOrdered[i] = Fname

    return FListOrdered, FId, refN

def GetRes(C,DataDir,Pts,Disp_Wall,Norm,Circ,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,ModelParChoice,ProfileChoice,ResChoice,ModelChoice):
    '''
    Function to calculate residual with a a given set of parameters
    
    Inputs:
    C - array of parameters used
    DataDir - Name of directory, e.g. 'tav02'
    Pts - Array of points
    Disp_Wall - Array of wall displacements for each point of mesh
    Norm - Array of Normal directions of longitudinal and circumferential arrays of points
    Circ - Array of circumferential directional arrays of points
    CellIds - Array of point Ids in each cell
    nCls - Number of cells
    STJ_Id - Point Ids for the STJ
    VAJ_Id - Point Ids for the VAJ
    FId - Array of the ordered frame Ids, starting the reference frame
    nF - Number of Frames
    CF - Id of the valve closing frame 
    PressureChoice - Choice to include pressure magnitude in least squares optimisation
    ModelParChoice - Choice to include model parameters in least squares optimisation
    ProfileChoice - Choice of pressure profile shape function
    RunLSChoice - Choice to run least quares optimisation
    ResChoice -  Choice of residual calculation method 
    ModelChoice - Choice of constitutive model
    '''
    #Define Timesteps
    TimeInterp = np.linspace(0,1,nF)
    
    # Use xml tree to update .feb file
    tree = et.parse('./FEB_Files/' + DataDir + '.feb')
    # Count the number of choices
    nC = 0
    # Define pressure magnitude
    if PressureChoice:
        tree.find('Loads').find('surface_load').find('pressure').text = str(C[0])
        nC += 1
      
    # Replace material properties with parameter estimates
    if ModelParChoice:
        if ModelChoice == 'MR':
            tree.find('Material').find('material').find('density').text = str(C[0+nC])
            tree.find('Material').find('material').find('c1').text = str(C[1+nC])
            tree.find('Material').find('material').find('c2').text = str(C[2+nC])
            tree.find('Material').find('material').find('k').text = str(C[3+nC])
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
            nC += 7
        elif ModelChoice == 'Ogden':
            tree.find('Material').find('material').find('k').text = str(C[0+nC])
            tree.find('Material').find('material').find('c1').text = str(C[1+nC])
            tree.find('Material').find('material').find('c2').text = str(C[2+nC])
            tree.find('Material').find('material').find('c3').text = str(C[3+nC])
            tree.find('Material').find('material').find('c4').text = str(C[4+nC])
            tree.find('Material').find('material').find('c5').text = str(C[5+nC])
            tree.find('Material').find('material').find('c6').text = str(C[6+nC])
            tree.find('Material').find('material').find('m1').text = str(C[7+nC])
            tree.find('Material').find('material').find('m2').text = str(C[8+nC])
            tree.find('Material').find('material').find('m3').text = str(C[9+nC])
            tree.find('Material').find('material').find('m4').text = str(C[10+nC])
            tree.find('Material').find('material').find('m5').text = str(C[11+nC])
            tree.find('Material').find('material').find('m6').text = str(C[12+nC])
            nC +=13
        elif ModelChoice == 'Fung':
            tree.find('Material').find('material').find('E1').text = str(C[0+nC])
            tree.find('Material').find('material').find('E2').text = str(C[1+nC])
            tree.find('Material').find('material').find('E3').text = str(C[2+nC])
            tree.find('Material').find('material').find('G12').text = str(C[3+nC])
            tree.find('Material').find('material').find('G23').text = str(C[4+nC])
            tree.find('Material').find('material').find('G31').text = str(C[5+nC])
            tree.find('Material').find('material').find('v12').text = str(C[6+nC])
            tree.find('Material').find('material').find('v23').text = str(C[7+nC])
            tree.find('Material').find('material').find('v31').text = str(C[8+nC])
            tree.find('Material').find('material').find('c').text = str(C[9+nC])
            tree.find('Material').find('material').find('k').text = str(C[10+nC])
            nC +=11
     
        
    # Define shape and value of pressure profile defined between 0 and 1. The pressure magnitude is defined separately 
    if ProfileChoice == 'Triangle':
        # A triangular profile that starts and ends with 0 and time of the peak is defined by parameter estimate
        Peak = C[0+nC]
        PressureInterp = np.zeros(nF)
        Line_A = (1/Peak)*TimeInterp
        Line_B = (-1/(1 - Peak))*TimeInterp + 1/(1-Peak)
        for i in range(nF):
            if TimeInterp[i] <=Peak:
                PressureInterp[i] = Line_A[i]
            else:
                PressureInterp[i] = Line_B[i]
        for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
            Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx])
    elif ProfileChoice == 'Step':
        # a box profile shape that is 1 in a range and 0 elsewhere
        PressureInterp = np.zeros(nF)
        for i in range(nF):
            if (TimeInterp[i] >=C[0+nC]) and (TimeInterp[i] <= C[1+nC]):
                PressureInterp[i] = 1.0
        for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
            Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx])
    elif ProfileChoice == 'SmoothStep':
        # A box shape with smoothers transitions with sigmoidal function
        PressureInterp = np.zeros(nF)
        for i in range(nF):
            PressureInterp[i] = (1/(1+e(-C[2+nC]*(TimeInterp[i]-C[0+nC]))))*(1-1/(1+e(-C[3+nC]*(TimeInterp[i]-C[1+nC]))))
        for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
            Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx])
    elif ProfileChoice == 'Fourier':
        # A profile shape defined by fourier series
        PressureInterp = np.zeros(nF)
        x=TimeInterp
        for i in range(nF):
            PressureInterp[i] = C[0+nC]*sin(pi*x[i])+C[1+nC]*sin(2*pi*x[i])+C[2+nC]*sin(4*pi*x[i])+C[3+nC]*sin(8*pi*x[i])
        for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
            Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx]) 
    elif ProfileChoice == 'Fitted':
        # A fourier series where every point is chosen by paramter estimates
        PressureInterp = np.zeros(nF)
        for i in range(nF-2):
            PressureInterp[i+1] = C[i+nC]
        for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
            Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx]) 
    elif ProfileChoice == 'Windkess':
        #Profile created with windkessel model
        qi_mag = C[0+nC]
        Rp = C[1+nC]
        Rd = C[2+nC]
        Cp  = C[3+nC]
        t = np.linspace(0, 100, 5001)
        qi = []
        CFt = CF/nF
        for ti in t:
            if ti%1>0 and ti%1<CFt:
                qi.append(qi_mag*sin(pi*(ti%1)/CFt))
                #qi.append(qi_mag)
            else:
                qi.append(0)
        qi = np.array(qi)
        qi_func = interp1d(t,qi, fill_value="extrapolate")
        def func(y,t):
            return (qi_func(t)-y)/Cp/Rd
        qo = odeint(func,0,t).flatten()
        # Get final profile once profile shape has settled
        last = (t>=99)*(t<100)
        Pressure = Rp*qi+Rd*qo
        PressureInterp = np.interp(np.linspace(99,100,nF),t[last],Pressure[last])
        PressureInterp = np.subtract(PressureInterp,PressureInterp[0])
        # PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
        for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
            Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx]) 
    elif ProfileChoice[0:3] == 'Set':
        #Choose to set the profile without paramter choices
        if ProfileChoice[3:] == 'Windkess':
            qi_mag = 0.01
            Rp     = 1.
            Rd     = 10.
            Cp     = 1.
            t = np.linspace(0, 100, 5001)
            qi = []
            CFt = CF/nF
            for ti in t:
                if ti%1>0 and ti%1<CFt:
                    qi.append(qi_mag*sin(pi*(ti%1)/CFt))
                    #qi.append(qi_mag)
                else:
                    qi.append(0)
            qi = np.array(qi)
            qi_func = interp1d(t,qi, fill_value="extrapolate")
            def func(y,t):
                return (qi_func(t)-y)/Cp/Rd
            qo = odeint(func,0,t).flatten()
            last = (t>=99)*(t<100)
            Pressure = Rp*qi+Rd*qo
            PressureInterp = np.interp(np.linspace(99,100,nF),t[last],Pressure[last])
            PressureInterp = np.subtract(PressureInterp,PressureInterp[0])
            # PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
            for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
                Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx]) 
       
    # Rewrite .feb file with paramter updates
    tree.write(DataDir + '.feb',xml_declaration=True,encoding="ISO-8859-1")
    
    #Run updated .feb file to create new .xplt file
    os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i '+ './FEB_Files/' + DataDir+'.feb >/dev/null 2>&1')
    # os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i '+ './FEB_Files/' + DataDir+'.feb')
    
    # Define file of newly created .xplt file
    XPLTfilename = './FEB_Files/' + DataDir + '.xplt'
    
    # Get data from .xplt 
    feb,file_size,nStates, mesh = GetFEB(XPLTfilename)
    
    nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
    
    #Collect the data, with all states
    displacement = GetData(feb,'displacement',nStates,nVar)
      
    #Gather stresses and displacements Note this is dependent on the Var Type
    Disp_X,Disp_Y,Disp_Z = np.zeros((nStates,nNodes)),np.zeros((nStates,nNodes)),np.zeros((nStates,nNodes))
    Disp_FEB, Res  = np.zeros((nStates,nNodes,3)), np.zeros((nStates,nNodes,3))  
    for i in range(nStates): 
        for j in range(nNodes):
            Disp_X[i,j]  = displacement[i][j*3+1]
            Disp_Y[i,j]  = displacement[i][j*3+2]
            Disp_Z[i,j]  = displacement[i][j*3+3]
        Disp_FEB[i,:,0]  = Disp_X[i,:]  
        Disp_FEB[i,:,1]  = Disp_Y[i,:]
        Disp_FEB[i,:,2]  = Disp_Z[i,:]
    
    # Interpolate data if the number of frames don't match
    if nF !=nStates:
        Disp_FEB_Interp = np.zeros(Disp_Wall.shape)
        for i in range(nNodes):
            for j in range(3):
                Disp_FEB_Interp[:,i,j] = np.interp(TimeInterp,StateTimes,Disp_FEB[:,i,j])
    else:
        Disp_FEB_Interp = Disp_FEB
            
    # Define method of calculating residuals
    if ResChoice == 'P2P':
        # Explicit point to point distances
        Res = np.subtract(Disp_Wall,Disp_FEB_Interp)
    elif ResChoice == 'CentreLine':
        # Define using difference in the distances from centre line
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
        Res = np.zeros((nStates,nNodes))
        for i in range(nStates):
            for j in range(nNodes):
                VTK_Pt = np.add(Pts[j],Disp_Wall[i,j])
                FEB_Pt = np.add(Pts[j],Disp_FEB[i,j])
                VTK2L = (np.linalg.norm(np.cross(np.subtract(VTK_Pt,VAJ_Cent),np.subtract(VTK_Pt,STJ_Cent))))/np.linalg.norm(np.subtract(STJ_Cent,VAJ_Cent))
                FEB2L = (np.linalg.norm(np.cross(np.subtract(FEB_Pt,VAJ_Cent),np.subtract(FEB_Pt,STJ_Cent))))/np.linalg.norm(np.subtract(STJ_Cent,VAJ_Cent))
            
                Res[i,j] = np.absolute(FEB2L - VTK2L)
    elif ResChoice == 'CellPlane':
        #Define points from VTK remeshed file
        VTK_Pts = np.zeros((nF,nNodes,3))
        for i in range(nF):
            for j in range(nNodes):
                VTK_Pts[i,j] = np.add(Pts[j],Disp_Wall[i,j])
         
        # Define pts of FEB simulation       
        FEB_Pts = np.zeros((nF,nNodes,3))
        for i in range(nF):
            for j in range(nNodes):
                FEB_Pts[i,j] = np.add(Pts[j],Disp_FEB_Interp[i,j])
           
        #Define empty arrays
        Res = np.zeros((nF,nNodes))
        CellPtSum    = [[] for _ in range(nNodes)]
        CellDistId   = [[] for _ in range(nNodes)]
        CellTotalSum = [[] for _ in range(nNodes)]
        
        # loop through  frames and points
        for i in range(nF):
            for j in range(nNodes): 
                # get closest points to using neighbouring cells, only for first frame for efficiency
                if i ==0:
                    #for each point get point ids of associated cells
                    for k in range(nCls):
                        if j in CellIds[k]:
                            CellPtSum[j].append([])
                            CellDistId[j].append(CellIds[k].astype(int))
                            
                    #For these ids of each get distance to the point with id=j
                    for k in range(len(CellDistId[j])):
                        for l in np.array(CellDistId[j][k]).astype(int):
                            CellPtSum[j][k].append(np.linalg.norm(np.subtract(VTK_Pts[i,j],FEB_Pts[i,l])))
                        CellTotalSum[j].append(np.sum(CellPtSum[j][k]))
                    
                
                MinCl = np.argmin(CellTotalSum[j])
                
                # Get closest three points
                Pt1, Pt2, Pt3 = FEB_Pts[i,CellDistId[j][MinCl][np.argpartition(CellPtSum[j][MinCl],3)[0]]], FEB_Pts[i,CellDistId[j][MinCl][np.argpartition(CellPtSum[j][MinCl],3)[1]]], FEB_Pts[i,CellDistId[j][MinCl][np.argpartition(CellPtSum[j][MinCl],3)[2]]]
                
                # Get normal of plane defined by the 3 closestpoints
                planeNorm = np.cross(np.subtract(Pt1,Pt2),np.subtract(Pt1,Pt3))
                planePt = Pt1
                
                # Get norm of point j
                PtNorm = Norm[i,j]
                Pt     = VTK_Pts[i,j]
                
                # Get point of intersection
                ndotu = planeNorm.dot(PtNorm)
                w = Pt - planePt
                si = -planeNorm.dot(w) / ndotu
                Intersect = np.add(np.add(w,np.multiply(si,PtNorm)),planePt)
                
                # Get distance between point and intersection
                Res[i,j] =  np.linalg.norm(np.subtract(Pt,Intersect))
            
    return Res.flatten()
    
def RunLS(DataDir,d,FListOrdered,FId,ref,CF,PressureChoice,ModelParChoice,ProfileChoice,RunLSChoice,ResChoice,ModelChoice):
    '''
    Function to run the least squares 
    
    Inputs:
    DataDir - Name of directory, e.g. 'tav02'
    FListOrdered - reordered list of files
    FId - Array of the ordered frame Ids, starting the reference frame
    ref - Reference remeshed file name
    CF - Id of the valve closing frame 
    PressureChoice - Choice to include pressure magnitude in least squares optimisation
    ModelParChoice - Choice to include model parameters in least squares optimisation
    ProfileChoice - Choice of pressure profile shape function
    RunLSChoice - Choice to run least quares optimisation
    ResChoice -  Choice of residual calculation method 
    ModelChoice - Choice of constitutive model
    '''
    # Read the source file
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(ref)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    polydata = reader.GetOutput()
    nPts     = polydata.GetNumberOfPoints()
    nCls     = polydata.GetNumberOfCells()
    
    STJ = vtk_to_numpy(polydata.GetPointData().GetArray('STJ'))
    VAJ = vtk_to_numpy(polydata.GetPointData().GetArray('VAJ'))
    Pts = vtk_to_numpy(polydata.GetPoints().GetData())
    
    #Get Number of frame
    nF = len(FId)
    
    #Create empty arrays
    STJ_Id = []
    VAJ_Id = []
    nSTJ = 0
    nVAJ = 0
    #Get STJ and VAJ ids and number of each type from remeshed files
    for i in range(nPts):
        if STJ[i] !=0.0:
            STJ_Id.append(i)
            nSTJ+=1
        if VAJ[i] !=0.0:
            VAJ_Id.append(i)
            nVAJ+=1
    
    CellIds = np.zeros((nCls,4))
    for i in range(nCls):
        for j in range(4):
            CellIds[i,j] = (polydata.GetCell(i).GetPointIds().GetId(j))
                
    Disp_Wall_STJ, Disp_Wall_VAJ, Circ, Long, Norm, Disp_Wall = np.zeros((nF,nSTJ,3)), np.zeros((nF,nVAJ,3)), np.zeros((nF,nPts,3)), np.zeros((nF,nPts,3)), np.zeros((nF,nPts,3)), np.zeros((nF,nPts,3))
    
    # Analyse each frame
    for X,Fname in enumerate(FListOrdered):
        # Read the source file
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(Fname)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        polydata = reader.GetOutput()
        
        STJ       = vtk_to_numpy(polydata.GetPointData().GetArray('STJ'))
        VAJ       = vtk_to_numpy(polydata.GetPointData().GetArray('VAJ'))
        Disp_Wall[X] = vtk_to_numpy(polydata.GetPointData().GetArray('Displacement_Wall'))
        Circ[X]      = vtk_to_numpy(polydata.GetPointData().GetArray('Circumferential'))
        Long[X]      = vtk_to_numpy(polydata.GetPointData().GetArray('Longtudinal'))
        Norm[X]      = vtk_to_numpy(polydata.GetPointData().GetArray('Normal'))
        
        Circ_Cls = np.zeros((nCls,3))
        
        for i in range(nCls):
            c = polydata.GetCell(i)
            Circ_Cls[i,0] = (Circ[0,c.GetPointId(0),0]+Circ[0,c.GetPointId(1),0]+Circ[0,c.GetPointId(2),0]+Circ[0,c.GetPointId(3),0])/4
            Circ_Cls[i,1] = (Circ[0,c.GetPointId(0),1]+Circ[0,c.GetPointId(1),1]+Circ[0,c.GetPointId(2),1]+Circ[0,c.GetPointId(3),1])/4
            Circ_Cls[i,2] = (Circ[0,c.GetPointId(0),2]+Circ[0,c.GetPointId(1),2]+Circ[0,c.GetPointId(2),2]+Circ[0,c.GetPointId(3),2])/4
            
        nPts = polydata.GetNumberOfPoints()
        j=0
        k=0
        for i in range(nPts):
            if STJ[i] !=0.0:
                Disp_Wall_STJ[X][j] = Disp_Wall[X,i]
                j+=1
            if VAJ[i] !=0.0:
                Disp_Wall_VAJ[X][k] = Disp_Wall[X,i]
                k+=1
    
    #Add choices to parameter array C and parameter ranges B_Min and B_Max
    if PressureChoice:
        C = [-0.005332]
        B_Min = [-1]
        B_Max = [0]
    else:
        C = []
        B_Min = []
        B_Max = []
        
    if ModelParChoice:
        if ModelChoice == 'MR':
            C = np.concatenate((C,[0.0,10,10,10]))
            B_Min = np.concatenate((B_Min,[0,0,0,0]))
            B_Max = np.concatenate((B_Max,[100,1000,1000,1000]))
        if ModelChoice == 'tiMR':
            C = np.concatenate((C,[10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0]))
            B_Min = np.concatenate((B_Min,[0,0,0,0,0,0,0,0]))
            B_Max = np.concatenate((B_Max,[1000,1000,1000,1000,1000,1000,1000,1000]))
        elif ModelChoice == 'Ogden':
            C = np.concatenate((C,[10,10,10,10,10,10,10,10,10,10,10,10,10]))
            B_Min = np.concatenate((B_Min,[0,0,0,0,0,0,0,0,0,0,0,0,0]))
            B_Max = np.concatenate((B_Max,[100,100,100,100,100,100,100,100,100,100,100,100,100]))
        elif ModelChoice == 'Fung':
            C = np.concatenate((C,[1,2,1.5,0.5,0.1,1,1.2,2,1.5,0.5,0.1]))
            B_Min = np.concatenate((B_Min,[0,0,0,0,0,0,0,0,0,0,0]))
            B_Max = np.concatenate((B_Max,[10,10,10,10,10,10,10,10,10,10,10]))
    else:
        C = np.concatenate((C,[]))
        B_Min = np.concatenate((B_Min,[]))
        B_Max = np.concatenate((B_Max,[]))
        
    if ProfileChoice == 'Triangle':
        C = np.concatenate((C,[1/2]))
        B_Min = np.concatenate((B_Min,[0]))
        B_Max = np.concatenate((B_Max,[1]))
    elif ProfileChoice == 'Step':
        C = np.concatenate((C,[2/nF,5/nF]))
        B_Min = np.concatenate((B_Min,[0,0]))
        B_Max = np.concatenate((B_Max,[1,1]))
    elif ProfileChoice == 'SmoothStep':
        C = np.concatenate((C,[2/nF, 5/nF, 50, 50]))
        B_Min = np.concatenate((B_Min,[0,0,0,0]))
        B_Max = np.concatenate((B_Max,[1,1,1000,1000]))
    elif ProfileChoice == 'Bio':
        C = np.concatenate((C,[]))
        B_Min = np.concatenate((B_Min,[]))
        B_Max = np.concatenate((B_Max,[]))
    elif ProfileChoice == 'Fourier':
        C = np.concatenate((C,[3.0, 1.0, 0.0, 0.05]))
        B_Min = np.concatenate((B_Min,[-100,-100,-100,-100]))
        B_Max = np.concatenate((B_Max,[100,100,100,100]))
    elif ProfileChoice == 'Fitted':
        PressureInterpStan = 0.0*np.ones(nF-2)
        C = np.concatenate((C,PressureInterpStan))
        B_Min = np.concatenate((B_Min,np.zeros(nF-2)))
        B_Max = np.concatenate((B_Max,np.ones(nF-2)))
    elif ProfileChoice == 'Windkess':
        C = np.concatenate((C,[0.01,1, 1 ,  0.5]))
        B_Min = np.concatenate((B_Min,[0,0,0,0]))
        B_Max = np.concatenate((B_Max,[1,100,100,100]))
    
    #Create .feb file of VTK remeshed case
    VTK2Feb_Func('./FEB_Files/' + DataDir,ref,nF,nCls,Disp_Wall_STJ,Disp_Wall_VAJ,STJ_Id,VAJ_Id,Circ_Cls,ProfileChoice,ModelChoice,CF)
    
    #Choose to run Least Squares optimisation or just run febio simulation
    if RunLSChoice:
        Out = least_squares(GetRes,C,bounds = [B_Min,B_Max],jac = '3-point', verbose=2,args=(DataDir,Pts,Disp_Wall,Norm,Circ,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,ModelParChoice,ProfileChoice,ResChoice,ModelChoice))
        Cs = Out.x
    else:
        os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i '+ './FEB_Files/' + DataDir+'.feb')
        Cs = C
     
    return Cs, Pts, Disp_Wall, Norm, CellIds, nCls, STJ_Id, VAJ_Id, FId, nF
    
if __name__=='__main__':
    
    #Get location of example script
    d = './Strains/medial_meshes/tav02'
    DataDir = 'tav02'
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
    #Choose if data needs remeshed
    PressureChoice = False           # Choose to vary pressure magnitude
    ModelParChoice = True            # Choose to vary modelparameters
    RunLSChoice    = False           # Choose to run least Squares (or default/initial guess)
    ProfileChoice  = ['Bio']         #['Triangle','Step','SmoothStep','Bio','Fourier','Fitted'] # Choose profile shapes, options are: 'Triangle','Step','SmoothStep','Bio', 'Fourier','Fitted'
    ResChoice      = ['CellPlane']   # Choose type of residual calculation method: 'P2P', 'CentreLine', 'CellPlane'
    ModelChoice    = ['MR']        #Choose model from 'MR','tiMR','Ogden' and 'Fung'
    
    #Create empty array for params
    Params = []
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
                    
                if MC:
                    if ModelChoice == 'MR':
                        Model_density = Out[0+nC]
                        Model_C1 = Out[1+nC]
                        Model_C2 = Out[2+nC]
                        Model_k  = Out[3+nC]
                        print('Model C1: ',Model_C1)
                        print('Model C2: ',Model_C2)
                        print('Model k:  ',Model_k)
                    if ModelChoice == 'tiMR':
                        Model_density = Out[0+nC]
                        Model_C1      = Out[1+nC]
                        Model_C2      = Out[2+nC]
                        Model_C3      = Out[3+nC]
                        Model_C4      = Out[4+nC]
                        Model_C5      = Out[5+nC]
                        Model_lam_max = Out[6+nC]
                        Model_k       = Out[7+nC]
                        print('Model density: ',Model_density)
                        print('Model C1: ',Model_C1)
                        print('Model C2: ',Model_C2)
                        print('Model C3: ',Model_C3)
                        print('Model C4: ',Model_C4)
                        print('Model C5: ',Model_C5)
                        print('Model lam_max:  ',Model_lam_max)
                        print('Model k:  ',Model_k)
                    elif ModelChoice == 'Ogden':
                        Model_bulk = Out[0+nC]
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
                    elif ModelChoice == 'Fung':
                        Model_E1   = Out[0+nC]
                        Model_E2   = Out[1+nC]
                        Model_E3   = Out[2+nC]
                        Model_G12  = Out[3+nC]
                        Model_G23  = Out[4+nC]
                        Model_G13  = Out[5+nC]
                        Model_v12  = Out[6+nC]
                        Model_v23  = Out[7+nC]
                        Model_v13  = Out[8+nC]
                        Model_c    = Out[9+nC]
                        Model_k    = Out[10+nC]
                        
                    nC+=3
                else:
                    print('Model C1: ',1)
                    print('Model C2: ',0)
                    print('Model k: ',10)
           
                if ModelChoice:
                    line = ''
                else:
                    line = '--'     
        
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
                elif PC == 'Fitted':
                    col = 'c'
                elif PC == 'Windkess':
                    col = 'y'
                    
                style = col +line
                
                XPLTfilename = DataDir + '.xplt'
                feb,file_size,nStates, mesh = GetFEB(XPLTfilename)
                nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
        
                TimeData = (np.genfromtxt('TimeProfile.csv', delimiter=','))
                PressureData = np.genfromtxt('PressureProfile.csv', delimiter=',')
                TimeInterp = np.linspace(0,1,1001)
                PressureInterp = np.interp(TimeInterp,TimeData,PressureData)
                PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
                
                Time = np.zeros(nStates)
                Pressure = np.zeros(nStates)
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
        
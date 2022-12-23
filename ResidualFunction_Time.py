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
from math import cos, sin, pi
from math import exp as e
from itertools import zip_longest

def SaveFiles(DataDir,FList,RemeshedFile,Pts,Disp_Wall,Norm,Circ_Cls,Long_Cls,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,ProfileChoice,RunLSChoice,ResChoice,ModelChoice,FibChoice,Out):
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
    
    
    #Order list to put reference frame first
    FListOrdered, FId, refN = OrderList(FList, nF-1, RemeshedFile )
    FListOrdered.append(FListOrdered[0])
    FId = np.concatenate([FId,[FId[0]]],axis=0)
    
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
        
        fname = './NewFiles/'+ DataDir+'/'+ModelChoice+'/'+RunDir + '/' + ProfileChoice+ '/' + Conditions+'/'+DataDir + '_' + str(fid+1)+'.vtk'
        print(fname)
        directory = os.path.dirname(fname)
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        writer = vtk.vtkDataSetWriter()
        writer.SetFileName(fname)
        writer.SetInputData(polydata)
        writer.Write()
        
        feb_name = './FEB_Files/'+ DataDir+'.feb'
        
        os.system('cp '+ feb_name + ' ./NewFiles/'+ DataDir+'/'+ModelChoice+'/'+RunDir + '/' + ProfileChoice+ '/' + Conditions+'/')
        
    return nStates, Residual
    
def OrderList(FList,N,ref):
    '''
    Reorders the list so that the first file is the reference file, and returns new list of filenames and their IDs 
    Keyword arguments:
    FList -- list of filenames
    N -- number of files
    ref -- reference file name
    For example: for a list [file1.vtk, file2.vtk, file3.vtk, file4.vtk, file7.vtk, file8.vtk] and ref = file3.vtk
    it will return [file3.vtk file4.vtk file7.vtk file8.vtk file1.vtk file2.vtk], [3 4 7 8 1 2], and 3
    '''
    # Order filenames so that reference frame goes first
    Fno = np.zeros(N)
    FId = np.zeros(N)

    FListOrdered = [None]*N
    common = os.path.commonprefix(FList)
    for i, Fname in enumerate(FList):
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
        for Fname in FList:
            X = Fname.replace(common,'')
            X = X.replace('.vtk','')
            X = np.fromstring(X, dtype=int, sep=' ')
            if X ==F:
                FListOrdered[i] = Fname

    return FListOrdered, FId, refN

def GetPressure(C,nC,nF,CF,ProfileChoice,TimeInterp):
    '''
    Function to calculate the pressure profile, for a given choice
    
    Inputs:
    C - array of parameters used
    nC - Count of parameter array
    nF - Number of Frames
    CF - Id of the valve closing frame 
    ProfileChoice - Choice of pressure profile shape function
    '''
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
    elif ProfileChoice == 'Step':
        # a box profile shape that is 1 in a range and 0 elsewhere
        PressureInterp = np.zeros(nF)
        for i in range(nF):
            if (TimeInterp[i] >=C[0+nC]) and (TimeInterp[i] <= C[1+nC]):
                PressureInterp[i] = 1.0
    elif ProfileChoice == 'SmoothStep':
        # A box shape with smoothers transitions with sigmoidal function
        PressureInterp = np.zeros(nF)
        for i in range(nF):
            PressureInterp[i] = (1/(1+e(-C[2+nC]*(TimeInterp[i]-C[0+nC]))))*(1-1/(1+e(-C[3+nC]*(TimeInterp[i]-C[1+nC]))))
    elif ProfileChoice == 'Virtual':
        TimeData = (np.genfromtxt('TimeProfile.csv', delimiter=','))
        PressureData = np.genfromtxt('PressureProfile.csv', delimiter=',')
        PressureInterp = np.interp(TimeInterp,TimeData,PressureData)
        PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
    elif ProfileChoice == 'Fourier':
        # A profile shape defined by fourier series
        PressureInterp = np.zeros(nF)
        x=TimeInterp
        for i in range(nF):
            PressureInterp[i] = C[0+nC]*sin(pi*x[i])+C[1+nC]*sin(2*pi*x[i])+C[2+nC]*sin(4*pi*x[i])+C[3+nC]*sin(8*pi*x[i])
    elif ProfileChoice == 'Fitted':
        # A fourier series where every point is chosen by paramter estimates
        PressureInterp = np.zeros(nF)
        for i in range(nF-2):
            PressureInterp[i+1] = C[i+nC]
    elif ProfileChoice == 'Windkessel':
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
    elif ProfileChoice[0:3] == 'Set':
        if ProfileChoice[3:] == 'Triangle':
            Peak = 0.5
            PressureInterp = np.zeros(nF)
            Line_A = (1/Peak)*TimeInterp
            Line_B = (-1/(1 - Peak))*TimeInterp + 1/(1-Peak)
            for i in range(nF):
                if TimeInterp[i] <=Peak:
                    PressureInterp[i] = Line_A[i]
                else:
                    PressureInterp[i] = Line_B[i]
            PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
        elif ProfileChoice[3:] == 'Step':
            PressureInterp = np.zeros(nF)
            for i in range(nF):
                if (TimeInterp[i] >=0.2) and (TimeInterp[i] <= 0.5):
                    PressureInterp[i] = 1.0
            PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
        elif ProfileChoice[3:] == 'SmoothStep':
            PressureInterp = np.zeros(nF)
            for i in range(nF):
                PressureInterp[i] = (1/(1+e(-50*(TimeInterp[i]-0.2))))*(1-1/(1+e(-50*(TimeInterp[i]-0.5))))
            PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
        elif ProfileChoice[3:] == 'Virtual':
            TimeData = (np.genfromtxt('TimeProfile.csv', delimiter=','))
            PressureData = np.genfromtxt('PressureProfile.csv', delimiter=',')
            PressureInterp = np.interp(TimeInterp,TimeData,PressureData)
            PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
        elif ProfileChoice[3:] == 'Fourier':
            PressureInterp = np.zeros(nF)
            x=TimeInterp
            for i in range(nF):
                PressureInterp[i] = 3.0*sin(pi*x[i])+1.0*sin(2*pi*x[i])+0.0*sin(4*pi*x[i])+0.05*sin(8*pi*x[i])
            PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
        elif ProfileChoice[3:] == 'Fitted':
            PressureInterp = 0.0*np.ones(nF)
        elif ProfileChoice[3:] == 'Windkessel':
            qi_mag = 11.0
            Rp     = 1.4
            Rd     = 14
            Cp     = 0.004
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
        
    PressureInterp = np.subtract(PressureInterp,min(PressureInterp))
    PressureInterp = np.divide(PressureInterp,max(PressureInterp))
        
    return PressureInterp
    
def GetRes(C,DataDir,Pts,Disp_Wall,Norm,Circ,Long,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,ProfileChoice,ResChoice,ModelChoice,FibChoice):
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
        # if ModelChoice == 'tiMR':
        tree.find('Loads').find('surface_load').find('pressure').text = str(-C[0])
        # else:
        #     tree.find('Loads').find('surface_load').find('pressure').text = str(C[0])
        nC += 1
    
    # Replace material properties with parameter estimates
    if ModelChoice[0:3] != 'Set':
        if ModelChoice == 'MR':
            tree.find('Material').find('material').find('density').text = str(C[0+nC])
            tree.find('Material').find('material').find('c1').text      = str(C[1+nC])
            tree.find('Material').find('material').find('c2').text      = str(C[2+nC])
            tree.find('Material').find('material').find('k').text       = str(C[3+nC])
            nC += 4
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
            tree.find('Material').find('material').find('density').text  = str(C[0+nC])
            tree.find('Material').find('material').find('k').text  = str(C[1+nC])
            tree.find('Material').find('material').find('c1').text = str(C[2+nC])
            tree.find('Material').find('material').find('c2').text = str(C[3+nC])
            tree.find('Material').find('material').find('c3').text = str(C[4+nC])
            tree.find('Material').find('material').find('c4').text = str(C[5+nC])
            tree.find('Material').find('material').find('c5').text = str(C[6+nC])
            tree.find('Material').find('material').find('c6').text = str(C[7+nC])
            tree.find('Material').find('material').find('m1').text = str(C[8+nC])
            tree.find('Material').find('material').find('m2').text = str(C[9+nC])
            tree.find('Material').find('material').find('m3').text = str(C[10+nC])
            tree.find('Material').find('material').find('m4').text = str(C[11+nC])
            tree.find('Material').find('material').find('m5').text = str(C[12+nC])
            tree.find('Material').find('material').find('m6').text = str(C[13+nC])
            nC +=14
        elif ModelChoice == 'Fung':
            tree.find('Material').find('material').find('density').text = str(C[0+nC])
            tree.find('Material').find('material').find('E1').text      = str(C[1+nC])
            tree.find('Material').find('material').find('E2').text      = str(C[2+nC])
            tree.find('Material').find('material').find('E3').text      = str(C[3+nC])
            tree.find('Material').find('material').find('G12').text     = str(C[4+nC])
            tree.find('Material').find('material').find('G23').text     = str(C[5+nC])
            tree.find('Material').find('material').find('G31').text     = str(C[6+nC])
            tree.find('Material').find('material').find('v12').text     = str(C[7+nC])
            tree.find('Material').find('material').find('v23').text     = str(C[8+nC])
            tree.find('Material').find('material').find('v31').text     = str(C[9+nC])
            tree.find('Material').find('material').find('c').text       = str(C[10+nC])
            tree.find('Material').find('material').find('k').text       = str(C[11+nC])
            nC +=12
        elif ModelChoice == 'HGO':
            TreeBranches = tree.find('Material').findall('material')
            for i in range(nCls):
                TreeBranches[i].find('c').text     = str(C[0+nC])
                TreeBranches[i].find('k1').text    = str(C[1+nC])
                TreeBranches[i].find('k2').text    = str(C[2+nC])
                TreeBranches[i].find('gamma').text = str(C[3+nC])
                # TreeBranches[i].find('kappa').text = str(C[4+nC])
                # TreeBranches[i].find('k').text     = str(C[5+nC])
            nC +=6
     
        
        if FibChoice and ModelChoice == 'tiMR':
            theta = C[0+nC]
            for i in range(nCls):
                C_tick = np.cross(np.cross(Circ[i],Long[i]),Circ[i])
                New = np.dot(cos(theta),Circ[i]) + np.dot(sin(theta),C_tick)
                TreeBranches[i].find('fiber').text = str(New[0])+str(',')+str(New[1])+str(',')+str(New[2])
            nC += 1
            
    #Get Pressure Profile
    if ProfileChoice[0:3] != 'Set':
        PressureInterp = GetPressure(C,nC,nF,CF,ProfileChoice,TimeInterp)
        for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
                Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx])
               
    # Rewrite .feb file with paramter updates
    tree.write('./FEB_Files/' + DataDir + '.feb',xml_declaration=True,encoding="ISO-8859-1")
    
    #Run updated .feb file to create new .xplt file
    os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i '+ './FEB_Files/' + DataDir+'.feb >/dev/null 2>&1')
    # os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i '+ './FEB_Files/' + DataDir+'.feb')
    
    # Define file of newly created .xplt file
    XPLTfilename = './FEB_Files/' + DataDir + '.xplt'
    
    if ModelChoice in ['MR']:
        nDoms = 1
    else:
        nDoms = nCls
    
    # Get data from .xplt 
    feb,file_size,nStates, mesh = GetFEB(XPLTfilename,nDoms,False)
    
    nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
    
    #Collect the data, with all states
    displacement = GetData(feb,'displacement',nStates)
      
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
                if abs(ndotu) < 1e-6:
                    Residual[i,j] =  0.0
                    raise RuntimeError("no intersection or line is within plane")
                else:
                    w = Pt - planePt
                    si = -planeNorm.dot(w) / ndotu
                    Intersect = np.add(np.add(w,np.multiply(si,PtNorm)),planePt)
                    
                    # Get distance between point and intersection
                    Res[i,j] =  np.linalg.norm(np.subtract(Pt,Intersect))
    
    return Res.flatten()
    
def RunLS(DataDir,FListOrdered,FId,ref,CF,C,PressureChoice,ProfileChoice,RunLSChoice,ResChoice,ModelChoice,FibChoice):
    '''
    Function to run the least squares 
    
    Inputs:
    DataDir - Name of directory, e.g. 'tav02'
    FListOrdered - reordered list of files
    FId - Array of the ordered frame Ids, starting the reference frame
    ref - Reference remeshed file name
    CF - Id of the valve closing frame 
    C - Initial Parameter estimation array
    PressureChoice - Choice to include pressure magnitude in least squares optimisation
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
        Circ[X]      = vtk_to_numpy(polydata.GetPointData().GetArray('Circumferential_Pt'))
        Long[X]      = vtk_to_numpy(polydata.GetPointData().GetArray('Longtudinal_Pt'))
        Norm[X]      = vtk_to_numpy(polydata.GetPointData().GetArray('Normal_Pt'))
        
        Circ_Cls = np.zeros((nCls,3))
        Long_Cls = np.zeros((nCls,3))
        
        for i in range(nCls):
            c = polydata.GetCell(i)
            Circ_Cls[i,0] = (Circ[0,c.GetPointId(0),0]+Circ[0,c.GetPointId(1),0]+Circ[0,c.GetPointId(2),0]+Circ[0,c.GetPointId(3),0])/4
            Circ_Cls[i,1] = (Circ[0,c.GetPointId(0),1]+Circ[0,c.GetPointId(1),1]+Circ[0,c.GetPointId(2),1]+Circ[0,c.GetPointId(3),1])/4
            Circ_Cls[i,2] = (Circ[0,c.GetPointId(0),2]+Circ[0,c.GetPointId(1),2]+Circ[0,c.GetPointId(2),2]+Circ[0,c.GetPointId(3),2])/4
            Long_Cls[i,0] = (Long[0,c.GetPointId(0),0]+Circ[0,c.GetPointId(1),0]+Long[0,c.GetPointId(2),0]+Long[0,c.GetPointId(3),0])/4
            Long_Cls[i,1] = (Long[0,c.GetPointId(0),1]+Circ[0,c.GetPointId(1),1]+Long[0,c.GetPointId(2),1]+Long[0,c.GetPointId(3),1])/4
            Long_Cls[i,2] = (Long[0,c.GetPointId(0),2]+Circ[0,c.GetPointId(1),2]+Long[0,c.GetPointId(2),2]+Long[0,c.GetPointId(3),2])/4
            
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
        B_Min = [0]
        B_Max = [5]
    else:
        B_Min = []
        B_Max = []
        
    if ModelChoice[0:3] != 'Set':
        if ModelChoice == 'MR':
            B_Min = np.concatenate((B_Min,[0,0,0,0]))
            B_Max = np.concatenate((B_Max,[100,1000,1000,2000]))
        if ModelChoice == 'tiMR':
            B_Min = np.concatenate((B_Min,[0,0,0,0,0,0,0,0]))
            B_Max = np.concatenate((B_Max,[1000,1000,1000,1000,1000,1000,1000,1000]))
        elif ModelChoice == 'Ogden':
            B_Min = np.concatenate((B_Min,[0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
            B_Max = np.concatenate((B_Max,[100,100,100,100,100,100,100,100,100,100,100,100,100,100]))
        elif ModelChoice == 'Fung':
            B_Min = np.concatenate((B_Min,[1,0,0,0,0,0,0,0,0,0,0,0]))
            B_Max = np.concatenate((B_Max,[10,10,10,10,10,10,10,10,10,10,10,10]))
        elif ModelChoice == 'HGO':
            B_Min = np.concatenate((B_Min,[ 0,    0,  0,  0]))
            B_Max = np.concatenate((B_Max,[100, 1000, 100, 50]))
    else:
        B_Min = np.concatenate((B_Min,[]))
        B_Max = np.concatenate((B_Max,[]))
    
    if FibChoice and ModelChoice == 'tiMR':
        B_Min = np.concatenate((B_Min,[-pi/4]))
        B_Max = np.concatenate((B_Max,[pi/4]))
    else:
        B_Min = np.concatenate((B_Min,[]))
        B_Max = np.concatenate((B_Max,[]))
        
    if ProfileChoice == 'Triangle':
        B_Min = np.concatenate((B_Min,[0]))
        B_Max = np.concatenate((B_Max,[1]))
    elif ProfileChoice == 'Step':
        B_Min = np.concatenate((B_Min,[0,0]))
        B_Max = np.concatenate((B_Max,[1,1]))
    elif ProfileChoice == 'SmoothStep':
        B_Min = np.concatenate((B_Min,[0,0,0,0]))
        B_Max = np.concatenate((B_Max,[1,1,1000,1000]))
    elif ProfileChoice == 'Virtual':
        B_Min = np.concatenate((B_Min,[]))
        B_Max = np.concatenate((B_Max,[]))
    elif ProfileChoice == 'Fourier':
        B_Min = np.concatenate((B_Min,[-100,-100,-100,-100]))
        B_Max = np.concatenate((B_Max,[100,100,100,100]))
    elif ProfileChoice == 'Fitted':
        B_Min = np.concatenate((B_Min,np.zeros(nF-2)))
        B_Max = np.concatenate((B_Max,np.ones(nF-2)))
    elif ProfileChoice == 'Windkessel':
        B_Min = np.concatenate((B_Min,[0,0,0,0]))
        B_Max = np.concatenate((B_Max,[100,100,100,100]))
    
    #Create .feb file of VTK remeshed case
    VTK2Feb_Func(DataDir,ref,C,nF,Disp_Wall_STJ,Disp_Wall_VAJ,STJ_Id,VAJ_Id,Circ_Cls,Long_Cls,Norm[0],PressureChoice,ProfileChoice,ModelChoice,FibChoice,CF)
    
    #Choose to run Least Squares optimisation or just run febio simulation
    if RunLSChoice:
        Out = least_squares(GetRes,C,bounds = [B_Min,B_Max],jac = '3-point', verbose=2,args=(DataDir,Pts,Disp_Wall,Norm,Circ_Cls,Long_Cls,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,ProfileChoice,ResChoice,ModelChoice,FibChoice))
        Cs = Out.x
    else:
        os.system('/Applications/FEBioStudio3.5/FEBioStudio.app/Contents/MacOS/febio3 -i ./FEB_Files/' + DataDir+'.feb')
        # os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i '+ './FEB_Files/' + DataDir+'.feb >/dev/null 2>&1')
        Cs = C
     
    return Cs, Pts, Disp_Wall, Norm, Circ_Cls, Long_Cls, CellIds, nCls, STJ_Id, VAJ_Id
    
if __name__=='__main__':
    
    #Get location of example script
    DataDir = 'tav02'
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
    FListOrdered.append(FListOrdered[0])
    FId = np.concatenate([FId,[FId[0]]],axis=0)
    nF+=1
    
    #Choose if data needs remeshed
    PressureChoice = False             # Choose to vary pressure magnitude
    RunLSChoice    = False             # Choose to run least Squares (or default/initial guess)
    FibChoice      = False             # Choose to vary fiber direction, as an angle from the circumferential direction
    ProfileChoice  = 'SetWindkessel'   # Choose profile shapes, options are: 'Triangle', 'Step', 'SmoothStep', 'Virtual', 'Fourier','Fitted', 'Windkessel'
    ResChoice      = 'P2P'             # Choose type of residual calculation method: 'P2P', 'CentreLine', 'CellPlane'
    ModelChoice    = 'HGO'             # Choose model from 'MR','tiMR','Ogden' and 'Fung',  'HGO'
    
    InitParams = [[],[],[],[]]
    
    # Choose initial estimate of pressure magnitude
    if PressureChoice:
        InitParams[0] = [1]
     
    # Choose initial model parameters
    if ModelChoice[0:3] != 'Set':
        if ModelChoice == 'MR':
            InitParams[1] = [1,10,10,10]
        if ModelChoice == 'tiMR':
            InitParams[1] = [10,14,0,2,60,600,100,10]
        elif ModelChoice == 'Ogden':
            InitParams[1] = [1,10,10,10,10,10,10,10,10,10,10,10,10,10]
        elif ModelChoice == 'Fung':
            InitParams[1] = [1,1,2,1.5,0.5,0.1,1,1.2,2,1.5,0.5,0.1]
        elif ModelChoice == 'HGO':     
            InitParams[1] = [100,25,90,25]#[ 5,  415,  5, 28]
      
    # Define initial fiber angle
    # Note if the chosen model is not 'tiMR' this parameter array will be made to be empty
    if FibChoice and ModelChoice =='tiMR':
        InitParams[2] = [0]

    #Choose initial paramters for the pressure profile
    if ProfileChoice == 'Triangle':
        InitParams[3] = [0.5]
    if ProfileChoice == 'Step':
        InitParams[3] = [0.2, 0.5]
    if ProfileChoice == 'SmoothStep':
        InitParams[3] = [0.2, 0.5,50,50]
    if ProfileChoice == 'Virtual':
        InitParams[3] = []
    if ProfileChoice == 'Fourier':
        InitParams[3] = [3.0, 1.0, 0.0, 0.05]
    if ProfileChoice == 'Fitted':
        InitParams[3] = np.zeros(nF-2)
    if ProfileChoice == 'Windkessel':
        InitParams[3] = [11,1.4,14,0.004]
    
    C = np.concatenate(InitParams, axis=0 )
    
    # If C is empty and RunLS is True, change to False
    if C.size == 0 and RunLS:
        print('Note: no parameters are being fitted, thus RunLSChoice is updated to be False')
        RunLSChoice == False
    
    #Create empty array for params
    Params = [['Condition'],['Choice'],['Model Parameter Name'],['Parameter Value'],['Pressure Parameter Name'],['Parameter Value'],['Model Time'],['Model Pressure'],['Interpolated Time'],['Interpolated Pressure']]

    #Run least squares script
    Out, Pts, Disp_Wall, Norm, Circ, Long, CellIds, nCls, STJ_Id, VAJ_Id = RunLS(DataDir,FListOrdered,FId,ref,CF,C,PressureChoice,ProfileChoice,RunLSChoice,ResChoice,ModelChoice,FibChoice)
    
    # Save new files
    nStates,Residual = SaveFiles(DataDir,FList,ref,Pts,Disp_Wall,Norm,Circ,Long,CellIds,nCls,STJ_Id,VAJ_Id,FId,nF,CF,PressureChoice,ProfileChoice,RunLSChoice,ResChoice,ModelChoice,FibChoice,Out)

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
        Params[5].append('N/A')
        print('Pressure Magnitude is Default')
        Conditions += 'PMagF_'
        
    if ModelChoice[0:3] != 'Set':
        Conditions += 'ModelT_'
        Params[0].append('Model_Parameters')
        Params[1].append('True')
        if ModelChoice == 'MR':
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
        elif ModelChoice == 'tiMR':
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
        elif ModelChoice == 'Ogden':
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
        elif ModelChoice == 'Fung':
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
        elif ModelChoice == 'HGO':
            Params[0].append('Model')
            Params[1].append('HGO')
            nP = 4
            ParamNames  = ['c','k1','k2','gamma']
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
        Params[1].append('FALSE')
        if ModelChoice == 'MR':
            Params[0].append('Model')
            Params[1].append('MR')
            ParamNames  = ['density','C1','C2','k']
            ParamValues = [1,10,10,10]
            for i in range(4):
                Params[2].append(ParamNames[i])
                Params[3].append(ParamValues[i])
        elif ModelChoice == 'tiMR':
            Params[0].append('Model')
            Params[1].append('tiMR')
            ParamNames  = ['density','C1','C2','C3','C4','C5','k','lam_max']
            ParamValues = [10,10,10,10,10,10,10,10]
            for i in range(8):
                Params[2].append(ParamNames[i])
                Params[3].append(ParamValues[i])
        elif ModelChoice == 'Ogden':
            Params[0].append('Model')
            Params[1].append('Ogden')
            ParamNames  = ['density','k','c1','c2','c3','c4','c5','c6','m1','m2','m3','m4','m5','m6']
            ParamValues = [1,10,10,10,10,10,10,10,10,10,10,10,10,10]
            for i in range(14):
                Params[2].append(ParamNames[i])
                Params[3].append(ParamValues[i])
        elif ModelChoice == 'Fung':
            Params[0].append('Model')
            Params[1].append('Fung')
            ParamNames  = ['density','E1','E2','E3','G12','G23','G31','v12','v23','v31','c','k']
            ParamValues = [1,1,2,1.5,0.5,0.1,1,1.2,2,1.5,0.5,0.1]
            for i in range(12):
                Params[2].append(ParamNames[i])
                Params[3].append(ParamValues[i])
        elif ModelChoice == 'HGO':
            Params[0].append('Model')
            Params[1].append('HGO')
            ParamNames  = ['c','k1','k2','gamma']
            ParamValues = [ 5,  415,  5, 28]
            for i in range(6):
                Params[2].append(ParamNames[i])
                Params[3].append(ParamValues[i])
    
    if RunLSChoice:
        Params[0].append('Run_Least_Squares')
        Params[1].append('TRUE')
    else:
        Params[0].append('Run_Least_Squares')
        Params[1].append('FALSE')
    
    if ProfileChoice == 'Triangle':
        nP = 1
        ParamNames  = ['Pressure Peak']
        ParamValues = Out[nC:nC+nP+1]
        for i in range(nP):
            print(ParamNames[i],': ',ParamValues[i])
            Params[4].append(ParamNames[i])
            Params[5].append(ParamValues[i])
        nC +=1
    elif ProfileChoice == 'Step':
        nP = 2
        ParamNames  = ['Pressure_Increase','Pressure_Decrease']
        ParamValues = Out[nC:nC+nP+1]
        for i in range(nP):
            print(ParamNames[i],': ',ParamValues[i])
            Params[4].append(ParamNames[i])
            Params[5].append(ParamValues[i])
        nC +=2
    elif ProfileChoice == 'SmoothStep':
        nP = 4
        ParamNames  = ['Pressure_Increase','Pressure_Decrease','Increase_Angle','Decrease_Angle']
        ParamValues = Out[nC:nC+nP+1]
        for i in range(nP):
            print(ParamNames[i],': ',ParamValues[i])
            Params[4].append(ParamNames[i])
            Params[5].append(ParamValues[i])
        nC +=4
    elif ProfileChoice == 'Fourier':
        nP = 4
        ParamNames  = ['Fourier_a1','Fourier_a2','Fourier_a3','Fourier_a4']
        ParamValues = Out[nC:nC+nP+1]
        for i in range(nP):
            print(ParamNames[i],': ',ParamValues[i])
            Params[4].append(ParamNames[i])
            Params[5].append(ParamValues[i])
        nC +=4
    elif ProfileChoice == 'Fitted':
        nP = nF-2
        ParamNames  = ['Pressure_P_'+ str(i+1) for i in range(nF)]
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
    elif ProfileChoice == 'SetTriangle':
        nP = 1
        ParamNames  = ['Pressure Peak']
        ParamValues = [0.5]
        for i in range(nP):
            print(ParamNames[i],': ',ParamValues[i])
            Params[4].append(ParamNames[i])
            Params[5].append(ParamValues[i])
    elif ProfileChoice == 'SetStep':
        nP = 2
        ParamNames  = ['Pressure_Increase','Pressure_Decrease']
        ParamValues = [0.2,0.5]
        for i in range(nP):
            print(ParamNames[i],': ',ParamValues[i])
            Params[4].append(ParamNames[i])
            Params[5].append(ParamValues[i])
    elif ProfileChoice == 'SetSmoothStep':
        nP = 4
        ParamNames  = ['Pressure_Increase','Pressure_Decrease','Increase_Angle','Decrease_Angle']
        ParamValues = [0.2,0.5,50,50]
        for i in range(nP):
            print(ParamNames[i],': ',ParamValues[i])
            Params[4].append(ParamNames[i])
            Params[5].append(ParamValues[i])
    elif ProfileChoice == 'SetFourier':
        nP = 4
        ParamNames  = ['Fourier_a1','Fourier_a2','Fourier_a3','Fourier_a4']
        ParamValues = [3.0,1.0,0.0,0.05]
        for i in range(nP):
            print(ParamNames[i],': ',ParamValues[i])
            Params[4].append(ParamNames[i])
            Params[5].append(ParamValues[i])
    elif ProfileChoice == 'SetFitted':
        nP = nF-2
        ParamNames  = ['Pressure_P_'+ str(i) for i in range(nF)]
        ParamValues = np.zeros(nF)
        for i in range(nP):
            print(ParamNames[i],': ',ParamValues[i])
            Params[4].append(ParamNames[i])
            Params[5].append(ParamValues[i])
    elif ProfileChoice == 'SetWindkessel':
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
    PressureData = GetPressure(ParamValues,0,nF,CF,ProfileChoice,TimeData)
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
    Params[1].append(ProfileChoice)
        
    if FibChoice and ModelChoice == 'tiMR':
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
        
    Conditions +=  ResChoice
     
    # Save Pressure Profile Choice
    Params[0].append('Residual_Choice')
    Params[1].append(ResChoice)
    
    if RunLSChoice:
        RunDir = 'RunLS'
    else:
        RunDir = 'RunDefault'
        
    with open('./NewFiles/'+ DataDir+'/'+ModelChoice+'/'+RunDir + '/' + ProfileChoice+ '/' + Conditions+'/Parameters.csv','w') as result_file:
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
    
    columns_data = zip_longest(*Params)
        
    with open('./NewFiles/'+ DataDir+'/'+ModelChoice+'/'+RunDir + '/' + ProfileChoice+ '/' + Conditions+'/Parameters.csv','w') as result_file:
        wr = csv.writer(result_file, dialect='excel')
        wr.writerows(columns_data)
         
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
        
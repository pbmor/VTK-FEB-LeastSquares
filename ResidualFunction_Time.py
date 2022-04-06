import numpy as np
import os
import vtk
import csv
import glob
from xml.etree import ElementTree as et
from scipy.optimize import least_squares
from vtk.util.numpy_support import vtk_to_numpy
from Read_XPLTfuncs import GetFEB, GetMeshInfo, GetData
from RemeshAll import Remesh_All
from vtk2feb_Disp import VTK2Feb_Func
import matplotlib.pyplot as plt
from math import exp as e

def SaveFiles(DataDir,RemeshedFile,Disp_Wall,FId):
    
    xpltName = DataDir + '.xplt'
    
    feb,file_size,nstates, mesh = GetFEB(xpltName)
    
    nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
    
    if len(StateTimes) != nstates:
        print('Warning: Noumber of Frames do n0t match.')
    
    #Collect the data
    displacement = GetData(feb,'displacement',nstates,nVar)
    stress = GetData(feb,'stress',nstates,nVar)
    
    #Gather stresses and displacements Note this is dependent on the Var Type
    Stress_X, Stress_Y, Stress_Z, Stress_XY, Stress_YZ, Stress_XZ = np.zeros((nstates,nElems)),np.zeros((nstates,nElems)),np.zeros((nstates,nElems)),np.zeros((nstates,nElems)),np.zeros((nstates,nElems)),np.zeros((nstates,nElems))
    Disp_X,Disp_Y,Disp_Z = np.zeros((nstates,nNodes)),np.zeros((nstates,nNodes)),np.zeros((nstates,nNodes))
    Disp_FEB  = np.zeros((nstates,nNodes,3))      
    
    for i in range(nstates): 
        for j in range(nNodes):
            Disp_X[i,j]  = displacement[i][j*3+1]
            Disp_Y[i,j]  = displacement[i][j*3+2]
            Disp_Z[i,j]  = displacement[i][j*3+3]
        Disp_FEB[i,:,0]  = Disp_X[i,:]  
        Disp_FEB[i,:,1]  = Disp_Y[i,:]  
        Disp_FEB[i,:,2]  = Disp_Z[i,:]
        
    for i in range(nstates): 
        for j in range(nElems):
            Stress_X[i,j]  = stress[i][j*6+1]
            Stress_Y[i,j]  = stress[i][j*6+2]
            Stress_Z[i,j]  = stress[i][j*6+3]
            Stress_XY[i,j] = stress[i][j*6+4]
            Stress_YZ[i,j] = stress[i][j*6+5]
            Stress_XZ[i,j] = stress[i][j*6+6]
     
    
    Residual = np.subtract(Disp_Wall,Disp_FEB)
        
    reader = vtk.vtkPolyDataReader()
    print(FId)
    
    #Save a new vtk file for each state
    for j, fid in enumerate((FId.astype(int)-1)):
        print(j)
        # Read the source file.
        reader.SetFileName(RemeshedFile)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        polydata = reader.GetOutput()
        
        # Add Cell Data on Points
        CellData = [Stress_X[j], Stress_Y[j], Stress_Z[j], Stress_XY[j], Stress_YZ[j], Stress_XZ[j]]
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
        PointData = [Stress_X_Pt,Stress_Y_Pt,Stress_Z_Pt,Stress_XY_Pt,Stress_YZ_Pt,Stress_XZ_Pt,Disp_X[j],Disp_Y[j],Disp_Z[j],Disp_FEB[j],Residual[j]]
        PointNames = ['Stress_X_Pt','Stress_Y_Pt','Stress_Z_Pt','Stress_XY_Pt','Stress_YZ_Pt','Stress_XZ_Pt','Disp_X','Disp_Y','Disp_Z','Disp_FEB','Residual']
        
        for i in range(len(PointNames)) :
            arrayPoint = vtk.util.numpy_support.numpy_to_vtk(PointData[i], deep=True)
            arrayPoint.SetName(PointNames[i])
            dataPoints = polydata.GetPointData()
            dataPoints.AddArray(arrayPoint)
            dataPoints.Modified()
            
        fname = './NewFiles/'+ xpltName[0:5] + '_' + str(fid+1)+'.vtk'
        directory = os.path.dirname(fname)
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        writer = vtk.vtkDataSetWriter()
        writer.SetFileName(fname)
        writer.SetInputData(polydata)
        print('Writing ',fname)
        writer.Write()
        
    return nstates
    
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

def GetRes(C,filename,DataDir,Disp_Wall,FId,nF,FT):
    
    #Define Timesteps
    TimeInterp = np.linspace(0,nF*float(FT)/1000,nF)
    # Use xml tree to update .feb file
    tree = et.parse(DataDir + '.feb')
    if len(C)==1:
        tree.find('Loads').find('surface_load').find('pressure').text = str(C[0])
    elif len(C)==2:
        tree.find('Loads').find('surface_load').find('pressure').text = str(C[0])
        TimeInterp = np.linspace(0,nF*float(FT)/1000,nF)
        PressureInterp = np.zeros(nF)
        Line_A = (1/C[1])*TimeInterp
        Line_B = (-1/(nF*float(FT)/1000 - C[1]))*TimeInterp + nF*float(FT)/1000/(nF*float(FT)/1000-C[1])
        for i in range(nF):
            if TimeInterp[i] <=C[1]:
                PressureInterp[i] = Line_A[i]
            else:
                PressureInterp[i] = Line_B[i]  
        for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
            Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx])
    elif len(C)==3:
        tree.find('Loads').find('surface_load').find('pressure').text = str(C[0])
        PressureInterp = np.zeros(nF)
        for i in range(nF):
            if (TimeInterp[i] >=C[1]) and (TimeInterp[i] <=C[2]):
                PressureInterp[i] = 1.0
        for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
            Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx])
    elif len(C)==4:
        tree.find('Material').find('material').find('c1').text = str(C[0])
        tree.find('Material').find('material').find('c2').text = str(C[1])
        tree.find('Material').find('material').find('k').text = str(C[2])
        tree.find('Loads').find('surface_load').find('pressure').text = str(C[3])
    # elif len(C)==5:
    #     tree.find('Loads').find('surface_load').find('pressure').text = str(C[0])
    #     Line_A = ( 1/(C[2]-C[1]))*(TimeInterp-C[1])
    #     Line_B = (-1/(C[4]-C[3]))*(TimeInterp-C[4])
    #     PressureInterp = np.zeros(nF)
    #     for i in range(nF):
    #         if (TimeInterp[i] >=C[1]) and (TimeInterp[i] <=C[2]):
    #             PressureInterp[i] = Line_A[i]
    #         elif (TimeInterp[i] >=C[2]) and (TimeInterp[i] <=C[3]):
    #             PressureInterp[i] = 1.0
    #         elif (TimeInterp[i] >=C[3]) and (TimeInterp[i] <=C[4]):
    #             PressureInterp[i] = Line_B[i]  
    #     for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
    #         Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx])
    # elif len(C)==8:
        # tree.find('Loads').find('surface_load').find('pressure').text = str(C[0])
        # # tree.find('Material').find('material').find('c1').text = str(C[1])
        # # tree.find('Material').find('material').find('c2').text = str(C[2])
        # # tree.find('Material').find('material').find('k').text = str(C[3])
        # PressureInterp = np.zeros(nF)
        # for i in range(nF):
        #     PressureInterp[i] = (1/(1+e(-C[6]*(TimeInterp[i]-C[4]))))*(1-1/(1+e(-C[7]*(TimeInterp[i]-C[5]))))
        # for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
        #     Pt.text = str(TimeInterp[idx])+','+str(PressureInterp[idx])
        
              
            
    tree.write(DataDir + '.feb',xml_declaration=True,encoding="ISO-8859-1")
    
    #Run updated .feb file to create new .xplt file
    os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i '+ DataDir+'.feb >/dev/null 2>&1')
    # os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i '+ DataDir+'.feb')
    
    XPLTfilename = DataDir + '.xplt'
    
    feb,file_size,nstates, mesh = GetFEB(XPLTfilename)
    
    nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
    
    #Collect the data, with all states
    displacement = GetData(feb,'displacement',nstates,nVar)
      
    #Gather stresses and displacements Note this is dependent on the Var Type
    Disp_X,Disp_Y,Disp_Z = np.zeros((nstates,nNodes)),np.zeros((nstates,nNodes)),np.zeros((nstates,nNodes))
    Disp_FEB, Res  = np.zeros((nstates,nNodes,3)), np.zeros((nstates,nNodes,3))  
    for i in range(nstates): 
        for j in range(nNodes):
            Disp_X[i,j]  = displacement[i][j*3+1]
            Disp_Y[i,j]  = displacement[i][j*3+2]
            Disp_Z[i,j]  = displacement[i][j*3+3]
        Disp_FEB[i,:,0]  = Disp_X[i,:]  
        Disp_FEB[i,:,1]  = Disp_Y[i,:]  
        Disp_FEB[i,:,2]  = Disp_Z[i,:]
    
    Disp_FEB_Interp = np.zeros(Disp_Wall.shape)
    for i in range(nNodes):
        for j in range(3):
            Disp_FEB_Interp[:,i,j] = np.interp(TimeInterp,StateTimes,Disp_FEB[:,i,j])
    
    Res = np.subtract(Disp_Wall,Disp_FEB_Interp)
    
    return Res.flatten()

def RemeshRef(DataDir,d,DataInfo):
    
    refN = int(DataInfo[4])-1 #Choose frame before valve opens
    
    #for 'medial_meshes' in List_of_subdirectories:
    fnames = sorted(glob.glob(d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/*.vtp'))

    common = os.path.commonprefix(fnames)
    for Fname in list(fnames):
        X = Fname.replace(common,'')
        X = X.replace('.vtp','')
        X = np.fromstring(X, dtype=int, sep=' ')
        X=X[0]
        if X==refN:
            ref=Fname
        else:
            Remesh_All(DataDir, Fname)
    
    Remesh_All(DataDir, ref)
    return 
    
def RunLS(DataDir,d,Remesh,FListOrdered,FId,ref):
    
    # Read the source file
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(ref)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    polydata = reader.GetOutput()
    nPts = polydata.GetNumberOfPoints()
    
    STJ = vtk_to_numpy(polydata.GetPointData().GetArray('STJ'))
    VAJ = vtk_to_numpy(polydata.GetPointData().GetArray('VAJ'))
    
    STJ_Id = []
    VAJ_Id = []
    nSTJ = 0
    nVAJ = 0
    for i in range(nPts):
        if STJ[i] !=0.0:
            STJ_Id.append(i)
            nSTJ+=1
        if VAJ[i] !=0.0:
            VAJ_Id.append(i)
            nVAJ+=1
                
    Disp_Wall_STJ, Disp_Wall_VAJ, Circ, Long, Disp_Wall = np.zeros((nF,nSTJ,3)), np.zeros((nF,nSTJ,3)), np.zeros((nF,nPts,3)), np.zeros((nF,nPts,3)),np.zeros((nF,nPts,3))
    
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
    
    BioPressure = False 
    
    VTK2Feb_Func(DataDir,ref,nF,FT,Disp_Wall_STJ,Disp_Wall_VAJ,STJ_Id,VAJ_Id,Circ,Long,BioPressure)

    # os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i '+ DataDir+'.feb')
    # Out = least_squares(GetRes,[1,0,10,-0.005332],jac = '3-point', verbose=2,args=(ref,DataDir,Disp_Wall,FId,nF,FT))
    # Out = least_squares(GetRes,[-0.005332],jac = '3-point', verbose=2,args=(ref,DataDir,Disp_Wall,FId,nF,FT))
    # Out = least_squares(GetRes,[-0.005332,float(FT)/1000,3*float(FT)/1000,5*float(FT)/1000,6*float(FT)/1000],jac = '3-point', verbose=2,args=(ref,DataDir,Disp_Wall,FId,nF,FT))
    Out = least_squares(GetRes,[-0.005332,1,0,10,0.2,0.5,50,50],bounds = [[-2,-100,-100,-100,0,0,0,0],[0,100,100,100,1,1,2000,2000]],jac = '3-point', verbose=2,args=(ref,DataDir,Disp_Wall,FId,nF,FT))
    # Out = least_squares(GetRes,[-0.005332,float(FT)*(float(nF)/(2*1000))],jac = '3-point', verbose=2,args=(ref,DataDir,Disp_Wall,FId,nF,FT))
    
    return Out, Disp_Wall
    
if __name__=='__main__':

    ''' 
    # Code to run through all data sets
    List_of_Subdirectories = sorted(glob.glob('../RunVTK/Strains/medial_meshes/*'))
    
    # List_of_Subdirectories = '../RunVTK/Strains/medial_meshes/bav02'
    ND = len(List_of_Subdirectories)
    
    CommonOfDir = os.path.commonprefix(List_of_Subdirectories)
    for d in List_of_Subdirectories:
        DataDir = d.replace(CommonOfDir,'')
        Out = RunLS(DataDir,d)
    '''
    
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
            
    #Define Frame Time Length
    FT = DataInfo[1]

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
    Remesh = False # Choice to remesh reference frame mesh is necessary
    
    #Run least squares script
    Out, Disp_Wall = RunLS(DataDir,d,Remesh,FListOrdered,FId,ref)
    
    RemeshedFile = glob.glob('./Remeshed/'+DataDir+'/*')[0]

    nstates = SaveFiles(DataDir,RemeshedFile,Disp_Wall,FId)
    
    # Time = np.zeros(nstates)
    # Pressure = np.zeros(nstates)
    
    # tree = et.parse(DataDir + '.feb')
    # Pressure_Mag = float(tree.find('Loads').find('surface_load').find('pressure').text)
    # for idx, Pt in enumerate(tree.find('LoadData').find('load_controller').find('points').findall("point")):
    #     for i in range(len(Pt.text)):
    #         if Pt.text[i] == ',':
    #             Time[idx] = float(Pt.text[0:i])
    #             Pressure[idx] = float(Pt.text[i+1:])
                
    # plt.plot(Time,np.multiply(Pressure,-Pressure_Mag))
    # plt.xlabel('Time')
    # plt.ylabel('Pressure')
    # plt.show()
        
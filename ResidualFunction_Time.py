import numpy as np
import os
#import sys
import vtk
import csv
import glob
from xml.etree import ElementTree as et
from scipy.optimize import least_squares
from vtk.util.numpy_support import vtk_to_numpy
from Read_XPLTfuncs import GetFEB, GetMeshInfo, GetData
from compute_strain import OrderList
from RemeshAll import Remesh_All
from vtk2feb_Func import VTK2Feb_Func

def GetRes(C,filename,DataDir,Disp_VTK):
    print('Model parameters: c1 = ',C[0],' c2 = ',C[1],'k = ',C[2] )
    
    '''
    Use xml tree to update .feb file
    '''
    tree = et.parse(DataDir + '.feb')
    tree.find('Material').find('material').find('c1').text = str(C[0])
    tree.find('Material').find('material').find('c2').text = str(C[1])
    tree.find('Material').find('material').find('k').text = str(C[2])
    tree.write(DataDir + '.feb',xml_declaration=True,encoding="ISO-8859-1")
    
    #Run updated .feb file to create new .xplt file
    os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i '+ DataDir+'.feb >/dev/null 2>&1')
    
    XPLTfilename = DataDir + '.xplt'
    
    feb,file_size,nstates = GetFEB(XPLTfilename)
    
    nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
    
    #Collect the data, with all states, note spaces in vairable names are replaced with _
    displacement = GetData(feb,'displacement',nstates,nVar)
      
    #Gather stresses and displacements Note this is dependent on the Var Type
    Disp_X    = np.zeros((nstates,nNodes))
    Disp_Y    = np.zeros((nstates,nNodes))
    Disp_Z    = np.zeros((nstates,nNodes))
    Disp_FEB  = np.zeros((nstates,nNodes,3))      
    for i in range(nstates): 
        for j in range(nNodes):
            Disp_X[i,j]  = displacement[i][j*3+1]
            Disp_Y[i,j]  = displacement[i][j*3+2]
            Disp_Z[i,j]  = displacement[i][j*3+3]
        Disp_FEB[i,:,0]  = Disp_X[i,:]  
        Disp_FEB[i,:,1]  = Disp_Y[i,:]  
        Disp_FEB[i,:,2]  = Disp_Z[i,:]
    
    Res = np.subtract(Disp_VTK,Disp_FEB)

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
    
    RemeshFilename = Remesh_All(DataDir, ref)
    return RemeshFilename
    
def RunLS(DataDir,d,Remesh):
    with open('echoframetime.csv') as csv_file:
        XLData = csv.reader(csv_file, delimiter=',')
        for row in XLData:
            if DataDir == row[0]:
                DataInfo = row
    refN = int(DataInfo[4])-1 #Choose frame before valve opens
    
    fnames = sorted(glob.glob(d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/*.vtp'))

    common = os.path.commonprefix(fnames)
    for Fname in list(fnames):
        X = Fname.replace(common,'')
        X = X.replace('.vtp','')
        X = np.fromstring(X, dtype=int, sep=' ')
        X=X[0]
        if X==refN:
            ref=Fname
            
    #Define Frame Time Length
    FT = DataInfo[1]
    
    if Remesh:
        RemeshFilename = RemeshRef(DataDir,d,DataInfo)
    else:
        RemeshFilename = os.path.split(os.path.splitext(ref)[0])[1] + '.vtk'

    VTKfilename = os.path.join('Remeshed/',DataDir,RemeshFilename)
    
    print('Reading Reference Frame:', VTKfilename)
    
    # Read the source file
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(VTKfilename)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    
    polydata = reader.GetOutput()
    Pts_re = vtk_to_numpy(polydata.GetPoints().GetData())

    STJ       = vtk_to_numpy(polydata.GetPointData().GetArray('STJ'))
    VAJ       = vtk_to_numpy(polydata.GetPointData().GetArray('VAJ'))
    NP = polydata.GetNumberOfPoints()
    nSTJ = np.count_nonzero(STJ == 1.)
    nVAJ = np.count_nonzero(VAJ == 1.)
    j=0
    k=0
    Pts_STJ_re = np.zeros((nSTJ,3))
    Pts_VAJ_re = np.zeros((nVAJ,3))
    for i in range(NP):
        if STJ[i] ==1:
            Pts_STJ_re[j] = Pts_re[i]
            j+=1
        if VAJ[i] ==1:
            Pts_VAJ_re[k] = Pts_re[i]
            k+=1
    
    #for 'medial_meshes' in List_of_subdirectories:
    fnames = sorted(glob.glob(d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/*.vtp'))
    
    nF = len(fnames)
    
    STJ_Motion = np.zeros((nF,nSTJ,3))
    VAJ_Motion = np.zeros((nF,nVAJ,3))
    Disp_Wall_re = np.zeros((nF,NP,3))    
    
    STJ_WallMotion = np.zeros((nF,nSTJ,3))
    VAJ_WallMotion = np.zeros((nF,nVAJ,3))
    STJ_Id = np.zeros((nF,nSTJ))
    VAJ_Id = np.zeros((nF,nVAJ))
    
    # Re-order filename list to start with reference frame
    FListOrdered, FId, refN = OrderList(fnames,nF,ref)
    # Analyse each frame
    for X,Fname in enumerate(FListOrdered):

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(Fname)
        reader.Update()
        polydata = reader.GetOutput()
        Pts = vtk_to_numpy(polydata.GetPoints().GetData())

        Disp_wall = vtk_to_numpy(polydata.GetPointData().GetArray('Displacement_Wall'))
        STJ       = vtk_to_numpy(polydata.GetPointData().GetArray('STJ'))
        VAJ       = vtk_to_numpy(polydata.GetPointData().GetArray('VAJ'))
            
        NP = polydata.GetNumberOfPoints()
        j=0
        k=0
        Pts_STJ = np.zeros((36,3))
        Pts_VAJ = np.zeros((36,3))
        Disp_wall_STJ = np.zeros((36,3))
        Disp_wall_VAJ = np.zeros((36,3))
        STJ_Id = np.zeros(36)
        VAJ_Id = np.zeros(36)
        for i in range(NP):
            if STJ[i] ==1:
                Pts_STJ[j] = Pts[i]
                STJ_Id[j] = i
                Disp_wall_STJ[j] = Disp_wall[i]
                j+=1
            if VAJ[i] ==1:
                Pts_VAJ[k] = Pts[i]
                VAJ_Id[k] = i
                Disp_wall_VAJ[k] = Disp_wall[i]
                k+=1
                
        STJ_Motion[X,:,0]   = np.interp(Pts_STJ_re[:,0],Pts_STJ[:,0],Disp_wall_STJ[:,0])
        VAJ_Motion[X,:,0]   = np.interp(Pts_VAJ_re[:,0],Pts_VAJ[:,0],Disp_wall_VAJ[:,0])
        Disp_Wall_re[X,:,0] = np.interp(Pts_re[:,0],Pts[:,0],Disp_wall[:,0])
        STJ_Motion[X,:,1]   = np.interp(Pts_STJ_re[:,1],Pts_STJ[:,1],Disp_wall_STJ[:,1])
        VAJ_Motion[X,:,1]   = np.interp(Pts_VAJ_re[:,1],Pts_VAJ[:,1],Disp_wall_VAJ[:,1])
        Disp_Wall_re[X,:,1] = np.interp(Pts_re[:,1],Pts[:,1],Disp_wall[:,1])
        STJ_Motion[X,:,2]   = np.interp(Pts_STJ_re[:,2],Pts_STJ[:,2],Disp_wall_STJ[:,2])
        VAJ_Motion[X,:,2]   = np.interp(Pts_VAJ_re[:,2],Pts_VAJ[:,2],Disp_wall_VAJ[:,2])
        Disp_Wall_re[X,:,2] = np.interp(Pts_re[:,2],Pts[:,2],Disp_wall[:,2])
    
    VTK2Feb_Func(DataDir,VTKfilename,nF,FT,STJ_WallMotion,VAJ_WallMotion,STJ_Id.astype(int),VAJ_Id.astype(int))

    Out = least_squares(GetRes,[1,0,10],jac = '3-point', verbose=2,args=(VTKfilename,DataDir,Disp_Wall_re))
    
    return Out
    
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
    Remesh = False # Choice to remesh reference frame mesh is necessary
    Out = RunLS(DataDir,d,Remesh)
    
    xpltName = DataDir + '.xplt'
    os.system('python Run_XPLT.py ' + xpltName + ' >/dev/null 2>&1')
    
    
    
import numpy as np
from AnalysisFuncs import GetRadii
from compute_strain import OrderList
import matplotlib.pyplot as plt
import csv
import os
import glob
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from matplotlib.pyplot import figure, Line2D
from mpl_toolkits import mplot3d


with open('echoframetime.csv') as csv_file:
    XLData = csv.reader(csv_file, delimiter=',')

List_of_Subdirectories = sorted(glob.glob('Strains/medial_meshes/*'))
ND = len(List_of_Subdirectories)

CommonOfDir = os.path.commonprefix(List_of_Subdirectories)
for figN, d in enumerate(List_of_Subdirectories):

    DataDir = d.replace(CommonOfDir,'')

    with open('echoframetime.csv') as csv_file:
        XLData = csv.reader(csv_file, delimiter=',')
        for row in XLData:
            if DataDir == row[0]:
                DataInfo = row

    #Define Frame Time Length
    FT = DataInfo[1]
    #Define Open Frame and Close Frame
    OF = DataInfo[4]
    CF = DataInfo[5]
    
    refN = int(DataInfo[4])-1 #Choose frame before valve opens

    #RootType = DataDir[0:3]
    
    #for 'medial_meshes' in List_of_subdirectories:
    fnames = sorted(glob.glob(d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/*.vtp'))

    common = os.path.commonprefix(fnames)
    for Fname in list(fnames):
        X = Fname.replace(common,'')
        X = X.replace('.vtk','')
        X = np.fromstring(X, dtype=int, sep=' ')
        X=X[0]
        if X==refN:
            ref=Fname
    Nf = len(fnames)

    STJ_WallMotion  = np.zeros((Nf,36,3))
    STJ_RootMotion  = np.zeros((Nf,36,3))
    STJ_TotalMotion = np.zeros((Nf,36,3))
    STJ_Pts = np.zeros((36,3))
    Time            = np.zeros(Nf)
    # Re-order filename list to start with reference frame
    FListOrdered, FId, refN = OrderList(fnames,Nf,ref)

    # Analyse each frame
    for X,Fname in enumerate(FListOrdered):
        if FId[0]<=FId[X]:
            Time[X] = float(FId[X])*float(FT)
        else :
            if FId[X]-FId[X-1]<0:
                Time[X] = Time[X-1]+1*float(FT)
            elif FId[X]-FId[X-1]>0:
                Time[X] = Time[X-1]+(FId[X]-FId[X-1])*float(FT)

        if Fname == ref:
            refID = X

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(Fname)
        reader.Update()
        polydata = reader.GetOutput()
        points = vtk_to_numpy(polydata.GetPoints().GetData())

        pointData = polydata.GetPointData()
        Disp_wall  = vtk_to_numpy(pointData.GetArray('Displacement_Wall'))
        Disp_root  = vtk_to_numpy(pointData.GetArray('Displacement_Root'))
        Disp_total = vtk_to_numpy(pointData.GetArray('Displacement_Total'))
        STJ        = vtk_to_numpy(pointData.GetArray('STJ'))
            
        NP = polydata.GetNumberOfPoints()
        j=0
        for i in range(NP):
            if STJ[i] == 1:
                STJ_WallMotion[X,j]  = Disp_wall[i]
                STJ_RootMotion[X,j]  = Disp_root[i]
                STJ_TotalMotion[X,j] = Disp_total[i]
                if X == refID:
                    STJ_Pts[j] = points[i]
                j += 1
   
    DataLocation = os.path.join('./STJMotionData/',DataDir+'_STJMotion.npz')
    np.savez(DataLocation,Time = Time ,STJ_Pts = STJ_Pts, STJ_WallMotion=STJ_WallMotion,STJ_RootMotion=STJ_RootMotion,STJ_TotalMotion=STJ_TotalMotion)


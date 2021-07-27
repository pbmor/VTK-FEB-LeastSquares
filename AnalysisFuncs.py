import numpy as np
import os
import sys
import glob
import csv
from matplotlib.pyplot import figure
from matplotlib import cm
from mpl_toolkits import mplot3d
import numpy as np


def FixAndRotate(Pts,Nf,NP):
    '''
    A function to return the points of the mesh, for all the times points, 
    that has been fixed to the point (0,0,0) and rotated so that the centre 
    line is set to be on the z axis.
    Inputs:
    Pts - Mesh points of all the time frames. dimensions (Nf, NP, 3)
    Nf - Number of Frames
    NP - Number of points 
    '''

    #Create empyt array
    Points = np.zeros((Nf,NP,3))

    for X in range(Nf):
        #Access points of a single frame
        points = Pts[X,:,:]
        ##############################
        # Redefine Points and Displacement vectors
        Mids = np.zeros((3))
        for i in range(3):
            Mids[i] = np.mean([np.amin(points[:,i]),np.amax(points[:,i])])
    
        #Move mesh centre to (0,0,0)
        for i in range(NP):
            points[i,:] = [points[i,j]-Mids[j] for j in range(3)]

        #################################
        # Rotate Mesh
        Pointsb = points.T
        Points2 = np.zeros((3,25,36))
        for i in range(3):
         for j in range(25):
          for k in range(36):
           Points2[i,j,k] = Pointsb[i,((j)%25+k*25)%900]
        
        #Find centre line
        centre_line = np.zeros((3,25))
        for i in range(3):
            for j in range(25):
                centre_line[i,j] = np.mean(Points2[i,j,:])

        #Rotate points and save to empty array
        a   = centre_line[:,13]
        mag = np.linalg.norm(a)
        a   = a/mag
        b   = [0,0,1]
        v   = np.cross(a,b)
        c   = np.dot(a,b)
        vx  = [[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]]
        vx2 = np.matmul(vx,vx)
        R = np.identity(3)+vx+vx2*(1/(1+c))
        RPoints = np.matmul(R,points.T)
        Points[X,:,:] = RPoints.T

    return Points



def GetRadii(Points,Nf,NP,Ndim):
    FramePoints = np.zeros((Nf,25,36,3))
    Radius = np.zeros((Nf,25,36))
    zAvg = np.zeros((Nf,25))
    for X in range(Nf):

        for i in range(25):
         for j in range(36):
          for k in range(3):
           FramePoints[X,i,j,k] = Points[X,((13+i)%25+j*25)%900,k]

        Ring = np.zeros((25,37,3))

        for i in range(25):
            for k in range(3):
                Ring[i,0:36,k] = FramePoints[X,i,:,k]
                Ring[i,36,k]   = Ring[i,0,k]

        for i in range(25):
            zAvg[X,i] = np.mean(FramePoints[X,i,0:35,2])
            for j in range(36):
                Radius[X,i,j] = np.sqrt(Ring[i,j,0]**2+Ring[i,j,1]**2)

    return Radius, FramePoints, zAvg


def TetraMesh(polydata,Points,NC,NP):
    import tetgen

    Points2 = np.zeros((3,25,36))
    TopIds = np.zeros(36)
    BottomIds = np.zeros(36)
    centreTop = np.zeros((1,3))
    centreBottom = np.zeros((1,3))
    Cells  = np.zeros((NC+36*2,3))

    for i in range(3):
     for j in range(25):
      for k in range(36):
       Points2[i,j,k] = Points[((j)%25+k*25)%900,i]

    for i in range(NC):
        #Assign point id
        Cells[i,:] = np.array([polydata.GetCell(i).GetPointId(j) for j in range(3)])


    for i in range(NP):
        Pt = np.array(polydata.GetPoint(i))
        for k in range(36):
            if Pt[0] == Points2[0,13,k] and Pt[1] == Points2[1,13,k] and Pt[2] == Points2[2,13,k]:
                TopIds[k] = i
            if Pt[0] == Points2[0,12,k] and Pt[1] == Points2[1,12,k] and Pt[2] == Points2[2,12,k]:
                BottomIds[k] = i

    for i in range(3):
     for j in range(25):
      for k in range(36):
       Points2[i,j,k] = Points[((j)%25+k*25)%900,i]


    for i in range(3):
        centreTop[0,i] = np.mean(Points2[i,13,:])
        centreBottom[0,i] = np.mean(Points2[i,12,:])

    Points = np.append(Points, centreTop, axis=0)
    Points = np.append(Points, centreBottom, axis=0)

    for i in range(36):
        Cells[NC+i,:] = [TopIds[i],TopIds[(i+1)%36], NP]

    for i in range(36):
        Cells[NC+36+i,:] = [BottomIds[i],BottomIds[(i+1)%36], NP+1]

    Cells=Cells.astype(int)

    tet = tetgen.TetGen(Points, Cells)
    tet.tetrahedralize(order=1, mindihedral=100, minratio=10)
    grid = tet.grid

    return grid

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
    # Sort frame labels
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
            if X[0] ==F:
                FListOrdered[i] = Fname

    return FListOrdered, FId, refN

def calDisp(polydata,Points,NP,RefPointsFixed):
    '''
    Calculate displacement with respect to the reference points split into total, wall and root displacements

    Keyword arguments:
    polydata -- vtk object of the current timeframe
    Points -- numpy array of the current coordinates
    NP -- number of points
    RefPointsFixed -- numpy array of the reference coordinates, assumed to be fixed with geometric center at (0,0,0)

    Returns three numpy arrays:
    first -- total displacement
    second -- displacement of the wall relative to the fixed root
    third -- displacement of the whole root without the movement of the wall
    first == second + third
    '''
    #######################################
    # Calculate Displacements and Find Points Fixed Root
    TotDisp  = np.zeros((NP,3))
    WallDisp = np.zeros((NP,3))
    RootDisp = np.zeros((NP,3))
    PointsFixed = np.zeros((NP,3))

    # Define total displacement, relative to reference frame
    for i in range(NP):
        TotDisp[i,:] = polydata.GetPoint(i) - RefPointsFixed[i,:]

    # Find mid points of current frameframe
    Ranges = np.array(polydata.GetPoints().GetBounds())
    Mids = np.zeros((3))
    Mids[0] = np.mean(Ranges[0:2])
    Mids[1] = np.mean(Ranges[2:4])
    Mids[2] = np.mean(Ranges[4:6])


    for i in range(NP):
        # Define points with the centre fixed at (0,0,0)
        PointsFixed[i,:] = Points[i,:] - Mids
        # Define wall displacement with a fixed root
        WallDisp[i,:]   = PointsFixed[i,:] - RefPointsFixed[i,:]
        # Define displacement of root, without wall movement
        RootDisp[i,:]   = TotDisp[i,:] - WallDisp[i,:]
    return TotDisp, WallDisp, RootDisp

def calAreaAndVol(polydata,ThickData,NC,NP):
    '''
    Calculate the wall area, wall volume, and the lumen volume

    Keyword arguments:
    polydata -- vtk object of the mesh
    ThickData -- numpy array of the thickness at the points
    NC -- number of cells
    NP -- number of points

    Returns three scalars: TotalWallArea, TotalWallVolume, TotalLumenVolume
    '''

    #######################################
    # Define Properties: Areas and Volumes

    TotalWallArea    = 0
    TotalWallVolume  = 0
    TotalLumenVolume = 0

    #Define Wall Area and Volume
    for i in range(NC):
        #Assign point id
        t = np.array([polydata.GetCell(i).GetPointId(j) for j in range(3)])
        #For each cell assign thickness for each point
        Thicks = [ThickData[t[j]] for j in range(3)]
        # Find and Sum Cell Areas and Volumes
        CA = polydata.GetCell(i).ComputeArea()
        CV = CA*(np.sum(Thicks))/3
        TotalWallArea   += CA
        TotalWallVolume += CV

    #Find Centre Points
    Ranges = np.array(polydata.GetPoints().GetBounds())
    Mids = np.zeros((3))
    Mids[0] = np.mean(Ranges[0:2])
    Mids[1] = np.mean(Ranges[2:4])
    Mids[2] = np.mean(Ranges[4:6])

    fedges = vtk.vtkFeatureEdges()
    fedges.BoundaryEdgesOn()
    fedges.FeatureEdgesOff()
    fedges.ManifoldEdgesOff()
    fedges.SetInputData(polydata)
    fedges.Update()
    ofedges = fedges.GetOutput()


    Connect = vtk.vtkPolyDataConnectivityFilter()
    Connect.SetInputData(ofedges)
    Connect.SetExtractionModeToAllRegions()
    Connect.ColorRegionsOn()
    Connect.Update()
    Connect = Connect.GetOutput()

    Ring1 = vtk.vtkThreshold()
    Ring1.ThresholdByUpper(0.5)
    Ring1.SetInputData(Connect)
    Ring1.Update()
    Ring1 = Ring1.GetOutput()

    Cap1 = vtk.vtkDelaunay2D()
    Cap1.SetProjectionPlaneMode(vtk.VTK_BEST_FITTING_PLANE)
    Cap1.SetInputData(Ring1)
    Cap1.Update()
    C1Ps = vtk_to_numpy(Cap1.GetOutput().GetPoints().GetData())

    Ring2 = vtk.vtkThreshold()
    Ring2.ThresholdByLower(0.5)
    Ring2.SetInputData(Connect)
    Ring2.Update()
    Ring2 = Ring2.GetOutput()

    Cap2 = vtk.vtkDelaunay2D()
    Cap2.SetProjectionPlaneMode(vtk.VTK_BEST_FITTING_PLANE)
    Cap2.SetInputData(Ring2)
    Cap2.Update()
    C2Ps = vtk_to_numpy(Cap2.GetOutput().GetPoints().GetData())

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(Cap1.GetOutputPort())
    normals.ComputePointNormalsOff()
    normals.ComputeCellNormalsOn()
    normals.Update()
    Norm1 = normals.GetOutput()

    vtkCenters=vtk.vtkCellCenters()
    vtkCenters.SetInputConnection(Cap1.GetOutputPort())
    vtkCenters.Update()
    centersOutput=vtkCenters.GetOutput()
    centers1=np.array([centersOutput.GetPoint(i) for i in range(Norm1.GetNumberOfCells())])

    Norms = vtk_to_numpy(Norm1.GetCellData().GetNormals())
    Ps    = vtk_to_numpy(Norm1.GetPoints().GetData())

    NAvg1 = np.array([np.mean(Norms[:,i]) for i in range(3)])
    PAvg1 = np.array([np.mean(Ps[:,i]) for i in range(3)])

    Cline1 = PAvg1-Mids
    u1 = NAvg1/np.linalg.norm(NAvg1)
    u2 = Cline1/np.linalg.norm(Cline1)
    dot_product = np.dot(u1, u2)
    angle = np.arccos(np.dot(u1,u2))

    if angle<= (np.pi/2):
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputConnection(Cap1.GetOutputPort())
        normals.ComputePointNormalsOff()
        normals.ComputeCellNormalsOn()
        normals.FlipNormalsOn()
        normals.Update()
        Cap1 = normals.GetOutput()
    else:
        Cap1 = Cap1.GetOutput()

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(Cap2.GetOutputPort())
    normals.ComputePointNormalsOff()
    normals.ComputeCellNormalsOn()
    normals.Update()
    Norm2 = normals.GetOutput()

    vtkCenters=vtk.vtkCellCenters()
    vtkCenters.SetInputConnection(Cap2.GetOutputPort())
    vtkCenters.Update()
    centersOutput=vtkCenters.GetOutput()
    centers2=np.array([centersOutput.GetPoint(i) for i in range(Norm2.GetNumberOfCells())])

    Norms = vtk_to_numpy(Norm2.GetCellData().GetNormals())
    Ps    = vtk_to_numpy(Norm2.GetPoints().GetData())

    NAvg2 = np.array([np.mean(Norms[:,i]) for i in range(3)])
    PAvg2 = np.array([np.mean(Ps[:,i]) for i in range(3)])

    Cline2 = PAvg2 - Mids
    u1 = NAvg2/np.linalg.norm(NAvg2)
    u2 = Cline2/np.linalg.norm(Cline2)
    dot_product = np.dot(u1, u2)
    angle = np.arccos(np.dot(u1,u2))

    if angle<= (np.pi/2):
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputConnection(Cap2.GetOutputPort())
        normals.ComputePointNormalsOff()
        normals.ComputeCellNormalsOn()
        normals.FlipNormalsOn()
        normals.Update()
        Cap2 = normals.GetOutput()
    else:
        Cap2 = Cap2.GetOutput()

    polydata.ShallowCopy(polydata)
    Cap1.ShallowCopy(Cap1)
    Cap2.ShallowCopy(Cap2)

    appendFilter = vtk.vtkAppendPolyData()
    appendFilter.AddInputData(Cap1)
    appendFilter.AddInputData(polydata)
    appendFilter.AddInputData(Cap2)
    appendFilter.Update()

    cleanFilter = vtk.vtkCleanPolyData()
    cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
    cleanFilter.Update()

    massUnion = vtk.vtkMassProperties()
    massUnion.SetInputConnection(cleanFilter.GetOutputPort())
    TotalLumenVolume = massUnion.GetVolume()

    return TotalWallArea, TotalWallVolume, TotalLumenVolume

def calStrains(polydata, RA, NC):
    '''
    Calculate strain invariants with respect to the reference basis

    Keyword arguments:
    polydata -- vtk object for the current time point
    RA -- numpy array of the reference basis vectors of size NC X 2 X 3
    NC -- number of cells

    Returns two numpy arrays of invariant J and I1 at the cells
    '''
    I1 = np.zeros(NC)
    J  = np.zeros(NC)

    for i in range(NC):
        Cell = np.array(polydata.GetCell(i).GetPoints().GetData())

        # Define Reference Vectors
        A1 = RA[i,0,:]
        A2 = RA[i,1,:]

        #Define Deformed Vectors
        a1 = Cell[2,:]-Cell[0,:]
        a2 = Cell[2,:]-Cell[1,:]

        G      = np.array([[np.dot(A1,A1),np.dot(A1,A2)],[np.dot(A2,A1),np.dot(A2,A2)]])
        invG   = np.linalg.inv(G)

        A1dual = invG[0,0]*A1 + invG[0,1]*A2
        A2dual = invG[1,0]*A1 + invG[1,1]*A2
        F   = np.outer(a1,A1dual) + np.outer(a2,A2dual)
        C   = np.dot(F.T,F)
        C2  = np.dot(C,C)

        # Define principle strains and J
        trC         = C.trace()
        trC2        = C2.trace()
        I1[i]    = trC
        J[i]     = np.sqrt((trC**2-trC2)/2)

    return J, I1

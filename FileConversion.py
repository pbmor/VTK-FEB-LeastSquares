import numpy as np
import os
import sys
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa
import tetgen

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



def ProcessData(flist,ref,FT,OF=None,CF=None,prefix='Strains/',FixAndRotate=True,opformat='vtp'):
    '''
    Process all the files in a given list, including calculation of strains with respect to a reference

    Keyword arguments:
    flist -- list of filenames
    ref -- name of the reference file
    FT -- time between frames (or framerate)
    OF -- file number of the frame when AV opens
    CF -- file number of the frame when AV closes
    prefix -- prefix to the directory where the resulting files are stored (default 'Strains/')
    FixAndRotate -- whether the points should be fixed and rotated or not (default True)
    opformat -- format of the output files, either vtp (default) or vtk

    The function creates new .vtp files with the following data added:
    calculated strains, curvatures, cell areas and volumes, and field data of total wall area, total wall volume, total lumen volume, AV open/closed status
    
    Then the function returns the following numpy arrays (which have been reordered using OrderList):
    WallArea, WallVol, LumenVol, Time, Pts, WallAreaRatio, WallVolRatio, LumenVolRatio
    and N (the number of frames)
    '''
    print('##############################')
    print('Starting ProcessData Function')
    print('##############################')

    reader = vtk.vtkPolyDataReader()
    #################################
    # Define Reference Vectors
    #################################
    print('Reading Reference Frame:', ref)

    # Read the source file.
    reader.SetFileName(ref)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()

    polydata = reader.GetOutput()
    dataset = polydata.GetPointData()
    RefPoints = vtk_to_numpy(polydata.GetPoints().GetData())
    # Get Number of Points and Cells
    NP = polydata.GetNumberOfPoints()
    NC = polydata.GetNumberOfCells()
    if not polydata.GetPointData().GetArray(0):
        print('Warning Wall Thickness data is missing')
        ThickData = np.zeros((NC,3))
    else:
        ThickData = np.zeros((NC,3))
        NumOfArr = polydata.GetPointData().GetNumberOfArrays()
        for i in range(NumOfArr):
           if polydata.GetPointData().GetArrayName(i) == 'Thickness':
               ThickData = vtk_to_numpy(polydata.GetPointData().GetArray(i))
    
    RefPointsFixed = np.zeros(RefPoints.shape)

    #Find Centre Points
    Ranges = np.array(polydata.GetPoints().GetBounds())
    Mids = np.zeros((3))
    Mids[0] = np.mean(Ranges[0:2])
    Mids[1] = np.mean(Ranges[2:4])
    Mids[2] = np.mean(Ranges[4:6])

    # Define points with mid point fixed at (0,0,0)
    for i in range(NP):
        RefPointsFixed[i,:] = RefPoints[i,:] - Mids


    RefWallArea, RefWallVoliume, RefLumenVolume = calAreaAndVol(polydata,ThickData,NC,NP)
    
    #Define initial frame of wall displacement for later use in motion calculation
    _, WallDispPF, _ = calDisp(polydata,RefPoints,NP,RefPointsFixed)

    # Define Empty Arrays for Reference Data
    RA = np.zeros((NC,2,3))

    for i in range(NC):
        Cells = vtk_to_numpy(polydata.GetCell(i).GetPoints().GetData())
        # Define refernce cell vectors
        RA[i,0,:] = Cells[2,:]-Cells[0,:]
        RA[i,1,:] = Cells[2,:]-Cells[1,:]


    ###########################################
    # Define Data At Every Time Frame
    ###########################################
    #Numer of frames
    N = len(flist)

    # Empty arrays for deformed frames

    Time          = np.zeros(N)
    ValvePosition = np.zeros(N)
    WallArea      = np.zeros(N)
    WallVol       = np.zeros(N)
    LumenVol      = np.zeros(N)
    LumenVolAlt   = np.zeros(N)
    WallAreaRatio = np.zeros(N)
    WallVolRatio  = np.zeros(N)
    LumenVolRatio = np.zeros(N)
    AvgJ          = np.zeros(N)
    AvgI1         = np.zeros(N)
    AvgJRatio     = np.zeros(N)
    AvgI1Ratio    = np.zeros(N)
    Pts           = np.zeros((N,NP,3))
    TotalMotion   = np.zeros(NP)
    
    # Re-order filename list to start with reference frame
    FListOrdered, FId, refN = OrderList(flist,N,ref)

    # Analyse each frame
    common = os.path.commonprefix(flist)
    for X,Fname in enumerate(FListOrdered):
        if OF is not None and CF is not None:
            OpenX  = float(OF)-float(refN)
            CloseX = float(CF)-float(refN)
            if float(X)>=OpenX and float(X)<CloseX:
                ValvePosition[X] = 1

        print('Reading Current Frame:', Fname)

        # Read the source file.
        reader.SetFileName(Fname)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()

        ######################################
        # Define File Data
        polydata  = reader.GetOutput()
        
        # Get time of frames
        if FId[0]<=FId[X]:
            Time[X] = float(FId[X])*float(FT) 
        else :
            if FId[X]-FId[X-1]<0:
                Time[X] = Time[X-1]+1*float(FT)
            elif FId[X]-FId[X-1]>0:
                Time[X] = Time[X-1]+(FId[X]-FId[X-1])*float(FT)

        if not polydata.GetPointData().GetArray(0): 
            print('Warning Wall Thickness data is missing')
            ThickData = np.zeros((NC,3))
        else:
            ThickData = np.zeros((NC,3))
            NumOfArr = polydata.GetPointData().GetNumberOfArrays()
            for i in range(NumOfArr):
               if polydata.GetPointData().GetArrayName(i) == 'Thickness':
                   ThickData = vtk_to_numpy(polydata.GetPointData().GetArray(i))

        #Check that the number of points and cells are consistent
        if NP != polydata.GetNumberOfPoints():
            raise ValueError('Error: the number of points between the reference and deformed frames differ')
        if NC != polydata.GetNumberOfCells():
            raise ValueError('Error: the number of cells between the reference and deformed frames differ')

        # Save Point coordinates to array to be saved
        Pts[X,:,:] = vtk_to_numpy(polydata.GetPoints().GetData())

        #########################################
        # Calculate Displacements
        TotDisp, WallDisp, RootDisp = calDisp(polydata,Pts[X,:,:],NP,RefPointsFixed)
        
        #########################################
        # Calculate Total Wall Area and Volume, and Lumen Volume
        TotalWallArea, TotalWallVolume, TotalLumenVolume = calAreaAndVol(polydata,ThickData,NC,NP)

        # Save areas and volumes for each frame
        WallArea[X]       = TotalWallArea
        WallVol[X]        = TotalWallVolume
        LumenVol[X]       = TotalLumenVolume

        #########################################
        # Calculate Motion Vector (difference between current and previous frame of wall displacement)
        MotionVector = WallDisp - WallDispPF
        Motion = np.zeros(900)

        #Define Motion magnitude
        for i in range(NP):
            Motion[i] = np.sqrt(sum(MotionVector[i,j]**2 for j in range (3)))

        #Update previous frame of wall disp for use in next frame
        WallDispPF = WallDisp
        
        #########################################
        #Mark Location of interatrial septum
        InterAtrialSeptum = np.zeros(900)
        for i in range(NP):
            if i<=24:
                InterAtrialSeptum[i] = 1

        #########################################
        #Define Properties: J and I1
        J, I1 = calStrains(polydata,RA,NC)

        #########################################
        # Define Curvature

        # Define Gaussian Curvature
        curvgauss = vtk.vtkCurvatures()
        curvgauss.SetInputConnection(reader.GetOutputPort())
        curvgauss.SetCurvatureTypeToGaussian()
        curvgauss.Update()
        CurvG = curvgauss.GetOutput()
        CurvG.GetPointData().GetScalars()
        NumOfArrs = CurvG.GetPointData().GetNumberOfArrays()
        # Get curvature data that has been saved to new array
        CDG = np.asarray(CurvG.GetPointData().GetArray(NumOfArrs-1))  

        # Define Mean Curvature
        curvmean = vtk.vtkCurvatures()
        curvmean.SetInputConnection(reader.GetOutputPort())
        curvmean.SetCurvatureTypeToMean()
        curvmean.Update()
        CurvM = curvmean.GetOutput()
        CurvM.GetPointData().GetScalars()
        NumOfArrs = CurvM.GetPointData().GetNumberOfArrays()
        # Get curvature data that has been saved to new array
        CDM = np.asarray(CurvM.GetPointData().GetArray(NumOfArrs-1))

        #####################################
        # Add Data to Files and Write Files

        # Add Cell Data
        CellData  = [I1,J]
        CellNames = ['I1','J'] 
        for i in range(len(CellNames)) :
            arrayCell = vtk.util.numpy_support.numpy_to_vtk(CellData[i], deep=True)
            arrayCell.SetName(CellNames[i])
            dataCells = polydata.GetCellData()
            dataCells.AddArray(arrayCell)

        NumArr = polydata.GetPointData().GetNumberOfArrays()

        # Convert Cell Data to Point Data
        c2p = vtk.vtkCellDataToPointData()
        c2p.AddInputData(polydata)
        c2p.Update()
        c2p.GetOutput()

        I1pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(NumArr))
        Jpt  = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(NumArr+1))

        # Add Point Data
        PointData = [CDG,CDM,I1pt,Jpt,Motion,InterAtrialSeptum]
        PointNames = ['CurvGaussian','CurvMean','I1_Pt','J_Pt','Motion','InterAtrialSeptum'] 
        for i in range(len(PointNames)) :
            arrayPoint = vtk.util.numpy_support.numpy_to_vtk(PointData[i], deep=True)
            arrayPoint.SetName(PointNames[i])
            dataPoints = polydata.GetPointData()
            dataPoints.AddArray(arrayPoint)
            dataPoints.Modified()

        # Add Vector Data
        VectorData = [TotDisp,WallDisp,RootDisp]
        VectorNames = ['Displacement_Total','Displacement_Wall','Displacement_Root'] 
        for i in range(len(VectorNames)) :
            arrayVector = vtk.util.numpy_support.numpy_to_vtk(VectorData[i], deep=True)
            arrayVector.SetName(VectorNames[i])
            dataVectors = polydata.GetPointData()
            dataVectors.AddArray(arrayVector)
            dataVectors.Modified()

        # Add Field Data
        if OF is not None and CF is not None:
            FieldData = [WallArea[X], WallVol[X],LumenVol[X],ValvePosition[X]]
            FieldNames = ['WallArea','WallVolume','LumenVolume','ValvePosition'] 
        else:
            FieldData = [WallArea[X], WallVol[X],LumenVol[X]]
            FieldNames = ['WallArea','WallVolume','LumenVolume'] 

        for i in range(len(FieldNames)) :
            arrayField = vtk.util.numpy_support.numpy_to_vtk(FieldData[i], deep=True)
            arrayField.SetName(FieldNames[i])
            dataFields = polydata.GetFieldData()
            dataFields.AddArray(arrayField)
            dataFields.Modified() 

        # Save points to be the reference frame with centre fixed at (0,0,0)
        if FixAndRotate == True:
            for i in range(NP):
                ptNew = RefPointsFixed[i,:]
                polydata.GetPoints().SetPoint(i, ptNew)

        #################################
        # Write data to vtp files
        if opformat == 'vtp':
            fname = prefix + os.path.splitext(Fname)[0] + '.vtp'
        elif opformat == 'vtk':
            fname = prefix + os.path.splitext(Fname)[0] + '.vtk'
        else:
            raise ValueError("Only vtp and vtk output formats are allowed")
        directory = os.path.dirname(fname)
        if not os.path.exists(directory):
            os.makedirs(directory)

        if opformat == 'vtp':
            writer = vtk.vtkXMLDataSetWriter()
        elif opformat == 'vtk':
            writer = vtk.vtkDataSetWriter()
        else:
            raise ValueError("Only vtp and vtk output formats are allowed")
        writer.SetFileName(fname)
        writer.SetInputData(polydata)
        print('Writing',fname)
        writer.Write()

    WallAreaRatio[:]  = WallArea[:]/WallArea[0]
    WallVolRatio[:]   = WallVol[:]/WallVol[0]
    LumenVolRatio[:]  = LumenVol[:]/LumenVol[0]

    return WallArea, WallVol, LumenVol, Time, Pts, WallAreaRatio, WallVolRatio, LumenVolRatio, AvgJ, AvgI1, AvgJRatio, AvgI1Ratio, TotalMotion, N

if __name__=='__main__':
    import glob
    FixAndRotate = True
    fnames = sorted(glob.glob('bav02/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/*.vtk'))
    print(fnames,len(fnames))
    fdir = os.path.dirname(fnames[0])
    # Check Directory
    if not os.path.exists(fdir):
        print('Error: Path does not exist:', fdir)
        sys.exit()

    OF,CF=3,10
    WallArea, WallVol, LumenVol, Time, Pts, WallAreaRatio, WallVolRatio, LumenVolRatio, AvgJ, AvgI1, AvgJRatio, AvgI1Ratio, TotalMotion, N = ProcessData(flist=fnames,ref=fnames[0],FT=100.,OF=OF,CF=CF,opformat='vtk')
    print('Total Wall Area =',WallArea)
    print('Total Wall Volume =',WallVol)
    print('Total Lumen Volume =',LumenVol)
   
    ###################################
    # Save data
    print(fdir)
    DataLocation = 'Strains/'+fdir+'/Data.npz'
    if OF is not None and CF is not None:
        np.savez(DataLocation,Time=Time,Pts=Pts,WallArea=WallArea,WallVol=WallVol, LumenVol=LumenVol, WallAreaRatio=WallAreaRatio, WallVolRatio=WallVolRatio, LumenVolRatio=LumenVolRatio, N=N, OF=OF, CF=CF)
    else:
        np.savez(DataLocation,Time=Time,Pts=Pts,WallArea=WallArea,WallVol=WallVol, LumenVol=LumenVol, WallAreaRatio=WallAreaRatio, WallVolRatio=WallVolRatio, LumenVolRatio=LumenVolRatio, N=N)


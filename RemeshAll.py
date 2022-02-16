import vtk
import numpy as np
import os
from vtk.util.numpy_support import vtk_to_numpy

def Remesh_All(DataDir, original_filename):
    # Read the original vtk file
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(original_filename)
    
    reader.Update()
    pdi = reader.GetOutput()
    writer = vtk.vtkSTLWriter()
    writer.SetInputData(pdi)
    writer.SetFileName('test2.stl')
    writer.Write()
    
    ####
    os.system('/Applications/Coreform-Cubit-2021.11.app/Contents/MacOS/python3/Python.framework/Versions/3.9/bin/python3 cubit-remesh.py')
    os.system('meshio-convert test-remesh2.e test-remesh2.vtk')
    
    # Read the original vtk file
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(original_filename)
    
    reader.Update()
    pdi = reader.GetOutput()
    
    num_pts = pdi.GetNumberOfPoints()
    num_cells = pdi.GetNumberOfCells()
    new_pts = pdi.GetPoints()
    in_pd = pdi.GetPointData()
    
    #calculate surface normals
    NORMALS = vtk.vtkPolyDataNormals()
    NORMALS.SetInputData(pdi)
    NORMALS.Update()
    NORMALS.GetOutput().GetPointData().GetArray('Normals')
    normals = NORMALS.GetOutput().GetPointData().GetArray('Normals')
    
    #correct the thickness at the ends
    thick = vtk_to_numpy(in_pd.GetArray('Thickness'))
    thick[thick<0.31] = np.mean(thick[thick>0.31])
    id_map = {}
    
    #create an unstructured grid object for the extruded version of original vtk file
    #extrusion is necessary to avoid floating point errors between two surfaces leading to points not belonging to any cell
    ugo = vtk.vtkUnstructuredGrid()
    ugo.Allocate()
    out_pd = ugo.GetPointData()
    out_pd.CopyAllocate(in_pd);
    scale = 2.
    
    for i in range(num_pts):
        point = np.array(pdi.GetPoint(i))
        normal = np.array(normals.GetTuple3(i))
        j = new_pts.InsertNextPoint(point+normal*scale)
        new_pts.SetPoint(i,point-normal*scale)
        out_pd.CopyData(in_pd, i, i)
        out_pd.CopyData(in_pd, i, j)
        id_map[i] = j
    
    ids = vtk.vtkIdList()
    for c_id in range(num_cells):
        pdi.GetCellPoints(c_id, ids)
        new_ids = ids
        for id in range(ids.GetNumberOfIds()):
            extruded_pt_id = id_map[ids.GetId(id)]
            new_ids.InsertNextId(extruded_pt_id)
    
        ugo.InsertNextCell(vtk.VTK_WEDGE, new_ids)
    
    ugo.SetPoints(new_pts)
    
    #read the remeshed version
    reader2 = vtk.vtkUnstructuredGridReader()
    reader2.SetFileName('test-remesh2.vtk')
    reader2.Update()
    
    reader3 = vtk.vtkDataSetSurfaceFilter()
    reader3.SetInputData(reader2.GetOutput())
    reader3.Update()
    polydata2 = reader3.GetOutput()
    
    #Using ResampleWithDataSet map the data from extruded original mesh onto remeshed mesh
    resample = vtk.vtkResampleWithDataSet()
    resample.AddInputData(polydata2)
    resample.SetComputeTolerance(False)
    resample.SetTolerance(0.01)
    resample.SetPassFieldArrays(True)
    resample.SetPassCellArrays(True)
    resample.SetPassPointArrays(True)
    resample.SetMarkBlankPointsAndCells(False)
    resample.SetSourceData(ugo)
    resample.Update()
    
    #write the output
    path = os.path.join('Remeshed/',DataDir)
    if not os.path.exists(path):
      # Create a new directory because it does not exist 
      os.mkdir(path)
      
    NewFilename = os.path.join('Remeshed/',DataDir,os.path.split(os.path.splitext(original_filename)[0])[1] + '.vtk')
    RemeshFilename = os.path.split(os.path.splitext(original_filename)[0])[1] + '.vtk'
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(resample.GetOutput())
    writer.SetFileName(NewFilename)
    writer.Write()
    
    #extrude the resampled data with thickness information
    pdi = resample.GetOutput()
    num_pts = pdi.GetNumberOfPoints()
    num_cells = pdi.GetNumberOfCells()
    new_pts = pdi.GetPoints()
    in_pd = pdi.GetPointData()
    
    NORMALS = vtk.vtkPolyDataNormals()
    NORMALS.SetInputData(pdi)
    NORMALS.Update()
    NORMALS.GetOutput().GetPointData().GetArray('Normals')
    normals = NORMALS.GetOutput().GetPointData().GetArray('Normals')
    thick = vtk_to_numpy(in_pd.GetArray('Thickness'))
    id_map = {}
    
    ugo = vtk.vtkUnstructuredGrid()
    ugo.Allocate()
    out_pd = ugo.GetPointData()
    out_pd.CopyAllocate(in_pd);
    
    for i in range(num_pts):
        point = np.array(pdi.GetPoint(i))
        normal = np.array(normals.GetTuple3(i))
        j = new_pts.InsertNextPoint(point+normal*thick[i]/2.)
        new_pts.SetPoint(i,point-normal*thick[i]/2.)
        out_pd.CopyData(in_pd, i, i)
        out_pd.CopyData(in_pd, i, j)
        id_map[i] = j
    
    ids = vtk.vtkIdList()
    for c_id in range(num_cells):
        pdi.GetCellPoints(c_id, ids)
        new_ids = ids
        for id in range(ids.GetNumberOfIds()):
            extruded_pt_id = id_map[ids.GetId(id)]
            new_ids.InsertNextId(extruded_pt_id)
        if new_ids.GetNumberOfIds()==6:
            ugo.InsertNextCell(vtk.VTK_WEDGE, new_ids)
        elif new_ids.GetNumberOfIds()==8:
            ugo.InsertNextCell(vtk.VTK_HEXAHEDRON, new_ids)
    
    ugo.SetPoints(new_pts)
    
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetInputData(ugo)
    writer.SetFileName('extrude.vtk')
    writer.Write()
    
    return RemeshFilename


import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

reader = vtk.vtkPolyDataReader()
reader.SetFileName('resample.vtk')
nlayers = 4 #number of element across the thickness

reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()
pdi = reader.GetOutput()
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
    new_point = point-normal*thick[i]/2.
    new_pts.SetPoint(i,new_point)
    out_pd.CopyData(in_pd, i, i)
    jold = i 
    for k in range(nlayers):
        j = new_pts.InsertNextPoint(new_point+(k+1)*normal*thick[i]/(nlayers))
        out_pd.CopyData(in_pd, i, j)
        id_map[jold] = j
        jold = j

ids = vtk.vtkIdList()
for c_id in range(num_cells):
    pdi.GetCellPoints(c_id, ids)
    for k in range(nlayers):
        if k>0:
            ids = next_ids
        new_ids = ids
        next_ids = vtk.vtkIdList()
        for id in range(ids.GetNumberOfIds()):
            extruded_pt_id = id_map[ids.GetId(id)]
            new_ids.InsertNextId(extruded_pt_id)
            next_ids.InsertNextId(extruded_pt_id)
        if new_ids.GetNumberOfIds()==6:
            ugo.InsertNextCell(vtk.VTK_WEDGE, new_ids)
        elif new_ids.GetNumberOfIds()==8:
            ugo.InsertNextCell(vtk.VTK_HEXAHEDRON, new_ids)

ugo.SetPoints(new_pts)

writer = vtk.vtkUnstructuredGridWriter()
writer.SetInputData(ugo)
writer.SetFileName('extrude2.vtk')
writer.Write()

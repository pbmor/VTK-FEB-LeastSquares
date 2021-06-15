import numpy as np
import os
import sys
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa
import vtkplotlib as vpl

fname = 'tav02_root_med-prop-recon_f4.vtk'

# reader for vtk data
reader = vtk.vtkPolyDataReader()
colors = vtk.vtkNamedColors()

# Read the source file.
reader.SetFileName(fname)
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()
polydata = reader.GetOutput()

print('Polydata Summary:',polydata)

#Mesh Information
NP = polydata.GetNumberOfPoints()
NC = polydata.GetNumberOfCells()
Ranges = np.array(polydata.GetPoints().GetBounds())
print('Number of Points = ',NP)
print('Number of Cells = ',NC)
print('Coordinate Ranges = ',Ranges)
print('[Xmin,Xmax] = ',Ranges[0:2])
print('[Ymin,Ymax] = ',Ranges[0:2])
print('[Zmin,Zmax] = ',Ranges[0:2])

#Get the array of points
Points = vtk_to_numpy(polydata.GetPoints().GetData())

#Get array data
if not polydata.GetPointData().GetArray(0):
    print('Warning: No saved data')
else:
    NumOfArr = polydata.GetPointData().GetNumberOfArrays()
    Data = np.zeros((NumOfArr,NP))
    print('List of Array Data:')
    for i in range(NumOfArr):
        ArrayName = polydata.GetPointData().GetArrayName(i)
        print(ArrayName)
        if i in range(1,11):
            Data[i,:] = vtk_to_numpy(polydata.GetPointData().GetArray(i))

#Ids of interest for Thickness, Mean Curvature, Jacobian, Motion
Ids = [1,6,8,9]

ScalarRange = polydata.GetScalarRange()



# Colour transfer function.
ctf = vtk.vtkColorTransferFunction()
ctf.SetColorSpaceToDiverging()
p1 = [0.0] + list(colors.GetColor3d('MidnightBlue'))
p2 = [1.0] + list(colors.GetColor3d('DarkOrange'))
ctf.AddRGBPoint(*p1)
ctf.AddRGBPoint(*p2)
cc = list()
for i in range(256):
    cc.append(ctf.GetColor(float(i) / 255.0))

# Lookup table
lut = vtk.vtkLookupTable()
lut.SetNumberOfColors(256)
for i, item in enumerate(cc):
    lut.SetTableValue(i, item[0], item[1], item[2], 1.0)
lut.SetRange(np.amin(Data[1,:]),np.amax(Data[1,:]))
lut.Build()


Cells = np.zeros((NC,3,3))
TD = np.zeros((NC,3))
for i in range(NC):
    Cells[i,:,:] = vtk_to_numpy(polydata.GetCell(i).GetPoints().GetData())
    #Assign point id
    t = np.array([polydata.GetCell(i).GetPointId(j) for j in range(3)])
    #For each cell assign thickness for each point
    TD[i,:] = [Data[1,t[j]] for j in range(3)]

plot = vpl.plots.MeshPlot.MeshPlot(Cells)
vpl.show()

plot = vpl.plots.MeshPlot.MeshPlot(Cells, tri_scalars=None, scalars=TD, color=None, opacity=None, cmap=lut, fig='gcf', label=None)

# Optionally the plot created by mesh_plot can be passed to color_bar
vpl.color_bar(plot, "Thickness")

vpl.show()

#Find Mid point of the mesh
Mids = np.zeros((3))
Mids[0] = np.mean(Ranges[0:2])
Mids[1] = np.mean(Ranges[2:4])
Mids[2] = np.mean(Ranges[4:6])

#Redefine the point locations
for i in range(NP):
    ptNew = np.array([Points[i,j] - Mids[j] for j in range(3)])
    polydata.GetPoints().SetPoint(i, ptNew)


#Create New Cell data 
NewCellData = np.random.randn(NC)
# Add Cell Data
CellData  = [NewCellData]
CellNames = ['NewCellData']
for i in range(len(CellNames)) :
    arrayCell = vtk.util.numpy_support.numpy_to_vtk(CellData[i], deep=True)
    arrayCell.SetName(CellNames[i])
    dataCells = polydata.GetCellData()
    dataCells.AddArray(arrayCell)

# Convert Cell Data to Point Data
c2p = vtk.vtkCellDataToPointData()
c2p.AddInputData(polydata)
c2p.Update()
c2p.GetOutput()

#Create new point data
NumArr = polydata.GetPointData().GetNumberOfArrays()
Cell2PointData = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(NumArr))
    
PointData = np.random.randn(NP)
# Add Point Data
PointData = [PointData, Cell2PointData]
PointNames = ['NewPointData','Cell2PointData']
for i in range(len(PointNames)) :
    arrayPoint = vtk.util.numpy_support.numpy_to_vtk(PointData[i], deep=True)
    arrayPoint.SetName(PointNames[i])
    dataPoints = polydata.GetPointData()
    dataPoints.AddArray(arrayPoint)
    dataPoints.Modified()

#Create New Cell data
VectorData = np.random.randn(NC,3)

# Add Vector Data
VectorData = [VectorData]
VectorNames = ['NewVectorData']
for i in range(len(VectorNames)) :
    arrayVector = vtk.util.numpy_support.numpy_to_vtk(VectorData[i], deep=True)
    arrayVector.SetName(VectorNames[i])
    dataVectors = polydata.GetPointData()
    dataVectors.AddArray(arrayVector)
    dataVectors.Modified()


#Save new data
#Choose format
opformat = 'vtp'

#define writer and new filename
if opformat == 'vtp':
    NewFname = 'NewVTKSummaryData.vtp'
    writer = vtk.vtkXMLDataSetWriter()
elif opformat == 'vtk':
    NewFname = 'NewVTKSummaryData.vtk'
    writer = vtk.vtkDataSetWriter()
        
writer.SetFileName(NewFname)
writer.SetInputData(polydata)
print('Writing',NewFname)
writer.Write()

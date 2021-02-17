import numpy as np
import pyevtk
import os
import sys
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

reader = vtk.vtkPolyDataReader()
map    = vtk.vtkPolyDataMapper()
actors = vtk.vtkActor()
colors = vtk.vtkNamedColors()
points = vtk.vtkPoints()
vertices = vtk.vtkCellArray()

os.chdir('bav02/medial meshes')

TF = 0

List = list(range(1,14))

for X in List :
	Filename = 'bav02_root_med_f' + str(X) + '.vtk'
	print(Filename)

	# Read the source file.
	reader.SetFileName(Filename)
	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.Update()

	output = reader.GetOutput()
	thickness = output.GetPointData().GetArray(1)

	dataset = output.GetPointData()


	data = np.asarray(thickness)
	#print(data.shape)
	#data = np.atleast_2d(data)
	nPoints = output.GetNumberOfPoints()
	#nRows, nCols = data.shape
	#assert(nPoints==nRows)
	array = vtk.vtkDoubleArray()
	array.SetName("Thick2")

	# For a one-dimensional array.
	array.SetNumberOfValues(900)
	for x in zip(range(900), data.T):
		array.SetValue(*x)
	
	dataset.AddArray(array)


	filename = 'bav02_New' + str(X) + '.vtp'	
	writer = vtk.vtkPolyDataWriter()
	writer = vtk.vtkXMLDataSetWriter()
	writer.SetFileName(filename)
	writer.SetInputData(output)
	writer.Write()





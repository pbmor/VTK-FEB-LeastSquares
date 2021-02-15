import numpy as np
import pyevtk
import os
import sys
import vtk
from vtk.util.numpy_support import vtk_to_numpy


#import matplotlib.pyplot as plt
#import matplotlib

#from matplotlib.pyplot import figure
#from mpl_toolkits.mplot3d import Axes3D

reader = vtk.vtkPolyDataReader()
map    = vtk.vtkPolyDataMapper()
actors = vtk.vtkActor()
colors = vtk.vtkNamedColors()
points = vtk.vtkPoints()
vertices = vtk.vtkCellArray()

os.chdir('bav02/medial meshes')

ctf = vtk.vtkColorTransferFunction()
ctf.SetColorSpaceToDiverging()
p1 = [0.0] + list(colors.GetColor3d("MidnightBlue"))
p2 = [1.0] + list(colors.GetColor3d("DarkOrange"))
ctf.AddRGBPoint(*p1)
ctf.AddRGBPoint(*p2)
cc = list()
for i in range(256):
   cc.append(ctf.GetColor(float(i) / 255.0))

# Lookup table.
lut = list()
lut = vtk.vtkLookupTable()
lut.SetNumberOfColors(256)
for i, item in enumerate(cc):
   lut.SetTableValue(i, item[0], item[1], item[2], 1.0)
lut.SetRange(-0.1, 0.1)
lut.Build()



TF = 0

List = list(range(1,2))

for X in List :
	Filename = 'bav02_root_med_f' + str(X) + '.vtk'
	print(Filename)

	reader.SetFileName(Filename)
	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.Update()

	polydata = reader.GetOutput()
	thickness = polydata.GetPointData().GetArray(1)

	npts = polydata.GetNumberOfPoints()
	ncls = polydata.GetNumberOfCells()
	#nptA = polydata.GetPointData().GetNumberOfPointArrays()

	print(npts)
	print(ncls)
	#print(nptA)

	points = polydata.GetPoints()
	triangles=  polydata.GetPolys().GetData()

	curv = vtk.vtkCurvatures()
	curv.SetInputConnection(reader.GetOutputPort())
	curv.SetCurvatureTypeToMean()
	curv.Update()
	#curvatures.SetCurvatureTypeToMean()
	pd = curv.GetOutput()

	#print(np.array(pd))	

	pd_port = reader.GetOutputPort()
	scalar_range = polydata.GetScalarRange()

	# Create the mapper that corresponds the objects of the vtk file
	# into graphics elements
	mapper = vtk.vtkDataSetMapper()
	mapper.SetInputConnection(curv.GetOutputPort())
	mapper.SetScalarRange(scalar_range)
	mapper.SetLookupTable(lut)
	mapper.SetUseLookupTableScalarRange(1)

	# Create the Actor
	actor = vtk.vtkActor()
	actor.SetMapper(mapper)
	actor.GetProperty().SetColor(1,0.2,0.2)

	# Create the Renderer
	renderer = vtk.vtkRenderer()
	renderer.AddActor(actor)
	renderer.SetBackground(1, 1, 1) # Set background to white

	# Create the RendererWindow
	renderer_window = vtk.vtkRenderWindow()
	renderer_window.AddRenderer(renderer)
	renderer_window.SetSize(1000,1000)


	# Create the RendererWindowInteractor and display the vtk_file
	interactor = vtk.vtkRenderWindowInteractor()
	interactor.SetRenderWindow(renderer_window)
	interactor.Initialize()
	interactor.Start()



	# Data info
	ThickData = np.array(thickness)
	
	#Check object attributes

	#Create empty arrays
	Cells = np.zeros((polydata.GetNumberOfCells(),3,3))
	RA = np.zeros((polydata.GetNumberOfCells(),3,2))
	DA = np.zeros((polydata.GetNumberOfCells(),3,3))
	G = np.zeros((polydata.GetNumberOfCells()))
	
	Thicks = np.zeros((polydata.GetNumberOfCells(),3))
	thi = np.zeros((polydata.GetNumberOfCells(),3))

	for i in range(polydata.GetNumberOfCells()):
		pts = polydata.GetCell(i).GetPoints()
		#thi[i,:] = np.array([pts.GetArray(1) for i in range(pts.GetNumberOfPoints())])

		#Assign three points to each cell
		Cells[i,:,:]  = np.array(np.array(pts.GetData()))
		#Assign point id
		t = np.array([polydata.GetCell(i).GetPointId(j) for j in range(3)])
		#For each cell assign thickness for each point
		for j in range(3) :
			Thicks[i,j] = ThickData[t[j]] 

	
	if X == 1:
		for i in range(polydata.GetNumberOfCells()):
				RA[i,:,0] = Cells[i,:,0]-Cells[i,:,1]
				RA[i,:,1] = Cells[i,:,0]-Cells[i,:,2]
				G[i] = np.dot(RA[i,:,0],RA[i,:,1])

	#Sum Cell Volumes
	Volume = 0
	for i in range(polydata.GetNumberOfCells()):
		CA = polydata.GetCell(i).ComputeArea()
		CV = CA*(np.sum(Thicks[i,:]))/3
		Volume = Volume +CV

	print('Volume = ' + str(Volume))
	


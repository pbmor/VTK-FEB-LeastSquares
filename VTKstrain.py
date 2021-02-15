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

	reader.SetFileName(Filename)
	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.Update()

	polydata = reader.GetOutput()
	thickness = polydata.GetPointData().GetArray(1)

	# Data info
	ThickData = np.array(thickness)
	
	#Create empty arrays
	Cells = np.zeros((polydata.GetNumberOfCells(),3,3))
	if X==1 :
		RA = np.zeros((polydata.GetNumberOfCells(),3,3))
		Ra = np.zeros((polydata.GetNumberOfCells(),3,3))
		I1 = np.zeros((13,polydata.GetNumberOfCells()))
		I2 = np.zeros((13,polydata.GetNumberOfCells()))
		I1alt = np.zeros((13,polydata.GetNumberOfCells()))
		J = np.zeros((13,polydata.GetNumberOfCells()))
		Jalt = np.zeros((13,polydata.GetNumberOfCells()))
		TotalVolume = np.zeros((13))
		TotalArea = np.zeros((13))
	
	Thicks = np.zeros((polydata.GetNumberOfCells(),3))

	for i in range(polydata.GetNumberOfCells()):
		pts = polydata.GetCell(i).GetPoints()
		#Assign three points to each cell
		Cells[i,:,:]  = np.array(np.array(pts.GetData()))
		#Assign point id
		t = np.array([polydata.GetCell(i).GetPointId(j) for j in range(3)])
		#For each cell assign thickness for each point
		for j in range(3) :
			Thicks[i,j] = ThickData[t[j]] 

	#Sum Cell Areas and Volumes
	for i in range(polydata.GetNumberOfCells()):
		CA = polydata.GetCell(i).ComputeArea()
		CV = CA*(np.sum(Thicks[i,:]))/3
		TotalArea[X-1]   += CA
		TotalVolume[X-1] += CV

	#Calculate J 
	for i in range(polydata.GetNumberOfCells()):
		if X == 1:
			RA[i,:,0] = Cells[i,:,0]-Cells[i,:,2]
			RA[i,:,1] = Cells[i,:,1]-Cells[i,:,2]
		
		Ra[i,:,0] = Cells[i,:,0]-Cells[i,:,2]
		Ra[i,:,1] = Cells[i,:,1]-Cells[i,:,2]	

		# Define refernce vectors
		RA1 = RA[i,:,0]
		RA2 = RA[i,:,1]
		
		# Define deformed vectors
		Ra1 = Ra[i,:,0]
		Ra2 = Ra[i,:,1]

		# Define reference and deformed normal vectors
		nA = np.cross(RA[i,:,0],RA[i,:,1])
		mag = np.linalg.norm(nA)
		na = np.cross(Ra[i,:,0],Ra[i,:,1])
		mag = np.linalg.norm(na)

		# Define stacked arrays and relevant inverse
		A = np.array([RA1,RA2,nA/mag])
		Ainv = np.linalg.inv(A)
		a = np.array([Ra1,Ra2,na/mag])

		
		F = np.dot(a,Ainv)
		FT = F.transpose()
		C = np.dot(F.T,F)
		C2 = np.dot(C,C)
		
		trC        = C.trace()
		trC2       = C2.trace()
		I1[X-1,i]  = trC
		I2[X-1,i]  = np.sqrt(0.5*(trC**2 - trC2)) 		
		J[X-1,i]   = np.linalg.det(F)

		#Define vectors for alt F
		A1=RA1
		A2=RA2

		a1=Ra1
		a2=Ra2
		G=np.array([[np.dot(A1,A1),np.dot(A1,A2)],[np.dot(A1,A2),np.dot(A2,A2)]])
		invG = np.linalg.inv(G)

		A1dual = invG[0,0]*A1 + invG[0,1]*A2
		A2dual = invG[1,0]*A1 + invG[1,1]*A2
		Falt = np.outer(a1,A1dual) + np.outer(a2,A2dual)
		Calt = np.dot(Falt.T,Falt)
		Calt2= np.dot(Calt,Calt)

		trC         = Calt.trace()
		trC2        = Calt2.trace()
		I1[X-1,i]   = Calt.trace()
		Jalt[X-1,i] = np.sqrt((trC**2-trC2)/2)


	Jmin = J[X-1,:].min()	
	Jmax = J[X-1,:].max()
	


	#Define array of colours
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
	lut.SetRange(Jmin, Jmax)
	lut.Build()

	
	#Set up to plot and write vtk files, not working
        # Create the mapper that corresponds the objects of the vtk file
        # into graphics elements
	#mapper = vtk.vtkDataSetMapper()
	#print(dir(mapper))
	#mapper.SetInputConnection(curv.GetOutputPort())
	##mapper.SetScalarRange(scalar_range)
	#mapper.SetLookupTable(lut)
	#mapper.SetUseLookupTableScalarRange(1)

        # Create the Actor
	#actor = vtk.vtkActor()
	#actor.SetMapper(mapper)
	#actor.GetProperty().SetColor(1,0.2,0.2)

        # Create the Renderer
	#renderer = vtk.vtkRenderer()
	#renderer.AddActor(actor)
	#renderer.SetBackground(1, 1, 1) # Set background to white

        # Create the RendererWindow
	#renderer_window = vtk.vtkRenderWindow()
	#renderer_window.AddRenderer(renderer)
	#renderer_window.SetSize(1000,1000)


	# Create the RendererWindowInteractor and display the vtk_file
	#interactor = vtk.vtkRenderWindowInteractor()
	#interactor.SetRenderWindow(renderer_window)
	#interactor.Initialize()
	#interactor.Start()

	# Add data set and write VTK file
	#meshNew = dsa.WrapDataObject(polydata)
	#print(dir(polydata.GetPointData()))
	#polydata.GetPointData().AddArray(J[X-1,:], )
	#FilenameNew = 'bav02_root_med_f' + str(X) + 'New.vtk'
	#print(FilenameNew)
	
	#writer = vtk.vtkPolyDataWriter()

	#writer.SetFileName(FilenameNew)
	#writer.SetInputData(meshNew.VTKObject)
	#writer.Write()

print(' ')
print(' ')
print('Total Volumes = ' + str(TotalVolume))
print('Total Areas = ' + str(TotalArea))

print(' ')
print(' ')
print('J = ' + str(J))

print(' ')
print(' ')
print('Alternative J = ' + str(Jalt))

print(' ')
print(' ')
print('Diff in Js = ' + str(J-Jalt))



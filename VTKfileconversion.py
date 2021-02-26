import numpy as np
import pyevtk
import os
import sys
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
	

reader = vtk.vtkPolyDataReader()

fdir = 'bav02/medial meshes/'


# Check Directory
if not os.path.exists(fdir):
    print('Error: Path does not exist:', fdir)
    sys.exit()


#Choose reference frame
RF = 1
#################################
# Define Reference Vectors
#################################
Filename = 'bav02_root_med_f' + str(RF) + '.vtk'
Fname = os.path.join(fdir, Filename)
print('Reading:', Fname)

# Read the source file.
reader.SetFileName(Fname)
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

polydata = reader.GetOutput()
dataset = polydata.GetPointData()
NC = polydata.GetNumberOfCells()
RA = np.zeros((NC,2,3))

for i in range(polydata.GetNumberOfCells()):
	Cells = np.array(polydata.GetCell(i).GetPoints().GetData())
	RA[i,0,:] = Cells[2,:]-Cells[0,:]
	RA[i,1,:] = Cells[2,:]-Cells[1,:]


###########################################
# Define Data At Every Time Frame
##########################################

for X in list(range(1,14)) :
	Filename = 'bav02_root_med_f' + str(X) + '.vtk'
	Fname = os.path.join(fdir, Filename)
	print('Reading:', Fname)

	# Read the source file.
	reader.SetFileName(Fname)
	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.Update()

	##############################
	# Define Properties: Areas, Volumes, J

	polydata = reader.GetOutput()
	thickness = polydata.GetPointData().GetArray(1)
	ThickData = np.array(thickness)

	#Create empty arrays
	NC = polydata.GetNumberOfCells()
	I1 = np.zeros((NC))
	I2 = np.zeros((NC))
	I1alt = np.zeros((NC))
	I2alt = np.zeros((NC))
	J = np.zeros((NC))
	Jalt = np.zeros((NC))
	TotalVolume = 0
	TotalArea = 0
	CA = np.zeros((NC))
	CV = np.zeros((NC))

	for i in range(NC):
		#Assign point id
		t = np.array([polydata.GetCell(i).GetPointId(j) for j in range(3)])
		#For each cell assign thickness for each point
		Thicks = [ThickData[t[j]] for j in range(3)]
		#Sum Cell Areas and Volumes
		CA[i] = polydata.GetCell(i).ComputeArea()
		CV[i] = CA[i]*(np.sum(Thicks))/3
		TotalArea   += CA[i]
		TotalVolume += CV[i]


	#Calculate J
	for i in range(NC):
		Cells = np.array(polydata.GetCell(i).GetPoints().GetData())

		# Define Reference Vectors
		A1 = RA[i,0,:]
		A2 = RA[i,1,:]
		
		#Define Deformed Vectors
		a1 = Cells[2,:]-Cells[0,:]
		a2 = Cells[2,:]-Cells[1,:]

		# Define reference and deformed normal vectors
		nA = np.cross(A1,A2)
		Mag = np.linalg.norm(nA)
		na = np.cross(a1,a2)
		mag = np.linalg.norm(na)

                # Define stacked arrays and relevant inverse
		A = np.array([A1,A2,nA/Mag])
		Ainv = np.linalg.inv(A)
		a = np.array([a1,a2,na/mag])

		F = np.dot(a,Ainv)
		FT = F.transpose()
		C = np.dot(F.T,F)
		C2 = np.dot(C,C)

		trC    = C.trace()
		trC2   = C2.trace()
		I1[i]  = trC
		I2[i]  = np.sqrt(0.5*(trC**2 - trC2))
		J[i]   = np.linalg.det(F)

		G      = np.array([[np.dot(A1,A1),np.dot(A1,A2)],[np.dot(A2,A1),np.dot(A2,A2)]])
		invG   = np.linalg.inv(G)
		
		A1dual = invG[0,0]*A1 + invG[0,1]*A2
		A2dual = invG[1,0]*A1 + invG[1,1]*A2
		Falt   = np.outer(a1,A1dual) + np.outer(a2,A2dual)
		Calt   = np.dot(Falt.T,Falt)
		Calt2  = np.dot(Calt,Calt)

		trC         = Calt.trace()
		trC2        = Calt2.trace()
		I1alt[i]    = trC
		Jalt[i]     = np.sqrt((trC**2-trC2)/2)

	#############################
	# Define Curvature

	# Define Gaussian Curvature
	curvgauss = vtk.vtkCurvatures()
	curvgauss.SetInputConnection(reader.GetOutputPort())
	curvgauss.SetCurvatureTypeToGaussian()
	curvgauss.Update()
	CurvG = curvgauss.GetOutput()
	CurvG.GetPointData().GetScalars()
	CDG = np.asarray(CurvG.GetPointData().GetArray(5))	


	# Define Mean Curvature
	curvmean = vtk.vtkCurvatures()
	curvmean.SetInputConnection(reader.GetOutputPort())
	curvmean.SetCurvatureTypeToMean()
	curvmean.Update()
	CurvM = curvmean.GetOutput()
	CurvM.GetPointData().GetScalars()
	CDM = np.asarray(CurvM.GetPointData().GetArray(5))
        
	##########################
	# Add Data to Files and Write Files

	# Add Point Data
	PointData = [CDG,CDM]
	PointNames = ['CurvGaussian','CurvMean'] 

	for i in range(len(PointNames)) :
		arrayPoint = vtk.util.numpy_support.numpy_to_vtk(PointData[i], deep=True)
		arrayPoint.SetName(PointNames[i])
		dataPoints = polydata.GetPointData()
		dataPoints.AddArray(arrayPoint)
		dataPoints.Modified()


	# Add Cell Data
	CellData = [CA,CV,I1,I2,J,Jalt]
	CellNames = ['CellArea','CellVolume','I1','I2','J','Jalt'] 

	for i in range(len(CellNames)) :
		arrayCell = vtk.util.numpy_support.numpy_to_vtk(CellData[i], deep=True)
		arrayCell.SetName(CellNames[i])
		dataCells = polydata.GetCellData()
		dataCells.AddArray(arrayCell)
		dataCells.Modified()


	# Add Field Data
	FieldData = [TotalArea, TotalVolume]
	FieldNames = ['TotalArea','TotalVolume'] 

	for i in range(len(FieldNames)) :
		arrayField = vtk.util.numpy_support.numpy_to_vtk(FieldData[i], deep=True)
		arrayField.SetName(FieldNames[i])
		dataFields = polydata.GetFieldData()
		dataFields.AddArray(arrayField)
		dataFields.Modified() 



	filename = 'bav02_New' + str(X) + '.vtp'	
	fname = os.path.join(fdir, filename)
	writer = vtk.vtkXMLDataSetWriter()
	writer.SetFileName(filename)
	writer.SetInputData(polydata)
	print('Writing',filename)
	writer.Write()

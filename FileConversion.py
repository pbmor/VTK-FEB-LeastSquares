import numpy as np
import pyevtk
import os
import sys
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from matplotlib.pyplot import figure
from matplotlib import cm
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt


reader = vtk.vtkPolyDataReader()
FixAndRotate = False

fdir = 'bav02/medial meshes/'

# Find rotation matrix from vector x to vector y
def rot_mat(x,y):
	xMag = np.linalg.norm(x)
	x = x/xMag
	yMag = np.linalg.norm(y)
	y = y/yMag
	b = x+y
	bMag = np.linalg.norm(b)
	b = b/bMag
	R = (np.identity(3)-2*(b*b.T))*(np.identity(3)-2*(x*x.T))
	return R

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
NP = polydata.GetNumberOfPoints()
RA = np.zeros((NC,2,3))

for i in range(NC):
	Cells = np.array(polydata.GetCell(i).GetPoints().GetData())
	RA[i,0,:] = Cells[2,:]-Cells[0,:]
	RA[i,1,:] = Cells[2,:]-Cells[1,:]

Points = np.array(polydata.GetPoints().GetData())

###########################################
# Define Data At Every Time Frame
###########################################
Time = np.zeros(13)
TotArea = np.zeros(13)
TotVol  = np.zeros(13)
Pts     = np.zeros((13,NP,3))

for X in range(13):
	Filename = 'bav02_root_med_f' + str(X+1) + '.vtk'
	Fname = os.path.join(fdir, Filename)
	print('Reading:', Fname)

	# Read the source file.
	reader.SetFileName(Fname)
	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.Update()

	##############################
	# Define Properties: Areas, Volumes, J

	polydata  = reader.GetOutput()
	Time[X] = polydata.GetFieldData().GetMTime()
	thickness = polydata.GetPointData().GetArray(1)
	ThickData = np.array(thickness)

	#Create empty arrays
	Disp = np.zeros((NC,3))
	I1 = np.zeros((NC))
	J = np.zeros((NC))
	TotalVolume = 0
	TotalArea = 0
	CA = np.zeros((NC))
	CV = np.zeros((NC))

	Pts[X,:,:] = np.array(polydata.GetPoints().GetData())
	##############################
	# Redefine Points and Displacement vectors
	if FixAndRotate == True:
		Ranges = np.array(polydata.GetPoints().GetBounds())
		Mids = np.zeros((3))
		Mids[0] = np.mean(Ranges[0:2])
		Mids[1] = np.mean(Ranges[2:4])
		Mids[2] = np.mean(Ranges[4:6])

		for i in range(NP):
			pt = polydata.GetPoint(i)
			ptNew = [pt[j]-Mids[j] for j in range(3)]
			polydata.GetPoints().SetPoint(i, ptNew)
			Disp[i,:] = pt-Points[i,:]

		##############################
		# Rotate Mesh
	
		Points = np.array(polydata.GetPoints().GetData())
		Pointsb = Points.T
		Points2 = np.zeros((3,25,36))
		for i in range(3):
		 for j in range(25):
		  for k in range(36):
		   Points2[i,j,k] = Pointsb[i,((j)%25+k*25)%900]
	
		centre_line = np.zeros((3,25))
		for i in range(3):
			for j in range(25):
				centre_line[i,j] = np.mean(Points2[i,j,:])

		a   = centre_line[:,13]
		mag = np.linalg.norm(a)
		a   = a/mag	
		b   = [0,0,1]
		v   = np.cross(a,b)
		c   = np.dot(a,b)

		vx  = [[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]]
		vx2 = np.matmul(vx,vx)
	
		R   = np.identity(3)+vx+vx2*(1/(1+c))
		RPoints = np.matmul(R,Points.T)
		Points = RPoints.T

		RPoints2 = np.zeros((3,25,36))
		for i in range(3):
		 for j in range(25):
		  for k in range(36):
		   RPoints2[i,j,k] = RPoints[i,((j)%25+k*25)%900]


		for i in range(NP):
			ptNew = Points[i,:]
			polydata.GetPoints().SetPoint(i, ptNew)

	#######################################
	# Define Properties: Areas, Volumes, J

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

	TotArea[X] = TotalArea
	TotVol[X]  = TotalVolume

	#Calculate J
	for i in range(NC):
		Cells = np.array(polydata.GetCell(i).GetPoints().GetData())

		# Define Reference Vectors
		A1 = RA[i,0,:]
		A2 = RA[i,1,:]
		
		#Define Deformed Vectors
		a1 = Cells[2,:]-Cells[0,:]
		a2 = Cells[2,:]-Cells[1,:]

		G      = np.array([[np.dot(A1,A1),np.dot(A1,A2)],[np.dot(A2,A1),np.dot(A2,A2)]])
		invG   = np.linalg.inv(G)
		
		A1dual = invG[0,0]*A1 + invG[0,1]*A2
		A2dual = invG[1,0]*A1 + invG[1,1]*A2
		F   = np.outer(a1,A1dual) + np.outer(a2,A2dual)
		C   = np.dot(F.T,F)
		C2  = np.dot(C,C)

		trC         = C.trace()
		trC2        = C2.trace()
		I1[i]    = trC
		J[i]     = np.sqrt((trC**2-trC2)/2)

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
        
	#####################################
	# Add Data to Files and Write Files

	# Add Cell Data
	CellData  = [CA,CV,I1,J]
	CellNames = ['CellArea','CellVolume','I1','J'] 
	for i in range(len(CellNames)) :
		arrayCell = vtk.util.numpy_support.numpy_to_vtk(CellData[i], deep=True)
		arrayCell.SetName(CellNames[i])
		dataCells = polydata.GetCellData()
		dataCells.AddArray(arrayCell)

	# Add Point Data
	PointData = [CDG,CDM]
	PointNames = ['CurvGaussian','CurvMean'] 
	for i in range(len(PointNames)) :
		arrayPoint = vtk.util.numpy_support.numpy_to_vtk(PointData[i], deep=True)
		arrayPoint.SetName(PointNames[i])
		dataPoints = polydata.GetPointData()
		dataPoints.AddArray(arrayPoint)
		dataPoints.Modified()

	# Add Vector Data
	VectorData = [Disp]
	VectorNames = ['Displacement'] 
	for i in range(len(VectorNames)) :
		arrayVector = vtk.util.numpy_support.numpy_to_vtk(VectorData[i], deep=True)
		arrayVector.SetName(VectorNames[i])
		dataVectors = polydata.GetPointData()
		#dataVectors.SetNumberOfTuples(3)
		#print(dir(dataVectors))
		dataVectors.SetVectors(arrayVector)
		dataVectors.Modified()

	# Add Field Data
	FieldData = [TotalArea, TotalVolume]
	FieldNames = ['TotalArea','TotalVolume'] 
	for i in range(len(FieldNames)) :
		arrayField = vtk.util.numpy_support.numpy_to_vtk(FieldData[i], deep=True)
		arrayField.SetName(FieldNames[i])
		dataFields = polydata.GetFieldData()
		dataFields.AddArray(arrayField)
		dataFields.Modified() 

	#################################
	# Write data to vtp files
	filename = 'bav02_New' + str(X+1) + '.vtp'	
	fname = os.path.join(fdir, filename)
	writer = vtk.vtkXMLDataSetWriter()
	writer.SetFileName(filename)
	writer.SetInputData(polydata)
	print('Writing',filename)
	writer.Write()


###################################
# Save data

os.chdir('Data')
np.save('Time.npy',Time)
np.save('TotArea.npy',TotArea)
np.save('TotVol.npy',TotVol)
np.save('Points.npy',Pts)
os.chdir('..')


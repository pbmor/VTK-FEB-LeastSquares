import numpy as np
import os
import sys
from matplotlib.pyplot import figure
from matplotlib import cm
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

os.chdir('Data')
Time    = np.load('Time.npy')
TotArea = np.load('TotArea.npy')
TotVol  = np.load('TotVol.npy')
Pts  = np.load('Points.npy')
os.chdir('..')

[Nf,NP,Ndim] = Pts.shape
VAJMax    = np.zeros((13))
SVMax     = np.zeros((13))
STJMax    = np.zeros((13))
VAJMin    = np.zeros((13))
SVMin     = np.zeros((13))
STJMin    = np.zeros((13))
CrossArea = np.zeros((13,25))
RootVol   = np.zeros((13))
xdata     = np.zeros((13,25,36))
ydata     = np.zeros((13,25,36))
zdata     = np.zeros((13,25,36))

for X in range(Nf) :
	Points = Pts[X,:,:]
	##############################
	# Redefine Points and Displacement vectors
	Mids = np.zeros((3))
	Mids[0] = np.mean([np.amin(Points[:,0]),np.amax(Points[:,0])])
	Mids[1] = np.mean([np.amin(Points[:,1]),np.amax(Points[:,1])])
	Mids[2] = np.mean([np.amin(Points[:,2]),np.amax(Points[:,2])])

	for i in range(NP):
		Points[i,:] = [Points[i,j]-Mids[j] for j in range(3)]

	
	##############################
	# Rotate Mesh
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

	R = np.identity(3)+vx+vx2*(1/(1+c))
	RPoints = np.matmul(R,Points.T)
	Points = RPoints.T

	###################################
	# Define Coordinates data to be saved to .npy file
	for i in range(25):
	 for j in range(36):
	  xdata[X,i,j] = Points[((13+i)%25+j*25)%900,0]
	  ydata[X,i,j] = Points[((13+i)%25+j*25)%900,1]
	  zdata[X,i,j] = Points[((13+i)%25+j*25)%900,2]

	xRing = np.zeros((25,37))
	yRing = np.zeros((25,37))

	for i in range(25):
		xRing[i,0:36] = xdata[X,i,:]
		yRing[i,0:36] = ydata[X,i,:]
		xRing[i,36]   = xRing[i,0]
		yRing[i,36]   = yRing[i,0]

	for i in range(25):
		for j in range(36):
			CrossArea[X,i] +=(xRing[i,j]*yRing[i,j-1] - xRing[i,j-1]*yRing[i,j])/2	

	for i in range(24):
		RootVol[X] += (np.mean(zdata[X,i,:])-np.mean(zdata[X,i+1,:]))*(CrossArea[X,i]+CrossArea[X,i+1])/2
	
	labelNames = ['VAJ','SV','STJ']
	Handles = [0,0,0]
	j=0
	cmap = cm.get_cmap('brg')
	Max_Colors = 25
	Color_N = 0
	figure(num=X)
	for i in range(25):	
		Color_N +=1
		Color = cmap(Color_N/Max_Colors)
		if i in [0,13,24]:
			plt.plot(xRing[24-i,:],yRing[24-i,:],color=Color,label=labelNames[j])
			j+=1	
	plt.legend()
	#plt.show()
	
	Max_Colors = 25
	Color_N = 0
	figure(num=X)
	for i in range(25):	
		Color_N +=1
		Color = cmap(Color_N/Max_Colors)
		plt.plot(xRing[24-i,:],yRing[24-i,:],color=Color)


	Radii = np.zeros((3,36))
	k=-1
	for j in [0,13,24]:
		k += 1
		for i in range(36):
			Radii[k,i] = np.sqrt(xdata[X,24-j,i]**2+ydata[X,24-j,i]**2)

	Title='Time ='+str(Time[X]) +'(ms)'
	plt.title(Title,size=20)
	plt.xlabel('x  (mm)',size=14)
	plt.ylabel('y (mm)',size=14)

		
	VAJMax[X]    = np.amax(Radii[0,:])
	VAJMin[X]    = np.amin(Radii[0,:])
	SVMax[X]     = np.amax(Radii[1,:])
	SVMin[X]     = np.amin(Radii[1,:])
	STJMax[X]    = np.amax(Radii[2,:])
	STJMin[X]    = np.amin(Radii[2,:])
	
	cmap = cm.get_cmap('winter')
	if X in [1]:
	 Max_Colors = 13
	 Color_N = 0

	Color = cmap(Color_N/Max_Colors)
	Color_N +=1

figure(num=14)
plt.plot(Time,VAJMax,'r')
plt.plot(Time,VAJMin,'g')
plt.ylim(10,20)
plt.title('VAJ',size=20)
plt.xlabel('Time (ms)',size=14)
plt.ylabel('Radius (mm)',size=14)
figure(num=15)
plt.plot(Time,SVMax,'r')
plt.plot(Time,SVMin,'g')
plt.ylim(10,20)
plt.title('SV',size=20)
plt.xlabel('Time (ms)',size=14)
plt.ylabel('Radius (mm)',size=14)
figure(num=16)
plt.plot(Time,STJMax,'r',label='Max Radius')
plt.plot(Time,STJMin,'g',label='Min Radius')
plt.ylim(10,20)
plt.title('STJ',size=20)
plt.xlabel('Time (ms)',size=14)
plt.ylabel('Radius (mm)',size=14)
plt.legend()
figure(num=17)
plt.plot(Time,TotArea)
plt.title('Root Surface Area',size=20)
plt.xlabel('Time (ms)',size=14)
plt.ylabel('Surface Area (mm$^2$)',size=14)
figure(num=18)
plt.plot(Time,RootVol)
plt.title('Root Volume',size=20)
plt.xlabel('Time (ms)',size=14)
plt.ylabel('Volume (mm$^3$)',size=14)

plt.show()

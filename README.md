# RunVTK
Run VTK files.

Currently three scripts are included.
  FileConversion.py
  GeometricAnalysis.py  
  bav02 data
  
FileConversion.py:
	Reads.vtk files and converts them to .vtp files. This script also adds
	Jacobians, I1, cells area and volumes, mean and gaussian curvatures, and 
	the total surface area and wall volume. The script also outputs the times 
	of each frame, the point coordiantes and the total wall area and volume values
	as .npy files. There is also the option to fix the mesh centre to be (0,0,0) 
	for each time frame and rotate the mesh to the z-axis.

GeometricAnalysis.py:
	Reads the .npy files saved in the 'Data' folder by the FileConversion.py script.
	This script fixes the mesh centre to be (0,0,0)	for each time frame 
	and rotate the mesh to the z-axis. This allows the circumferential cross-sections
	to be plotted, as well as calculate the radius of these cross-sections and 
	estimate the lumen volume. This script also plots the the max and min radii 
	of the STJ, SV and VAJ cross-sections against time (although this can be done 
	for all cross-sections), as well as  the total surface areas and volumes against time


The bav02 data is also included again, so that the scripts should run as written.
 
 

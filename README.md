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
	as .npy files. There is also the option to fix the mesh centre to be the centre
    of the reference frame.

process-all.py:
    This script runs through all the data, filtering out the excluded data sets and
    saving the new vtp files in the Strains directory.


WiggerAreaVol.py:
    Creates a Wigger plot of the root surface area, root wall volume, and lumen volume.

AnalysisFuncs.py:
    A script of functions that are useful, some of which are no longer used. 

VTKSummary:
    a python script and pdf that offer an intro to using vtk data

echoframetime.csv:
    csv files of information about the data sets, importantly, the time of opening and closing and if the data set is included

The bav02 data is also included, which is used when FileConversion.py is run

 
 

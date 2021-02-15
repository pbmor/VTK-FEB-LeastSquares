# RunVTK
Run VTK files.

Currently three scripts are included.
  VTKcurv.py
  VTKcurvmean.py
  VTKstrain.py
  bav02 data
  
VTKcurv.py and VTKcurvmean.py are in essence identical, 
      with VTKcurv.py finding and plotting the gaussian curvature, 
      and VTKmean.py uses the mean curvature. Note the different 
      "SetRange" in the lookup table
      
VTKstrain.py finds the Jacobian of each cell at each time frame with
      both methods, which are stored in the arrays J and Jalt.
      VTKstrain also find the total area and volume at each time frame
      which are stored in the Array TotalArea and TotalVolume, respectively
      
 The bav02 data is also included again, so that the scripts should run as written.
 
 

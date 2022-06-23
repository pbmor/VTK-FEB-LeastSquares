
This repository contains the scripts:
    ResidualFunction_Time.py
    RunLeastSquares.py
    Read_XPLTfuncs.py
    vtk2feb_Disp.py
    WindKessel.py
    PressureProfile.csv
    TimeProfile.csv
    echoframetime.csv
    RemeshAll.py
    
Also included is the directory:
    Remeshed: the .vtp files of the remeshed of all the tav and bav data
    FEB_Files: The .feb and .xplt files of the Febio simulations
Note: that a new directory is NewFiles is created by SaveFiles, a function
    that is utilisied by ResidualFunction_Time and RunLeastSquares
    

ResidualFunction_Time.py:
    This script contains several functions:
        SaveFiles - A Function that saves the xplt output data to the remeshed 
                    vtk files. It also saves the output of the least squares 
                    optimisation, or the initial parameters if the 
                    optimisation isn't chosen. This function relies on GetFEB, 
                    GetMeshInfo
                    and GetData (from Read_XPLTfuncs).
        OrderList - A function to order the list of frames files, to start 
                    with a specific frame and end with the frames that were 
                    before the given frame.
        GetRes    - A function to find the residual of the remeshed file and 
                    the output from Febio, for a given parameter set
                    This function returns the flattened array of the residual
                    so that it can be used with scipy.optimize.least_squares. 
                    This function updates exisiting .feb files and relies on 
                    GetFEB, GetMeshInfo and GetData (from Read_XPLTfuncs).
        RunLS     - A function to create a .feb file for the chosen properties 
                    and initial paramter estimation, then set off the least 
                    squares optimisation. Note there is an option to only run
                    the Febio simulation and get the output for the initial 
                    parameter estimation. This function relies on VTK2FEB_Func
                    (from vtk2feb_Disp) and GetRes.
    If the ResidualFunction_Time is run directly it will run a default case
    with the tav02 case. This relies of GetFEB and GetMeshInfo and SaveFiles. 
    This produces an array Param that collects the fitted parameters.
    Several choices can be made:
        PressureChoice - Choose to optimise the pressure magnitude.
        ModelParChoice - Choose to vary modelparameters
        RunLSChoice    - Choose to run least Squares (or default/initial guess)
        ProfileChoice  - Choose profile shapes, options are: 'Triangle','Step',
                         'SmoothStep','Bio', 'Fourier','Fitted'.
                         Note: 'Bio' relies on the PressureProfile.csv and
                         TimeProfile.csv files
        ResChoice      - Choose type of residual calculation method: 'P2P', 
                         'CentreLine', 'CellPlane'
        ModelChoice    - Choose constitutive model from Mooney-Rivlin ('MR'), 
                         trans-isotropic Mooney Rivilin 'tiMR', Ogden ('Ogden') 
                         and Fung ('Fung')
                        
                        
RunLeastSquares.py:
    This script is set up to run through multiple cases, with the same options
    as the default choices in ResidualFunction_Time script.
    

Read_XPLTfuncs.py
    This script reads the .xplt files (the output of the febio simulations) 
    and retrieves the new data, specificially the displacement and stresses. 
    The three key functions are: 
        GetFeb which gets the data tree, 
        GetMeshData which gets fundamental output info, such as the number of 
        frames, elements, nodes, and also the frame times and the number of 
        variables, their names and types (the types dictate the format)
        GetData retrieves the physical data, e.g. the displacements and stresses
    
    
vtk2feb_Disp.py
    This script contains the function VTK2Feb_Func that constructs the .feb 
    file for Febio, determined by the choices made, specifically the pressure 
    profile and constitutive model.
    
    
vtk2feb_Test.py
    Create a test .feb file with a simple cylindrical mesh. You can choose the
    resolution of the mesh, the number of frames and time between frames.


WindKessel.py
    This script solves the windkessel equations to find the pressure profile
    for a given set of parameters
    
    
PressureProfile.csv and TimeProfile.csv
    The physiological data for the aortic root, found in the repository 
    http://haemod.uk/original
    
    
echoframetime.csv
    .csv file of the info for the different aortic root data sets, e.g. opening
    and closing frame, time per frame, number of frames
    

RemeshAll.py
    Contains the function Remesh_All that creates a remeshed version of the 
    medial meshes from the pipeline.
    Note: the remeshed files are included in the git and does not need to be 
          repeated.
    
    
    
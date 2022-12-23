import vtk
import numpy as np
from numpy import genfromtxt
from vtk.util.numpy_support import vtk_to_numpy
from math import exp as e
from math import cos, sin, pi
from scipy.integrate import odeint
from scipy.interpolate import interp1d

def VTK2Feb_Func(DataDir,Fname,C,nF,STJ_WallMotion,VAJ_WallMotion,STJ_Id,VAJ_Id,Circ,Long,Norm,PressureChoice,ProfileChoice,ModelChoice,FibChoice,CF):
    '''
    This script contains the function VTK2Feb_Func that constructs the .feb 
    file for Febio, determined by the choices made, specifically the pressure 
    profile and constitutive model.


    DataDir        - The data set name.
    Fname          - The filename of the frame that will define the mesh, 
                     usually chosen to be the reference frame.
    C              - Array of parameters
    nF             - The number of frames.
    STJ_WallMotion - The displacement of the STJ points relative to the root 
                     centre (i.e. overall root movement is not included).
    VAJ_WallMotion - The displacement of the VAJ points relative to the root 
                     centre (i.e. overall root movement is not included).
    STJ_Id         - The point ids of the STJ points.
    VAJ_Id         - The point ids of the VAJ points.
    Circ           - An array of the circumferential vectors of each node.
    Long           - An array of the longitudinal vectors of each node.
    Norm           - An array of the normal vectors of each node.
    ProfileChoice  - The choice for pressure profile.
    ModelChoice    - The choice for model.
    FiberChoice    - The choice of whether the fibers are included.
    CF             - The frame in which the valve closes.

    '''
    # Count number of parameters
    nC = 0    
    
    # Choose initial estimate of pressure magnitude
    if PressureChoice:
        PressureMag = C[0]
        nC =+ 1
    else:
        PressureMag = 0.126
     
    # Choose initial model parameters
    if ModelChoice[0:3] != 'Set':
        if ModelChoice == 'MR':
            ModelPar = C[nC:nC+4]
            nC += 4
        if ModelChoice == 'tiMR':
            ModelPar = C[nC:nC+8]
            nC += 8
        elif ModelChoice == 'Ogden':
            ModelPar = C[nC:nC+14]
            nC += 14
        elif ModelChoice == 'Fung':
            ModelPar = C[nC:nC+12]
            nC += 12
        elif ModelChoice == 'HGO':     
            ModelPar = C[nC:nC+6]
    else:   
        if ModelChoice == 'MR':
            ModelPar = [1,10,10,10]
            nC +=4
        if ModelChoice == 'tiMR':
            ModelPar = [10,10,10,10,10,10,10,10]
            nC +=8
        elif ModelChoice == 'Ogden':
            ModelPar = [1,10,10,10,10,10,10,10,10,10,10,10,10,10]
            nC +=14
        elif ModelChoice == 'Fung':
            ModelPar = [1,1,2,1.5,0.5,0.1,1,1.2,2,1.5,0.5,0.1]
            nC+=12
        elif ModelChoice == 'HGO':     
            ModelPar = [ 5,  415,  5, 2]
            # ModelPar = [ 5,  415,  5, 28, 0.3, 10]
            nC+6
            
      
    # Define initial fiber angle
    # Note if the chosen model is not 'tiMR' this parameter array will be made to be empty
    if FibChoice and ModelChoice =='tiMR':
        FiberAng = C[nC]
        nC +=1
    else:
        FiberAng = [0]
        
    #Choose initial paramters for the pressure profile
    if ProfileChoice == 'Triangle':
        PressurePar = [0.5]
    if ProfileChoice == 'Step':
        PressurePar = [0.2, 0.5]
    if ProfileChoice == 'SmoothStep':
        PressurePar = [0.2, 0.5,50,50]
    if ProfileChoice == 'Virtual':
        PressurePar = []
    if ProfileChoice == 'Fourier':
        PressurePar = [3.0, 1.0, 0.0, 0.05]
    if ProfileChoice == 'Fitted':
        PressurePar = np.zeros(nF-2)
    if ProfileChoice == 'Windkessel':
        PressurePar = [11,1.4,14,0.004]
    
    
    TimeInterp = np.linspace(0,1,nF)
    
    if ProfileChoice == 'Triangle':
        Peak = PressurePar[0]
        PressureInterp = np.zeros(nF)
        Line_A = (1/Peak)*TimeInterp
        Line_B = (-1/(1 - Peak))*TimeInterp + 1/(1-Peak)
        for i in range(nF):
            if TimeInterp[i] <=Peak:
                PressureInterp[i] = Line_A[i]
            else:
                PressureInterp[i] = Line_B[i]
        PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
    elif ProfileChoice == 'Step':
        PressureInterp = np.zeros(nF)
        for i in range(nF):
            if (TimeInterp[i] >=PressurePar[0]) and (TimeInterp[i] <= PressurePar[1]):
                PressureInterp[i] = 1.0
        PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
    elif ProfileChoice == 'SmoothStep':
        PressureInterp = np.zeros(nF)
        for i in range(nF):
            PressureInterp[i] = (1/(1+e(-PressurePar[2]*(TimeInterp[i]-PressurePar[0]))))*(1-1/(1+e(-PressurePar[3]*(TimeInterp[i]-PressurePar[1]))))
        PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
    elif ProfileChoice == 'Virtual':
        TimeData = (genfromtxt('TimeProfile.csv', delimiter=','))
        PressureData = genfromtxt('PressureProfile.csv', delimiter=',')
        PressureInterp = np.interp(TimeInterp,TimeData,PressureData)
        PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
    elif ProfileChoice == 'Fourier':
        PressureInterp = np.zeros(nF)
        x=TimeInterp
        for i in range(nF):
            PressureInterp[i] = PressurePar[0]*sin(pi*x[i])+PressurePar[1]*sin(2*pi*x[i])+PressurePar[2]*sin(4*pi*x[i])+PressurePar[3]*sin(8*pi*x[i])
        PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
    elif ProfileChoice == 'Fitted':
        PressureInterp = 0.0*np.ones(nF)
    elif ProfileChoice == 'Windkessel':
        qi_mag = PressurePar[0]
        Rp = PressurePar[1]
        Rd = PressurePar[2]
        Cp  = PressurePar[3]
        t = np.linspace(0, 100, 5001)
        qi = []
        CFt = CF/nF
        for ti in t:
            if ti%1>0 and ti%1<=CFt:
                qi.append(qi_mag*sin(pi*(ti%1)/CFt))
            else:
                qi.append(0)
        qi = np.array(qi)
        qi_func = interp1d(t,qi, fill_value="extrapolate")
        def func(y,t):
            return (qi_func(t)-y)/Cp/Rd
        qo = odeint(func,0,t).flatten()
        last = (t>=99)*(t<100)
        Pressure = Rp*qi+Rd*qo
        PressureInterp = np.interp(np.linspace(99,100,nF),t[last],Pressure[last])
        
    elif ProfileChoice[0:3] == 'Set':
        if ProfileChoice[3:] == 'Triangle':
            Peak = 0.5
            PressureInterp = np.zeros(nF)
            Line_A = (1/Peak)*TimeInterp
            Line_B = (-1/(1 - Peak))*TimeInterp + 1/(1-Peak)
            for i in range(nF):
                if TimeInterp[i] <=Peak:
                    PressureInterp[i] = Line_A[i]
                else:
                    PressureInterp[i] = Line_B[i]
            PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
        elif ProfileChoice[3:] == 'Step':
            PressureInterp = np.zeros(nF)
            for i in range(nF):
                if (TimeInterp[i] >=0.2) and (TimeInterp[i] <= 0.5):
                    PressureInterp[i] = 1.0
            PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
        elif ProfileChoice[3:] == 'SmoothStep':
            PressureInterp = np.zeros(nF)
            for i in range(nF):
                PressureInterp[i] = (1/(1+e(-50*(TimeInterp[i]-0.2))))*(1-1/(1+e(-50*(TimeInterp[i]-0.5))))
            PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
        elif ProfileChoice[3:] == 'Virtual':
            TimeData = (genfromtxt('TimeProfile.csv', delimiter=','))
            PressureData = genfromtxt('PressureProfile.csv', delimiter=',')
            PressureInterp = np.interp(TimeInterp,TimeData,PressureData)
            PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
        elif ProfileChoice[3:] == 'Fourier':
            PressureInterp = np.zeros(nF)
            x=TimeInterp
            for i in range(nF):
                PressureInterp[i] = 3.0*sin(pi*x[i])+1.0*sin(2*pi*x[i])+0.0*sin(4*pi*x[i])+0.05*sin(8*pi*x[i])
            PressureInterp = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
        elif ProfileChoice[3:] == 'Fitted':
            PressureInterp = 0.0*np.ones(nF)
        elif ProfileChoice[3:] == 'Windkessel':
            qi_mag = 11.0
            Rp     = 1.4
            Rd     = 14
            Cp     = 0.004
            t = np.linspace(0, 100, 5001)
            qi = []
            CFt = CF/nF
            for ti in t:
                if ti%1>0 and ti%1<CFt:
                    qi.append(qi_mag*sin(pi*(ti%1)/CFt))
                    #qi.append(qi_mag)
                else:
                    qi.append(0)
            qi = np.array(qi)
            qi_func = interp1d(t,qi, fill_value="extrapolate")
            def func(y,t):
                return (qi_func(t)-y)/Cp/Rd
            qo = odeint(func,0,t).flatten()
            last = (t>=99)*(t<100)
            Pressure = Rp*qi+Rd*qo
            PressureInterp = np.interp(np.linspace(99,100,nF),t[last],Pressure[last])
    
    PressureInterp = np.subtract(PressureInterp,min(PressureInterp))
    PressureInterp = np.divide(PressureInterp,max(PressureInterp))
    
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(Fname)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    pdi = reader.GetOutput()
    num_pts = pdi.GetNumberOfPoints()
    num_cells = pdi.GetNumberOfCells()
    pts = pdi.GetPoints()
    in_pd = pdi.GetPointData()
    X = vtk_to_numpy(pts.GetData())
    thick = vtk_to_numpy(in_pd.GetArray('Thickness'))
    
    NewFebName = './FEB_Files/' + DataDir+'.feb'
    f=open(NewFebName,'w')
    
    STJ_Cent = [0,0,0]
    nSTJ = len(STJ_Id)
    for i in STJ_Id:
        STJ_Cent[0] += X[i,0]
        STJ_Cent[1] += X[i,1]
        STJ_Cent[2] += X[i,2]
    STJ_Cent[0] /= nSTJ
    STJ_Cent[1] /= nSTJ
    STJ_Cent[2] /= nSTJ
    
    VAJ_Cent = [0,0,0]
    nVAJ = len(VAJ_Id)
    for i in STJ_Id:
        VAJ_Cent[0] += X[i,0]
        VAJ_Cent[1] += X[i,1]
        VAJ_Cent[2] += X[i,2]
    VAJ_Cent[0] /= nVAJ
    VAJ_Cent[1] /= nVAJ
    VAJ_Cent[2] /= nVAJ

    Root_axis = np.subtract(STJ_Cent,VAJ_Cent)
    
    #Generic starting information
    f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
    f.write('<febio_spec version="3.0">\n')
    f.write('\t<Module type="solid"/>\n')
    
    #Solver settings
    f.write('\t<Control>\n')
    f.write('\t\t<analysis>STATIC</analysis>\n')
    f.write('\t\t<time_steps>'+str(nF)+'</time_steps>\n')
    f.write('\t\t<step_size>'+str(1/(nF-1))+'</step_size>\n')
    f.write('\t\t<plot_level>PLOT_MUST_POINTS</plot_level>\n')
    f.write('\t\t<solver>\n')
    f.write('\t\t\t<max_refs>15</max_refs>\n')
    f.write('\t\t\t<max_ups>10</max_ups>\n')
    f.write('\t\t\t<diverge_reform>1</diverge_reform>\n')
    f.write('\t\t\t<reform_each_time_step>1</reform_each_time_step>\n')
    f.write('\t\t\t<dtol>0.001</dtol>\n')
    f.write('\t\t\t<etol>0.01</etol>\n')
    f.write('\t\t\t<rtol>0</rtol>\n')
    f.write('\t\t\t<lstol>0.9</lstol>\n')
    f.write('\t\t\t<min_residual>1e-20</min_residual>\n')
    f.write('\t\t\t<qnmethod>BFGS</qnmethod>\n')
    f.write('\t\t\t<rhoi>0</rhoi>\n')
    f.write('\t\t</solver>\n')

    f.write('\t\t<time_stepper>\n')
    f.write('\t\t\t<dtmin>'+str(0.9/(nF-1))+'</dtmin>\n')
    # f.write('\t\t\t<dtmax>'+str(1.1/nF)+'</dtmax>\n')
    f.write('\t\t\t<dtmax lc="1">1</dtmax>\n')
    f.write('\t\t\t<max_retries>5</max_retries>\n')
    f.write('\t\t\t<opt_iter>10</opt_iter>\n')
    f.write('\t\t\t<aggressiveness>1</aggressiveness>\n')
    f.write('\t\t</time_stepper>\n')
    f.write('\t</Control>\n')
    
    f.write('\t<Globals>\n')
    f.write('\t\t<Constants>\n')
    f.write('\t\t\t<T>0</T>\n')
    f.write('\t\t\t<R>0</R>\n')
    f.write('\t\t\t<Fc>0</Fc>\n')
    f.write('\t\t</Constants>\n')
    f.write('\t</Globals>\n')
    
    #Material model
    f.write('\t<Material>\n')
    if ModelChoice == 'MR':
        f.write('\t\t<material id="1" name="Material1" type="Mooney-Rivlin">\n')
        f.write(f'\t\t\t<density>{ModelPar[0]}</density>\n')
        f.write(f'\t\t\t<c1>{ModelPar[1]}</c1>\n')
        f.write(f'\t\t\t<c2>{ModelPar[2]}</c2>\n')
        f.write(f'\t\t\t<k>{ModelPar[3]}</k>\n')
        f.write('\t\t\t<mat_axis type="angles">\n') 
        f.write('\t\t\t\t<theta>20</theta>\n')
        f.write('\t\t\t\t<phi>60</phi>\n')
        # f.write(f'\t\t\t\t<vector>{Circ[0,0]},{Circ[0,1]},{Circ[0,2]}</vector>\n')
        f.write('\t\t\t</mat_axis>\n') 
        f.write('\t\t</material>\n') 
        # for i in range(num_cells):
        #     f.write(f'\t\t<material id="{i+1}" name="Material{i+1}" type="Mooney-Rivlin">\n')
        #     f.write(f'\t\t\t<density>{ModelPar[0]}</density>\n')
        #     f.write(f'\t\t\t<c1>{ModelPar[1]}</c1>\n')
        #     f.write(f'\t\t\t<c2>{ModelPar[2]}</c2>\n')
        #     f.write(f'\t\t\t<k>{ModelPar[3]}</k>\n')
        #     f.write('\t\t\t<mat_axis type="vector">\n') 
        #     f.write(f'\t\t\t\t<a>{Circ[i,0]},{Circ[i,1]},{Circ[i,2]}</a>\n')
        #     f.write(f'\t\t\t\t<d>{Long[i,0]},{Long[i,1]},{Long[i,2]}</d>\n')
        #     f.write('\t\t\t</mat_axis>\n') 
        #     f.write('\t\t</material>\n') 
    elif ModelChoice =='tiMR':
        for i in range(num_cells):
            f.write(f'\t\t<material id="{i+1}" name="Material{i+1}" type="2D trans iso Mooney-Rivlin">\n') 
            f.write(f'\t\t\t<density>{ModelPar[0]}</density>\n') 
            f.write(f'\t\t\t<c1>{ModelPar[1]}</c1>\n') 
            f.write(f'\t\t\t<c2>{ModelPar[2]}</c2>\n') 
            f.write(f'\t\t\t<c3>{ModelPar[3]}</c3>\n') 
            f.write(f'\t\t\t<c4>{ModelPar[4]}</c4>\n') 
            f.write(f'\t\t\t<c5>{ModelPar[5]}</c5>\n') 
            f.write(f'\t\t\t<k>{ModelPar[6]}</k>\n') 
            f.write('\t\t\t<lam_max>1.0</lam_max>\n') 
            if FibChoice:
                theta = FiberAng
                C_prime = np.cross(np.cross(Circ[i],Long[i]),Circ[i])
                New = np.dot(cos(theta),Circ[i]) + np.dot(sin(theta),C_prime)
                f.write(f'\t\t\t<fiber type="vector">{New[0]},{New[1]},{New[2]}</fiber>\n') 
            else:
                Fib_X = Circ[i,0]
                Fib_Y = Circ[i,1]
                Fib_Z = Circ[i,2]
                f.write(f'\t\t\t<fiber type="vector">{Fib_X},{Fib_Y},{Fib_Z}</fiber>\n') 
            f.write('\t\t</material>\n') 
    elif ModelChoice =='Ogden':
        for i in range(num_cells):
            f.write(f'\t\t<material id="{i+1}" name="Material{i+1}" type="Ogden">\n')
            f.write(f'\t\t\t<density>{ModelPar[0]}<</density>\n')
            f.write(f'\t\t\t<k>{ModelPar[1]}</k>\n')
            f.write(f'\t\t\t<c1>{ModelPar[2]}</c1>\n')
            f.write(f'\t\t\t<c2>{ModelPar[3]}</c2>\n')
            f.write(f'\t\t\t<c3>{ModelPar[4]}</c3>\n')
            f.write(f'\t\t\t<c4>{ModelPar[5]}</c4>\n')
            f.write(f'\t\t\t<c5>{ModelPar[6]}</c5>\n')
            f.write(f'\t\t\t<c6>{ModelPar[7]}</c6>\n')
            f.write(f'\t\t\t<m1>{ModelPar[8]}</m1>\n')
            f.write(f'\t\t\t<m2>{ModelPar[9]}</m2>\n')
            f.write(f'\t\t\t<m3>{ModelPar[10]}</m3>\n')
            f.write(f'\t\t\t<m4>{ModelPar[11]}</m4>\n')
            f.write(f'\t\t\t<m5>{ModelPar[12]}</m5>\n')
            f.write(f'\t\t\t<m6>{ModelPar[13]}</m6>\n')
            f.write('\t\t</material>\n')
    elif ModelChoice =='Fung':
        for i in range(num_cells):
            f.write(f'\t\t<material id="{i+1}" name="Material{i+1}" ype="Fung orthotropic">\n')
            f.write(f'\t\t\t<density>{ModelPar[0]}</density>\n')
            f.write(f'\t\t\t<E1>{ModelPar[1]}</E1>\n')
            f.write(f'\t\t\t<E2>{ModelPar[2]}</E2>\n')
            f.write(f'\t\t\t<E3>{ModelPar[3]}</E3>\n')
            f.write(f'\t\t\t<G12>{ModelPar[4]}</G12>\n')
            f.write(f'\t\t\t<G23>{ModelPar[5]}</G23>\n')
            f.write(f'\t\t\t<G31>{ModelPar[6]}</G31>\n')
            f.write(f'\t\t\t<v12>{ModelPar[7]}</v12>\n')
            f.write(f'\t\t\t<v23>{ModelPar[8]}</v23>\n')
            f.write(f'\t\t\t<v31>{ModelPar[9]}</v31>\n')
            f.write(f'\t\t\t<c>{ModelPar[10]}</c>\n')
            f.write(f'\t\t\t<k>{ModelPar[11]}</k>\n')
            f.write('\t\t</material>\n')
    elif ModelChoice == 'HGO':
        for i in range(num_cells):
            f.write(f'\t\t<material id="{i+1}" name="Material{i+1}" type="Holzapfel-Gasser-Ogden">\n')
            f.write(f'\t\t\t<c>{ModelPar[0]}</c>\n')
            f.write(f'\t\t\t<k1>{ModelPar[1]}</k1>\n')
            f.write(f'\t\t\t<k2>{ModelPar[2]}</k2>\n')
            f.write(f'\t\t\t<gamma>{ModelPar[3]}</gamma>\n')
            f.write(f'\t\t\t<kappa>{1/3}</kappa>\n')
            f.write('\t\t\t<k>100</k>\n')
            f.write('\t\t\t<mat_axis type="vector">\n') 
            f.write(f'\t\t\t\t<a>{Circ[i,0]},{Circ[i,1]},{Circ[i,2]}</a>\n')
            f.write(f'\t\t\t\t<d>{Long[i,0]},{Long[i,1]},{Long[i,2]}</d>\n')
            f.write('\t\t\t</mat_axis>\n') 
            f.write('\t\t</material>\n')
                
    f.write('\t</Material>\n')
    #Mesh and geometry
    f.write('\t<Mesh>\n')
    f.write('\t\t<Nodes name="Object01">\n')
    for i in range(num_pts):
        f.write(f'\t\t\t<node id="{i+1}">{X[i,0]},{X[i,1]},{X[i,2]}</node>\n')        
        
    f.write('\t\t</Nodes>\n')   
    
    if ModelChoice in ['MR','tiMR','HGO']:
        for i in range(num_cells):
            f.write(f'\t\t<Elements type="quad4" name="Part{i+1}">\n')
            c = pdi.GetCell(i)
            if c.GetNumberOfPoints()!=4:
                raise ValueError()
            f.write(f'\t\t\t<elem id="{i+1}">{1+c.GetPointId(0)}')
            for j in range(1,4):
                f.write(f',{1+c.GetPointId(j)}')
            f.write('</elem>\n')
            f.write('\t\t</Elements>\n')
    else:
        f.write('\t\t<Elements type="quad4" name="Part1">\n')
        for i in range(num_cells):
            c = pdi.GetCell(i)
            if c.GetNumberOfPoints()!=4:
                raise ValueError()
            f.write(f'\t\t\t<elem id="{i+1}">{1+c.GetPointId(0)}')
            for j in range(1,4):
                f.write(f',{1+c.GetPointId(j)}')
            f.write('</elem>\n')
        f.write('\t\t</Elements>\n')
            
    #Create node sets to apply boundary conditions
    xd = vtk_to_numpy(in_pd.GetArray('VAJ')).astype(int)
    xd = np.where(xd==1)
    
    for i,j in enumerate(STJ_Id):
        f.write('\t\t<NodeSet name="STJ_'+str(i+1)+'">\n')
        f.write(f'\t\t\t<n id="{j+1}"/>\n')
        f.write('\t\t</NodeSet>\n')
    for i,j in enumerate(VAJ_Id):
        f.write('\t\t<NodeSet name="VAJ_'+str(i+1)+'">\n')
        f.write(f'\t\t\t<n id="{j+1}"/>\n')
        f.write('\t\t</NodeSet>\n')
    
    #Loading is applied on the same surface as the mesh
    f.write('\t\t<Surface name="PressureLoad1">\n')
    for i in range(0,num_cells):
        c = pdi.GetCell(i)
        if c.GetNumberOfPoints()!=4:
            raise ValueError()
        f.write(f'\t\t\t<quad4 id="{i+1}">{1+c.GetPointId(0)},{1+c.GetPointId(1)},{1+c.GetPointId(2)},{1+c.GetPointId(3)}')
        f.write('</quad4>\n')
    f.write('\t\t</Surface>\n')
    f.write('\t</Mesh>\n')
    
    #Mesh domain assigns material to the surface
    f.write('\t<MeshDomains>\n')
    if ModelChoice in ['MR','tiMR','HGO']:
        for i  in range(num_cells):
            f.write(f'\t\t<ShellDomain name="Part{i+1}" mat="Material{i+1}">\n')
            f.write('\t\t\t<shell_normal_nodal>1</shell_normal_nodal>\n')
            f.write('\t\t</ShellDomain>\n')
    else:
        f.write('\t\t<ShellDomain name="Part1" mat="Material1">\n')
        f.write('\t\t\t<shell_normal_nodal>1</shell_normal_nodal>\n')
        f.write('\t\t</ShellDomain>\n')
    f.write('\t</MeshDomains>\n')
    
    #Mesh data where we assign thickness (and possibly fiber directions)
    f.write('\t<MeshData>\n')
    if ModelChoice in ['MR','tiMR','HGO']:
        for i in range(num_cells):
            f.write(f'\t\t<ElementData var="shell thickness" elem_set="Part{i+1}">\n')
            c = pdi.GetCell(i)
            f.write(f'\t\t\t<e lid="1">{thick[c.GetPointId(0)]}')
            for j in range(1,4):
                f.write(f',{thick[c.GetPointId(j)]}')
            f.write('</e>\n')
            f.write('\t\t</ElementData>\n') 
    else:
        f.write('\t\t<ElementData var="shell thickness" elem_set="Part1">\n')
        for i in range(num_cells):
            c = pdi.GetCell(i)
            f.write(f'\t\t\t<elem lid="{i+1}">{thick[c.GetPointId(0)]}')
            for j in range(1,4):
                f.write(f',{thick[c.GetPointId(j)]}')
            f.write('</elem>\n')
        f.write('\t\t</ElementData>\n')  
    f.write('\t</MeshData>\n')
    
    #Boundary conditions
    f.write('\t<Boundary>\n')
    for k, j in enumerate(['x','y','z']):
        for i, Id in enumerate(STJ_Id):
            f.write('\t\t<bc name="STJ_'+str(i+1)+'" type="prescribe" node_set="STJ_'+str(i+1)+'">\n')
            f.write('\t\t\t<dof>'+j+'</dof>\n')
            f.write('\t\t\t<scale lc="'+str(1+2+i + len(STJ_Id)*k)+'">1</scale>\n')
            f.write('\t\t\t<relative>0</relative>\n')
            f.write('\t\t</bc>\n')
    for k, j in enumerate(['x','y','z']):
        for i, Id in enumerate(VAJ_Id):
            f.write('\t\t<bc name="VAJ_'+str(i+1)+'" type="prescribe" node_set="VAJ_'+str(i+1)+'">\n')
            f.write('\t\t\t<dof>'+j+'</dof>\n')
            f.write('\t\t\t<scale lc="'+str(1+2+3*len(STJ_Id)+i + len(VAJ_Id)*k)+'">1</scale>\n')
            f.write('\t\t\t<relative>0</relative>\n')
            f.write('\t\t</bc>\n')
    f.write('\t</Boundary>\n')
    
    #Apply load due to pressure
    f.write('\t<Loads>\n')
    f.write('\t\t<surface_load name="PressureLoad1" type="pressure" surface="PressureLoad1">\n')
    # if ModelChoice == 'tiMR':
    f.write(f'\t\t\t<pressure lc="1">{-PressureMag}</pressure>\n')
    # else:
    #     f.write(f'\t\t\t<pressure lc="1">{PressureMag}</pressure>\n')
    f.write('\t\t\t<linear>0</linear>\n')
    f.write('\t\t\t<symmetric_stiffness>2</symmetric_stiffness>\n')
    f.write('\t\t</surface_load>\n')
    f.write('\t</Loads>\n')
    
    f.write('\t<LoadData>\n')
    f.write('\t\t<load_controller id="1" type="loadcurve">\n')
    f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
    f.write('\t\t\t<points>\n')
    for i in range(nF):
        f.write('\t\t\t\t<point>' + str(i/(nF-1)) + ',' + str(1/(nF-1))+'</point>\n')
    f.write('\t\t\t</points>\n')
    f.write('\t\t</load_controller>\n')
    f.write('\t\t<load_controller id="2" type="loadcurve">\n')
    f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
    f.write('\t\t\t<points>\n')
    for i in range(nF):
        f.write('\t\t\t\t<point>' + str(i/(nF-1)) + ',' + str(PressureInterp[i])+'</point>\n')
    f.write('\t\t\t</points>\n')
    f.write('\t\t</load_controller>\n')
    for k, j in enumerate(['x','y','z']):
        for i in range(len(STJ_Id)):
            f.write('\t\t<load_controller id="'+str(1+2+i + len(STJ_Id)*k)+'" type="loadcurve">\n')
            f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
            f.write('\t\t\t<points>\n')
            for l in range(nF):
                f.write('\t\t\t\t<point>' + str(l/(nF-1)) + ',' + str(STJ_WallMotion[l,i,k])+'</point>\n')
            f.write('\t\t\t</points>\n')
            f.write('\t\t</load_controller>\n')
    for k, j in enumerate(['x','y','z']):
        for i in range(len(VAJ_Id)):
            f.write('\t\t<load_controller id="'+str(1+2+3*len(STJ_Id)+i + len(VAJ_Id)*k)+'" type="loadcurve">\n')
            f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
            f.write('\t\t\t<points>\n')
            for l in range(nF):
                f.write('\t\t\t\t<point>' + str(l/(nF-1)) + ',' + str(VAJ_WallMotion[l,i,k])+'</point>\n')
            f.write('\t\t\t</points>\n')
            f.write('\t\t</load_controller>\n')
    f.write('\t</LoadData>\n')
    
    #Define the outputs needed in the resulting output file
    f.write('\t<Output>\n')
    f.write('\t\t<plotfile type="febio">\n')
    f.write('\t\t\t<var type="displacement"/>\n')
    # f.write('\t\t\t<var type="relative volume"/>\n')
    f.write('\t\t\t<var type="stress"/>\n')
    f.write('\t\t</plotfile>\n')
    f.write('\t</Output>\n')
    f.write('</febio_spec>')
    f.close()

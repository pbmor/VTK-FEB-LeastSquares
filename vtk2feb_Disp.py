import vtk
import numpy as np
from numpy import genfromtxt
from vtk.util.numpy_support import vtk_to_numpy
from math import exp as e
from math import cos, sin, pi
from scipy.integrate import odeint
from scipy.interpolate import interp1d

def VTK2Feb_Func(DataDir,Fname,nF,nCells,STJ_WallMotion,VAJ_WallMotion,STJ_Id,VAJ_Id,Circ,Long,ProfileChoice,ModelChoice,FiberChoice,CF):
    
    TimeInterp = np.linspace(0,1,nF)
    
    if ProfileChoice == 'Triangle':
        Peak = 1/2
        PressureInterp = np.zeros(nF)
        Line_A = (1/Peak)*TimeInterp
        Line_B = (-1/(1 - Peak))*TimeInterp + 1/(1-Peak)
        for i in range(nF):
            if TimeInterp[i] <=Peak:
                PressureInterp[i] = Line_A[i]
            else:
                PressureInterp[i] = Line_B[i]
        PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
    elif ProfileChoice == 'Step':
        PressureInterp = np.zeros(nF)
        for i in range(nF):
            if (TimeInterp[i] >=2/nF) and (TimeInterp[i] <= 5/nF):
                PressureInterp[i] = 1.0
        PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
    elif ProfileChoice == 'SmoothStep':
        PressureInterp = np.zeros(nF)
        for i in range(nF):
            PressureInterp[i] = (1/(1+e(-50*(TimeInterp[i]-2/nF))))*(1-1/(1+e(-50*(TimeInterp[i]-5/nF))))
        PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
    elif ProfileChoice == 'Bio':
        TimeData = (genfromtxt('TimeProfile.csv', delimiter=','))
        PressureData = genfromtxt('PressureProfile.csv', delimiter=',')
        PressureInterp = np.interp(TimeInterp,TimeData,PressureData)
        PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
    elif ProfileChoice == 'Fourier':
        PressureInterp = np.zeros(nF)
        x=TimeInterp
        for i in range(nF):
            PressureInterp[i] = 3.0*sin(pi*x[i])+1.0*sin(2*pi*x[i])+0.0*sin(4*pi*x[i])+0.05*sin(8*pi*x[i])
        PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))
    elif ProfileChoice == 'Fitted':
        PressureInterpStan = 0.0*np.ones(nF)
    elif ProfileChoice == 'Windkess':
        qi_mag = 1
        Rp = 1
        Rd = 1
        Cp  = 1
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
        PressureInterpStan = np.subtract(PressureInterp,PressureInterp[0])
        
    elif ProfileChoice[0:3] == 'Set':
        if ProfileChoice[3:] == 'Windkess':
            qi_mag = 1.08833313e+01
            Rp     = 1.39329372e+00
            Rd     = 1.43243521e+01
            Cp     = 4.19180868e-03
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
            PressureInterpStan = np.subtract(PressureInterp,PressureInterp[0])
    
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
    
    NewFebName = './FEB_Files/' + DataDir+'.feb'
    f=open(NewFebName,'w')
    
    #Generic starting information
    f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
    f.write('<febio_spec version="3.0">\n')
    f.write('\t<Module type="solid"/>\n')
    
    #Solver settings
    f.write('\t<Control>\n')
    f.write('\t\t<analysis>STATIC</analysis>\n')
    f.write('\t\t<time_steps>'+str(nF)+'</time_steps>\n')
    f.write('\t\t<step_size>'+str(1/nF)+'</step_size>\n')
    f.write('\t\t<plot_level>PLOT_MUST_POINTS</plot_level>\n')
    f.write('\t\t<time_stepper>\n')
    f.write('\t\t\t<dtmax lc="1">'+str(1/nF)+'</dtmax>\n')
    f.write('\t\t\t<max_retries>5</max_retries>\n')
    f.write('\t\t\t<opt_iter>10</opt_iter>\n')
    f.write('\t\t</time_stepper>\n')
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
        f.write('\t\t\t<density>0.0</density>\n')
        f.write('\t\t\t<c1>10</c1>\n')
        f.write('\t\t\t<c2>10</c2>\n')
        f.write('\t\t\t<k>10</k>\n')
        f.write('\t\t</material>\n') 
    elif ModelChoice =='tiMR':
        for i in range(nCells):
            f.write(f'\t\t<material id="{i+1}" name="Material{i+1}" type="2D trans iso Mooney-Rivlin">\n') 
            f.write('\t\t\t<density>10.0</density>\n') 
            f.write('\t\t\t<c1>14.0</c1>\n') 
            f.write('\t\t\t<c2>0.0</c2>\n') 
            f.write('\t\t\t<c3>2.0</c3>\n') 
            f.write('\t\t\t<c4>60.0</c4>\n') 
            f.write('\t\t\t<c5>600.0</c5>\n') 
            f.write('\t\t\t<k>100.0</k>\n') 
            f.write('\t\t\t<lam_max>1.0</lam_max>\n') 
            if FiberChoice:
                theta = 0
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
        f.write('\t\t<material id="1" name="Material1" type="Ogden">\n')
        f.write('\t\t\t<density>10</density>\n')
        f.write('\t\t\t<k>10</k>\n')
        f.write('\t\t\t<c1>10</c1>\n')
        f.write('\t\t\t<c2>10</c2>\n')
        f.write('\t\t\t<c3>10</c3>\n')
        f.write('\t\t\t<c4>10</c4>\n')
        f.write('\t\t\t<c5>10</c5>\n')
        f.write('\t\t\t<c6>10</c6>\n')
        f.write('\t\t\t<m1>10</m1>\n')
        f.write('\t\t\t<m2>10</m2>\n')
        f.write('\t\t\t<m3>10</m3>\n')
        f.write('\t\t\t<m4>10</m4>\n')
        f.write('\t\t\t<m5>10</m5>\n')
        f.write('\t\t\t<m6>10</m6>\n')
        f.write('\t\t</material>\n')
    elif ModelChoice =='Fung':
        f.write('\t\t<material id="1" name="Material1" type="Fung orthotropic">\n')
        f.write('\t\t\t<density>1</density>\n')
        f.write('\t\t\t<E1>1</E1>\n')
        f.write('\t\t\t<E2>1</E2>\n')
        f.write('\t\t\t<E3>1</E3>\n')
        f.write('\t\t\t<G12>1</G12>\n')
        f.write('\t\t\t<G23>1</G23>\n')
        f.write('\t\t\t<G31>1</G31>\n')
        f.write('\t\t\t<v12>1</v12>\n')
        f.write('\t\t\t<v23>1</v23>\n')
        f.write('\t\t\t<v31>1</v31>\n')
        f.write('\t\t\t<c>1</c>\n')
        f.write('\t\t\t<k>1</k>\n')
        f.write('\t\t</material>\n')
    elif ModelChoice == 'HGO':
        f.write('\t\t<material id="1" name="Material1" type="Holzapfel-Gasser-Ogden">\n')
        f.write('\t\t\t<c>7.6</c>\n')
        f.write('\t\t\t<k1>1000.0</k1>\n')
        f.write('\t\t\t<k2>520.0</k2>\n')
        f.write('\t\t\t<gamma>50.0</gamma>\n')
        f.write('\t\t\t<kappa>0.20</kappa>\n')
        f.write('\t\t\t<k>100000.0</k>\n')
        f.write('\t\t</material>\n')
                
    f.write('\t</Material>\n')
    #Mesh and geometry
    f.write('\t<Mesh>\n')
    f.write('\t\t<Nodes name="Object01">\n')
    for i in range(num_pts):
        f.write(f'\t\t\t<node id="{i+1}">{X[i,0]},{X[i,1]},{X[i,2]}</node>\n')
    f.write('\t\t</Nodes>\n')
    
    
    if ModelChoice == 'tiMR' or ModelChoice == 'HGO_unc':
        for i in range(nCells):
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
    # f.write('\t\t<NodeSet name="VAJ">\n')
    # for i in range(len(xd[0])):
    #     f.write(f'\t\t\t<n id="{1+xd[0][i]}"/>\n')
    # f.write('\t\t</NodeSet>\n')
    # xd = vtk_to_numpy(in_pd.GetArray('STJ')).astype(int)
    # xd = np.where(xd==1)
    # f.write('\t\t<NodeSet name="STJ">\n')
    # for i in range(len(xd[0])):
    #     f.write(f'\t\t\t<n id="{1+xd[0][i]}"/>\n')
    # f.write('\t\t</NodeSet>\n')
    for i in STJ_Id:
        f.write('\t\t<NodeSet name="STJ_'+str(i+1)+'">\n')
        f.write(f'\t\t\t<n id="{i+1}"/>\n')
        f.write('\t\t</NodeSet>\n')
    for i in VAJ_Id:
        f.write('\t\t<NodeSet name="VAJ_'+str(i+1)+'">\n')
        f.write(f'\t\t\t<n id="{i+1}"/>\n')
        f.write('\t\t</NodeSet>\n')
    
    #Loading is applied on the same surface as the mesh
    f.write('\t\t<Surface name="PressureLoad1">\n')
    for i in range(0,num_cells):
        c = pdi.GetCell(i)
        if c.GetNumberOfPoints()!=4:
            raise ValueError()
        f.write(f'\t\t\t<quad4 id="{i+1}">{1+c.GetPointId(0)}')
        for j in range(1,4):
            f.write(f',{1+c.GetPointId(j)}')
        f.write('</quad4>\n')
    f.write('\t\t</Surface>\n')
    f.write('\t</Mesh>\n')
    
    #Mesh domain assigns material to the surface
    f.write('\t<MeshDomains>\n')
    if ModelChoice == 'tiMR' or ModelChoice == 'HGO_unc':
        for i  in range(nCells):
            f.write(f'\t\t<ShellDomain name="Part{i+1}" mat="Material{i+1}">\n')
            f.write('\t\t\t<shell_normal_nodal>1</shell_normal_nodal>\n')
            f.write('\t\t</ShellDomain>\n')
    else:
        f.write('\t\t<ShellDomain name="Part1" mat="Material1">\n')
        f.write('\t\t\t<shell_normal_nodal>1</shell_normal_nodal>\n')
        f.write('\t\t</ShellDomain>\n')
    f.write('\t</MeshDomains>\n')
    
    #Mesh data where we assign thickness (and possibly fiber directions)
    thick = vtk_to_numpy(in_pd.GetArray('Thickness'))
    f.write('\t<MeshData>\n')
    if ModelChoice == 'tiMR' or ModelChoice == 'HGO_unc':
        for i in range(nCells):
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
        if ModelChoice == 'HGO':
            f.write('\t\t<ElementData var="mat_axis" elem_set="Part1">\n')
            for i in range(num_cells):
                f.write(f'\t\t\t<elem lid="{i+1}">\n')
                f.write(f'\t\t\t\t<a>{Circ[i,0]},{Circ[i,1]},{Circ[i,2]}</a>\n')
                f.write(f'\t\t\t\t<d>{-Long[i,0]},{-Long[i,1]},{-Long[i,2]}</d>\n')
                f.write('\t\t\t</elem>\n')
            f.write('\t\t</ElementData>\n')   
    f.write('\t</MeshData>\n')
    
    #Boundary conditions
    f.write('\t<Boundary>\n')
    # f.write('\t\t<bc name="FixedDisplacement1" type="fix" node_set="STJ">\n')
    # f.write('\t\t\t<dofs>x,y,z</dofs>\n')
    # f.write('\t\t</bc>\n')
    # f.write('\t\t<bc name="FixedDisplacement1" type="fix" node_set="VAJ">\n')
    # f.write('\t\t\t<dofs>x,y,z</dofs>\n')
    # f.write('\t\t</bc>\n')
    for k, j in enumerate(['x','y','z']):
        for i, Id in enumerate(STJ_Id):
            f.write('\t\t<bc name="STJ_'+str(Id+1)+'" type="prescribe" node_set="STJ_'+str(Id+1)+'">\n')
            f.write('\t\t\t<dof>'+j+'</dof>\n')
            f.write('\t\t\t<scale lc="'+str(1+1+i + len(STJ_Id)*k)+'">1</scale>\n')
            f.write('\t\t\t<relative>0</relative>\n')
            f.write('\t\t</bc>\n')
    for k, j in enumerate(['x','y','z']):
        for i, Id in enumerate(VAJ_Id):
            f.write('\t\t<bc name="VAJ_'+str(Id+1)+'" type="prescribe" node_set="VAJ_'+str(Id+1)+'">\n')
            f.write('\t\t\t<dof>'+j+'</dof>\n')
            f.write('\t\t\t<scale lc="'+str(1+1+3*len(STJ_Id)+i + len(VAJ_Id)*k)+'">1</scale>\n')
            f.write('\t\t\t<relative>0</relative>\n')
            f.write('\t\t</bc>\n')
    f.write('\t</Boundary>\n')
    
    #Apply load due to pressure
    f.write('\t<Loads>\n')
    f.write('\t\t<surface_load name="PressureLoad1" type="pressure" surface="PressureLoad1">\n')
    f.write('\t\t\t<pressure lc="1">-1.0</pressure>\n')
    f.write('\t\t\t<linear>0</linear>\n')
    f.write('\t\t\t<symmetric_stiffness>1</symmetric_stiffness>\n')
    f.write('\t\t</surface_load>\n')
    f.write('\t</Loads>\n')
    
    f.write('\t<LoadData>\n')
    f.write('\t\t<load_controller id="1" type="loadcurve">\n')
    f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
    f.write('\t\t\t<points>\n')
    for i in range(nF):
        f.write('\t\t\t\t<point>' + str(TimeInterp[i]) + ',' + str(1/nF)+'</point>\n')
    f.write('\t\t\t</points>\n')
    f.write('\t\t</load_controller>\n')
    f.write('\t\t<load_controller id="2" type="loadcurve">\n')
    f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
    f.write('\t\t\t<points>\n')
    for i in range(nF):
        f.write('\t\t\t\t<point>' + str(TimeInterp[i]) + ',' + str(PressureInterpStan[i])+'</point>\n')
    f.write('\t\t\t</points>\n')
    f.write('\t\t</load_controller>\n')
    for k, j in enumerate(['x','y','z']):
        for i in range(len(STJ_Id)):
            f.write('\t\t<load_controller id="'+str(2+1+i + len(STJ_Id)*k)+'" type="loadcurve">\n')
            f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
            f.write('\t\t\t<points>\n')
            for l in range(nF):
                f.write('\t\t\t\t<point>' + str(TimeInterp[l]) + ',' + str(STJ_WallMotion[l,i,k])+'</point>\n')
            f.write('\t\t\t</points>\n')
            f.write('\t\t</load_controller>\n')
    for k, j in enumerate(['x','y','z']):
        for i in range(len(VAJ_Id)):
            f.write('\t\t<load_controller id="'+str(2+1+3*len(STJ_Id)+i + len(VAJ_Id)*k)+'" type="loadcurve">\n')
            f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
            f.write('\t\t\t<points>\n')
            for l in range(nF):
                f.write('\t\t\t\t<point>' + str(TimeInterp[l]) + ',' + str(VAJ_WallMotion[l,i,k])+'</point>\n')
            f.write('\t\t\t</points>\n')
            f.write('\t\t</load_controller>\n')
    f.write('\t</LoadData>\n')
    
    #Define the outputs needed in the resulting output file
    f.write('\t<Output>\n')
    f.write('\t\t<plotfile type="febio">\n')
    f.write('\t\t\t<var type="displacement"/>\n')
    f.write('\t\t\t<var type="relative volume"/>\n')
    f.write('\t\t\t<var type="stress"/>\n')
    f.write('\t\t</plotfile>\n')
    f.write('\t</Output>\n')
    f.write('</febio_spec>')
    f.close()

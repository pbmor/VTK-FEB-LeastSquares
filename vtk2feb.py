import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
reader = vtk.vtkPolyDataReader()
reader.SetFileName('resample.vtk')
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()
pdi = reader.GetOutput()
num_pts = pdi.GetNumberOfPoints()
num_cells = pdi.GetNumberOfCells()
pts = pdi.GetPoints()
in_pd = pdi.GetPointData()
X = vtk_to_numpy(pts.GetData())

f=open('test.feb','w')

#Generic starting information
f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
f.write('<febio_spec version="3.0">\n')
f.write('\t<Module type="solid"/>\n')

#Solver settings
f.write('\t<Control>\n')
f.write('\t\t<analysis>STATIC</analysis>\n')
f.write('\t\t<time_steps>10</time_steps>\n')
f.write('\t\t<step_size>0.1</step_size>\n')
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
f.write('\t\t\t<dtmin>0.01</dtmin>\n')
f.write('\t\t\t<dtmax>0.1</dtmax>\n')
f.write('\t\t\t<max_retries>5</max_retries>\n')
f.write('\t\t\t<opt_iter>10</opt_iter>\n')
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
f.write('\t\t<material id="1" name="Material1" type="Mooney-Rivlin">\n')
f.write('\t\t\t<density>1</density>\n')
f.write('\t\t\t<c1>1</c1>\n')
f.write('\t\t\t<c2>0</c2>\n')
f.write('\t\t\t<k>10</k>\n')
f.write('\t\t</material>\n')
f.write('\t</Material>\n')

#Mesh and geometry
f.write('\t<Mesh>\n')
f.write('\t\t<Nodes name="part1">\n')
for i in range(0,num_pts):
    f.write(f'\t\t\t<node id="{i+1}">{X[i,0]},{X[i,1]},{X[i,2]}</node>\n')
f.write('\t\t</Nodes>\n')
f.write('\t\t<Elements type="quad4" name="Part1">\n')
for i in range(0,num_cells):
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
f.write('\t\t<NodeSet name="VAJ">\n')
for i in range(len(xd[0])):
    f.write(f'\t\t\t<n id="{1+xd[0][i]}"/>\n')
f.write('\t\t</NodeSet>\n')
xd = vtk_to_numpy(in_pd.GetArray('STJ')).astype(int)
xd = np.where(xd==1)
f.write('\t\t<NodeSet name="STJ">\n')
for i in range(len(xd[0])):
    f.write(f'\t\t\t<n id="{1+xd[0][i]}"/>\n')
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
f.write('\t<MeshDomains>\n\t\t<ShellDomain name="Part1" mat="Material1">\n\t\t\t<shell_normal_nodal>1</shell_normal_nodal>\n\t\t</ShellDomain>\n\t</MeshDomains>\n')

#Mesh data where we assign thickness (and possibly fiber directions)
thick = vtk_to_numpy(in_pd.GetArray('Thickness'))
f.write('\t<MeshData>\n')
f.write('\t\t<ElementData var="shell thickness" elem_set="Part1">\n')
for i in range(0,num_cells):
    c = pdi.GetCell(i)
    f.write(f'\t\t\t<e lid="{i+1}">{thick[c.GetPointId(0)]}')
    for j in range(1,4):
        f.write(f',{thick[c.GetPointId(j)]}')
    f.write('</e>\n')
    #For constant thickness
    #f.write(f'\t\t\t<e lid="{i+1}">1,1,1,1</e>\n')
f.write('\t\t</ElementData>\n')
f.write('\t</MeshData>\n')

#Boundary conditions
f.write('\t<Boundary>\n')
f.write('\t\t<bc name="FixedDisplacement1" type="fix" node_set="VAJ">\n')
f.write('\t\t\t<dofs>x,y,z</dofs>\n')
f.write('\t\t</bc>\n')
f.write('\t\t<bc name="FixedDisplacement2" type="fix" node_set="STJ">\n')
f.write('\t\t\t<dofs>x,y,z</dofs>\n')
f.write('\t\t</bc>\n')
f.write('\t</Boundary>\n')

#Apply load due to pressure
f.write('\t<Loads>\n')
f.write('\t\t<surface_load name="PressureLoad1" type="pressure" surface="PressureLoad1">\n')
f.write('\t\t\t<pressure lc="1">-0.005332</pressure>\n')
f.write('\t\t\t<linear>0</linear>\n')
f.write('\t\t\t<symmetric_stiffness>1</symmetric_stiffness>\n')
f.write('\t\t</surface_load>\n')
f.write('\t</Loads>\n')

f.write('\t<LoadData>\n')
f.write('\t\t<load_controller id="1" type="loadcurve">\n')
f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
f.write('\t\t\t<points>\n')
f.write('\t\t\t\t<point>0,0</point>\n')
f.write('\t\t\t\t<point>1,1</point>\n')
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

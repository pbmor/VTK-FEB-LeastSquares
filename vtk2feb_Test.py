
import numpy as np
from numpy import genfromtxt
    

Circ_Pts = 8
Height_Pts = 5
Nf = 10
FT = 100

num_pts = Circ_Pts*Height_Pts
num_cells = Circ_Pts*(Height_Pts-1)
TimeData = (genfromtxt('TimeProfile.csv', delimiter=','))
    
PressureData = genfromtxt('PressureProfile.csv', delimiter=',')

TimeInterp = np.linspace(0,Nf*float(FT)/1000,Nf)
PressureInterp  = np.interp(TimeInterp,TimeData,PressureData)
PressureInterpStan = np.subtract(PressureInterp,np.min(PressureInterp))/np.max(np.subtract(PressureInterp,np.min(PressureInterp)))

X = np.zeros((num_pts,3))
r = np.linspace(0,10,Height_Pts)
k=0
for i in range(Height_Pts):
    for j in range(Circ_Pts):
        X[k,:] = [np.power((r[i]+1),1/3)*np.cos(j*2*np.pi/Circ_Pts),np.power((r[i]+1),1/3)*np.sin(j*2*np.pi/Circ_Pts),(r[i]-5)/2]
        k+=1

X_Cell = np.zeros((num_cells,4))
k=0
for i in range(Height_Pts-1):
    for j in range(Circ_Pts):
        X_Cell[k,:] = [(j%Circ_Pts)+Circ_Pts*i+1,((j+1)%Circ_Pts)+Circ_Pts*i+1,((j+1)%Circ_Pts)+Circ_Pts*(i+1)+1,(j%Circ_Pts)+Circ_Pts*(i+1)+1]
        k+=1

VAJ_Cent = [np.sum(X[0:Circ_Pts,0])/Circ_Pts,np.sum(X[0:Circ_Pts,1])/Circ_Pts,np.sum(X[0:Circ_Pts,2])/Circ_Pts]
STJ_Cent = [np.sum(X[Circ_Pts*(Height_Pts-1):Circ_Pts*(Height_Pts),0])/Circ_Pts,np.sum(X[Circ_Pts*(Height_Pts-1):Circ_Pts*(Height_Pts),1])/Circ_Pts,np.sum(X[Circ_Pts*(Height_Pts-1):Circ_Pts*(Height_Pts),2])/Circ_Pts]

VAJ_Motion = np.zeros((Circ_Pts,3))
STJ_Motion = np.zeros((Circ_Pts,3))
for i in range(Circ_Pts):
    VAJ_Motion[i] = [(VAJ_Cent[0]-X[i,0])/np.sqrt(np.sum(np.fromiter(((X[i,j]-VAJ_Cent[j])**2 for j in range(3)),float))),
                     (VAJ_Cent[1]-X[i,1])/np.sqrt(np.sum(np.fromiter(((X[i,j]-VAJ_Cent[j])**2 for j in range(3)),float))),
                     (VAJ_Cent[2]-X[i,2])/np.sqrt(np.sum(np.fromiter(((X[i,j]-VAJ_Cent[j])**2 for j in range(3)),float)))]
    STJ_Motion[i] = [(STJ_Cent[0]-X[Circ_Pts*(Height_Pts-1)+i,0])/np.sqrt(np.sum(np.fromiter(((X[Circ_Pts*(Height_Pts-1)+i,j]-STJ_Cent[j])**2 for j in range(3)),float))),
                     (STJ_Cent[1]-X[Circ_Pts*(Height_Pts-1)+i,1])/np.sqrt(np.sum(np.fromiter(((X[Circ_Pts*(Height_Pts-1)+i,j]-STJ_Cent[j])**2 for j in range(3)),float))),
                     (STJ_Cent[2]-X[Circ_Pts*(Height_Pts-1)+i,2])/np.sqrt(np.sum(np.fromiter(((X[Circ_Pts*(Height_Pts-1)+i,j]-STJ_Cent[j])**2 for j in range(3)),float)))]

NewFebName = 'Test.feb'
print(NewFebName)
f=open(NewFebName,'w')

#Generic starting information
f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
f.write('<febio_spec version="3.0">\n')
f.write('\t<Module type="solid"/>\n')

#Solver settings
f.write('\t<Control>\n')
f.write('\t\t<analysis>STATIC</analysis>\n')
f.write('\t\t<time_steps>'+str(Nf)+'</time_steps>\n')
f.write('\t\t<step_size>'+str(float(FT)/1000)+'</step_size>\n')
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
f.write('\t\t\t<dtmin>'+str(0.9*float(FT)/1000)+'</dtmin>\n')
f.write('\t\t\t<dtmax>'+str(1.1*float(FT)/1000)+'</dtmax>\n')
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
    f.write(f'\t\t\t<elem id="{i+1}">'+str(int(X_Cell[i,0])))
    for j in range(1,4):
        f.write(','+str(int(X_Cell[i,j])))
    f.write('</elem>\n')
f.write('\t\t</Elements>\n')

#Create node sets to apply boundary conditions
f.write('\t\t<NodeSet name="VAJ">\n')
for i in range(Circ_Pts):
    f.write(f'\t\t\t<n id="{i+1}"/>\n')
f.write('\t\t</NodeSet>\n')
f.write('\t\t<NodeSet name="STJ">\n')
for i in range(Circ_Pts):
    f.write('\t\t\t<n id="'+str(Circ_Pts*(Height_Pts-1)+i+1)+'"/>\n')
f.write('\t\t</NodeSet>\n')
f.write('\t\t<NodeSet name="ANCHOR">\n')
f.write('\t\t\t<n id="'+str(int(Circ_Pts*(Height_Pts-1)/2+1))+'"/>\n')
f.write('\t\t</NodeSet>\n')
for i in range(Circ_Pts):
    f.write(f'\t\t<NodeSet name="VAJ_{i+1}">\n')
    f.write(f'\t\t\t<n id="{i+1}"/>\n')
    f.write('\t\t</NodeSet>\n')
for i in range(Circ_Pts):
    f.write(f'\t\t<NodeSet name="STJ_{i+1}">\n')
    f.write('\t\t\t<n id="'+str(Circ_Pts*(Height_Pts-1)+i+1)+'"/>\n')
    f.write('\t\t</NodeSet>\n')

#Loading is applied on the same surface as the mesh
f.write('\t\t<Surface name="PressureLoad1">\n')
for i in range(0,num_cells):
    f.write(f'\t\t\t<quad4 id="{i+1}">'+str(int(X_Cell[i,0])))
    for j in range(1,4):
        f.write(','+str(int(X_Cell[i,j])))
    f.write('</quad4>\n')
f.write('\t\t</Surface>\n')
f.write('\t</Mesh>\n')

#Mesh domain assigns material to the surface
f.write('\t<MeshDomains>\n\t\t<ShellDomain name="Part1" mat="Material1">\n\t\t\t<shell_normal_nodal>1</shell_normal_nodal>\n\t\t</ShellDomain>\n\t</MeshDomains>\n')

#Mesh data where we assign thickness (and possibly fiber directions)
f.write('\t<MeshData>\n')
f.write('\t\t<ElementData var="shell thickness" elem_set="Part1">\n')
for i in range(0,num_cells):
    f.write(f'\t\t\t<e lid="{i+1}">0.1')
    for j in range(1,4):
        f.write(',0.1')
    f.write('</e>\n')
f.write('\t\t</ElementData>\n')
f.write('\t\t<ElementData var="material_dir" elem_set="Part1">\n')
f.write('\t\t\t<mat_axis type="vector">\n')
f.write('\t\t\t\t<a>1,0,0</a>\n')
f.write('\t\t\t\t<d>0,1,0</d>\n')
f.write('\t\t\t</mat_axis>\n')
f.write('\t\t</ElementData>\n')
# f.write('\t\t<ElementData var="fibs" elem_set="Part1">\n')
# for i in range(0,num_cells):
#     f.write(f'\t\t\t<e lid="{i+1}">0.1')
#     for j in range(1,4):
#         f.write(',0.1')
#     f.write('</e>\n')
# f.write('\t\t</ElementData>\n')
f.write('\t</MeshData>\n')

#Boundary conditions
f.write('\t<Boundary>\n')
# f.write('\t\t<bc name="FixedDisplacement1" type="fix" node_set="VAJ">\n')
# f.write('\t\t\t<dofs>x,y,z</dofs>\n')
# f.write('\t\t</bc>\n')
# f.write('\t\t<bc name="FixedDisplacement1" type="fix" node_set="STJ">\n')
# f.write('\t\t\t<dofs>x,y,z</dofs>\n')
# f.write('\t\t</bc>\n')

# f.write('\t\t<bc name="FixedDisplacement1" type="fix" node_set="ANCHOR">\n')
# f.write('\t\t\t<dofs>x,y,z</dofs>\n')
# f.write('\t\t</bc>\n')

# f.write('\t\t<bc type="prescribe" node_set="VAJ">\n')
# f.write('\t\t\t<dof>x</dof>\n')
# f.write('\t\t\t<scale lc="1">2.0</scale>\n')
# f.write('\t\t\t<relative>0</relative>\n')
# f.write('\t\t</bc>\n')
# f.write('\t\t<bc type="prescribe" node_set="STJ">\n')
# f.write('\t\t\t<dof>x</dof>\n')
# f.write('\t\t\t<scale lc="1">2.0</scale>\n')
# f.write('\t\t\t<relative>0</relative>\n')
# f.write('\t\t</bc>\n')

# f.write('\t\t<bc type="prescribe" node_set="STJ">\n')
# f.write('\t\t\t<dof>y</dof>\n')
# f.write('\t\t\t<scale lc="1">-2.0</scale>\n')
# f.write('\t\t\t<relative>0</relative>\n')
# f.write('\t\t</bc>\n')
        
for i in range(Circ_Pts):
    for k, j in enumerate(['x','y']):
        f.write(f'\t\t<bc type="prescribe" node_set="VAJ_{i+1}">\n')
        f.write('\t\t\t<dof>'+j+'</dof>\n')
        # f.write('\t\t\t<scale lc="1">1.0</scale>\n')
        # f.write('\t\t\t<scale lc="1">'+str(0.1*(np.random.normal(1)-1)**2)+'</scale>\n')
        if k<=1:
            f.write('\t\t\t<scale lc="1">'+str(VAJ_Motion[i,k])+'</scale>\n')
        else:
            f.write('\t\t\t<scale lc="1">'+str(-1)+'</scale>\n')
        f.write('\t\t\t<relative>0</relative>\n')
        f.write('\t\t</bc>\n')
        
for i in range(Circ_Pts):
    for  j in ['z']:
        f.write(f'\t\t<bc type="prescribe" node_set="VAJ_{i+1}">\n')
        f.write('\t\t\t<dof>'+j+'</dof>\n')
        # f.write('\t\t\t<scale lc="1">1.0</scale>\n')
        # f.write('\t\t\t<scale lc="1">'+str(0.1*(np.random.normal(1)-1)**2)+'</scale>\n')
        f.write('\t\t\t<scale lc="2">'+str(-1)+'</scale>\n')
        f.write('\t\t\t<relative>0</relative>\n')
        f.write('\t\t</bc>\n')

# for i in range(Circ_Pts):
#     for k, j in enumerate(['x','y','z']):
#         f.write(f'\t\t<bc type="prescribe" node_set="STJ_{i+1}">\n')
#         f.write('\t\t\t<dof>'+j+'</dof>\n')
#         # f.write('\t\t\t<scale lc="1">1.0</scale>\n')
#         # f.write('\t\t\t<scale lc="1">'+str(0.1*(np.random.normal(1)-1)**2)+'</scale>\n')
#         if k<=1:
#             f.write('\t\t\t<scale lc="1">'+str(STJ_Motion[i,k])+'</scale>\n')
#         else:
#             f.write('\t\t\t<scale lc="1">'+str(1)+'</scale>\n')
#         f.write('\t\t\t<relative>0</relative>\n')
#         f.write('\t\t</bc>\n')

for i in range(Circ_Pts):
    for k, j in enumerate(['x','y','z']):
        f.write(f'\t\t<bc type="prescribe" node_set="STJ_{i+1}">\n')
        f.write('\t\t\t<dof>'+j+'</dof>\n')
        f.write('\t\t\t<scale lc="'+str(3+i)+'">'+str(1)+'</scale>\n')
        f.write('\t\t\t<relative>0</relative>\n')
        f.write('\t\t</bc>\n')
'''
for i in range(Circ_Pts):
    f.write(f'\t\t<bc type="prescribe" node_set="STJ_{i+1}">\n')
    f.write('\t\t\t<dof>z</dof>\n')
    # f.write('\t\t\t<scale lc="1">1.0</scale>\n')
    f.write('\t\t\t<scale lc="1">'+str(0.1*(np.random.normal(1)-1)**2)+'</scale>\n')
    f.write('\t\t\t<relative>0</relative>\n')
    f.write('\t\t</bc>\n')
    
for i in range(Circ_Pts):
    f.write(f'\t\t<bc type="prescribe" node_set="STJ_{i+1}">\n')
    f.write('\t\t\t<dof>x</dof>\n')
    # f.write('\t\t\t<scale lc="1">1.0</scale>\n')
    f.write('\t\t\t<scale lc="1">'+str(0.1*(np.random.normal(1)-1)**2)+'</scale>\n')
    f.write('\t\t\t<relative>0</relative>\n')
    f.write('\t\t</bc>\n')
    
    
for i in range(Circ_Pts):
    f.write(f'\t\t<bc type="prescribe" node_set="STJ_{i+1}">\n')
    f.write('\t\t\t<dof>y</dof>\n')
    # f.write('\t\t\t<scale lc="1">1.0</scale>\n')
    f.write('\t\t\t<scale lc="1">'+str(0.1*(np.random.normal(1)-1)**2)+'</scale>\n')
    f.write('\t\t\t<relative>0</relative>\n')
    f.write('\t\t</bc>\n')
    
for i in range(Circ_Pts):
    f.write(f'\t\t<bc type="prescribe" node_set="VAJ_{i+1}">\n')
    f.write('\t\t\t<dof>z</dof>\n')
    # f.write('\t\t\t<scale lc="1">1.0</scale>\n')
    f.write('\t\t\t<scale lc="1">'+str(0.1*(np.random.normal(1)-1)**2)+'</scale>\n')
    f.write('\t\t\t<relative>0</relative>\n')
    f.write('\t\t</bc>\n')
    
for i in range(Circ_Pts):
    f.write(f'\t\t<bc type="prescribe" node_set="VAJ_{i+1}">\n')
    f.write('\t\t\t<dof>x</dof>\n')
    # f.write('\t\t\t<scale lc="1">1.0</scale>\n')
    f.write('\t\t\t<scale lc="1">'+str(0.1*(np.random.normal(1)-1)**2)+'</scale>\n')
    f.write('\t\t\t<relative>0</relative>\n')
    f.write('\t\t</bc>\n')
    
    
for i in range(Circ_Pts):
    f.write(f'\t\t<bc type="prescribe" node_set="VAJ_{i+1}">\n')
    f.write('\t\t\t<dof>y</dof>\n')
    # f.write('\t\t\t<scale lc="1">1.0</scale>\n')
    f.write('\t\t\t<scale lc="1">'+str(0.1*(np.random.normal(1)-1)**2)+'</scale>\n')
    f.write('\t\t\t<relative>0</relative>\n')
    f.write('\t\t</bc>\n')

'''
# for i in range(1,Circ_Pts):
#     f.write(f'\t\t<bc type="prescribe" node_set="VAJ_{i+1}">\n')
#     f.write('\t\t\t<dof>z</dof>\n')
#     # f.write('\t\t\t<scale lc="1">-10.0</scale>\n')
#     f.write('\t\t\t<scale lc="1">'+str((np.random.normal(1)-1)**2)+'</scale>\n')
#     f.write('\t\t\t<relative>0</relative>\n')
#     f.write('\t\t</bc>\n')
f.write('\t</Boundary>\n')

# f.write('\t<MeshData>\n')
# f.write('\t\t<mesh_data param="fem.bc_prescribed[0]">\n')
# f.write('\t\t\t<node lid="1">1.0</node>\n')
# # for i in range(Circ_Pts):
#     # f.write('\t\t\t<node lid="'+str(Circ_Pts*(Height_Pts-1)+i+1)+'">1.0</node>\n')
#     # f.write('\t\t\t<node lid="'+str(i+1)+'">'+str(i+1)+'</node>\n')
# f.write('\t\t</mesh_data>\n')
# f.write('\t</MeshData>\n')
        
#Apply load due to pressure
f.write('\t<Loads>\n')
f.write('\t\t<surface_load name="PressureLoad1" type="pressure" surface="PressureLoad1">\n')
f.write('\t\t\t<pressure lc="1">-0.05332</pressure>\n')
# f.write('\t\t\t<pressure lc="1">-0.02</pressure>\n')
# f.write('\t\t\t<pressure lc="1">-0.03</pressure>\n')
f.write('\t\t\t<linear>0</linear>\n')
f.write('\t\t\t<symmetric_stiffness>1</symmetric_stiffness>\n')
f.write('\t\t</surface_load>\n')
f.write('\t</Loads>\n')

f.write('\t<LoadData>\n')
f.write('\t\t<load_controller id="1" type="loadcurve">\n')
f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
f.write('\t\t\t<points>\n')
for i in range(Nf):
    f.write('\t\t\t\t<point>' + str(TimeInterp[i]) + ',' + str(PressureInterpStan[i])+'</point>\n')
f.write('\t\t\t</points>\n')
f.write('\t\t</load_controller>\n')
f.write('\t\t<load_controller id="2" type="loadcurve">\n')
f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
f.write('\t\t\t<points>\n')
for i in range(Nf):
    f.write('\t\t\t\t<point>' + str(TimeInterp[i]) + ',' + str(np.sin(TimeInterp*4*np.pi)[i])+'</point>\n')
f.write('\t\t\t</points>\n')
f.write('\t\t</load_controller>\n')
for i in range(Circ_Pts):
    f.write('\t\t<load_controller id="'+str(3+i)+'" type="loadcurve">\n')
    f.write('\t\t\t<interpolate>SMOOTH</interpolate>\n')
    f.write('\t\t\t<points>\n')
    for i in range(Nf):
        f.write('\t\t\t\t<point>' + str(TimeInterp[i]) + ',' + str(0.1*np.random.normal(1))+'</point>\n')
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

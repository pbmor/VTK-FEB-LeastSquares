import numpy as np
import os
from xml.etree import ElementTree as et
from scipy.optimize import least_squares
from Read_XPLTfuncs import GetFEB, GetMeshInfo, GetData


def GetRes(C,Disp_2):
    
    '''
    Use xml tree to update .feb file
    '''
    tree = et.parse('test.feb')
    tree.find('Material').find('material').find('c1').text = str(C[0])
    if len(C)==2:
        tree.find('Material').find('material').find('c2').text = str(C[1])
    elif len(C)==3:
        tree.find('Material').find('material').find('k').text = str(C[2])
    tree.write('test.feb',xml_declaration=True,encoding="ISO-8859-1")
    
    #Run updated .feb file to create new .xplt file
    os.system('/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio3 -i test.feb >/dev/null 2>&1')
    
    XPLTfilename = 'test.xplt'
    
    feb,file_size,nstates = GetFEB(XPLTfilename)
    
    nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
    
    #Collect the data, with all states, note spaces in vairable names are replaced with _
    displacement = GetData(feb,'displacement',nstates,nVar)
      
    #Gather stresses and displacements Note this is dependent on the Var Type
    Disp_X    = np.zeros((nstates,nNodes))
    Disp_Y    = np.zeros((nstates,nNodes))
    Disp_Z    = np.zeros((nstates,nNodes))
    Disp_1    = np.zeros((nstates,nNodes,3))      
    for i in range(nstates): 
        for j in range(nNodes):
            Disp_X[i,j]  = displacement[i][j*3+1]
            Disp_Y[i,j]  = displacement[i][j*3+2]
            Disp_Z[i,j]  = displacement[i][j*3+3]
        Disp_1[i,:,0]  = Disp_X[i,:]  
        Disp_1[i,:,1]  = Disp_Y[i,:]  
        Disp_1[i,:,2]  = Disp_Z[i,:]
        
    Res = np.subtract(Disp_1,Disp_2)

    return Res.flatten()
    
    
if __name__=='__main__':
    # Standard case with c1=1, c2=0 and k=10
    feb,file_size,nstates = GetFEB('test2.xplt')
    
    nNodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)
    
    #Collect the data, with all states, note spaces in vairable names are replaced with _
    displacement = GetData(feb,'displacement',nstates,nVar)
      
    #Gather stresses and displacements Note this is dependent on the Var Type
    Disp_X    = np.zeros((nstates,nNodes))
    Disp_Y    = np.zeros((nstates,nNodes))
    Disp_Z    = np.zeros((nstates,nNodes))
    Disp_2    = np.zeros((nstates,nNodes,3))      
    for i in range(nstates): 
        for j in range(nNodes):
            Disp_X[i,j]  = displacement[i][j*3+1]
            Disp_Y[i,j]  = displacement[i][j*3+2]
            Disp_Z[i,j]  = displacement[i][j*3+3]
        Disp_2[i,:,0]  = Disp_X[i,:]  
        Disp_2[i,:,1]  = Disp_Y[i,:]  
        Disp_2[i,:,2]  = Disp_Z[i,:]
        
    
    Out = least_squares(GetRes,[1.1,0.1,10.1], verbose = 2, jac = '3-point', args = [Disp_2])
    
    print(Out)
    
    
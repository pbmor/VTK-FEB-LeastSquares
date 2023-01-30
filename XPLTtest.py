from Read_XPLTfuncs import GetFEB,GetMeshInfo,GetData

    
# xpltname = '/Users/petermortensen/Downloads/jobs/tav02.xplt'
xpltname = '/Users/petermortensen/Documents/RAposition/Code/semi-automated segmentations - root.nosync/VTK-FEB-LeastSquares/FEB_Files/bav07.xplt'
xpltName = './FEB_Files/tav12.xplt'
# xpltname = 'NewFiles/tav02/HGO/RunLS/SetWindkessel/PMagF_ModelT_FibF_CellPlane/jobs/tav02.xplt'
# xpltname = 'NewFiles/tav02/MR/RunDefault/SetWindkessel/PMagF_ModelF_FibF_CellPlane/jobs/tav02.xplt'
# Run GetFeb, using .xplt filename, the number of domains and your to output the data tree using True or False
feb, _,nStates, _ = GetFEB(xpltname,1086,True)
Nodes, nElems, nVar, StateTimes, VarNames, VarType = GetMeshInfo(feb)


displacement = GetData(feb,'displacement',nStates)
stress = GetData(feb,'stress',nStates)
# RV = GetData(feb,'relative volume',nStates)

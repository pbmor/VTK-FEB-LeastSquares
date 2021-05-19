import csv
import glob
import numpy as np
from natsort import natsorted # pip install natsort
from FileConversion import *

with open('echoframetime.csv') as csv_file:
    XLData = csv.reader(csv_file, delimiter=',')

FixAndRotate=False
List_of_Subdirectories = sorted(glob.glob('./medial_meshes/*'))
#List_of_Subdirectories = sorted(glob.glob('./medial_meshesOld/*'))
TotalMotionAvgBAV = np.zeros(900)
TotalMotionAvgTAV = np.zeros(900)
bavN = 0
tavN = 0
CommonOfDir = os.path.commonprefix(List_of_Subdirectories)
#print(List_of_Subdirectories)
for d in List_of_Subdirectories:
        
    DataDir = d.replace(CommonOfDir,'')

    with open('echoframetime.csv') as csv_file:
        XLData = csv.reader(csv_file, delimiter=',')
        for row in XLData:
            if DataDir == row[0]:
                DataInfo = row
    
    #Define Frame Time Length
    FT = DataInfo[1]
    #Define Open Frame and Close Frame
    OF = DataInfo[4]
    CF = DataInfo[5] 


    #for 'medial_meshes' in List_of_subdirectories:    
    fnames = sorted(glob.glob(d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/*.vtk')) 
#    fnames = sorted(glob.glob(d + '/medial meshes - propagated from reference/*.vtk'))
#    fnames = sorted(glob.glob(d + '/medial meshes/*.vtk'))

    
    #Choose reference frame
    refN = int(DataInfo[4])-1 #Choose frame before valve opens
    #refN = 0 #Choose first saved frame

    common = os.path.commonprefix(fnames)
    for Fname in list(fnames):
        X = Fname.replace(common,'')
        X = X.replace('.vtk','')
        X = np.fromstring(X, dtype=int, sep=' ')
        X=X[0]
        if X==refN:
            ref=Fname

    if not fnames:
        print(DataDir," is empty")
    else:
        fdir = os.path.dirname(fnames[0])
        # Check Directory
        if not os.path.exists(fdir):
            print('Error: Path does not exist:', fdir)
            sys.exit()
        if DataInfo[9] == 'n':
            print(DataDir,' is excluded')
        else:
            print(DataDir,' is included')
            WallArea, WallVol, LumenVol, Time, Pts, WallAreaRatio, WallVolRatio, LumenVolRatio, AvgJ, AvgI1, AvgJRatio, AvgI1Ratio, TotalMotion, N = ProcessData(flist=fnames,ref=fnames[0],FT=100.,OF=OF,CF=CF,opformat='vtp')
    
            
            if DataDir[0] == 'b':
                bavN += 1
                TotalMotionAvgBAV += TotalMotion
            elif DataDir[0] == 't':
                tavN += 1
                TotalMotionAvgTAV += TotalMotion

            print('Total Wall Area =',WallArea)
            print('Total Wall Volume =',WallVol)
            print('Total Lumen Volume =',LumenVol)
        
            ###################################
            # Save data
            DataLocation = 'Strains/' + d + '/medial meshes - propagated from reference/Data.npz'
            np.savez(DataLocation,Time=Time,Pts=Pts,WallArea=WallArea,WallVol=WallVol, LumenVol=LumenVol, WallAreaRatio=WallAreaRatio, WallVolRatio=WallVolRatio, LumenVolRatio=LumenVolRatio, AvgJ=AvgJ,AvgI1=AvgI1, AvgJRatio=AvgJRatio, AvgI1Ratio=AvgI1Ratio, N=N, OF=OF, CF=CF,refN = refN)


TotalMotionAvgBAV = TotalMotionAvgBAV/bavN
TotalMotionAvgTAV = TotalMotionAvgTAV/tavN

np.save('./TotalMotionAvgBAV',TotalMotionAvgBAV)
np.save('./TotalMotionAvgTAV',TotalMotionAvgTAV)


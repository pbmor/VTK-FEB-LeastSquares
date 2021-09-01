import csv
import os
import glob
import numpy as np
from natsort import natsorted 
from FileConversion import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib import cm


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
    NX = len(fnames)

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
            WallArea, WallVol, LumenVol, Time, Pts, WallAreaRatio, WallVolRatio, LumenVolRatio, AvgJ, AvgI1, AvgJRatio, AvgI1Ratio, TotalMotion, N = ProcessData(flist=fnames,ref=ref,FT=100.,OF=OF,CF=CF,opformat='vtp')

    
            
            if DataDir[0] == 'b':
                bavN += 1
                TotalMotionAvgBAV += TotalMotion[NX-1,:]
            elif DataDir[0] == 't':
                tavN += 1
                TotalMotionAvgTAV += TotalMotion[NX-1,:]

            print('Total Wall Area =',WallArea)
            print('Total Wall Volume =',WallVol)
            print('Total Lumen Volume =',LumenVol)
        
            ###################################
            # Save data
            DataLocation = 'Strains/' + d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/Data.npz'
            np.savez(DataLocation,Time=Time,Pts=Pts,WallArea=WallArea,WallVol=WallVol, LumenVol=LumenVol, WallAreaRatio=WallAreaRatio, WallVolRatio=WallVolRatio, LumenVolRatio=LumenVolRatio, AvgJ=AvgJ,AvgI1=AvgI1, AvgJRatio=AvgJRatio, AvgI1Ratio=AvgI1Ratio, N=N, OF=OF, CF=CF,refN = refN)
            
            L = len(Time)

            CSVDataRaw = np.zeros((L,15))
            CSVDataStandardA = np.zeros((33,15))
            CSVDataStandardB = np.zeros((68,15))
            CSVDataStandard = np.zeros((101,15))

            CSVDataRaw[:,0]  = Time
            CSVDataRaw[:,1]  = WallArea
            CSVDataRaw[:,2]  = WallVol
            CSVDataRaw[:,3]  = LumenVol
            CSVDataRaw[:,4]  = WallAreaRatio
            CSVDataRaw[:,5]  = WallVolRatio
            CSVDataRaw[:,6]  = LumenVolRatio
            CSVDataRaw[:,7]  = AvgJ
            CSVDataRaw[:,8]  = AvgI1
            CSVDataRaw[:,9]  = AvgJRatio
            CSVDataRaw[:,10] = AvgI1Ratio
            CSVDataRaw[:,11] = np.full(L,N)
            CSVDataRaw[:,12] = np.full(L,OF)
            CSVDataRaw[:,13] = np.full(L,CF)
            CSVDataRaw[:,14] = np.full(L,refN)

            TimeOpen  = np.linspace(0,320,33)
            TimeClose = np.linspace(330,1000,68)
            TimeStandard  = np.concatenate((TimeOpen,TimeClose))
            print(TimeStandard)
            
            Oid = int(OF) - int(refN)
            Cid = int(CF) - int(refN)
            
            CSVDataStandard[:,0] = TimeStandard
            for i in range(1,15):
                CSVDataStandardA[:,i] = np.interp(np.linspace(Time[0],Time[Cid],33),Time[0:Cid+1],CSVDataRaw[0:Cid+1,i])
                CSVDataStandardB[:,i] = np.interp(np.linspace(Time[Cid],Time[L-1],68),Time[Cid:L],CSVDataRaw[Cid:L,i])
                CSVDataStandard[:,i]  = np.concatenate((CSVDataStandardA[:,i],CSVDataStandardB[:,i]))

            CSVDataRawFile = 'Strains/' + d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/CSVDataRaw.csv'
            CSVDataStandardFile = 'Strains/' + d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/CSVDataStandard.csv'

            np.savetxt(CSVDataRawFile, CSVDataRaw, delimiter=",",header="Time,WallArea,WallVol, LumenVol, WallAreaRatio, WallVolRatio, LumenVolRatio, AvgJ,AvgI1, AvgJRatio, AvgI1Ratio, N, OF, CF,refN")
            np.savetxt(CSVDataStandardFile, CSVDataStandard, delimiter=",",header="Time,WallArea,WallVol, LumenVol, WallAreaRatio, WallVolRatio, LumenVolRatio, AvgJ,AvgI1, AvgJRatio, AvgI1Ratio, N, OF, CF,refN")
            

            plt.figure(num=1)
            plt.plot(Time, WallAreaRatio,'r', Time,WallVolRatio,'b',Time,LumenVolRatio,'k')
            plt.set_title=('Raw Time Data')
            plt.xlabel('Time (ms)',size=12)
            plt.ylabel('Ratio ',size=12)
            plt.legend(['Wall Area Ratio', 'Wall Volume Ratio','Lumen Volume Ratio'])
            plt.savefig('RatioCompRaw.png', format='png')
    

            plt.figure(num=2)
            plt.plot(CSVDataStandard[:,0], CSVDataStandard[:,4],'r', CSVDataStandard[:,0],CSVDataStandard[:,5],'b',CSVDataStandard[:,0],CSVDataStandard[:,6],'k')
            plt.set_title=('Standardised Time Data')
            plt.xlabel('Time (ms)',size=12)
            plt.ylabel('Ratio ',size=12)
            plt.legend(['Wall Area Ratio', 'Wall Volume Ratio','Lumen Volume Ratio'])
            plt.savefig('RatioCompStandard.png', format='png')

            os.system("pdflatex Summary.tex")
            pdfname = 'Strains/' + d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/Summary.pdf'
            os.rename("Summary.pdf", pdfname )
            Fig1name = 'Strains/' + d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/RatioCompRaw.png'
            Fig2name = 'Strains/' + d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/RatioCompStandard.png'
            os.rename("RatioCompRaw.png", Fig1name)
            os.rename("RatioCompStandard.png", Fig2name)

TotalMotionAvgBAV = TotalMotionAvgBAV/bavN
TotalMotionAvgTAV = TotalMotionAvgTAV/tavN

np.save('./TotalMotionAvgBAV',TotalMotionAvgBAV)
np.save('./TotalMotionAvgTAV',TotalMotionAvgTAV)



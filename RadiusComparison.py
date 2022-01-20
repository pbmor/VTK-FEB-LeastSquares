import numpy as np
from AnalysisFuncs import GetRadii
import matplotlib.pyplot as plt
import csv
import os
import glob
from matplotlib.pyplot import figure, Line2D
from mpl_toolkits import mplot3d

RootRadAvg = np.zeros((2,22))
with open('echoframetime.csv') as csv_file:
    XLData = csv.reader(csv_file, delimiter=',')

List_of_Subdirectories = sorted(glob.glob('Strains/medial_meshes/*'))
ND = len(List_of_Subdirectories)

common = os.path.commonprefix(List_of_Subdirectories)
for figN, d in enumerate(List_of_Subdirectories):
    DataDir = d.replace(common,'')
    print(DataDir,figN)
    RootType = DataDir[0:3]
    DataLocation = d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/Data.npz'

    with np.load(DataLocation) as Data:
        Time           = Data['Time']
        WallArea       = Data['WallArea']
        WallVol        = Data['WallVol']
        LumenVol       = Data['LumenVol']
        WallAreaRatio  = Data['WallAreaRatio']
        WallVolRatio   = Data['WallVolRatio']
        LumenVolRatio  = Data['LumenVolRatio']
        Pts            = Data['Pts']
        N              = Data['N']
        OF             = Data['OF']
        CF             = Data['CF']
        refN           = Data['refN']
        
    Oid = int(OF) - int(refN)
    Cid = int(CF) - int(refN)

    [Nf,NP,Ndim]     = Pts.shape

    Radius, _, _, Ring = GetRadii(Pts,Nf,NP,1)

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.scatter3D(Pts[0,:,0], Pts[0,:,1], Pts[0,:,2])
    # ax.plot3D(Ring[0,10,:,0],Ring[0,10,:,1],Ring[0,10,:,2])
    # plt.show()
    
    MeanRadius = np.zeros(Nf)
    for i in range(Nf):
        MeanRadius[i] = np.mean(Radius[i,10,:])

    Time = np.linspace(0,1,Nf)

    custom_lines = [Line2D([0], [0], color='r', lw=4),
                Line2D([0], [0], color='g', lw=4)]
    figure(num = 1)
    if RootType == 'bav':
        plt.plot(Time,MeanRadius,'r')
    elif RootType == 'tav':
        plt.plot(Time,MeanRadius,'g')
    plt.ylabel('Radius (mm)')
    plt.xlabel('Time (ms)')
    plt.legend(custom_lines, ['BAVs', 'TAVs'], loc='upper right')

    figure(num = 2)
    if RootType == 'bav':
        plt.plot(Time,MeanRadius/MeanRadius[0],'r')
    elif RootType == 'tav':
        plt.plot(Time,MeanRadius/MeanRadius[0],'g')
    plt.ylabel('Radius')
    plt.xlabel('Time (ms)')
    plt.legend(custom_lines, ['BAVs', 'TAVs'], loc='upper right')
    
    TimeOpen  = np.linspace(0,1,11)
    TimeClose = np.linspace(1,3,22)
    TimeNorm  = np.concatenate((TimeOpen,TimeClose))
    RadiusAvg = np.zeros((Nf,36))
    RadiusSD  = np.zeros((Nf,36))

    MeanRadiusNorma = np.zeros(11)
    MeanRadiusNormb = np.zeros(22)
    MeanRadiusNorm  = np.zeros(33)

    MeanRadiusNorma  = np.interp(np.linspace(Time[0],Time[Cid],11),Time[0:Cid+1],MeanRadius[0:Cid+1])
    MeanRadiusNormb  = np.interp(np.linspace(Time[Cid],Time[N-1],22),Time[Cid:N],MeanRadius[Cid:N])
    MeanRadiusNorm   = np.concatenate((MeanRadiusNorma,MeanRadiusNormb))
    
    
    figure(num = 3)
    if RootType == 'bav':
        plt.plot(TimeNorm,MeanRadiusNorm,'r')
    elif RootType == 'tav':
        plt.plot(TimeNorm,MeanRadiusNorm,'g')
    plt.ylabel('Radius (mm)')
    plt.xlabel('Time (ms)')
    plt.legend(custom_lines, ['BAVs', 'TAVs'], loc='upper right')

    figure(num = 4)
    if RootType == 'bav':
        plt.plot(TimeNorm,MeanRadiusNorm/MeanRadiusNorm[0],'r')
    elif RootType == 'tav':
        plt.plot(TimeNorm,MeanRadiusNorm/MeanRadiusNorm[0],'g')
    plt.ylabel('Radius')
    plt.xlabel('Time (ms)')
    plt.legend(custom_lines, ['BAVs', 'TAVs'], loc='upper right')
    
    
    id_min = np.argmin(MeanRadius)
    MeanRadius_fromMin = np.concatenate((MeanRadius[id_min:Nf],MeanRadius[0:id_min]))
    
    figure(num = 5)
    if RootType == 'bav':
        plt.plot(Time,MeanRadius_fromMin,'r')
    elif RootType == 'tav':
        plt.plot(Time,MeanRadius_fromMin,'g')
    plt.ylabel('Radius (mm)')
    plt.xlabel('Time (ms)')
    plt.title('Starting from minimum mean radius')
    plt.legend(custom_lines, ['BAVs', 'TAVs'], loc='upper right')

    figure(num = 6)
    if RootType == 'bav':
        plt.plot(Time,MeanRadius_fromMin/MeanRadius[0],'r')
    elif RootType == 'tav':
        plt.plot(Time,MeanRadius_fromMin/MeanRadius[0],'g')
    plt.ylabel('Radius')
    plt.xlabel('Time (ms)')
    plt.title('Starting from minimum mean radius, relative to the reference frame')
    plt.legend(custom_lines, ['BAVs', 'TAVs'], loc='upper right')

    DataLocation = os.path.join('./RadiusData/',DataDir+'_RadiusData.npz')
    np.savez(DataLocation,Time = Time ,MeanRadius_fromMin=MeanRadius_fromMin,MeanRadiusRatio_fromMin=(MeanRadius_fromMin/MeanRadius[0]))

plt.show()


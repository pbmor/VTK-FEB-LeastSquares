import numpy as np
import os
import sys
import glob
import csv
from matplotlib.pyplot import figure
from matplotlib import cm
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

BAVN = 0
TAVN = 0
WAavg      = np.zeros((2,22))
WARavg     = np.zeros((2,22))
WVavg      = np.zeros((2,22))
WVRavg     = np.zeros((2,22))
LVavg      = np.zeros((2,22))
LVRavg     = np.zeros((2,22))

WASD       = np.zeros((2,22))
WARSD      = np.zeros((2,22))
WVSD       = np.zeros((2,22))
WVRSD      = np.zeros((2,22))
LVSD       = np.zeros((2,22))
LVRSD      = np.zeros((2,22))

with open('echoframetime.csv') as csv_file:
    XLData = csv.reader(csv_file, delimiter=',')

List_of_Subdirectories = sorted(glob.glob('Strains/medial_meshes/*'))

common = os.path.commonprefix(List_of_Subdirectories)
for d in List_of_Subdirectories:
    DataDir = d.replace(common,'')
    print(DataDir)
    RootType = DataDir[0:3]

    with open('echoframetime.csv') as csv_file:
        XLData = csv.reader(csv_file, delimiter=',')
        for row in XLData:
            if DataDir == row[0]:
                DataInfo = row

    #Choose reference frame
    refN = int(DataInfo[4])-1 #Choose frame before valve opens

    DataLocation = d + '/medial meshes/Data.npz'
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

    TimeOpen  = np.linspace(0,1,11)
    TimeClose = np.linspace(1,2,11)
    TimeNorm  = np.concatenate((TimeOpen,TimeClose))

    WAa  = np.interp(np.linspace(Time[0],Time[Cid],11),Time[0:Cid+1],WallArea[0:Cid+1])
    WAb  = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],WallArea[Cid:N])
    WARa = np.interp(np.linspace(Time[0],Time[Cid],11),Time[0:Cid+1],WallAreaRatio[0:Cid+1])
    WARb = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],WallAreaRatio[Cid:N])
    WVa  = np.interp(np.linspace(Time[0],Time[Cid],11),Time[0:Cid+1],WallVol[0:Cid+1])
    WVb  = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],WallVol[Cid:N])
    WVRa = np.interp(np.linspace(Time[0],Time[Cid],11),Time[0:Cid+1],WallVolRatio[0:Cid+1])
    WVRb = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],WallVolRatio[Cid:N])
    LVa  = np.interp(np.linspace(Time[0],Time[Cid],11),Time[0:Cid+1],LumenVol[0:Cid+1])
    LVb  = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],LumenVol[Cid:N])
    LVRa = np.interp(np.linspace(Time[0],Time[Cid],11),Time[0:Cid+1],LumenVolRatio[0:Cid+1])
    LVRb = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],LumenVolRatio[Cid:N])

    WallAreaNorm      = np.concatenate((WAa,WAb))
    WallAreaRatioNorm = np.concatenate((WARa,WARb))
    WallVolNorm       = np.concatenate((WVa,WVb))
    WallVolRatioNorm  = np.concatenate((WVRa,WVRb))
    LumenVolNorm      = np.concatenate((LVa,LVb))
    LumenVolRatioNorm = np.concatenate((LVRa,LVRb))

    if RootType == 'bav':
        FN = 0        
        Title = 'BAV'
        BAVN += 1
        WAavg[0,:]  += WallAreaNorm
        WARavg[0,:] += WallAreaRatioNorm
        WVavg[0,:]  += WallVolNorm
        WVRavg[0,:] += WallVolRatioNorm
        LVavg[0,:]  += LumenVolNorm
        LVRavg[0,:] += LumenVolRatioNorm
    else :
        FN = 10
        Title = 'TAV'
        TAVN += 1
        WAavg[1,:]  += WallAreaNorm
        WARavg[1,:] += WallAreaRatioNorm
        WVavg[1,:]  += WallVolNorm
        WVRavg[1,:] += WallVolRatioNorm
        LVavg[1,:]  += LumenVolNorm
        LVRavg[1,:] += LumenVolRatioNorm


    figure(num=FN+1)
    plt.suptitle(Title)
    plt.subplot(2,3,1)  
    plt.plot(Time,WallArea)
    plt.plot(Time[Oid],WallArea[Oid],'.',color='red')
    plt.plot(Time[Cid],WallArea[Cid],'.',color='red')
    plt.title('Wall Area',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Area (mm$^2$)',size=14)

    plt.subplot(2,3,4)
    plt.plot(Time,WallAreaRatio)
    plt.plot(Time[Oid],WallAreaRatio[Oid],'.',color='red')
    plt.plot(Time[Cid],WallAreaRatio[Cid],'.',color='red')
    plt.title('Wall Area Ratio',size=14)
    plt.xlabel('Time (ms)',size=14)
    plt.ylabel('Area Ratio',size=14)

    plt.subplot(2,3,2)
    plt.plot(Time,WallVol)
    plt.plot(Time[Oid],WallVol[Oid],'.',color='red')
    plt.plot(Time[Cid],WallVol[Cid],'.',color='red')
    plt.title('Wall Volume',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Volume (mm$^3$)',size=14)

    plt.subplot(2,3,5)
    plt.plot(Time,WallVolRatio)
    plt.plot(Time[Oid],WallVolRatio[Oid],'.',color='red')
    plt.plot(Time[Cid],WallVolRatio[Cid],'.',color='red')
    plt.title('Wall Volume Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Volume Ratio',size=14)

    plt.subplot(2,3,3)
    plt.plot(Time,LumenVol)
    plt.plot(Time[Oid],LumenVol[Oid],'.',color='red')
    plt.plot(Time[Cid],LumenVol[Cid],'.',color='red')
    plt.title('Lumen Volume',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Volume (mm$^3$)',size=14)

    plt.subplot(2,3,6)
    plt.plot(Time,LumenVolRatio)
    plt.plot(Time[Oid],LumenVolRatio[Oid],'.',color='red')
    plt.plot(Time[Cid],LumenVolRatio[Cid],'.',color='red')
    plt.title('Lumen Volume Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Volume Ratio',size=14)


    figure(num=FN+2)
    plt.suptitle(Title)
    plt.subplot(2,3,1)
    plt.plot(TimeNorm,WallAreaNorm)
    plt.title('Wall Area',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Area (mm$^2$)',size=14)

    plt.subplot(2,3,4)
    plt.plot(TimeNorm,WallAreaRatioNorm)
    plt.title('Wall Area Ratio',size=14)
    plt.xlabel('Time (ms)',size=14)
    plt.ylabel('Area Ratio',size=14)

    plt.subplot(2,3,2)
    plt.plot(TimeNorm,WallVolNorm)
    plt.title('Wall Volume',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Volume (mm$^3$)',size=14)

    plt.subplot(2,3,5)
    plt.plot(TimeNorm,WallVolRatioNorm)
    plt.title('Wall Volume Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Volume Ratio',size=14)

    plt.subplot(2,3,3)
    plt.plot(TimeNorm,LumenVolNorm)
    plt.title('Lumen Volume',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Volume (mm$^3$)',size=14)

    plt.subplot(2,3,6)
    plt.plot(TimeNorm,LumenVolRatioNorm)
    plt.title('Lumen Volume Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Volume Ratio',size=14)




common = os.path.commonprefix(List_of_Subdirectories)
for d in List_of_Subdirectories:
    DataDir = d.replace(common,'')
    print(DataDir)
    RootType = DataDir[0:3]

    with open('echoframetime.csv') as csv_file:
        XLData = csv.reader(csv_file, delimiter=',')
        for row in XLData:
            if DataDir == row[0]:
                DataInfo = row

    #Choose reference frame
    refN = int(DataInfo[4])-1 #Choose frame before valve opens

    DataLocation = d + '/medial meshes/Data.npz'
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

    TimeOpen  = np.linspace(0,1,11)
    TimeClose = np.linspace(1,2,11)
    TimeNorm  = np.concatenate((TimeOpen,TimeClose))

    WAa  = np.interp(np.linspace(0,Time[Cid],11),Time[0:Cid+1],WallArea[0:Cid+1])
    WAb  = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],WallArea[Cid:N])
    WARa = np.interp(np.linspace(0,Time[Cid],11),Time[0:Cid+1],WallAreaRatio[0:Cid+1])
    WARb = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],WallAreaRatio[Cid:N])
    WVa  = np.interp(np.linspace(0,Time[Cid],11),Time[0:Cid+1],WallVol[0:Cid+1])
    WVb  = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],WallVol[Cid:N])
    WVRa = np.interp(np.linspace(0,Time[Cid],11),Time[0:Cid+1],WallVolRatio[0:Cid+1])
    WVRb = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],WallVolRatio[Cid:N])
    LVa  = np.interp(np.linspace(0,Time[Cid],11),Time[0:Cid+1],LumenVol[0:Cid+1])
    LVb  = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],LumenVol[Cid:N])
    LVRa = np.interp(np.linspace(0,Time[Cid],11),Time[0:Cid+1],LumenVolRatio[0:Cid+1])
    LVRb = np.interp(np.linspace(Time[Cid],Time[N-1],11),Time[Cid:N],LumenVolRatio[Cid:N])

    WallAreaNorm      = np.concatenate((WAa,WAb))
    WallAreaRatioNorm = np.concatenate((WARa,WARb))
    WallVolNorm       = np.concatenate((WVa,WVb))
    WallVolRatioNorm  = np.concatenate((WVRa,WVRb))
    LumenVolNorm      = np.concatenate((LVa,LVb))
    LumenVolRatioNorm = np.concatenate((LVRa,LVRb))

    if RootType == 'bav':
        RN = BAVN
        WASD[0,:]  += (WallAreaNorm - WAavg[0,:]/RN)**2
        WARSD[0,:] += (WallAreaRatioNorm - WARavg[0,:]/RN)**2
        WVSD[0,:]  += (WallVolNorm - WVavg[0,:]/RN)**2
        WVRSD[0,:] += (WallVolRatioNorm - WVRavg[0,:]/RN)**2
        LVSD[0,:]  += (LumenVolNorm - LVavg[0,:]/RN)**2
        LVRSD[0,:] += (LumenVolRatioNorm - LVRavg[0,:]/RN)**2
    if RootType == 'tav':
        RN = TAVN
        WASD[1,:]  += (WallAreaNorm - WAavg[1,:]/RN)**2
        WARSD[1,:] += (WallAreaRatioNorm - WARavg[1,:]/RN)**2
        WVSD[1,:]  += (WallVolNorm - WVavg[1,:]/RN)**2
        WVRSD[1,:] += (WallVolRatioNorm - WVRavg[1,:]/RN)**2
        LVSD[1,:]  += (LumenVolNorm - LVavg[1,:]/RN)**2
        LVRSD[1,:] += (LumenVolRatioNorm - LVRavg[1,:]/RN)**2


for i in range(2):
    if i == 0:
        RN = BAVN
        FN = 0
        Title = 'BAV'
        Colour = 'firebrick'
    elif i== 1:
        RN = TAVN
        FN = 10
        Title = 'TAV'
        Colour = 'forestgreen'

    WASD[i,:]  = np.sqrt(WASD[i,:]/RN)
    WARSD[i,:] = np.sqrt(WARSD[i,:]/RN)
    WVSD[i,:]  = np.sqrt(WVSD[i,:]/RN)
    WVRSD[i,:] = np.sqrt(WVRSD[i,:]/RN)
    LVSD[i,:]  = np.sqrt(LVSD[i,:]/RN)
    LVRSD[i,:] = np.sqrt(LVRSD[i,:]/RN)


    figure(num=FN+2)
    plt.suptitle(Title)
    plt.subplot(2,3,1)
    plt.plot(TimeNorm,WAavg[i,:]/RN,color='k')
    plt.fill_between(TimeNorm,WAavg[i,:]/RN+WASD[i,:],WAavg[i,:]/RN-WASD[i,:],color='b',alpha=0.3)
    plt.title('Wall Area',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Area (mm$^2$)',size=14)

    plt.subplot(2,3,4)
    plt.plot(TimeNorm,WARavg[i,:]/RN,color='k')
    plt.fill_between(TimeNorm,WARavg[i,:]/RN+WARSD[i,:],WARavg[i,:]/RN-WARSD[i,:],color='b',alpha=0.3)
    plt.title('Wall Area Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Area Ratio',size=14)

    plt.subplot(2,3,2)
    plt.plot(TimeNorm,WVavg[i,:]/RN,color='k')
    plt.fill_between(TimeNorm,WVavg[i,:]/RN+WVSD[i,:],WVavg[i,:]/RN-WVSD[i,:],color='b',alpha=0.3)
    plt.title('Wall Volume',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Volume (mm$^2$)',size=14)

    plt.subplot(2,3,5)
    plt.plot(TimeNorm,WVRavg[i,:]/RN,color='k')
    plt.fill_between(TimeNorm,WVRavg[i,:]/RN+WVRSD[i,:],WVRavg[i,:]/RN-WVRSD[i,:],color='b',alpha=0.3)
    plt.title('Wall Volume Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Volume Ratio',size=14)

    plt.subplot(2,3,3)
    plt.plot(TimeNorm,LVavg[i,:]/RN,color='k')
    plt.fill_between(TimeNorm,LVavg[i,:]/RN+LVSD[i,:],LVavg[i,:]/RN-LVSD[i,:],color='b',alpha=0.3)
    plt.title('Lumen Volume',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Volume (mm$^3$)',size=14)

    plt.subplot(2,3,6)
    plt.plot(TimeNorm,LVRavg[i,:]/RN,color='k')
    plt.fill_between(TimeNorm,LVRavg[i,:]/RN+LVRSD[i,:],LVRavg[i,:]/RN-LVRSD[i,:],color='b',alpha=0.3)
    plt.title('Lumen Volume Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Volume Ratio',size=14)




    figure(num=100)
    plt.subplot(2,3,1)
   # plt.errorbar(TimeNorm,WAavg[i,:]/RN,WASD[i,:],color=Colour)
    plt.plot(TimeNorm,WAavg[i,:]/RN,color=Colour,label=Title)
    plt.fill_between(TimeNorm,WAavg[i,:]/RN+WASD[i,:],WAavg[i,:]/RN-WASD[i,:],color=Colour,alpha=0.3)
    plt.title('Wall Area',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Area (mm$^2$)',size=14)
    plt.legend()

    plt.subplot(2,3,4)
    #plt.errorbar(TimeNorm,WARavg[i,:]/RN,WARSD[i,:],color=Colour)    
    plt.plot(TimeNorm,WARavg[i,:]/RN,color=Colour)
    plt.fill_between(TimeNorm,WARavg[i,:]/RN+WARSD[i,:],WARavg[i,:]/RN-WARSD[i,:],color=Colour,alpha=0.3)
    plt.title('Wall Area Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Area Ratio',size=14)

    plt.subplot(2,3,2)
#    plt.errorbar(TimeNorm,WVavg[i,:]/RN,WVSD[i,:],color=Colour)
    plt.plot(TimeNorm,WVavg[i,:]/RN,color=Colour)
    plt.fill_between(TimeNorm,WVavg[i,:]/RN+WVSD[i,:],WVavg[i,:]/RN-WVSD[i,:],color=Colour,alpha=0.3)
    plt.title('Wall Volume',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Volume (mm$^2$)',size=14)

    plt.subplot(2,3,5)
    #plt.errorbar(TimeNorm,WVRavg[i,:]/RN,WVRSD[i,:],color=Colour)
    plt.plot(TimeNorm,WVRavg[i,:]/RN,color=Colour)
    plt.fill_between(TimeNorm,WVRavg[i,:]/RN+WVRSD[i,:],WVRavg[i,:]/RN-WVRSD[i,:],color=Colour,alpha=0.3)
    plt.title('Wall Volume Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Volume Ratio',size=14)

    plt.subplot(2,3,3)
    #plt.errorbar(TimeNorm,LVavg[i,:]/RN,LVSD[i,:],color=Colour)
    plt.plot(TimeNorm,LVavg[i,:]/RN,color=Colour)
    plt.fill_between(TimeNorm,LVavg[i,:]/RN+LVSD[i,:],LVavg[i,:]/RN-LVSD[i,:],color=Colour,alpha=0.3)
    plt.title('Lumen Volume',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Volume (mm$^3$)',size=14)

    plt.subplot(2,3,6)
    #plt.errorbar(TimeNorm,LVRavg[i,:]/RN,LVRSD[i,:],color=Colour)
    plt.plot(TimeNorm,LVRavg[i,:]/RN,color=Colour)
    plt.fill_between(TimeNorm,LVRavg[i,:]/RN+LVRSD[i,:],LVRavg[i,:]/RN-LVRSD[i,:],color=Colour,alpha=0.3)
    plt.title('Lumen Volume Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Volume Ratio',size=14)


plt.show()

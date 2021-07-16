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
WAavg      = np.zeros((2,202))
WARavg     = np.zeros((2,202))
WVavg      = np.zeros((2,202))
WVRavg     = np.zeros((2,202))
LVavg      = np.zeros((2,202))
LVRavg     = np.zeros((2,202))

WASD       = np.zeros((2,202))
WARSD      = np.zeros((2,202))
WVSD       = np.zeros((2,202))
WVRSD      = np.zeros((2,202))
LVSD       = np.zeros((2,202))
LVRSD      = np.zeros((2,202))

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
    
    if DataInfo[9] == 'n':
        print(DataDir,' is excluded')
    else:
        print(DataDir,' is included')
        DataLocation = d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/Data.npz'
        print(d)
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

        TimeOpen  = np.linspace(0,1,101)
        TimeClose = np.linspace(1,3,101)
        TimeNorm  = np.concatenate((TimeOpen,TimeClose))

        WAa  = np.interp(np.linspace(Time[0],Time[Cid],101),Time[0:Cid+1],WallArea[0:Cid+1])
        WAb  = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],WallArea[Cid:N])
        WARa = np.interp(np.linspace(Time[0],Time[Cid],101),Time[0:Cid+1],WallAreaRatio[0:Cid+1])
        WARb = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],WallAreaRatio[Cid:N])
        WVa  = np.interp(np.linspace(Time[0],Time[Cid],101),Time[0:Cid+1],WallVol[0:Cid+1])
        WVb  = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],WallVol[Cid:N])
        WVRa = np.interp(np.linspace(Time[0],Time[Cid],101),Time[0:Cid+1],WallVolRatio[0:Cid+1])
        WVRb = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],WallVolRatio[Cid:N])
        LVa  = np.interp(np.linspace(Time[0],Time[Cid],101),Time[0:Cid+1],LumenVol[0:Cid+1])
        LVb  = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],LumenVol[Cid:N])
        LVRa = np.interp(np.linspace(Time[0],Time[Cid],101),Time[0:Cid+1],LumenVolRatio[0:Cid+1])
        LVRb = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],LumenVolRatio[Cid:N])

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

    if DataInfo[9] == 'n':
        print(DataDir,' is excluded')
    else:
        print(DataDir,' is included')
        DataLocation = d + '/medial meshes - propagated from reference/with point data - recon from propagated boundary meshes/Data.npz'
        print(d)
        #DataLocation = d + '/medial meshes/Data.npz'
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

        TimeOpen  = np.linspace(0,1,101)
        TimeClose = np.linspace(1,3,101)
        TimeNorm  = np.concatenate((TimeOpen,TimeClose))
    
        WAa  = np.interp(np.linspace(0,Time[Cid],101),Time[0:Cid+1],WallArea[0:Cid+1])
        WAb  = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],WallArea[Cid:N])
        WARa = np.interp(np.linspace(0,Time[Cid],101),Time[0:Cid+1],WallAreaRatio[0:Cid+1])
        WARb = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],WallAreaRatio[Cid:N])
        WVa  = np.interp(np.linspace(0,Time[Cid],101),Time[0:Cid+1],WallVol[0:Cid+1])
        WVb  = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],WallVol[Cid:N])
        WVRa = np.interp(np.linspace(0,Time[Cid],101),Time[0:Cid+1],WallVolRatio[0:Cid+1])
        WVRb = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],WallVolRatio[Cid:N])
        LVa  = np.interp(np.linspace(0,Time[Cid],101),Time[0:Cid+1],LumenVol[0:Cid+1])
        LVb  = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],LumenVol[Cid:N])
        LVRa = np.interp(np.linspace(0,Time[Cid],101),Time[0:Cid+1],LumenVolRatio[0:Cid+1])
        LVRb = np.interp(np.linspace(Time[Cid],Time[N-1],101),Time[Cid:N],LumenVolRatio[Cid:N])

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
    



    fig = plt.figure(num=1)
    ax1 = fig.add_subplot(311)
    ax1.plot(TimeNorm, WARavg[i,:]/RN, c=Colour,zorder=1,label=Title)
    ax1.fill_between(TimeNorm,WARavg[i,:]/RN+WARSD[i,:],WARavg[i,:]/RN-WARSD[i,:],color=Colour,alpha=0.3)
    ax1.axes.get_xaxis().set_ticks([])
    ax1.set_ylabel('Wall Area Ratio',size=12)
    ax1.axes.get_xaxis().set_visible(False)
    top_side = ax1.spines["top"]
    top_side.set_visible(False)
    bottom_side = ax1.spines["bottom"]
    bottom_side.set_visible(False)
    right_side = ax1.spines["right"]
    right_side.set_visible(False)
    ax1.legend()

    ax2 = fig.add_subplot(312)
    ax2.plot(TimeNorm, WVRavg[i,:]/RN, c=Colour,zorder=1)
    ax2.fill_between(TimeNorm,WVRavg[i,:]/RN+WVRSD[i,:],WVRavg[i,:]/RN-WVRSD[i,:],color=Colour,alpha=0.3)
    ax2.axes.get_xaxis().set_ticks([])
    top_side = ax2.spines["top"]
    top_side.set_visible(False)
    bottom_side = ax2.spines["bottom"]
    bottom_side.set_visible(False)
    right_side = ax2.spines["right"]
    right_side.set_visible(False)
    ax2.set_ylabel('Wall Volume Ratio',size=12)
    
    ax3 = fig.add_subplot(313)
    ax3.plot(TimeNorm, LVRavg[i,:]/RN, c=Colour,zorder=1)
    ax3.fill_between(TimeNorm,LVRavg[i,:]/RN+LVRSD[i,:],LVRavg[i,:]/RN-LVRSD[i,:],color=Colour,alpha=0.3)
    ax3.axes.get_xaxis().set_ticks([])
    ax3.set_ylabel('Lumen Volume Ratio',size=12)
    top_side = ax3.spines["top"]
    top_side.set_visible(False)
    bottom_side = ax3.spines["bottom"]
    bottom_side.set_visible(False)
    right_side = ax3.spines["right"]
    right_side.set_visible(False)
    str_list = ['Open Valve','Closed Valve']
    ax3.set_xticks([0.5,2])
    ax3.tick_params(axis = "x", which = "both", bottom = False, top = False)
    ax3.set_xticklabels(str_list)
    for x in range(2):
        print(x)
        ax1.axvline(x=x,ymin=-1.2,ymax=1,c="black",linestyle = '--',linewidth=1,zorder=0, clip_on=False)
        ax2.axvline(x=x,ymin=0,ymax=1.2,c="black",linestyle = '--',linewidth=1, zorder=0,clip_on=False)
        ax3.axvline(x=x,ymin=0,ymax=1.2,c="black",linestyle = '--',linewidth=1, zorder=0,clip_on=False)
    plt.draw()





    figure(num=100)
    plt.subplot(2,3,1)
    plt.plot(TimeNorm,WAavg[i,:]/RN,color=Colour,label=Title)
    plt.fill_between(TimeNorm,WAavg[i,:]/RN+WASD[i,:],WAavg[i,:]/RN-WASD[i,:],color=Colour,alpha=0.3)
    plt.title('Wall Area',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Area (mm$^2$)',size=14)
    plt.legend()

    plt.subplot(2,3,4)    
    plt.plot(TimeNorm,WARavg[i,:]/RN,color=Colour)
    plt.fill_between(TimeNorm,WARavg[i,:]/RN+WARSD[i,:],WARavg[i,:]/RN-WARSD[i,:],color=Colour,alpha=0.3)
    plt.title('Wall Area Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Area Ratio',size=14)

    plt.subplot(2,3,2)
    plt.plot(TimeNorm,WVavg[i,:]/RN,color=Colour)
    plt.fill_between(TimeNorm,WVavg[i,:]/RN+WVSD[i,:],WVavg[i,:]/RN-WVSD[i,:],color=Colour,alpha=0.3)
    plt.title('Wall Volume',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Volume (mm$^3$)',size=14)

    plt.subplot(2,3,5)
    plt.plot(TimeNorm,WVRavg[i,:]/RN,color=Colour)
    plt.fill_between(TimeNorm,WVRavg[i,:]/RN+WVRSD[i,:],WVRavg[i,:]/RN-WVRSD[i,:],color=Colour,alpha=0.3)
    plt.title('Wall Volume Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Volume Ratio',size=14)

    plt.subplot(2,3,3)
    plt.plot(TimeNorm,LVavg[i,:]/RN,color=Colour)
    plt.fill_between(TimeNorm,LVavg[i,:]/RN+LVSD[i,:],LVavg[i,:]/RN-LVSD[i,:],color=Colour,alpha=0.3)
    plt.title('Lumen Volume',size=14)
    plt.xlabel('',size=14)
    plt.ylabel('Volume (mm$^3$)',size=14)

    plt.subplot(2,3,6)
    plt.plot(TimeNorm,LVRavg[i,:]/RN,color=Colour)
    plt.fill_between(TimeNorm,LVRavg[i,:]/RN+LVRSD[i,:],LVRavg[i,:]/RN-LVRSD[i,:],color=Colour,alpha=0.3)
    plt.title('Lumen Volume Ratio',size=14)
    plt.xlabel('Time  (ms)',size=14)
    plt.ylabel('Volume Ratio',size=14)
    plt.tight_layout()

plt.show()

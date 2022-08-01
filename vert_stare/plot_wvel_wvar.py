
import numpy as np
import glob
import os
from copy import deepcopy
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as colors     
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


filename = [20210531, 20210611, 20210709]
snr_th = ['0p998', '0p994', '1p005']

fig = plt.figure(figsize=(6*2, 3*3))
gs = gridspec.GridSpec(3,2, width_ratios=[1,1] ) 

folderpath = '/'

# Plot vertical velocity
lab = ['(a)', '(b)', '(c)']
for nf in range(len(filename)):
    
    #=============================================================================================
    #   open doppler lidar vertical stare
    #=============================================================================================
    inputfile = folderpath+'/fval_fmi_dlidST00_l1_any_v00_'+str(int(filename[nf]))+'.nc'
    ax1 = fig.add_subplot(gs[nf,0])
    xrs = xr.open_dataset(inputfile)
    data = deepcopy(xrs.dv.values)
    date = xrs.time
    beta = xrs.beta.values  
    height = xrs.range.values
    intensity = xrs.intensity.values
    
    data[xrs.intensity.values<1.005]=np.nan
    data[:,:4]=np.nan
    
    beta[xrs.intensity.values<1.005]=np.nan
    beta[:,:4]=np.nan
    
    intensity[xrs.intensity.values<1.005]=np.nan
    intensity[:,:4]=np.nan

    cloud_layer = []
    cloud_layer = deepcopy(beta)
    cloud_index = np.where(np.log10(beta)<-4)
    cloud_layer[cloud_index] = np.nan 
    
    sh = 5
    eh = 19
    st = 60*5
    et = 60*19
    ids = np.where(date.dt.hour ==sh)[0][0] #10
    ide = np.where(date.dt.hour ==eh)[0][0]#14
    xrs2 = xr.open_dataset(folderpath+'/fval_fmi_dlidST00_l2_any_v00_'+str(int(filename[nf]))+'.nc')
    x = xrs2.time2[st:et]
    y = xrs2.height.values
    z = xrs2.w.values[st:et,:]
    mlh_ = xrs2.mld.values
    x3 = xrs2.time
    z3 = xrs2.w2.values

    v1=[1e-7, 1e-6, 1e-5, 1e-4, 1e-3]#, 1e1]  
    v= (-5,-2,-1,-0.2,0.2, 1,2,5)
    interval = np.arange(0,24,4)

    #ax = fig.add_subplot(gs[nf])


    cs = ax1.contourf(x, y, z.T, v, cmap='seismic',extend='both')#,norm = colors.LogNorm()) #"blue","lightblue", "gold", 
    mlh, = ax1.plot(x3[7:41], mlh_[7:41], '--k', lw = 2, label = 'Estimated MLH')
    cs2 = ax1.contourf(date[ids:ide], y, cloud_layer[ids:ide,:].T, v1, cmap=colors.ListedColormap(["darkgreen"]), extend='both',norm = colors.LogNorm())
    ax1.set_ylabel('Height [m]', fontsize = 14) #set y label
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M')) #set format tick label
    ax1.xaxis.set_major_locator(mdates.HourLocator(interval))
    ax1.xaxis.set_minor_locator(mdates.MinuteLocator([0]))
    j=0
    n = 6
    ax1.tick_params(axis='x', which='major', length=6)
    ax1.set_ylim(0,3000)
    ax1.tick_params(labelsize=16)
    plt.text(0.01, 0.85, lab[nf], transform=ax1.transAxes, fontsize = 18)
    ax1.set_xlabel('Time [UTC]', fontsize = 14) #set x label 

    
import matplotlib.patches as mpatches
ax1 = fig.add_subplot(gs[0,0])
#cloud_legend = mpatches.Patch(color='navy', label='The red data')
ceilo = mpatches.Patch(color='darkgreen', label='The red data')
ax1.legend([mlh, ceilo], ['Estimated MLH', 'Cloud'], fontsize = 10, loc='upper right')
    
cbar_ax = fig.add_axes([0.1, -0.01, .4, 0.02]) 
cbar=fig.colorbar(cs, cax=cbar_ax, pad = 0.08, orientation='horizontal') 
cbar.ax.set_title("w [m s$^{-1}]$", loc = 'right', )
ax1.set_xlabel('Time [UTC]', fontsize = 14)  

# Plot variance
lab = ['(d)', '(e)', '(f)']
for nf in range(len(filename)):
    #=============================================================================================
    #   open doppler lidar vertical stare
    #=============================================================================================
    print(filename[nf])
    ax2 = fig.add_subplot(gs[nf,1])
    #load precipitation file
    inputfile_3 = folderpath+"sups_rao_mets00_l1_precip_v00_"+str(int(filename[nf]))+".nc"
    rainfile = xr.open_dataset(inputfile_3)
    t_file3 = rainfile.time
    precip = rainfile.precip.values

    xrs2 = xr.open_dataset(folderpath+'/fval_fmi_dlidST00_l2_any_v00_'+str(int(filename[nf]))+'.nc')
    x = xrs2.time2[st:et]
    y = xrs2.height.values
    z = xrs2.w.values[st:et,:]
    mlh_ = xrs2.mlh.values
    x3 = xrs2.time
    z3 = xrs2.w2.values

    
    
    #averaging x- minute average  
    t_rain = []
    for i in range(0,24):
        for j in [0,30]:
            try:
                time = np.where((t_file3.dt.day==int(str(filename[nf])[-2:])) & (t_file3.dt.hour==i) & (t_file3.dt.minute==j))[0][0]
            except IndexError:
                continue

            t_rain.append(time)
    t_rain = np.append(0, t_rain)

    rain_idx = np.zeros(48)
    rain_idx [:]= np.nan
    for i in range(1,48):
        rain_idx[i-1] = np.nanmean(precip[t_rain[i-1]:t_rain[i]+1])
    rain = np.where(rain_idx >0 )
        
    v1=[1e-7, 1e-6, 1e-5, 1e-4, 1e-3]#, 1e1]  

    v = np.arange(0,3.1,0.1)
    interval = np.arange(0,24,4)
    print(rain)
    try:
        z3[np.array(rain[0]), :] = np.nan
    except IndexError:
        pass
        
    cs = ax2.contourf(x3[sh*2:eh*2], y, z3[sh*2:eh*2, :].T, v, cmap='viridis',extend='both')#,norm = colors.LogNorm()) #"blue","lightblue", "gold", 
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M')) #set format tick label
    ax2.xaxis.set_major_locator(mdates.HourLocator(interval))
    ax2.xaxis.set_minor_locator(mdates.MinuteLocator([0]))
    j=0
    n = 6
    ax2.tick_params(axis='x', which='major', length=6)
    ax2.set_ylim(0,3000)
    ax2.tick_params(labelsize=16)
    ax2.set_yticklabels([])
    import matplotlib.patches as mpatches
    cloud_legend = mpatches.Patch(color='navy', label='The red data')
    ceilo = mpatches.Patch(color='cyan', label='The red data')
    plt.text(0.01, 0.85, lab[nf], transform=ax2.transAxes, fontsize = 18)

       
plt.figtext(0.48, 0.955+0.005, '31 May 2021', fontsize=16)
plt.figtext(0.48, 0.65-0.01, '11 June 2021', fontsize=16)   
plt.figtext(0.48, 0.34-0.015, '09 July 2021', fontsize=16)   
cbar_ax = fig.add_axes([0.57, -0.01, .4, 0.02]) 
cbar=fig.colorbar(cs, cax=cbar_ax, pad = 0.08, orientation='horizontal') 
cbar.ax.set_title(" $\sigma$$^{2}_w$ [m$^{-2}$ s$^{-2}$]", loc = 'right', )
    
ax2.set_xlabel('Time [UTC]', fontsize = 14) #set x label 
plt.subplots_adjust(wspace=0.05)
fig.tight_layout(pad=3)

fig.savefig('./example_case.pdf', format='pdf', dpi = 300, bbox_inches='tight')
plt.show() 





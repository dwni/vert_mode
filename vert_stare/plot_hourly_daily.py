
import matplotlib.cm as cm
import warnings
warnings.filterwarnings("ignore")
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from matplotlib import gridspec
import numpy as np
import glob
import os
import xarray as xr
import pandas as pd
from scipy import interpolate
from metpy.units import units
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

xrs = xr.open_dataset('var_merged.nc')
df3 = pd.read_feather('meteorological_param.feather')
filename = xrs.filename.values
var_new = xrs.var_m.values
y_new = np.arange(0.05, 1.05, 0.05)  

fig = plt.figure(figsize=(6*2, 5))
gs = gridspec.GridSpec(1, 2, width_ratios=[1,1] ) 

cs_cldtp = [32,33,31,34,35,36,67,68,69,70,71,0,1,2,3,23,4,5,24,25,26,6,22,27,7,8,9,10,11,12,28,29,13,14,30,15,16,17,18,19,20,
         21,54,53,44,38,45,51,59,56,60,72,47,65,52,46,62,49,48,57,58,66,63,41,61,104,39,37,50,43,]#,39,37,50,43,]

n=11
colors = plt.cm.brg(np.linspace(0,1,n))

#analytical profile Sorbjan (1989)
Sb = 1.17*(y_new*(1-y_new))**(2/3)

#analytical profile Lenschow (1980)
Le = 1.8*(y_new)**(2/3)*(1-(0.8*y_new))**2


# Hourly plot
st1 = np.arange(15, 35, 2)
et1 = np.arange(17, 37, 2)

ax2 = fig.add_subplot(gs[0])

var_ = var_new[cs_cldtp]
var_hourly_ = np.nanmean(var_,0)
for loop in range(len(st1)-1):
    st = st1[loop]
    et = et1[loop]
    
    var_fin = np.nanmean(var_hourly_[st:et,:],0)

    #Plotting     
    ax2.plot(var_fin, y_new,  color=colors[loop], label = str(int((st1[loop]+1)/2))+'-'+str(int((st1[loop]+1)/2)+1))#, label = str(xrs2.date.dt.hour[st].values)+'-'+str(xrs2.date.dt.hour[et].values)) 
    ax2.vlines(0,0,1,colors='grey',linestyles='dashed')
    ax2.set_ylabel('z/$z_i$', fontsize=16)
    ax2.tick_params(labelsize=16)
    ax2.set_xlabel('$\sigma$$^{2}_w$ / $w^2_*$', fontsize=16)  
    ax2.set_xlim(-0.1,1.)         
    ax2.set_ylim(0,1)
ax2.legend(loc="upper right", fontsize=12) 
ax2.grid(linestyle = 'dashed')
ax2.text(0.9, 0.05, '(a)', fontsize=18) 



# Daily plot
colors = plt.cm.brg(np.linspace(0,1,n))
scaling =''

y_new = np.arange(0.05, 1.05, 0.05)  
st = 19 
et = 29 
    
#analytical profile Lenschow (1980)
Le = 1.8*(y_new)**(2/3)*(1-(0.8*y_new))**2

ax1 = fig.add_subplot(gs[1])

var_fin_mean = []
var_fin = np.nanmean(var_[:,st:et+1,:],1)

for hr in range(len(cs_cldtp)):
    ax1.plot(var_fin[hr,:], y_new, color='grey'  , alpha=.9)
    var_fin_mean.append(var_fin[hr,:])
    
    ax1.vlines(0,0,1,colors='grey',linestyles='dashed')
    ax1.set_ylabel('z/$z_i$', fontsize=16)
    ax1.tick_params(labelsize=16)

    if scaling =='':
        ax1.set_xlabel('$\sigma$$^{2}_w$ / $w^2_*$', fontsize=16)  
    else:
        ax1.set_xlabel('$\sigma$$^{2}_w$ / $w^2_m$', fontsize=16)  
    ax1.set_xlim(-0.1, 1)        
    ax1.set_ylim(0,1)
    ax1.minorticks_on()
    ax1.tick_params(axis='x', which='minor', direction='out')   

ax1.plot(Le,y_new, '--k', lw=3, label ='Ref:L80')
ax1.plot(np.nanmean(var_fin_mean,0), y_new, '-r', lw=3, label = 'mean')    
ax1.legend(loc="upper right", fontsize =10)
ax1.grid(linestyle = 'dashed')
ax1.text(0.9, 0.05, '(b)', fontsize=18) 
plt.tight_layout(pad=2)

#save
fig.savefig('./hourly_daily.pdf', format='pdf', dpi = 300, bbox_inches='tight')
plt.show() 


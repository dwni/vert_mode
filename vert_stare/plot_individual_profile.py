
import matplotlib.cm as cm
import warnings
warnings.filterwarnings("ignore")
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from matplotlib import gridspec
import numpy as np
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

fig = plt.figure(figsize=(11, 3*3))
gs = gridspec.GridSpec(4, 2, height_ratios=[4, 1, 1, 1]) 
medianprops = dict(linestyle='-.', linewidth=3.5, color='black')

scaling =''

st = 19
et = 29 

y_new = np.arange(0.05, 1.05, 0.05)  

#analytical profile Lenschow (1980)
Le = 1.8*(y_new)**(2/3)*(1-(0.8*y_new))**2

datalow = lowstab
datahigh = highstab

var_fin_meanlow = []
var_fin_meanmed = []
var_fin_meanhigh = []    

ax1 = fig.add_subplot(gs[0,0])
var_fin_ = np.nanmean(var_new,0)
var_fin = np.nanmean(var_new[:,st:et+1,:],1)

for hr in datalow: 
    var_fin_meanlow.append(var_fin[hr,:])


for hr in datahigh:         
    var_fin_meanhigh.append(var_fin[hr,:])    

ax1.vlines(0,0,1,colors='grey',linestyles='dashed')
ax1.set_ylabel('z/$z_i$', fontsize=16)
ax1.legend()
ax1.tick_params(labelsize=16)
if scaling == '':
    ax1.set_xlabel('$\sigma$$^{2}_w$ / $w^2_*$', fontsize=16)  
else:
    ax1.set_xlabel('$\sigma$$^{2}_w$ / $w^2_m$', fontsize=16)  
ax1.set_xlim(-0.1, 1)        
ax1.set_ylim(0,1)
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='out')        
ax1.grid(linestyle = 'dashed')
ax1.plot(Le,y_new, '--k', lw=3, label ='Ref:L80')
ax1.plot(np.nanmean(var_fin_meanlow,0), y_new, '-b', lw=2,  label = 'Low (n='+str(len(datalow))+')')  
ax1.plot(np.nanmean(var_fin_meanhigh,0), y_new, '-r', lw=2, label = 'High (n='+str(len(datahigh))+')')
ax1.legend(loc="upper right", fontsize =9)
ax1.text(-0.09, 0.9, '(a)', fontsize=18) 
ax1.set_title("Bulk stability", fontsize = 16)

#####
datalow = lowu
datamed = medu
datahigh = highu

var_fin_meanlow = []
var_fin_meanmed = []
var_fin_meanhigh = []    

ax1 = fig.add_subplot(gs[0,1])
var_fin_ = np.nanmean(var_new,0)
var_fin = np.nanmean(var_new[:,st:et+1,:],1)

for hr in datalow: 
    var_fin_meanlow.append(var_fin[hr,:])

for hr in datamed: 
    var_fin_meanmed.append(var_fin[hr,:])
      
for hr in datahigh:         
    var_fin_meanhigh.append(var_fin[hr,:])    

ax1.vlines(0,0,1,colors='grey',linestyles='dashed')
ax1.set_ylabel('z/$z_i$', fontsize=16)
ax1.legend()
ax1.tick_params(labelsize=16)
if scaling == '':
    ax1.set_xlabel('$\sigma$$^{2}_w$ / $w^2_*$', fontsize=16)  
else:
    ax1.set_xlabel('$\sigma$$^{2}_w$ / $w^2_m$', fontsize=16)  
ax1.set_xlim(-0.1, 1)        
ax1.set_ylim(0,1)
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='out')        
ax1.grid(linestyle = 'dashed')
ax1.plot(Le,y_new, '--k', lw=3, label ='Ref:L80')
ax1.plot(np.nanmean(var_fin_meanlow,0), y_new, '-b', lw=2,  label = 'Low (n='+str(len(datalow))+')')  
ax1.plot(np.nanmean(var_fin_meanmed,0), y_new, '-', color ='darkorange', lw=2, label = 'Med (n='+str(len(datamed))+')')    
ax1.plot(np.nanmean(var_fin_meanhigh,0), y_new, '-r', lw=2, label = 'High (n='+str(len(datahigh))+')')
ax1.legend(loc="upper right", fontsize =9)
ax1.text(-0.09, 0.9, '(b)', fontsize=18) 
ax1.set_title("Friction velocity", fontsize = 16)

plt.tight_layout()
plt.savefig('./stab_fric_.pdf', format='pdf', dpi = 300, bbox_inches='tight')    
plt.show()


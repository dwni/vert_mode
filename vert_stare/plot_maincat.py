
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

datalow = rainy
datamed = cldtp
datahigh = csky

xrs = xr.open_dataset('var_merged.nc')
df3 = pd.read_feather('meteorological_param.feather')
filename = xrs.filename.values
var_new = xrs.var_m.values

y_new = np.arange(0.05, 1.05, 0.05)  

st1 = [19]
et1 = [29] 

scaling =''
data_title = ''

#analytical profile Lenschow (1980)
Le = 1.8*(y_new)**(2/3)*(1-(0.8*y_new))**2

fig = plt.figure(figsize=(5*3,4))
gs = gridspec.GridSpec(2, 3, height_ratios=[1, 1]) 
medianprops = dict(linestyle='-.', linewidth=2.5, color='black')

 
var_fin_mean = []
var_fin_meanlow = []
var_fin_meanmed = []
var_fin_meanhigh = []    
var_fin_meancs = []        

ax1 = fig.add_subplot(gs[:,0])

### Rainy days
i=0 
for hr in datalow: 
    print(filename[hr])
    if str(filename[hr])[:4]=='2021':
        if int(filename[hr]<20210807):
            inputfile_3 = "/Mast/sups_rao_mets00_l1_precip_v00_"+str(int(filename[hr]))+".nc"
        else:
             inputfile_3 = "/Mast/sups_rao_mets00_l1_precip_v01_"+str(int(filename[hr]))+".nc"           
    else:
        inputfile_3 = "/Mast/sups_rao_mets00_l1_precip_v00_"+str(int(filename[hr]))+".nc"

    rainfile = xr.open_dataset(inputfile_3)
    t_file3 = rainfile.time
    precip = rainfile.precip.values

    #averaging x- minute average  
    t_rain = []
    for i in range(0,24):
        for j in [0,30]:
            try:
                time = np.where((t_file3.dt.day==int(str(int(filename[hr]))[-2:])) & (t_file3.dt.hour==i) & (t_file3.dt.minute==j))[0][0]
            except IndexError:
                continue

            t_rain.append(time)
    t_rain = np.append(0, t_rain)
    t_rain = np.append(t_rain, 143)

    rain_idx = np.zeros(48)
    rain_idx [:]= np.nan
    for i in range(0,48):
        rain_idx[i] = np.nanmean(precip[t_rain[i]:t_rain[i+1]])


    rain = np.where(rain_idx >0 )
    v1=[1e-7, 1e-6, 1e-5, 1e-4, 1e-3]#, 1e1]  
    v = np.arange(0,3.1,0.1)
    interval = np.arange(0,24,4)

    try:
        var_new[hr, np.array(rain[0]), :] = np.nan
        if int(filename[hr]) == 20210709:
            var_new[hr, 27,:] = np.nan
            print(20210709)
    except IndexError:
        pass



    st = st1[loop]
    et = et1[loop]
    var_fin = np.nanmean(var_new[hr,st:et+1,:],0)  

    if np.isnan(np.nanmean(var_fin)) == True:
        continue
    var_fin_meanlow.append(var_fin)

# Cloudy days
st = st1[loop]
et = et1[loop]
var_fin_ = np.nanmean(var_new,0)
var_fin = np.nanmean(var_new[:,st:et+1,:],1)

for hr in datamed: 
    var_fin_meanmed.append(var_fin[hr,:])
    
    
# Clear-sky days
for hr in datahigh:    
    var_fin_meanhigh.append(var_fin[hr,:])
    
# Plot
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
ax1.set_title(data_title, fontsize=18)


### Add boxplot 
xlab='Latent heat flux [$Wm^{-2}$]'
axs = fig.add_subplot(gs[0,1])
box = axs.boxplot([df3.LE[datahigh], df3.LE[datamed],df3.LE[datalow] ], vert=False,showbox=True, patch_artist=True, medianprops =medianprops )
colors = ['blue', 'darkorange', 'red']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
axs.set_xlabel(xlab, fontsize=14) 
axs.set_xlim(0,250)
axs.tick_params(labelsize=13)
axs.set_yticklabels([])
axs.text(4, 3, '(b)', fontsize=14)  

xlab='Relative humidity [%]'
axs = fig.add_subplot(gs[1,1])
box = axs.boxplot([df3.RH[datahigh]*100, df3.RH[datamed]*100,df3.RH[datalow]*100 ], vert=False,showbox=True, patch_artist=True, medianprops =medianprops )
colors = ['blue', 'darkorange', 'red']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
axs.set_xlabel(xlab, fontsize=14) 
axs.set_xlim(0,100)
axs.tick_params(labelsize=13)
axs.set_yticklabels([])
axs.text(1.8, 2.9, '(c)', fontsize=14) 

ax1.grid(linestyle = 'dashed')
ax1.plot(Le,y_new, '--k', lw=3, label ='Ref:L80')
ax1.plot(np.nanmean(var_fin_meanhigh,0), y_new, '-b', lw=2,  label = 'Clear-sky(n='+str(len(datahigh))+')')  
ax1.plot(np.nanmean(var_fin_meanmed,0), y_new, '-', color ='darkorange', lw=2, label = 'Cloud-topped(n='+str(len(datamed))+')')    
ax1.plot(np.nanmean(var_fin_meanlow,0), y_new, '-r', lw=2, label = 'Rainy(n='+str(len(datalow))+')')
ax1.text(-0.09, 0.93, '(a)', fontsize=14)                         
ax1.legend(loc="upper right", fontsize =10)#
ax1.grid(linestyle = 'dashed')

#### Add plot from radiosonde 
ax2 = fig.add_subplot(gs[:,2])
ax2.set_xlim(0,100)  
ax2.plot(np.nanmean(relh_interp[:11],0), y/1000, 'blue', label = 'clear-sky')  #11   24
ax2.plot(np.nanmean(relh_interp[11:69],0), y/1000, 'darkorange', label = 'cloud-topped')   #11-69    24-39
ax2.plot(np.nanmean(relh_interp[69:],0), y/1000, 'red', label = 'rainy')     #69   39
ax2.set_xlabel('Relative humidity [%]', fontsize=15)
ax2.set_ylabel('Height [km]', fontsize=15)    
ax2.grid(linestyle = 'dashed')
ax2.tick_params(labelsize=13)  
ax2.text(2., 2.7, '(d)', fontsize=14) 
plt.tight_layout()

#### save
plt.savefig('./threemain_.pdf', format='pdf', dpi = 300, bbox_inches='tight')    
plt.show()


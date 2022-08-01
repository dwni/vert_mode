
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

hrr = 0
n =  len(filename)-hrr
colors = plt.cm.brg(np.linspace(0,1,n))
scaling =''

y_new = np.arange(0.05, 1.05, 0.05)  

#analytical profile Lenschow (1980)
Le = 1.8*(y_new)**(2/3)*(1-(0.8*y_new))**2


st = 19 
et = 29

data2020 = [32,33,31,34,35,36,]
data2021 = [67,68,69,70,71,]
datacs = [32,33,31,34,35,36,67,68,69,70,71,]

fig = plt.figure(figsize=(6*2, 5))
gs = gridspec.GridSpec(1, 2, width_ratios=[1,1] ) 

var_fin_mean = []
var_fin_2020 = []
var_fin_2021 = []

var_fin_ = np.nanmean(var_new,0)
var_fin = np.nanmean(var_new[:,st:et+1,:],1)

ax1 = fig.add_subplot(gs[0])
nn = []
colors = plt.cm.brg(np.linspace(0,1,n))

i=0 
for hr in data2020: 
    if np.isnan(np.nanmean(var_fin[hr,:])) == True:
        continue
    hhri = hhri+1
    print('data2020',filename[hr],hr)       
    ax1.plot(var_fin[hr,:], y_new, '--', color='darkorange', alpha=.7)
    var_fin_2020.append(var_fin[hr,:])



for hr in data2021: 
    hhri = hhri+1
    print('data2021',filename[hr],hr)   
    ax1.plot(var_fin[hr,:], y_new, '--', color='b', alpha=.7)
    var_fin_2021.append(var_fin[hr,:])


       
for hr in datacs: 
    hhri = hhri+1
    print('datacs',filename[hr],hr)           
    var_fin_mean.append(var_fin[hr,:])    

#plotting    
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

ax1.plot(Le,y_new, '--k', lw=3, label ='Ref:L80')
ax1.plot(var_fin[32,:], y_new, '--', color='darkorange', alpha=1, label = '2020-daily case') 
ax1.plot(var_fin[67,:], y_new, '--', color='b', alpha=1,  label = '2021-daily case') 
ax1.plot(np.nanmean(var_fin_2020,0), y_new, '-', color ='darkorange', lw=2, label = '2020 mean(n='+str(len(data2020))+')')    
ax1.plot(np.nanmean(var_fin_2021,0), y_new, '-b', lw=2,  label = '2021 mean (n='+str(len(data2021))+')')  
ax1.plot(np.nanmean(var_fin_mean,0), y_new, '-r', lw=2, label = 'Clear-sky mean (n='+str(len(datacs))+')')
ax1.legend(fontsize =10, bbox_to_anchor=(0.6, 0.595), )# loc="upper right"
ax1.grid(linestyle = 'dashed')
ax1.text(-0.08, 0.92, '(a)', fontsize=18) 


##############
# meteo param
##############
x3 = [ df3.RH, df3.LE,  df3.H, 
      df3.windspeed,  df3['w*'].values , df3['T'],
      df3.ustar, df3.zeta, df3.soilmoist , 
      df3.Bowen_ratio, df3.var_mean_n, (df3.H+df3.LE),  
      df3.Theta_s,  df3.Theta_e_marq, df3.Theta_l]
y3 = [df3.Bowen_ratio,  df3.Cloud_fraction ,]

ylab = ['Bowen ratio', 'Cloud Fraction'] 
xlab = ['RH', 'LE',  'H', 
        'Windspeed', '$w_*$', 'T',         
        '$u_*$','-zi/L', 'soil moisture',
        'Bowen ratio', '$\sigma$$^{2}_w$/$w^2_*$','[H+LE]',
        'Theta_e', 'Theta_e_marq', 'Theta_l']

xi=[]
yi=[]

data2020=np.array([32,33,31,34,35,36,])
data2021=np.array([67,68,69,70,71,])
print(df3.Date[data2020])
print(df3.Date[data2021])

ax2 = fig.add_subplot(gs[1])
ylimmin = 0
ylimmax = 3
ax2.set_ylim(ylimmin, ylimmax)  
ax2.set_xlim(0, 300)     
cs = ax2.scatter(np.array(df3.LE[data2020]), np.array(df3.Bowen_ratio[data2020]), s=80, c='darkorange', label ='LHF - 2020')#np.arange(0,len(np.array(x[i])[CF03]),1),  cmap= cm.rainbow) #np.array(f.var_mean_n)[CF03],
cs = ax2.scatter(np.array(df3.LE[data2021]), np.array(df3.Bowen_ratio[data2021]), s=80, c= 'darkblue', label ='LHF - 2021')#np.arange(0,len(np.array(x2[i])[CF03_2021]),1), cmap= cm.rainbow, label ='2021') #np.array(data2021.var_mean_n)[CF03_2021]

ax2.set_ylim(ylimmin, ylimmax)    
cs = ax2.scatter(np.array(df3.H[data2020]), np.array(df3.Bowen_ratio[data2020]), s=80, c='darkorange', marker="v", label ='SHF - 2020')#np.arange(0,len(np.array(x[i])[CF03]),1),  cmap= cm.rainbow) #np.array(f.var_mean_n)[CF03],
cs = ax2.scatter(np.array(df3.H[data2021]), np.array(df3.Bowen_ratio[data2021]), s=80, c= 'darkblue', marker="v",  label ='SHF - 2021')#np.arange(0,len(np.array(x2[i])[CF03_2021]),1), cmap= cm.rainbow, label ='2021') #np.array(data2021.var_mean_n)[CF03_2021]

ax2.set_xlabel("$W m^{-2}$" , fontsize=16) 
ax2.set_ylabel("Bowen ratio" , fontsize=16) 
ax2.tick_params(labelsize=16) 
ax2.legend( fontsize =10, bbox_to_anchor=(0.7, 1.01),)
ax2.grid(linestyle = 'dashed')   
ax2.text(1, 2.75, '(b)', fontsize=18) 
plt.tight_layout(pad=3)
plt.show
fig.savefig('./Clear-sky.pdf', format='pdf', dpi = 300, bbox_inches='tight')





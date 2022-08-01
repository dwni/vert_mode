
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
import matplotlib.cm as cm


xrs = xr.open_dataset('var_merged.nc')
df3 = pd.read_feather('meteorological_param.feather')
filename = xrs.filename.values
var_new = xrs.var_m.values
y_new = np.arange(0.05, 1.05, 0.05)  

case = cldtp

x3 = [ df3.RH*100, df3.LE,  df3.H+df3.LE, 
      df3.windspeed, df3.Bowen_ratio, df3['T'],
      df3.ustar, df3.zeta, df3.Cloud_fraction ,
      df3['w*'].values , df3.soilmoist ]

y3 = [df3.var_max_n]

ylab = ['$\sigma$$^{2}_w$/$w^2_*$'] 
xlab = ['Relative humidity [%]', 'LHF [W$\,$m$^{-2}$]',  'SHF+LHF [W$\,$m$^{-2}$]', 
        'Wind speed [m$\,$s$^{-1}$]', 'Bowen ratio', 'Temperature [K]',         
        '$u_*$','-zi/L', 'Cloud fraction']
xi=[]
yi=[]

ylimmin = -0.05
ylimmax = 1.05

fig, axs = plt.subplots(3,3 , figsize=(9, 10))
for i in range(0,9):
    if  i == 0:
        xi = 0 
        yi = 0
        axs[xi,yi].set_xlim(0.,100)
        axs[xi,yi].set_ylim(ylimmin, ylimmax)
        axs[xi,yi].set_ylabel(ylab[0], fontsize=16)              
    elif i ==1: #LE
        xi = 0
        yi = 1            
        axs[xi,yi].set_xlim(0, 250)
        axs[xi,yi].set_ylim(ylimmin, ylimmax)      

    elif i ==2: #H
        xi = 0
        yi = 2            
        axs[xi,yi].set_xlim(0, 400)
        axs[xi,yi].set_ylim(ylimmin, ylimmax)
    ##
    elif i ==3: #wspd
        xi = 1
        yi = 0            
        axs[xi,yi].set_xlim(0, 7)            
        axs[xi,yi].set_ylim(ylimmin, ylimmax)

    elif i ==4:#w*
        xi = 1
        yi = 1
        axs[xi,yi].set_xlim(0, 5)#(0, 250)           
        axs[xi,yi].set_ylim(ylimmin, ylimmax)

    elif i ==5:#T
        xi = 1
        yi = 2            
        axs[xi,yi].set_xlim(285, 310)
        axs[xi,yi].set_ylim(ylimmin, ylimmax) 
    ##
    elif i ==6:#U*
        xi = 2
        yi = 0            
        axs[xi,yi].set_xlim(0., 0.5) 
        axs[xi,yi].set_ylim(ylimmin, ylimmax)   

    elif i ==7:#stab
        xi = 2
        yi = 1            
        axs[xi,yi].set_xlim(0, 250) #(80, 340)
        axs[xi,yi].set_ylim(ylimmin, ylimmax)

    elif i ==8:#qsoil
        xi = 2
        yi = 2            
        axs[xi,yi].set_xlim(-0.1, 1.1)        
        axs[xi,yi].set_ylim(ylimmin, ylimmax) 

    cs = axs[xi,yi].scatter(x3[i][case], y3[0][case], c=df3.RH[case]*100, vmin = 30, vmax = 60, cmap= cm.viridis, label ='merge') #, vmax=1
    axs[xi,yi].set_xlabel(xlab[i], fontsize=12) 
    axs[xi,0].set_ylabel(ylab[0], fontsize=12) 

    coefficient_of_dermination = np.corrcoef(x3[i][case], y3[0][case])[0,1]
    axs[xi,yi].text(0.5, 0.85, 'R ='+str(np.round(coefficient_of_dermination,2)), transform=axs[xi,yi].transAxes, fontsize=14) 
    axs[xi,yi].tick_params(labelsize=12)                
    axs[xi,yi].grid(linestyle = 'dashed')            


plt.tight_layout(pad=3,w_pad=2, h_pad=4)
cbar =fig.colorbar(cs, ax=axs, pad = 0.08, shrink=0.6, orientation='horizontal')
cbar.ax.tick_params(labelsize=14)
cbar.set_ticks([30,40,50,60])
cbar.set_ticklabels([30,40,50,60])
cbar.set_label('Relative humidity [%]', fontsize=13)

plt.savefig('./scatterplot_.pdf', format='pdf', dpi = 300, bbox_inches='tight')    
plt.show()


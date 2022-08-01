
######################################################
#  Calculate normalized variance and 
#  daily averaged meteorological parameters
#  Required input file :
#  1. L2 vertical stare
#  2. Temperature
#  3. Pressure
#  4. Surface heat flux
######################################################

import numpy as np
import glob
import os
import xarray as xr
Import pandas as pd
from scipy import interpolate
from utils.eq_func import time2slice
from metpy.units import units
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Obukhov length
def obukovL(ustar, theta, Q0):
    L = - (ustar**3)*theta/(0.4*9.81*Q0)
    return L

g = 9.8
rho = 1.21
cp = 1005
##
scaling = ''

## clear sky 
filename = [20210626,]
snr_th = ['0p994',]

ntime = 48
y_new = np.arange(0.05, 1.05, 0.05)            
var_new = np.zeros((len(filename), ntime, len(y_new)))
var_new[:,:,:] = np.nan
Var_ = np.zeros((len(filename), ntime,100))
Var_[:,:,:] = np.nan  
var_new_n = np.zeros((len(filename), ntime, len(y_new)))
var_new_n[:,:,:] = np.nan
Var_scale_n = np.zeros((len(filename), ntime,100))
Var_scale_n[:,:,:] = np.nan  

for nf in range(len(filename)):
    print(nf, filename[nf])
    xrs2 = xr.open_dataset( './202105/fig/hom_SNR'+str(snr_th[nf])+'_'+str(filename[nf])+'_hom_30.nc')
    z1 = xrs2.mld.values

    # convective velocity scale
    hfss = np.zeros(48)
    hfls = np.zeros(48)    
    ustar = np.zeros(48)
    hfss[:] = np.nan
    hfls[:] = np.nan
    ustar [:] = np.nan

    hf_data = xr.open_dataset('/Mast/sups_rao_turb00_l2_any_v00_'+str(filename[nf])+'.nc')
    hfss[:-1] = hf_data.hfss[1:].values
    hfls[:-1] = hf_data.hfls[1:].values    
    ustar[:-1] = hf_data.ustar[1:].values     
        
    kinematik_heat_flux = hfss / (1.21 * 1005)

    try:
        T_data = xr.open_dataset('/Mast/sups_rao_mett00_l1_any_v00_'+str(filename[nf])+'.nc')
    except FileNotFoundError:
        T_data = xr.open_dataset('/Mast/sups_rao_mett00_l1_any_v01_'+str(filename[nf])+'.nc')      
    T = T_data.ta.values[:, 0]
    Rhum = T_data.hur.values[:,0]
    
    try:
        P_data = xr.open_dataset('/Mast/sups_rao_mets00_l1_pa_v00_'+str(filename[nf])+'.nc')
    except FileNotFoundError:
        P_data = xr.open_dataset('/Mast/sups_rao_mets00_l1_pa_v01_'+str(filename[nf])+'.nc')

    p = P_data.pa.values[:]/100
    #
    mixing_ratio = mpcalc.mixing_ratio_from_relative_humidity(p*units.hPa, T_data.ta.values[:, 0]*units.degK, Rhum)
    Tv = np.array(mpcalc.virtual_temperature(T_data.ta.values[:, 0]*units.degK, mixing_ratio))
            
    #
    Tv30min = []    
    rh30min = []  
    #average from 10minutes to 30minutes
    for i in range(0,24):
        # minutes 10 to 30
        try:
            id1 = np.where((T_data.time.dt.hour==i) & (T_data.time.dt.minute==10))[0]
            id2 = np.where((T_data.time.dt.hour==i) & (T_data.time.dt.minute==30))[0]
        except IndexError:
            id1 = np.where((T_data.time.dt.hour==i) & (T_data.time.dt.minute==10))[0]

        Tv30min.append(np.nanmean(Tv[id1[0]:id2[0]+1]))
        rh30min.append(np.nanmean(Rhum[id1[0]:id2[0]+1])) 
        # minutes 40 to 0
        id3 = np.where((T_data.time.dt.hour==i) & (T_data.time.dt.minute==40))[0]
        if i+1 !=24:
            i=i
        else:
            i=-1
        id4 = np.where((T_data.time.dt.hour==i+1) & (T_data.time.dt.minute==0))[0]  

        Tv30min.append(np.nanmean(Tv[id3[0]:id4[0]+1]))        
        rh30min.append(np.nanmean(Rhum[id3[0]:id4[0]+1])) 


    buoy = (hfss+(0.07*hfls))/(rho*cp)
    w_star = ((9.81/np.array(Tv30min))*buoy*z1)**(1/3)

    if scaling == 'MS':
        print('MS Scaling')
        for i in range(48):
            Var_scale[nf,i,:] = xrs2.w2.values[i,:]/(((w_star[i]**3)**(2/3))+((5*ustar[i]**3)**(2/3))) # (w_star[i]**2)+(4*ustar[i]**2)) # 
    else:
        print('w* scaling')
        for i in range(48):
            Var_scale [nf, i, :] = xrs2.w2.values[i,:]/(w_star[i]**2)   


    Var_[nf,:,:] = xrs2.w2[:,:]
    
    hh = np.tile(height, (48,1)).T
    z_zi = (hh/z1).T
    y_old = z_zi

    for it in range(48):        
        f_var = interpolate.interp1d( y_old[it, :], Var_[nf,it, :] )
        f_var_n = interpolate.interp1d( y_old[it, :], Var_scale[nf,it, :] )        
        try:
            var_new[nf, it, :] = f_var(y_new)
            var_new_n[nf, it, :] = f_var_n(y_new)            

        except ValueError:
            var_new[nf, it, :] = 0
            var_new_n[nf, it, :] = 0            
     
                
    var_mean[nf, 0] = np.nanmean(var_new[nf, st:et+1, 4:12])   
    var_mean_n[nf, 0] = np.nanmean(var_new_n[nf, st:et+1, 4:12]) 
    
    #======================================
    # Cloud fraction
    #======================================
    # averaging the data to 60 minute average  
    t     = time2slice(date, 60)        
    cdate = date[t[1:]]

    # open variable
    cloud_1h     = np.zeros((len(t)-1))
    cloud_1h [:] = np.nan

    n_data       = np.zeros((len(t)-1))
    n_data [:]   = np.nan

    cloud_layer  = []
    
    #
    cloud_layer = deepcopy(beta)
    cloud_index = np.where(np.log10(beta)<-4) # find the cloud
    cloud_layer[cloud_index] = np.nan # set other than cloud to nan
    cloud_layer_nan = np.nanmean(cloud_layer,1) # average over height
    
    # count how many cloud within 1 hour
    for i in range(1,len(t)):
        cloud_1h[i-1] = np.count_nonzero(~np.isnan(cloud_layer_nan [t[i-1]:t[i]]))
        n_data[i-1] = len(cloud_layer_nan[t[i-1]:t[i]])

    # cloud fraction 
    cf = cloud_1h/n_data          
    
    
    #=======================================================
    # precipitation
    #=======================================================    
    if str(filename[nf])[:4]=='2021':
        if int(filename[nf]<20210807):
            #load precipitation file
            inputfile_3 = "/Mast/sups_rao_mets00_l1_precip_v00_"+str(int(filename[nf]))+".nc"
        else:
             inputfile_3 = "/Mast/sups_rao_mets00_l1_precip_v01_"+str(int(filename[nf]))+".nc"           
    else:
        inputfile_3 = "/Mast/sups_rao_mets00_l1_precip_v00_"+str(int(filename[nf]))+".nc"

    rainfile = xr.open_dataset(inputfile_3)
    t_file3 = rainfile.time
    precip = rainfile.precip.values
 
    #averaging x- minute average  
    t_rain = []
    for i in range(0,24):
        for j in [0,30]:
            try:
                time = np.where((t_file3.dt.day==int(str(int(filename[nf]))[-2:])) & (t_file3.dt.hour==i) & (t_file3.dt.minute==j))[0][0]
            except IndexError:
                time = np.nan
                continue

            t_rain.append(time)
    t_rain = np.append(0, t_rain)
    t_rain = np.append(t_rain, len(t_file3)-1)  
    
    rain_m = np.zeros(48)
    rain_m [:]= np.nan
    for i in range(0,48):
        rain_m[i] = np.nanmean(precip[t_rain[i]:t_rain[i+1]])
 
    rain_idx = np.where(rain_m >0 )[0]
    

    #precipitation masked
    try:
        ustar[rain_idx] = np.nan
        Tv30min[rain_idx] = np.nan 
        kinematik_heat_flux[rain_idx] = np.nan    
        hfls[rain_idx] = np.nan   
        hfss[rain_idx] = np.nan   
        T30min[rain_idx] = np.nan   
        z1[rain_idx] = np.nan     
        wstar[rain_idx] = np.nan   
        ustar[rain_idx] = np.nan   
        rh30min[rain_idx] = np.nan   
        windspeed[:precip.shape[0]][precip>0] = np.nan

    except TypeError:
        pass
              
    # Calculate Obukov L
    L = obukovL(ustar, Tv30min, kinematik_heat_flux)      
    
    #=======================================================            
    # Calculate daily mean value
    #=======================================================  
    st = 19
    et = 29
    # -z/L
    LObukov[nf]= np.nanmean(L[st:et+1])
    zetamean[nf] = np.nanmean(z1[st:et+1])/-LObukov[nf] #-z1/L
    
    # mixed layer height
    zimean[nf] = np.nanmean(z1[st:et+1])
    
    # w* mean
    wstarmean[nf] = np.nanmean(w_star[st:et+1])
    
    # heat flux 
    HF[nf] = np.nanmean(hfss[st:et+1])
    LE[nf] = np.nanmean(hfls[st:et+1])
    bowen[nf] = np.nanmean(hfss[st:et+1]/hfls[st:et+1])
    
    # u*
    ustarmean[nf] = np.nanmean(ustar[st:et+1])
    
    # rh_mean
    rhmean[nf] = np.nanmean(rh30min[st:et+1])
    
    if (st+1)%2==0:
        ids_ws = np.where(T_data.time.dt.hour.values ==(st+1)/2)[0][0] #11##
    else:
        ids_ws = np.where((T_data.time.dt.hour.values ==(st)/2 )& (T_data.time.dt.minute.values == 30))[0][0] #14##)[0][0] #11##

    if (et+1)%2==0:      
        ide_ws = np.where(T_data.time.dt.hour.values ==(et+1)/2)[0][0] 
    else:
        ide_ws = np.where((T_data.time.dt.hour.values ==(et)/2) & (T_data.time.dt.minute.values == 30))[0][0] #14##

   
    cst = int((st+1)/2)
    cet = int((et+1)/2)

    # wind speed
    wsmean[nf] = np.nanmean(windspeed[ids_ws:ide_ws+1],0)
    cloud_frac[nf] = np.nanmean(cf[cst:cet+1],0)     
    
print("save to nc file....") 
ds = xr.Dataset(data_vars={"filename":(filename[:]),
                           "norm_w2":(["idx","timestep","height"],var_new[:,:,:]),
                          },
            coords={"idx": np.arange(0,len(filename),1),
                    "timestep": np.arange(0,48,1),
                    "height": ( y_new)
                    })


#=======================================================            
# Output to the ncfile (normalized variance)
#=======================================================   
# add variable attribute metadata
ds['norm_w2'].attrs={'units':'m2/s2', 'long_name':'all hours'}
ds['height'].attrs={'units':'m', 'long_name':'height'}
outfile = 'var_2021_allcase_v1.nc'
if os.path.isfile(outfile):
    os.remove(outfile)   
ds.to_netcdf(outfile)


#=======================================================            
# Output to the table (average values)
#=======================================================        
df3 = pd.DataFrame(
              {"Date" : filename,
               "zeta" : np.round(zetamean,2),
               "zi" : zimean.astype(int),
               "windspeed" : np.round(wsmean,2),
               "w*" : np.round(wstarmean,2),
               "H" : np.round(HF,2),
               "LE" : np.round(LE,2),
               "L" : np.round(LObukov,2),
               "ustar" : np.round(ustarmean,2),
               "Cloud_fraction" : np.round(cloud_frac,4),
               "Bowen_ratio" : np.round(bowen,2),
               "RH" : np.round(rhmean,2),
               "var_mean_n" : np.round(var_mean_n[:,0],5),
               "var_mean" : np.round(var_mean[:,0],5),   
               })

df3.to_feather(‘meteo.feather')

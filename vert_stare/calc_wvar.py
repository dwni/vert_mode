#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import glob
import os
import xarray as xr
import pandas as pd
import warnings
warnings.filterwarnings(action='ignore', message='Mean of empty slice')
from netCDF4 import Dataset

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from copy import deepcopy

from utils.eq_func import time2slice
from utils.histat import histats
#=========================================================================

# Input variable

folderpath = ''
filename=[]
for fn in sorted(glob.glob(os.path.join('/L1_2020/', 'sups_rao*'))):
    filename.append(str(fn)[-11:-3])
    #filename.append(fn)
filename= filename[21:]
snr_th = 1.005

for nf in range(len(filename)):
    #=============================================================================================
    #   open doppler lidar vertical stare
    #=============================================================================================
    if str(filename[nf])[:4] == '2020':
        if int(filename[nf]) < 20200722:
            inputfile = '/L1_2020/sups_rao_dlidST161_l1_any_v04_'+str(int(filename[nf]))+'.nc'
        else:
            inputfile = '/L1_2020/sups_kit_dlidST161_l1_any_v00_'+str(int(filename[nf]))+'.nc'    
    else:
        inputfile = '/L1/fval_fmi_dlidST00_l1_any_v00_'+str(int(filename[nf]))+'.nc'

    
    xrs = xr.open_dataset(inputfile)
    radvel = deepcopy(xrs.dv.values)
    beta = xrs.beta.values    
    intensity = xrs.intensity.values
    date = xrs.time
    height = xrs.range.values   

    if str(filename[nf])[:4] == '2020':    
        radvel[xrs.intensity.values<snr_th] = np.nan
        radvel[:,:4]=np.nan

        beta[xrs.intensity.values<snr_th] = np.nan
        beta[:,:4]=np.nan

        intensity[xrs.intensity.values<snr_th] = np.nan
        intensity[:,:4]=np.nan

        for test in [5]:#def 5
            beta_snr = deepcopy(beta)
            radvel_data = deepcopy(radvel)
            inten_snr = deepcopy(intensity)
            d_raw =deepcopy(radvel)
            n = 4 #def 4
            for i in range(radvel_data.shape[0]):
                for j in range(radvel_data.shape[1]):
                    if j < n:
                        sj=0
                    else:
                        sj=j-n 

                    if i <n:
                        si=0
                    else:
                        si=i-n

                    if i+n>(radvel_data.shape[0]):
                        ei = radvel_data.shape[0]
                    else:
                        ei =  i+n

                    if j+n>(radvel_data.shape[0]):
                        ej = radvel_data.shape[0]
                    else:
                        ej =  j+n    


                    d=d_raw[si:ei, sj:ej]

                    c = np.abs(d-d_raw[i,j])
                    c_idx = np.where(c>test)[0]

                    if len(c_idx)>0:
                        radvel_data[i,j]=np.nan
                        beta_snr[i,j]=np.nan  
                        inten_snr[i,j]=np.nan     

    else:
        radvel[xrs.intensity.values<snr_th] = np.nan
        radvel[:,:2]=np.nan

        beta[xrs.intensity.values<snr_th] = np.nan
        beta[:,:2]=np.nan

        intensity[xrs.intensity.values<snr_th] = np.nan
        intensity[:,:2]=np.nan        
        
        radvel_data = radvel
    variance, t_skew,  = histats(date,height,radvel_data,30,40)

    varnan = np.full((1, len(height)), np.nan)
    if t_skew.shape[0]<48:
        if t_skew.dt.hour[0]!=0:
            for i in range (0, 48-t_skew.shape[0]):
                t_skew = np.append(t_skew[0]-np.timedelta64(30,'m'),t_skew)
                variance = np.append(varnan, variance, axis=0)
                
        else:
            print('BREAK!')
            break
            
            
    t=time2slice(date, 1)        
    data1min=deepcopy(radvel_data[:,:])
    wdl=np.zeros((len(t)-1, data1min.shape[1]))
    wdl[:,:]=np.nan
    
    for j in range(data1min.shape[1]):
        for i in range (len(t)-1):
            st=t[i]
            et=t[i+1]
            wdl[i,j]=np.nanmean(data1min[st:et,j])
 
    w=wdl[:,:]
    time2=date.values[t[1:]]          
    
    # mixing height layer
    zi = np.zeros((48))
    
    # set threshold
    if str(filename[nf])[:4] == '2020':
        mlh_th = [0.2**2]
        # find the first layer which variance below the threshold
        for hh in range(len(t_skew)):
            a = np.where(variance[hh,:] < mlh_th)
            if a[0].size == 0: #prevent during night time, if can't detect the first layer, the previous value will be taken
                zi[hh] = zi[hh-1]
            else:
                zi[hh] = height[int(a[0][0])]
        mlh = zi[:] 
    else:
        mlh_th = [0.3**2]
        # find the first layer which variance below the threshold
        for hh in range(len(t_skew)):
            a = np.where((variable_var[hh,2:] < mlh_th)| (np.isnan(variable_var[hh,2:])==True))
            if a[0].size == 0: #prevent during night time, if can't detect the first layer, the previous value will be taken
                zi[hh] = zi[hh-1]
            else:
                zi[hh] = height[int(a[0][0])]
        mlh = zi[:] 

    now = datetime.now()

    print("save to nc file...."+str(filename[nf]))
    ds = xr.Dataset(data_vars={"w2":(["time","height"],variance),
                               "mld":(["time"],mlh),
                               "w":(["time2","height"],w),
                               "lat":(np.float(xrs.lat.values)),
                               "lon":(np.float(xrs.lon.values))},
                            
                coords={"time": (["time"], t_skew),
                        "height": (["height"], height),
                        "time2": (["time2"], time2),                        
                        
                        })
    
    # add variable attribute metadata
    if str(filename[nf])[:4] == '2020'  :  
        outfile = 'fval_kit_dlidST00_l2_any_v00_'+str(filename[nf])+'.nc'
        ds.attrs['history'] = str(now)
        ds.attrs['title'] = 'L2 products from vertical stare Doppler lidar 161'
        ds.attrs['description'] = ''   
        ds.attrs['institution'] = ''
        ds.attrs['source'] = 'HALO Photonics XR Doppler lidar (production number: 161)'
        ds.attrs['Author'] = 'Noviana Dewani (dewani@iau.uni-frankfurt.de)'    
        ds.attrs['latitude'] = '52.167206'
        ds.attrs['longitude'] = '14.122907'     
        ds['mld'].attrs={'standard_name':'atmosphere_boundary_layer_thickness', 'units':'m', 'long_name':'mixing layer height',  'comments':'mixing layer height is estimated at the first level where the variance below 0.04 m2/s2'}      
    else:
        outfile = 'fval_fmi_dlidST00_l2_any_v00_'+str(filename[nf])+'.nc'
        ds.attrs['history'] = str(now)
        ds.attrs['title'] = 'L2 products from vertical stare Doppler lidar 146'
        ds.attrs['description'] = ''   
        ds.attrs['institution'] = ''
        ds.attrs['source'] = 'HALO Photonics XR Doppler lidar (production number: 146)'
        ds.attrs['Author'] = 'Noviana Dewani (dewani@iau.uni-frankfurt.de)'    
        ds.attrs['latitude'] = '52.167152'
        ds.attrs['longitude'] = '14.122753'
        ds['mld'].attrs={'standard_name':'atmosphere_boundary_layer_thickness', 'units':'m', 'long_name':'mixing layer height', 'comments':'mixing layer height is estimated at the first level where the variance below 0.09 m2/s2'}    

        ds['w'].attrs={'standard_name':'upward_air_velocity', 'units':'m/s', 'long_name':'vertical velocity ','comments':'one-minute averaged of vertical velocity '}    
    ds['w2'].attrs={'standard_name':'variance of w', 'units':'m2/s2', 'long_name':'variance of vertical velocity','comments':'variance of vertical velocity over 30-minutes resolution corrected using Lenschow et al. (2000) method'}    
    ds['lat'].attrs={'standard_name' : 'latitude',  'comments' : 'latitude of sensor' , 'long_name' : 'latitude', 'units': 'degrees_north'}
    ds['lon'].attrs={'standard_name' : 'longitude','comments' : 'longitude of sensor', 'long_name' : 'longitude', 'units' :'degrees_east'}
    ds['time'].attrs={'long_name': 'time','comments':'thirty minutes interval'}
    ds['time2'].attrs={'long_name': 'time','comments':'one minute interval'}
    ds['height'].attrs={'units':'m', 'long_name':'height'}
    ds.to_netcdf(folderpath+outfile)


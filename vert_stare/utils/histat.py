import numpy as np
import xarray as xr
from scipy import interpolate, signal
from scipy.optimize import curve_fit
from VS.eq_func import time2slice, M_star,dtrend
import os


def histats(date, height, dv_filt, nt, num_lag):
    t_idx = time2slice (date,nt)
    nr = len(height)
    ndata_avg = np.zeros((len(t_idx))) 

    
    for nd in range(len(t_idx)-1): #number of data in one-time window
        ndata_avg[nd] = (t_idx[nd+1]-t_idx[nd])
        
    ndata= np.ceil(np.nanmax(ndata_avg)) 
    ndata = int(ndata)+1

    M11= np.zeros([len(t_idx)-1, ndata, 50, nr])
    M11[:,:,:,:]=np.nan

    M21= np.zeros([len(t_idx)-1, ndata, 50, nr])
    M21[:,:,:,:]=np.nan

    M22= np.zeros([len(t_idx)-1, ndata, 50, nr])
    M22[:,:,:,:]=np.nan

    M31= np.zeros([len(t_idx)-1, ndata, 50, nr])
    M31[:,:,:,:]=np.nan

    rand_error = np.zeros([len(t_idx)-1,nr])
    rand_error[:,:] = np.nan    

    integral_tscale = np.zeros([len(t_idx)-1,nr])
    integral_tscale[:,:] = np.nan        
    
    wwind = np.zeros([len(t_idx)-1,nr])
    wwind[:,:] = np.nan

    w_prime =np.zeros([len(dv_filt), nr])
    w_prime[:,:] = np.nan

    dv_filt2 =np.zeros([len(dv_filt), nr])
    dv_filt2[:,:] = np.nan

    Vr_f= np.zeros([len(t_idx)-1,nr])
    Vr_f[:,:]=np.nan

    dv_filt2 = dv_filt
    

    for tt in range(1,len(t_idx)):
        stime = t_idx[tt-1] # start time
        etime = t_idx[tt] # end time

        wwind[tt-1,:] = np.nanmean(dv_filt2[stime:etime, :], axis=0) 
        w_prime[stime:etime,:] = dv_filt2[stime:etime, :] - wwind[tt-1,:] 
        Vr_f[tt-1,:] = np.nanmean(w_prime[stime:etime,:]**2, axis=0) 
        

        for se in range(0, (etime-stime)):#-num_lag):
            sstime = stime+se
            for nlag in [0]:
                M11[tt-1, se, nlag, :] = (w_prime[sstime, :])*(w_prime[sstime, :]) #M11(t)
                M21[tt-1, se, nlag, :] = (w_prime[sstime, :]**2)*w_prime[sstime, :]
                M22[tt-1, se, nlag, :] = (w_prime[sstime, :]**2)*(w_prime[sstime, :]**2)
                M31[tt-1, se, nlag, :] = (w_prime[sstime, :]**3)*w_prime[sstime, :]
                
        for se in range(0, (etime-stime)-50):#-num_lag):
            sstime = stime+se        
            for nlag in range(1,50): #lag
                M11[tt-1, se, nlag, :] = (w_prime[sstime, :])*(w_prime[sstime+nlag, :]) #M11(t)
                M21[tt-1, se, nlag, :] = (w_prime[sstime, :]**2)*w_prime[sstime+nlag, :]
                M22[tt-1, se, nlag, :] = (w_prime[sstime, :]**2)*(w_prime[sstime+nlag, :]**2)
                M31[tt-1, se, nlag, :] = (w_prime[sstime, :]**3)*w_prime[sstime+nlag, :]


    M11_0 = np.nanmean(M11,axis = 1) #tt, nlag, range
    M21_0 = np.nanmean(M21,axis = 1) #tt, nlag, range
    M22_0 = np.nanmean(M22,axis = 1) #tt, nlag, range
    M31_0 = np.nanmean(M31,axis = 1) #tt, nlag, range
    

    interp_result_11 = np.zeros([len(t_idx)-1,num_lag,nr])
    interp_result_21 = np.zeros([len(t_idx)-1,num_lag,nr])
    interp_result_22 = np.zeros([len(t_idx)-1,num_lag,nr])
    interp_result_31 = np.zeros([len(t_idx)-1,num_lag,nr])

    M11_to_0 = np.zeros([len(t_idx)-1, nr])
    M21_to_0 = np.zeros([len(t_idx)-1, nr])
    M22_to_0 = np.zeros([len(t_idx)-1, nr])
    M31_to_0 = np.zeros([len(t_idx)-1, nr])

    nn_ind = np.arange(0,100,1) #number of lag

    w_11 = np.zeros([len(t_idx)-1, nr]) #for M11(->0)
    w_21 = np.zeros([len(t_idx)-1, nr]) #for M21(->0)
    w_22 = np.zeros([len(t_idx)-1, nr]) #for M22(->0)
    w_31 = np.zeros([len(t_idx)-1, nr]) #for M22(->0)

    C_11 = np.zeros([len(t_idx)-1, nr])
    C_21 = np.zeros([len(t_idx)-1, nr])
    C_22 = np.zeros([len(t_idx)-1, nr])
    C_31 = np.zeros([len(t_idx)-1, nr])

    print("calculate high order moments")
    for j in range(nr):
        for tt in range(1,len(t_idx)):
            stime = t_idx[tt-1]
            etime = t_idx[tt]
            m_input_11 = M11_0[tt-1, :num_lag ,j]
            m_input_21 = M21_0[tt-1, :num_lag ,j]
            m_input_22 = M22_0[tt-1, :num_lag ,j]
            m_input_31 = M31_0[tt-1, :num_lag ,j]

            good_ind_11 = np.where(~np.isnan(m_input_11)) #searching good data index
            good_ind_21 = np.where(~np.isnan(m_input_21))
            good_ind_22 = np.where(~np.isnan(m_input_22))
            good_ind_31 = np.where(~np.isnan(m_input_31))

            m_input_good_11 = m_input_11[good_ind_11]
            m_input_good_21 = m_input_21[good_ind_21]
            m_input_good_22 = m_input_22[good_ind_22]
            m_input_good_31 = m_input_31[good_ind_31]

            nn_ind_good_11 = nn_ind[good_ind_11]
            nn_ind_good_21 = nn_ind[good_ind_21]
            nn_ind_good_22 = nn_ind[good_ind_22]
            nn_ind_good_31 = nn_ind[good_ind_31]

            #----------------------------------------------------
            #using linear least square method
            #----------------------------------------------------
            interp_data_11 = interpolate.interp1d(nn_ind[1:num_lag], M11_0[tt-1, 1:num_lag ,j], kind='linear', fill_value = 'extrapolate' ) #create function
            interp_data_21 = interpolate.interp1d(nn_ind[1:num_lag], M21_0[tt-1, 1:num_lag ,j], kind='linear', fill_value = 'extrapolate' ) #create function
            interp_data_22 = interpolate.interp1d(nn_ind[1:num_lag], M22_0[tt-1, 1:num_lag ,j], kind='linear', fill_value = 'extrapolate' ) #create function
            interp_data_31 = interpolate.interp1d(nn_ind[1:num_lag], M31_0[tt-1, 1:num_lag ,j], kind='linear', fill_value = 'extrapolate' ) #create function

            interp_result_11[tt-1,:,j] = interp_data_11(nn_ind[:num_lag])
            interp_result_21[tt-1,:,j] = interp_data_21(nn_ind[:num_lag])
            interp_result_22[tt-1,:,j] = interp_data_22(nn_ind[:num_lag])
            interp_result_31[tt-1,:,j] = interp_data_31(nn_ind[:num_lag])

            M11_to_0[tt-1,j] = interp_result_11[tt-1,0,j]
            M21_to_0[tt-1,j] = interp_result_21[tt-1,0,j]
            M22_to_0[tt-1,j] = interp_result_22[tt-1,0,j]
            M31_to_0[tt-1,j] = interp_result_31[tt-1,0,j]

           # plt.plot(nn_ind_good[:14], interp_result[tt,:14,j], '-*', nn_ind_good[:14], m_input_good,'-o', nn_ind_good[:14], M_star(nn_ind_good[:14],*param),'--')
            if M11_to_0[tt-1,j] >= M11_0[tt-1,0,j]:
                M11_to_0[tt-1,j] = M11_0[tt-1,0,j]
            if M11_to_0[tt-1,j] < 0:
                M11_to_0[tt-1,j] = 0 # M11_0[tt-1,0,j]


            if M21_to_0 [tt-1,j] >= M21_0[tt-1,0,j]:
                M21_to_0[tt-1,j] = M21_0[tt-1,0,j]
  

            if M22_to_0[tt-1,j]>= M22_0[tt-1,0,j]:
                M22_to_0[tt-1,j]=M22_0[tt-1,0,j]


            if M31_to_0[tt-1,j]>= M31_0[tt-1,0,j]:
                M31_to_0[tt-1,j]=M31_0[tt-1,0,j]
            if M31_to_0[tt-1,j]<= 0:
                M31_to_0[tt-1,j]=0


    #
    dM11 = M11_0[:,0,:] - M11_to_0

    return(M11_to_0, date[t_idx[1:]])
##==============================


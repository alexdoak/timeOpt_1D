import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr
from scipy.signal import hilbert
from scipy.fft import fftfreq, fft, ifft
import pickle
import xarray as xr
import pandas as pd
from filters import apply_high_bandpass, smooth, apply_low_bandpass
from filters import taner_filter



def timeOpt(data_x,x,targetE,targetP,sedmin=0.00005,sedmax=0.0005,numsed=100, 
            flow=0.035,fhigh=0.065,roll=1000,padfac=2,genplot=False):
    """
    Code written by Alex Doak.
    
    This timeOpt code is a replica of the R code in the astrochron package. We
    rebuilt the code in python to compare with new methodologies, which we are 
    building in Python. We have not included many of the checks that the 
    original code provides. Hence, for interested readers, we suggest using this
    code with caution, and instead direct you to the astrochron package for
    which this code has been more significantly more carefully implemented. 
    
    https://cran.r-project.org/web/packages/astrochron/index.html
    
    Credit to Stephen Meyers et al for their work in developing the timeOpt 
    algorithm which has inspired this new research. 
    
    Some notation has been taken from the astrochron package, others we have
    given different names according to our preference. A subtle difference is
    we perform operations on the data transformed to time, rather than simply
    scaling spatially processed data (since sedrate is assumed constant). This
    choice is related to our work where we consider non-constant sedrates.
    
    VARIABLES:
        
    data_x: 1D array - the chemical data as a function of x,
    x: 1D array - the corresponding x values assossiated with data_x
    targetE/P: 1D array - target eccentricity/precession periods
    sedmin/sedmax - minimum and maxiumum sedemintation rate to consider
    numsed - number of sedrates to test inbetween sedmin and sedmax
    flow/fhigh/roll - parameters related to the Taner filter, using the same 
                      notation as astrochron
    padfac - how many zeros added for zero-padding (padfac=2 doubles data size)
    """
    # ------------------------------------------------------------------------ #
    # Functions to build the matrix X
    # ------------------------------------------------------------------------ #
    # The matrix X containing all the sine and cosines of target frequencies
    def genCycles(target_periods, Npnts, time):
        X = np.zeros((Npnts, 2*target_periods.size))
        for i, period in enumerate(target_periods):
            freq = 2*np.pi/ period
            X[:,2*i] = np.cos(freq*time)
            X[:,2*i+1] = np.sin(freq*time)
        return X

    # -------------------------------------------------------------------------- #
    # Cycle through sedimentation rates, buidling the interpolater X(t) at each
    # rate
    # -------------------------------------------------------------------------- #
    sedrate_test = np.linspace(sedmin,sedmax,numsed)
    
    ### r^2 for each model    
    r2_model1 = np.array([])
    r2_model2 = np.array([])
    
    for sedrate in sedrate_test:
        ### Get time as a function of x
        T = x / sedrate
        dt = T[2]-T[1]
        ### Express data_x(x) as data_t(T(x))
        data_t = xr.DataArray(data_x, dims=("time"), coords={"time": T})
     

        # ------------------------------------------------------------------- #
        # Model 2: fitting full data to eccentricity and precession.
        # ------------------------------------------------------------------- #
        ### Step 1: Fit model and compare fitted model to data
        target_periods = np.concatenate((targetE,targetP))  
        X = genCycles(target_periods, np.size(data_t), T)
        model2 = LinearRegression().fit(X, data_t)
        predicted2 = model2.predict(X)
        r_val, _ = pearsonr(data_t, predicted2)
        r2_model2 = np.append(r2_model2, [r_val**2])
        
        # ----------------------------------------------------------------------- #
        # Model 1: fitting amplitude of precession to eccentricity (modulation)
        # ----------------------------------------------------------------------- #
      
        ### Step 1: Pad data with zeros
        data_padded = np.concatenate((data_t,np.zeros(round(data_t.size*(padfac-1)))))
     
        ### Step 2: Filter spectrum
        ft = fft(data_padded)
        [data_filtered,ft_filtered,taper] = taner_filter(ft, dt, flow, fhigh, roll, padfac)
        
        ### Step 3: Recover amplitude using hilbert transform
        data_prec_amp = np.abs(hilbert(data_filtered))
        
        ### Step 4: Fit model and compare fitted model to data
        target_periods = targetE
        X = genCycles(target_periods, np.size(data_t), T)
        model1 = LinearRegression().fit(X, data_prec_amp)
        predicted1 = model1.predict(X)
        r_val, _ = pearsonr(data_prec_amp, predicted1)
        r2_model1 = np.append(r2_model1, [r_val**2])
       
        
        if genplot:
            fig, axs = plt.subplots(2)
            fig.suptitle(f'sedrate = {sedrate}')
            axs[0].plot(T,data_t,'k')
            axs[0].plot(T,predicted2,'r')
            
            axs[1].plot(T,data_prec_amp,'k')
            axs[1].plot(T,predicted1,'r')
            

        
    return r2_model1, r2_model2





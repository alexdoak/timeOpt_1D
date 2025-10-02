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
from filters import taner_filter, plot_spectrum
from timeOpt import timeOpt
from spectrum_methods import global_spectrum_method
import pywt

# --------------------------------------------------------------------------- #
# Load data for TimeOpt. Here, we are loading the data from the Josso et al
# 2021 paper. We are reproducing figure 5, so the sedmin and sedmax are chosen
# as the values they use in the R code. Likewise, we use their targetE and 
# targetP.
# --------------------------------------------------------------------------- #
chemdata = pd.read_csv('Data_Depth5000_Detrended_Interpolated.csv', header=None)
chemdata = chemdata.values
x = chemdata[:,0]
x = x-x[0]
dx = x[1]-x[0]
data_x = chemdata[:,1]
data_x = data_x-np.mean(data_x)
data_x = data_x / np.std(data_x) 
del(chemdata)

### Store data in xarrays: a package with a useful set of operations
data_x = xr.DataArray(data_x, dims=("space"), coords={"space": x})
    

# --------------------------------------------------------------------------- #
# Choose target periods for constructing the matrix X in timeOpt
# --------------------------------------------------------------------------- #
targetE = np.array([405.6795,130.719,123.839,98.86307,94.87666])
targetP = np.array([23.62069,22.31868,19.06768,18.91979])
# Express as frequencies
targetE_freq = 1/targetE
targetP_freq = 1/targetP

# --------------------------------------------------------------------------- #
# Choose parameters in Taner filter to pick out precession amplitude 
# --------------------------------------------------------------------------- #
flow = 0.035
fhigh = 0.065
roll = 1000
padfac = 2

# --------------------------------------------------------------------------- #
# Choose sedementation rates to test
# --------------------------------------------------------------------------- #
sedmin = 0.00005/100
sedmax = 0.0005/100
numsed = 100
sedrate_test = np.linspace(sedmin,sedmax,numsed)


# *************************************************************************** #
# --------------------------------------------------------------------------- #
# TIMEOPT
# --------------------------------------------------------------------------- #
# *************************************************************************** #

[r2_model1,r2_model2] = timeOpt(data_x,x,targetE,targetP,sedmin,sedmax,numsed,
                                flow,fhigh,roll,padfac,genplot=False)


fig, axs = plt.subplots(2)
fig.suptitle('Fig 5E and 5F from Josso et al')
axs[0].plot(sedrate_test,r2_model1,'ro')
axs[0].plot(sedrate_test,r2_model2,'k')
axs[0].set(xlabel='Sedementation rate', ylabel='Pearson squared correlation')


axs[1].plot(sedrate_test,r2_model2*r2_model1,'k')
axs[1].set(xlabel='Sedementation rate', ylabel='product of correlation')

















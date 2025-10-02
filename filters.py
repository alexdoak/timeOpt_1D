import numpy as np
from scipy.fft import fftfreq, fft, ifft
import matplotlib.pyplot as plt

def plot_spectrum(data,dt,padfac=2,absolute=False,normalise=True):
    ### Step 1: Pad data with zeros
    data_padded = np.concatenate((data,np.zeros(round(data.size*(padfac-1)))))
   
    ### Step 2: FFT
    ft = fft(data_padded)
    if normalise:
        ft = ft/np.max(np.abs(ft))
    if absolute:
        ft = np.abs(ft)
    ### Step 3: get periods and plot
    nf = ft.size
    freq = np.fft.fftfreq(nf, d=dt)
    freq = np.fft.fftshift(freq)
    ft = np.fft.fftshift(ft)
    
#    plt.plot(freq,ft,'k')
    plt.plot(freq,np.abs(ft),'k')

    return 

    


def taner_filter(ft, dt, flow=0.035, fhigh=0.065, roll=1000, padfac=2):
    """
    Apply a Taner filter to padded fourier transform data ft, identical to the
    code in the R timeOpt code. 
    
    ft - zero padded fft data
    flow - low frequency filter
    fhigh - high frequency filter
    roll - parameter determing the slope of the cutoo
    dt - spacing in the time sampling
    padfac = 2 - how much zero padding was applied to data. Best to choose an 
                 integer
    """
    
    # Frequencies related to sample size and sample spacing
    nf = len(ft)
    df = 1/(nf*dt)
    nyqfreq = nf // 2
    negfreq = nyqfreq -1
    
    # As in timeOpt, the filter is parameterised in terms of angular frequencies
    fcent = (flow + fhigh) / 2
    wl = 2 * np.pi * flow
    wc = 2 * np.pi * fcent
    wh = 2 * np.pi * fhigh
    bw = wh - wl

    # Shape parameter
    amp2 = 1/2**.5
    arg1 = 1 - (roll * np.log(10)) / (20 * np.log(amp2))
    arg1 = np.log(arg1)
    arg2 = ((bw + 2) / bw) ** 2
    arg2 = np.log(arg2)
    twod = 2 * arg1 / arg2

    # Positive frequencies
    dw = 2 * np.pi / (nf * dt)
    w_pos = np.arange(0, nyqfreq+1) * dw
    arg = (2 * np.abs(w_pos - wc) / bw) ** twod
    filter_pos = 1 / np.sqrt(2) * np.exp(-arg)
    filter_pos /= np.max(filter_pos)

    # Negative frequencies
    w_neg = -w_pos[-2:0:-1]
    aw = np.abs(w_neg)
    arg = (2 * np.abs(aw - wc) / bw) ** twod
    filter_neg = 1 / np.sqrt(2) * np.exp(-arg)
    filter_neg /= np.max(filter_neg)

    # Full taper
    taper = np.concatenate([filter_pos, filter_neg])
    
    # Apply the taper filter and apply ifft
    ft_filtered = ft * taper
    data_filtered = np.real(ifft(ft_filtered))
    
    # Remove padded region of data
    ndata = round(nf/padfac)
    data_filtered = data_filtered[0:ndata]
    
    return data_filtered, ft_filtered, taper



def smooth(x, window_len=11):
    window = np.ones(window_len) / window_len
    smoothened=np.convolve(x, window, mode='same')
    smoothened[1] = x[1]
    return np.convolve(x, window, mode='same')


# --- Butterworth Bandpass Filter ---
def high_bandpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    cut = cutoff / nyq
    b, a = butter(order, cut, 'high')
    return b, a

def apply_high_bandpass(data, cutoff, fs,padfac=2, order=5):
    N = data.size
    data_pad = np.concatenate((data,np.zeros(round(data.size*(padfac-1)))))
    b, a = high_bandpass(cutoff, fs, order=order)
    return filtfilt(b, a, data_pad)[0:N]

# --- Butterworth Bandpass Filter ---
def low_bandpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    cut = cutoff / nyq
    b, a = butter(order, cut, 'low')
    return b, a

def apply_low_bandpass(data, cutoff, fs, padfac=2, order=5):
    data_pad = np.concatenate((data,np.zeros(round(data.size*(padfac-1)))))
    b, a = low_bandpass(cutoff, fs, order=order)
    return filtfilt(b, a, data_pad)





def butter_bandpass(lowcut, highcut, fs, order=5):
        from scipy.signal import butter, sosfilt, sosfreqz 
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

def appply_butter_bandpass(data, lowcut, highcut, fs, order=5):
        from scipy.signal import butter, sosfilt, sosfreqz 
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfilt(sos, data)
        return y



from scipy.signal import butter, filtfilt
from scipy.fft import fftfreq, fft, ifft
import numpy as np


def zero_pad_sym(x, padfac=2):
    """
    Padding a vector x with numel(x)*(padfac-1) zeros, symmtrically about the 
    vector x. Only needed if performing convolutions. For taner filter, we can 
    add zeros to the end
    """
    n = len(x)
    if padfac <= 1:
        return np.asarray(x)

    total_len = n * padfac
    pad_len = total_len - n

    # Split zeros as evenly as possible
    left = pad_len // 2
    right = pad_len - left

    return np.pad(x, (left, right), mode="constant", constant_values=0)



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



















# chemdata = pd.read_csv('Data_Depth5000_Detrended_Interpolated.csv', header=None)
# chemdata = chemdata.values
# x_raw = chemdata[:,0]
# x_raw = x_raw-x_raw[0]
# data = chemdata[:,1]
# data = data-np.mean(data)

# sedrate = 0.00027/100
# t = x_raw / sedrate
# dt=t[1]-t[0]
# flow=0.035
# fhigh=0.065
# roll = 1000



# data_padded = np.concatenate((data,np.zeros(round(data.size*(padfac-1)))))
# ft = fft(data_padded)
# data_padded2 = np.real(ifft(ft))

# [data_filtered, ft_filtered, taper] = taner_filter(ft, flow, fhigh, roll, dt)

# df = 1/(nf*dt)
# freq = np.fft.fftfreq(nf, d=dt)


# plt.plot(ft_filtered,'x')
# plt.plot(ft)
# plt.plot(taper*2000,'o')
# plt.xlim((0,500))

# #ft2 = xrft.fft(data)
# #data2 = np.real(xrft.ifft(ft))



# freq_test = 0.06
# data_test = 0.9*np.cos(2*np.pi*freq_test*t)

# data=data_test

# data_padded = np.concatenate((data,np.zeros(round(data.size*(padfac-1)))))
# ft = fft(data_padded)
# data_padded2 = np.real(ifft(ft))

# [data_filtered, ft_filtered, taper] = taner_filter(ft, flow, fhigh, roll, dt)
# plt.plot(data_filtered)
# plt.show()
# plt.plot(taper)
# plt.show()


# flow = 0.05
# fhigh = 0.07



# [data_filtered, ft_filtered, taper] = taner_filter(ft, flow, fhigh, roll, dt)
# plt.plot(data_filtered)
# plt.show()
# #plt.plot(taper)
# #plt.show()
# plt.plot(freq,ft_filtered,'x')
# plt.plot(freq,ft)
# plt.plot(freq,taper*2000,'o')
# plt.xlim((0,.1))








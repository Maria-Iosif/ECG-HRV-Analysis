import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from ecg_hrv_analysis.time_domain import *

def frequencydomain(rr_interpolated):
    results, fxx, pxx = frequency_domain(rr_interpolated)
    print("Frequency domain metrics:")
    for k, v in results.items():
        print("- %s: %.2f" % (k, v))

    return


def frequency_domain(rr_interpolated, rrs,fs=4, sampling_rate = 360):
    period = 1 / sampling_rate
    Rrss = [rrs[i] * period * 1000 for i in range(len(rrs))]
    rr_manual = nn_intervals(Rrss)

    x = np.cumsum(rr_manual) / 1000
    f = interp1d(x, rr_manual, kind='cubic')

    # sample rate for interpolation
    steps = 1 / sampling_rate

    # now we can sample from interpolation function
    xx = np.arange(1, np.max(x), steps)
    rr_interpolated = f(xx)

    # Estimate the spectral density using Welch's method
    fxx, pxx = signal.welch(x=rr_interpolated, fs=fs)

    '''
    Segement found frequencies in the bands 
     - Very Low Frequency (VLF): 0-0.04Hz 
     - Low Frequency (LF): 0.04-0.15Hz 
     - High Frequency (HF): 0.15-0.4Hz
    '''
    cond_vlf = (fxx >= 0) & (fxx < 0.04)
    cond_lf = (fxx >= 0.04) & (fxx < 0.15)
    cond_hf = (fxx >= 0.15) & (fxx < 0.4)

    # calculate power in each band by integrating the spectral density
    vlf = trapz(pxx[cond_vlf], fxx[cond_vlf])
    lf = trapz(pxx[cond_lf], fxx[cond_lf])
    hf = trapz(pxx[cond_hf], fxx[cond_hf])
    # sum these up to get total power
    total_power = vlf + lf + hf

    # find which frequency has the most power in each band
    peak_vlf = fxx[cond_vlf][np.argmax(pxx[cond_vlf])]
    peak_lf = fxx[cond_lf][np.argmax(pxx[cond_lf])]
    peak_hf = fxx[cond_hf][np.argmax(pxx[cond_hf])]

    # fraction of lf and hf
    lf_nu = 100 * lf / (lf + hf)
    hf_nu = 100 * hf / (lf + hf)

    results = {}
    results['Power VLF (ms2)'] = vlf
    results['Power LF (ms2)'] = lf
    results['Power HF (ms2)'] = hf
    results['Power Total (ms2)'] = total_power

    results['LF/HF'] = (lf / hf)
    results['Peak VLF (Hz)'] = peak_vlf
    results['Peak LF (Hz)'] = peak_lf
    results['Peak HF (Hz)'] = peak_hf

    results['Fraction LF (nu)'] = lf_nu
    results['Fraction HF (nu)'] = hf_nu
    return rr_interpolated, results, fxx, pxx


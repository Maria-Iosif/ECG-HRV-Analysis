import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from ecg_hrv_analysis.signal_processing import *
from ecg_hrv_analysis.time_domain import avg_fpfn


def ecg_signals(i_pat, patient,chan, f_low, f_high ):
    ecg = patient[i_pat].get('p_signal')[:, chan]
    fs = patient[i_pat].get('fs')
    filtered_ecg = bandpass_filt(ecg, fs, f_low, f_high)
    dif_ecg = differentiation(filtered_ecg)
    square_ecg = squaring(dif_ecg)
    mwi_ecg = mwi(square_ecg, fs)

    fig, axs = plt.subplots(5, 1, figsize=(20, 10), sharex=True)
    x = np.linspace(0, 1800, num=650000, endpoint=True)
    axs[0].plot(x[:5400], ecg[:5400])
    axs[1].plot(x[:5400], filtered_ecg[:5400])
    axs[2].plot(x[:5400], dif_ecg[:5400])
    axs[3].plot(x[:5400], square_ecg[:5400])
    axs[4].plot(x[:5400], mwi_ecg[:5400])
    axs[4].set_xlabel('time [s]')
    axs[0].set(ylabel='Unfiltered ECG')
    axs[1].set(ylabel='Bandpass ECG')
    axs[2].set(ylabel='Dif ECG')
    axs[3].set(ylabel='Squared ECG ')
    axs[4].set(ylabel='MWI ECG')

def r_total(pat, an, paper, method):
    rpeaks_dict, fpfn_R_an_dict, fpfn_R_R_dict = rpeaks_extr(pat, an, 5, 12)
    x = [pat[i].get('record_name') for i in range(len(pat))]
    paper = [2273, 1865, 2187, 2230, 2572, 2027, 2137, 1763, 2532, 2124, 2539, 1795, 1879, 1953, 2412, 1535, 2275, 1987,
             1863, 2476, 1518, 1619, 2601, 1963, 2136, 2982, 2656, 1862, 2956, 3004,
             2647, 2748, 3251, 2262, 3363, 2208, 2154, 2048, 2427, 2484, 2605, 2053, 2256, 1886, 1780, 3079, 2753, 2084]
    method = [len(rpeaks_dict[i]) for i in range(len(rpeaks_dict))]

    fig = plt.figure(figsize=(10, 7))
    plt.scatter(x, paper, color="purple", label="total number from paper")
    plt.scatter(x, method, color="purple", alpha=0.5, label="total number from method")
    plt.grid("on")
    plt.xticks(rotation=45)
    plt.legend()
    plt.title("Total number of detected R-peaks")
    plt.ylabel("number of R-peaks")
    plt.xlabel("Record number of subject")


def tp_fp_fn(pat, an, low, high):
    rpeaks_dict, fpfn_R_an_dict, fpfn_R_R_dict = rpeaks_extr(pat, an, low, high)
    avg_dict = avg_fpfn(fpfn_R_an_dict)

    data = [avg_dict["tp"], avg_dict["fp"], avg_dict["fn"]]

    fig = plt.figure(figsize=(10, 7))
    x_ax = ["TP", "FP", "FN"]
    ax = sns.boxplot(data, palette="husl")
    ticks = [0, 1, 2]
    ax.set_xticklabels(x_ax)

def count_beats(pat, an, num=True):
    symbols = {'!', '"', '+', '/', 'A', 'E', 'F', 'J', 'L', 'N', 'Q', 'R', 'S', 'V', '[', ']', 'a', 'e', 'f', 'j', 'x',
               '|', '~'}
    beats_dict = {s_: list(np.zeros(len(pat))) for s_ in symbols}
    list_set = list(symbols)

    for i in range(len(list_set)):

        ss = list_set[i]

        for p in range(len(an)):

            AA = an[p].get('symbol')
            total_beats = len(AA)
            count_ = AA.count(ss)

            if num == True:
                beats_dict[ss][p] = count_
            else:
                beats_dict[ss][p] = count_ / total_beats

    return beats_dict

def healthy_not(an, percentage):
    # Create a list of indices corresponding to healthy individuals
    healthy_idx = [i for i in range(len(an)) if
                   (an[i].get('symbol').count('N')) >= percentage * len(an[i].get('symbol'))]

    # Create a list of indices corresponding to unhealthy individuals
    unhealthy_idx = [i for i in range(len(an)) if i not in healthy_idx]

    return healthy_idx, unhealthy_idx
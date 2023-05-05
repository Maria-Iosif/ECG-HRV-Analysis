import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from ecg_hrv_analysis.signal_processing import *
from ecg_hrv_analysis.time_domain import *
from ecg_hrv_analysis.frequency_domain import *
from scipy.interpolate import interp1d
from scipy.stats import norm


def ecg_signals(i_pat, patient, chan, f_low, f_high):
    ecg = patient[i_pat].get("p_signal")[:, chan]
    fs = patient[i_pat].get("fs")
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
    axs[4].set_xlabel("time [s]")
    axs[0].set(ylabel="Unfiltered ECG")
    axs[1].set(ylabel="Bandpass ECG")
    axs[2].set(ylabel="Dif ECG")
    axs[3].set(ylabel="Squared ECG ")
    axs[4].set(ylabel="MWI ECG")


def rpeaks_plot(r_peaks, mwi_ecg, annotations, ind):
    x = np.linspace(0, 1800, num=650000, endpoint=True)

    check = [p for p in r_peaks if p < 5400]
    k = len(check)
    xx = [rp / 360 for rp in r_peaks[:k]]
    yy = mwi_ecg[r_peaks[:k]]

    xx2 = [rp / 360 for rp in annotations[ind].get("sample")[1 : k + 1]]
    yy2 = mwi_ecg[annotations[ind].get("sample")[1 : k + 1]]

    fig, axs = plt.subplots(1, 1, figsize=(12, 7), sharex=True)
    axs.plot(x[:5400], mwi_ecg[:5400])
    axs.scatter(xx, yy, c="r", marker="x", label="detected R peaks")
    axs.scatter(xx2, yy2, c="g", marker="x", label="Annotations")
    plt.legend()
    plt.xlabel("time [s]")
    plt.ylabel("MWI ECG")


def histogram(hr1, hr2):
    fig, axs = plt.subplots(
        3, 1, sharex=True, figsize=(12, 7), gridspec_kw={"height_ratios": [7, 1, 1]}
    )

    min_val = min(min(hr1), min(hr2)) - 5
    max_val = max(max(hr1), max(hr2[1:])) + 5
    # bins = np.linspace(min_val, max_val, num=150)

    mu1, std1 = norm.fit(hr1)
    mu2, std2 = norm.fit(hr2)

    axs[0].hist(hr1, bins=100, histtype="step", color="r")
    axs[0].hist(
        hr1,
        bins=100,
        alpha=0.3,
        label="Heart rate from detected R, \n mean HR %.2f" % mu1,
        color="r",
    )
    axs[0].hist(hr2[1:], bins=100, histtype="step", color="blue")
    axs[0].hist(
        hr2[1:],
        bins=100,
        alpha=0.3,
        label="Heart rate from annotated R, \n mean HR %.2f" % mu2,
        color="blue",
    )
    sns.boxplot(hr1, ax=axs[1], color="r", boxprops=dict(alpha=0.3), orient="h")
    sns.boxplot(hr2[1:], ax=axs[2], color="blue", boxprops=dict(alpha=0.3), orient="h")

    axs[1].axis("off")
    axs[2].axis("off")

    axs[0].tick_params(reset=True)
    plt.xlim(min_val, max_val)
    axs[0].legend()


def r_total(rpeaks_dict, pat, an):
    # rpeaks_dict, fpfn_R_an_dict, fpfn_R_R_dict = rpeaks_extr(pat, an, 5, 12)

    x = [pat[i].get("record_name") for i in range(len(pat))]
    paper = [
        2273,
        1865,
        2187,
        2230,
        2572,
        2027,
        2137,
        1763,
        2532,
        2124,
        2539,
        1795,
        1879,
        1953,
        2412,
        1535,
        2275,
        1987,
        1863,
        2476,
        1518,
        1619,
        2601,
        1963,
        2136,
        2982,
        2656,
        1862,
        2956,
        3004,
        2647,
        2748,
        3251,
        2262,
        3363,
        2208,
        2154,
        2048,
        2427,
        2484,
        2605,
        2053,
        2256,
        1886,
        1780,
        3079,
        2753,
        2084,
    ]
    method = [len(rpeaks_dict[i][0]) for i in range(len(rpeaks_dict))]

    fig = plt.figure(figsize=(10, 7))
    plt.scatter(x, paper, color="purple", label="total number from paper")
    plt.scatter(x, method, color="purple", alpha=0.5, label="total number from method")
    plt.grid("on")
    plt.xticks(rotation=45)
    plt.legend()
    plt.title("Total number of detected R-peaks")
    plt.ylabel("number of R-peaks")
    plt.xlabel("Record number of subject")


def tpfpfn_plot(avg_dict):
    data = [avg_dict["tp"], avg_dict["fp"], avg_dict["fn"]]

    fig = plt.figure(figsize=(10, 7))
    x_ax = ["TP", "FP", "FN"]
    ax = sns.boxplot(data, palette="husl")
    ticks = [0, 1, 2]

    ax.set_xticklabels(x_ax)
    ax.set_title(
        "The TP, FP, FN among R-peaks from code implementation and annotations"
    )


def heart_rate_plot(pat, annotations, mean_hr):
    means = mean_hr_an(annotations)

    fig, axs = plt.subplots(1, 1, figsize=(15, 10))
    axs.plot([i for i in range(len(pat))], mean_hr, color="b")
    axs.plot([i for i in range(len(pat))], means, color="r")
    axs.scatter(
        [i for i in range(len(pat))], mean_hr, color="b", label="Mean heart rate"
    )
    axs.scatter(
        [i for i in range(len(pat))], means, color="r", label="Mean heart rate from an."
    )
    plt.legend()
    plt.xlabel("Individuals")
    plt.ylabel("Mean heart rate")


def bar_plot(count_beats, pat):
    data = count_beats
    df = pd.DataFrame(data, index=[i for i in range(len(pat))])

    # create stacked bar plot
    ax = df.plot(kind="bar", stacked=True, figsize=(16, 9), cmap="PiYG")

    # add labels and legend
    ax.set_xlabel("Individuals")
    ax.set_ylabel("Number of Labels")
    ax.legend(title="Labels")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    plt.show()


def heatmap_plot(count_beats):
    df = pd.DataFrame.from_dict(count_beats)
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))
    ax = sns.heatmap(df, linewidth=0.0, cmap="PiYG")


def sdnn_plot(h_sdnn, n_sdnn):
    fig, axs = plt.subplots(1, 1, figsize=(10, 7))
    axs.boxplot([h_sdnn, n_sdnn], labels=["Healthy", "Non-healthy"])

    x = [1 for i in range(len(h_sdnn))]
    plt.scatter(x, h_sdnn, s=15, color="r", alpha=0.3)

    x2 = [2 for i in range(len(n_sdnn))]
    plt.scatter(x2, n_sdnn, s=15, color="green", alpha=0.3)

    plt.ylabel("SDNN")


def rmssd_plot(h_rmssd, n_rmssd):
    fig, axs = plt.subplots(1, 1, figsize=(10, 7))
    axs.boxplot([h_rmssd, n_rmssd], labels=["Healthy", "Non-healthy"])

    x = [1 for i in range(len(h_rmssd))]
    plt.scatter(x, h_rmssd, s=15, color="r", alpha=0.3)

    x2 = [2 for i in range(len(n_rmssd))]
    plt.scatter(x2, n_rmssd, s=15, color="green", alpha=0.3)

    plt.ylabel("RMSSD")


def frequency_plot(x, rr_manual, xx, rr_interpolated):
    plt.figure(figsize=(20, 15))
    plt.subplot(211)
    plt.title("RR intervals")
    plt.plot(
        x,
        rr_manual,
        color="k",
        markerfacecolor="#A651D8",
        markeredgewidth=0,
        marker="o",
        markersize=8,
    )
    plt.xlabel("Time (s)")
    plt.ylabel("RR-interval (ms)")
    plt.title("Interpolated")
    plt.gca().set_xlim(0, 20)

    plt.subplot(212)
    plt.title("RR-Intervals (cubic interpolation)")
    plt.plot(
        xx,
        rr_interpolated,
        color="k",
        markerfacecolor="#51A6D8",
        markeredgewidth=0,
        marker="o",
        markersize=8,
    )
    # plt.gca().set_xlim(0, 20)
    plt.xlabel("Time (s)")
    plt.ylabel("RR-interval (ms)")
    plt.show()


def fft(fxx, pxx):
    plt.figure(figsize=(20, 7))
    plt.plot(fxx, pxx, color="k", linewidth=0.3)
    plt.title("FFT Spectrum (Welch's periodogram)")

    # create interpolation function for plotting frequency bands
    psd_f = interp1d(fxx, pxx)

    # setup frequency bands for plotting
    x_vlf = np.linspace(0, 0.04, 100)
    x_lf = np.linspace(0.04, 0.15, 100)
    x_hf = np.linspace(0.15, 0.4, 100)

    plt.gca().fill_between(x_vlf, psd_f(x_vlf), alpha=0.2, color="#A651D8", label="VLF")
    plt.gca().fill_between(x_lf, psd_f(x_lf), alpha=0.2, color="#51A6D8", label="LF")
    plt.gca().fill_between(x_hf, psd_f(x_hf), alpha=0.2, color="#D8A651", label="HF")

    plt.gca().set_xlim(0, 0.5)
    plt.gca().set_ylim(0)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Density")
    plt.legend()
    plt.show()


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
    symbols = {
        "!",
        '"',
        "+",
        "/",
        "A",
        "E",
        "F",
        "J",
        "L",
        "N",
        "Q",
        "R",
        "S",
        "V",
        "[",
        "]",
        "a",
        "e",
        "f",
        "j",
        "x",
        "|",
        "~",
    }
    beats_dict = {s_: list(np.zeros(len(pat))) for s_ in symbols}
    list_set = list(symbols)

    for i in range(len(list_set)):
        ss = list_set[i]

        for p in range(len(an)):
            AA = an[p].get("symbol")
            total_beats = len(AA)
            count_ = AA.count(ss)

            if num == True:
                beats_dict[ss][p] = count_
            else:
                beats_dict[ss][p] = count_ / total_beats

    return beats_dict

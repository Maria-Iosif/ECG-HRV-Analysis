import glob
import pathlib
import sklearn.metrics
import numpy as np
import wfdb
from scipy import signal


def data_dict(data_path):
    patients = [
        pathlib.Path(p).name.replace(".hea", "")
        for p in glob.glob(str(data_path / "*.hea"))
    ]

    patients_dict = {num_p: {} for num_p in range(len(patients))}
    annotation_dict = {num_p: {} for num_p in range(len(patients))}

    for i in range(len(patients)):
        patient = patients[i]
        annotation = wfdb.rdann(str(data_path / patient), extension="atr")
        record = wfdb.rdrecord(data_path / patient)

        rec_keys = list(record.__dict__.keys())
        an_keys = list(annotation.__dict__.keys())

        for key in rec_keys:
            patients_dict[i][key] = record.__dict__.get(key)

        for key in an_keys:
            annotation_dict[i][key] = annotation.__dict__.get(key)

    return patients_dict, annotation_dict


def bandpass_filt(ecg, fs, f_low, f_high):
    """
    Apply bandpass filter to the ECG signal

    Input parameters:
      ecg: The ECG signal
      f_low:  The low cutoff frequency
      f_high: The high cutoff frequency

    Output:
      The filtered ECG signal

    """
    #   b, a = butter(2, [5, 15], btype='bandpass', fs=fs)
    # filtered_signal = filtfilt(b, a, ecg_signal)

    # Define filter parameters
    order = 1  # Filter order

    # Create bandpass filter
    b, a = signal.butter(order, [f_low, f_high], btype="bandpass", fs=fs)

    filtered_ecg = signal.lfilter(b, a, ecg)

    return filtered_ecg


def differentiation(ecg):
    """
    Find the derivative of the ECG signal
    d(fy)
    """
    diff = np.diff(ecg)

    return diff


def squaring(ecg):
    """
    Find the square of the ECG signal
    (dy)^2
    """

    squared = ecg**2

    return squared


def mwi(ecg, fs):
    """
    Apply moving-window integration to the given ECG signal
    """
    max_time = 0.083
    window_size = int(max_time * fs)
    ones = np.ones(window_size)

    # integrated_signal = np.convolve(squared_signal, ones, mode='valid')
    ecg_mwi = np.convolve(ecg, ones, mode="valid")

    return ecg_mwi


def r_peak_detection(ecg_detection, fs):
    min_distance = int(0.15 * fs)

    signal_peaks = [0]
    noise_peaks = []

    SPKI = 0.0
    NPKI = 0.0

    threshold_I1 = 0.0
    threshold_I2 = 0.0

    RR_missed = 0
    index = 0
    indexes = []

    missed_peaks = []
    peaks = []

    for i in range(1, len(ecg_detection) - 1):
        if (
            ecg_detection[i - 1] < ecg_detection[i]
            and ecg_detection[i + 1] < ecg_detection[i]
        ):
            peak = i
            peaks.append(i)

            if (
                ecg_detection[peak] > threshold_I1
                and (peak - signal_peaks[-1]) > 0.3 * fs
            ):
                signal_peaks.append(peak)
                indexes.append(index)
                SPKI = 0.125 * ecg_detection[signal_peaks[-1]] + 0.875 * SPKI
                if RR_missed != 0:
                    if signal_peaks[-1] - signal_peaks[-2] > RR_missed:
                        missed_section_peaks = peaks[indexes[-2] + 1 : indexes[-1]]
                        missed_section_peaks2 = []
                        for missed_peak in missed_section_peaks:
                            if (
                                missed_peak - signal_peaks[-2] > min_distance
                                and signal_peaks[-1] - missed_peak > min_distance
                                and ecg_detection[missed_peak] > threshold_I2
                            ):
                                missed_section_peaks2.append(missed_peak)

                        if len(missed_section_peaks2) > 0:
                            signal_missed = [
                                ecg_detection[i] for i in missed_section_peaks2
                            ]
                            index_max = np.argmax(signal_missed)
                            missed_peak = missed_section_peaks2[index_max]
                            missed_peaks.append(missed_peak)
                            signal_peaks.append(signal_peaks[-1])
                            signal_peaks[-2] = missed_peak

            else:
                noise_peaks.append(peak)
                NPKI = 0.125 * ecg_detection[noise_peaks[-1]] + 0.875 * NPKI

            threshold_I1 = NPKI + 0.25 * (SPKI - NPKI)
            threshold_I2 = 0.5 * threshold_I1

            if len(signal_peaks) > 8:
                RR = np.diff(signal_peaks[-9:])
                RR_ave = int(np.mean(RR))
                RR_missed = int(1.66 * RR_ave)

            index = index + 1

    signal_peaks.pop(0)

    return signal_peaks


def Pan_Tompkins(ecg, fs, f_low, f_high):
    """ """

    filtered_ecg = bandpass_filt(ecg, fs, f_low, f_high)
    dif_ecg = differentiation(filtered_ecg)
    square_ecg = squaring(dif_ecg)
    mwi_ecg = mwi(square_ecg, fs)

    r_peaks = r_peak_detection(mwi_ecg, fs)

    return mwi_ecg, r_peaks


def rpeaks_extr(pat, an, low, high):
    def avg_fpfn(fpfn, pat):
        avg_dict = {"tp": [], "fp": [], "fn": []}

        for i in range(len(pat)):
            total = np.sum(list(fpfn[i].values()))
            for j in ["tp", "fp", "fn"]:
                avg_val = (fpfn[i][j][0]) / total
                avg_dict[j].append(avg_val)

        return avg_dict

    rpeaks_dict = {i: [] for i in pat}
    rpeaks_signal_dict = {i: [] for i in pat}
    fpfn_R_an_dict = {i: {"tp": [], "fp": [], "fn": []} for i in pat}
    fpfn_R_R_dict = {i: {"tp": [], "fp": [], "fn": []} for i in pat}

    for i in range(len(pat)):
        ecg_1 = pat[i].get("p_signal")[:, 0]
        fs = pat[i].get("fs")
        ecg_filt_1, rpeaks_1 = Pan_Tompkins(ecg_1, fs, low, high)

        rpeaks_dict[i].append(rpeaks_1)

        peaks_1 = signal.find_peaks(ecg_filt_1, distance=150)[0]

        # Calculate false positives and false negatives
        tp = 0
        fp = 0
        fn = 0

        length_an = len(an[i].get("sample"))

        for j in range(length_an):
            if an[i].get("sample")[j] in rpeaks_1:
                tp += 1
            elif np.any(
                np.logical_and(
                    rpeaks_1 > an[i].get("sample")[j] - int(0.083 * fs),
                    rpeaks_1 < an[i].get("sample")[j] + int(0.083 * fs),
                )
            ):
                fp += 1
            else:
                fn += 1

        fpfn_R_an_dict[i]["tp"].append(tp)
        fpfn_R_an_dict[i]["fp"].append(fp)
        fpfn_R_an_dict[i]["fn"].append(fn)

        # Calculate false positives and false negatives
        tp = 0
        fp = 0
        fn = 0

        length_an = len(peaks_1)

        for j in range(length_an):
            if peaks_1[j] in rpeaks_1:
                tp += 1
            elif np.any(
                np.logical_and(
                    rpeaks_1 > peaks_1[j] - int(0.083 * fs),
                    rpeaks_1 < peaks_1[j] + int(0.083 * fs),
                )
            ):
                fp += 1
            else:
                fn += 1

        fpfn_R_R_dict[i]["tp"].append(tp)
        fpfn_R_R_dict[i]["fp"].append(fp)
        fpfn_R_R_dict[i]["fn"].append(fn)

    return rpeaks_dict, fpfn_R_an_dict, fpfn_R_R_dict

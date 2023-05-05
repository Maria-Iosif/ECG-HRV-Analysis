import numpy as np


def timedomain(peaks, sampling_rate = 360):

    rr = nn_intervals(peaks)
    hr = calculate_heart_rate(peaks,sampling_rate = 360)

    results = {}

    results["Mean NN (ms)"] = np.mean(rr)
    results["STD NN/SDNN (ms)"] = np.std(rr)
    results["Mean HR (beats/min)"] = np.mean(hr)
    results["STD HR (beats/min)"] = np.std(hr)
    results["Min HR (beats/min)"] = np.min(hr)
    results["Max HR (beats/min)"] = np.max(hr)
    results["RMSSD (ms)"] = rmssd_measure(peaks)
    results["SDNN"] = sdnn_measure(peaks)
    # results["SDANN"] = sdann_measure(peaks)
    # results['NNxx'] = np.sum(np.abs(np.diff(rr)) > 50)*1
    # results['pNNxx (%)'] = 100 * np.sum((np.abs(np.diff(rr)) > 50)*1) / len(rr)

    print("Time domain metrics - for given NN-intervals:")
    for k, v in results.items():
        print("- %s: %.2f" % (k, v))

    return

def calculate_heart_rate(peaks, sampling_rate):
  """
  Calculate heart rate from R-R intervals in seconds.

  Args:
  rr_intervals (list of floats): list of R-R intervals in seconds
  sampling_rate (float): sampling rate used to acquire ECG signal (in Hz)

  Returns:
  float: the average heart rate in beats per minute
  """

  period = 1/sampling_rate
  rs = [peaks[i]*period*1000 for i in range(len(peaks))]
  nn_intervals = np.diff(rs)
  hr = 60000/nn_intervals

  return hr

def mean_heart_rate(r_p_dict ,sampling_rate=360):
    """
    Input: A dictionary with r-peaks
    Sampling rate (Hz) (default 360Hz)

    Output: the mean Heart rate for each individual
    """
    :param r_p_dict:
    :param sampling_rate:
    :return:

    means = np.zeros(len(r_p_dict))
    for i in range(len(r_p_dict)):

    hr = calculate_heart_rate(r_p_dict[i][0], sampling_rate)
    hr_avg = np.mean(hr)

    means[i] = hr_avg

    return means


def nn_intervals(peaks):
    peak_s = [peaks[i] for i in range(len(peaks))]

    p_size = len(peak_s)
    nn_inter = np.zeros(p_size - 1)

    for i in range(p_size - 1):
        nn_inter[i] = peak_s[i + 1] - peak_s[i]

    return nn_inter


def sdnn_measure(peaks):
    nn_inter = nn_intervals(peaks)

    sdnn = np.std(nn_inter, ddof=1)

    return sdnn


def peaks_per_inter(peaks, time_range=300):
    """
    time_range : Time range to split in seconds (Default: 300[s] = 5[mins])

    """

    nn_inter = nn_intervals(peaks).tolist()

    tot_duration = np.sum(nn_inter)

    number_of_splits = int(np.ceil(tot_duration / time_range))

    split_dict = {str(i) + "_inter": [] for i in range(number_of_splits)}

    peak_int = [int(peak_ / time_range) for peak_ in peaks]

    for i in range(number_of_splits):
        k = peak_int[i]

        split_dict[str(k) + "_inter"].append(nn_inter[i])

    return split_dict


def sdann_measure(peaks, time_range=300):
    nn_inter = nn_intervals(peaks)

    split_dict = peaks_per_inter(peaks)

    keys_ = list(split_dict.keys())

    avg_split_dict = {key: [] for key in keys_}

    for key in keys_:
        avg_split_dict[key] = max(0, np.mean(peaks_per_inter(peaks)[key]))

    avg_list = list(avg_split_dict.values())
    # flat_avg_list =  [num for sublist in avg_list for num in sublist]

    sdann = np.std(avg_list, ddof=1)

    return sdann


def rmssd_measure(peaks):
    nn_inter = nn_intervals(peaks)
    nn_dif = [nn_inter[i + 1] - nn_inter[i] for i in range(len(nn_inter) - 1)]

    rmssd = np.sum([x_**2 for x_ in nn_dif])
    rmssd = np.sqrt(rmssd / len(nn_dif))

    return rmssd


def mean_for_an():
    means = []
    sampling_rate = 360
    for i in range(len(an)):
        hr = calculate_heart_rate(an[i].get('sample'), sampling_rate)
        mean_hr = np.mean(hr)

        means.append(mean_hr)
    return means

def avg_fpfn(fpfn):
  avg_dict = {'tp':[], 'fp':[], 'fn':[]}

  for i in range(len(pat)):

    total = np.sum(list(fpfn[i].values()))
    for j in ['tp', 'fp', 'fn']:
      avg_val = (fpfn[i][j][0])/total
      avg_dict[j].append(avg_val)

  return avg_dict
import pathlib

from ecg_hrv_analysis.frequency_domain import *
from ecg_hrv_analysis.plot import *
from ecg_hrv_analysis.signal_processing import *
from ecg_hrv_analysis.time_domain import *

data_path = pathlib.Path("./ecg_data/physionet.org/files/mitdb/1.0.0")

patiens, annotations = data_dict(data_path)

rpeaks_dict, fpfn_R_an_dict, fpfn_R_R_dict = rpeaks_extr(pat, an, low, high)



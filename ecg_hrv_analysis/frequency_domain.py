from scipy.interpolate import interp1d


def frequencydomain():
    results, fxx, pxx = frequency_domain(rr_interpolated)
    print("Frequency domain metrics:")
    for k, v in results.items():
        print("- %s: %.2f" % (k, v))

    return

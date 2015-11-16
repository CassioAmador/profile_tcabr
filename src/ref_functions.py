"""
Function: Common functions and constants associated with Reflectometry.
TODO: ?
"""

from numpy import sqrt

# k = (electron mass)/(electron charge)
k = 8.9786629998410721
# (1/k)^2
k_2 = 0.012404425565580237


def freq2den(freq):
    """converts density (m^-3) in frequency (Hz)"""

    return freq * freq * k_2


def den2freq(den):
    """converts frequency (Hz) in density (m^-3)"""

    return k * sqrt(den)


def refrac(w_p, w):
    """returns refraction index in O mode, given the plasma
    oscillation frequency (w_p), and the probing frequency (w)"""

    return sqrt(1 - pow(w_p / w, 2))


def smooth_signal(self, signal, window=7, order=3):
    """ Smooth signal (low-pass filter). Work for scipy version >= 0.14.0."""

    from scipy.signal import savgol_filter
    try:
        signal_smoothed = savgol_filter(signal, window_length=window, polyorder=order)
    except:
        return signal
    return signal_smoothed

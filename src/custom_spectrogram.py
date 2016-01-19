"""
Custom spectrogram function. Returned values are different than scipy.
(c) Cassio Amador 2015
TODO: DECIDE IF CUTTING FREQUENCY SHOULD RETURN A FULL MATRIX WITH NaNS, 
OR JUST A PARTIAL ONE.
"""
import numpy as np
from scipy import signal, fftpack


def spectrogram(sig, window_size=256, step_scale=4, zer_pad=2, time_array=None,fft_type='fft', log=False, normalize=0, dc_cut=0, fft_shift=0,filtered=0, freq_mask=[]):
    """Evaluate and plot spectrogram (SFFT) of beating signal.

    Input:

    window_size: size of each segment, in number of points; defaults to 256
    step_scale: step for next segment, proportional to window_size (step=window_size/step_scale); defaults to 4
    zer_pad: size of each segment, including zero padding, defaults to 2
    time_array: array with time for each
    fft_type: 'fft' for standard fft (from scipy.fftpack); 'welch' for scipy.signal.welch transform
    log: True for log spectrogram; defaults to False
    normalize: set True to normalize each spectrum; defaults to False
    dc_cut: set True to cut the zero frequency; defaults to False
    filtered: set True to apply a filfilt filter to the spectrum. Not recommended unless you are sure about the output. Defaults to False
    freqs_window: [f_min,f_max] array with maximum and minimum beating frequency, in terms of index. The 'forbidden'frequencies are set to None. Defaults to (0,window_size)

    Output:

    matrix: spectrogram as a matrix with n spectrums and m frequencies. [n x m]
    time_spec: if time_array is given, returns an array with time of the sliding window center for each spectrum
    beat_freq: if time_array is given, returns the beating frequency in (1/time_array) unity.
    """

    # alias for sig length
    sig=np.concatenate((np.zeros(window_size/4),sig,np.zeros(window_size/4)))
    N = len(sig)

    # SFFT step size,
    step = int(window_size / step_scale)

    if time_array is not None:
        # beating frequency
        if len(time_array) == len(sig):
            beat_freq=eval_beat_freq(time_array,window_size=window_size,zer_pad=zer_pad,fft_shift=fft_shift)
        else:
            raise ValueError('length of time array and signal are different to evaluate spectrogram')

        # time array for spectrogram
        time_spec = np.linspace(time_array[window_size], time_array[-window_size], num=(N - window_size) / step)

    # creates the window function that will 'slide' trough the signal
    # to evaluate each FFT. Kaiser seems to be the cleanest one
    #window_func = np.hanning(window_size)
    #window_func = signal.tukey(window_size, 0.25)
    window_func = signal.kaiser(window_size, 10)

    # if not shifting, treats as if real signal
    factor=2
    if fft_shift:
        factor=1

    # create a matrix to receive the spectra
    mat_Y=window_size*zer_pad/factor
    #if len(freq_mask)!=0:
    #    mat_Y=len(np.where(freq_mask)[0])
    matrix = np.empty(((N - window_size) / step, mat_Y))

    if filtered == 1:
        b, a, zi = _init_filter()

    # slide window trough signal, and evaluates the FFT.
    for i in range(int((N - window_size) / step)):
        t = i * step
        new_sig = sig[t:t + window_size]
        #print(len(new_sig))
        try:
            new_sig = np.multiply(new_sig, window_func)
        except ValueError:
            print(len(new_sig), i, t)
            if t < window_size:
                new_sig = np.multiply(new_sig, window_func[:len(new - sig)])
            elif t > window_size:
                new_sig = np.multiply(new_sig, window_func[-len(new - sig):])
        if fft_type == 'fft':
            fft_sig = fftpack.fft(new_sig, n=zer_pad * window_size)[:window_size*zer_pad]
            #fft_sig = fftpack.rfft(new_sig, zer_pad * window_size)[:window_size]
        elif fft_type == 'welch':
            freqs, fft_sig = signal.welch(new_sig, nfft=zer_pad * window_size*zer_pad)
            fft_sig = fft_sig[1:window_size]
        if dc_cut == True:
            fft_sig = np.concatenate(
                ([fft_sig[1]], fft_sig[1:-1], [fft_sig[1]]))
        if fft_shift == 1:
            fft_sig = np.fft.fftshift(fft_sig)
        else:
            fft_sig=fft_sig[:len(fft_sig)/factor]
        fft_sig=abs(fft_sig)
        # if len(freq_mask)!=0:
        #     fft_sig[freq_mask == False] = np.nan
        #     fft_sig=fft_sig[freq_mask]
        # if normalize == True:
        #     fft_sig *= (1. / fft_sig.max())
        # if log == True:
        #     fft_sig = np.log(fft_sig)
        if filtered == 1:
            fft_sig = _butter_filter(fft_sig, b, a, zi)
        if 0:
            import matplotlib.pyplot as plt
            plt.figure('sfft')
            plt.clf()
            print(i, t, t + window_size, len(sig[t:t + window_size]))
            plt.plot(sig[t:t + window_size], 'b',label='signal')
            plt.plot(window_func, 'k',label='window')
            plt.plot(new_sig, 'r',label='signal w/ window')
            plt.legend(loc='best')
            plt.twinx()
            plt.plot(fft_sig, 'c')
            plt.draw()
            input('')

        matrix[i] = fft_sig

    if len(freq_mask)!=0:
        matrix=matrix[:,freq_mask]

    if normalize == True:
        matrix /= matrix.max(axis=1)[:, None]

    if log == True:
        matrix = np.log(matrix)

    if time_array is not None:
        return matrix.transpose(), time_spec, beat_freq
    else:
        return matrix.transpose()

def eval_beat_freq(time_array,window_size,step_scale=4, zer_pad=1,fft_shift=0):
    # evaluates acquisition rate
    acq_rate = 1 / (time_array[1] - time_array[0])
    # create beating frequency
    if fft_shift == 0:
        beat_freq = np.linspace(0, acq_rate / 2, num=window_size*zer_pad/2)
    elif fft_shift == 1:
        beat_freq = np.linspace(-acq_rate / 2, acq_rate / 2, num=window_size*zer_pad)
            # time array for spectrogram
        # SFFT step size,

    time_spec = np.linspace(time_array[window_size*0.25], time_array[-window_size*0.25], num=((len(time_array)-window_size/2)*step_scale / window_size))

    return time_spec,beat_freq

def eval_mask(beat_freq,window_size,freq_min,freq_max,zer_pad=1,fft_shift=0):
    # find the index position for min and max frequency
    fmin=abs(beat_freq-freq_min).argmin()
    if freq_max is None:
        fmax=len(beat_freq)
    else:
        fmax=abs(beat_freq-freq_max).argmin()
    factor=2
    if fft_shift:
        #fmin-=window_size/2
        #fmax-=window_size/2
        factor=1

    # creates mask array for frequency window:
    mask_array = np.arange(window_size*zer_pad/factor)
    mask = np.logical_and(mask_array > fmin, mask_array < fmax)
    #print(mask,mask_array)
    return fmin,fmax,mask

def _init_filter():
    from signal import butter, lfilter_zi
    # Create an order 3 lowpass butterworth filter.
    b, a = butter(3, 0.05)
    # Apply the filter to xn.  Use lfilter_zi to choose the initial condition
    # of the filter.
    zi = lfilter_zi(b, a)
    return b, a, zi


def _butter_filter(sig, b, a, zi):
    from signal import lfilter, filtfilt
    z, _ = lfilter(b, a, sig, zi=zi * sig[0])

    # Apply the filter again, to have a result filtered at an order
    # the same as filtfilt.
    z2, _ = lfilter(b, a, z, zi=zi * z[0])

    # Use filtfilt to apply the filter.
    return filtfilt(b, a, sig)

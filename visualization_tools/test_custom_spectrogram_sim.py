"""Compare custom spectrogram with scipy's implementation"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import specgram
from scipy import signal
import time
import sys
sys.path.insert(0, './../src/')

import custom_spectrogram as cs


tem = np.arange(1024 * 4)
w = tem * 0.0001
sinal = np.sin(tem * w * np.pi)
win = 32
zer_pad = 2
step_scale = 4

print('\n time cost:')
# measure custom spectrogram time
# for some unkown reason, the step must be scaled down to compare to scipy.
time0 = time.time()
for i in range(10):
    time_spec, beat_freq = cs.eval_beat_freq(tem, window_size=win, step_scale=step_scale / 4, zer_pad=zer_pad)
    freq_min, freq_max = 0.1, 0.4
    # mask=[]
    fmin, fmax, mask = cs.eval_mask(beat_freq, win, freq_min, freq_max, zer_pad=zer_pad)
    # print(len(mask),fmin,fmax,len(beat_freq))

    mat_cs = cs.spectrogram(
        sinal, window_size=win, zer_pad=zer_pad, step_scale=step_scale / 4, freq_mask=mask)
time1 = time.time()
print('\nCUSTOM: {0} ms'.format(100 * (time1 - time0)))

plt.subplot(4, 1, 1)
plt.plot(tem, sinal)
plt.title('signal')
plt.subplot(4, 1, 2)
plt.pcolormesh(time_spec, beat_freq[mask], mat_cs)
plt.plot(tem, w, 'k', lw=2)
plt.ylim(freq_min, freq_max)
plt.title('custom')

# measure scipy spectrogram time
time0 = time.time()
for i in range(10):
    freqs, tempo, Sxx_scipy = signal.spectrogram(sinal, fs=1, nperseg=win * zer_pad, noverlap=(1 - 1 / step_scale) * win)
time1 = time.time()
print('\nSCIPY: {0} ms'.format(100 * (time1 - time0)))

plt.subplot(4, 1, 3)
plt.pcolormesh(tempo, freqs, Sxx_scipy)
plt.ylim(freq_min, freq_max)
plt.plot(tem, w, 'k', lw=2)
plt.ylim(freq_min, freq_max)
plt.title('scipy')


# measure matplotlib spectrogram time
plt.subplot(4, 1, 4)
time0 = time.time()
for i in range(10):
    Sxx_mlab, freqs, bins = specgram(sinal, Fs=1, NFFT=win * zer_pad, noverlap=(1 - 1. / step_scale) * win)
time1 = time.time()
print('\nMATPLOTLIB: {0} ms'.format(100 * (time1 - time0)))
print('\ntime resolution:\n CUSTOM: {} \t SCIPY: {} \t MLAB: {}'.format(mat_cs.shape, Sxx_scipy.shape, Sxx_mlab.shape))
print('\ndark line is simulated frequency\n')
print(len(tempo))

plt.plot(tem, w, 'k', lw=2)
plt.pcolormesh(bins, freqs, Sxx_mlab)
plt.ylim(freq_min, freq_max)
plt.title('mlab')
plt.tight_layout()
plt.show()

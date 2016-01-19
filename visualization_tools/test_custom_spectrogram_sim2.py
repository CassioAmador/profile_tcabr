"""Compare custom spectrogram with scipy's implementation"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import specgram
from scipy import signal
import time
import sys
sys.path.insert(0, './../src/')

import custom_spectrogram as cs


fs = 100  # 100 MHz
tem = np.arange(int(8 * fs)) / fs
f = 12 - 4 * np.sqrt(np.arange(tem.size) / tem.size)
sinal = np.sin(2 * np.pi * f * tem) + 0.8 * np.sin(2 * np.pi * (f + 6) * tem) + 0.8 * np.sin(2 * np.pi * (f + 12) * tem)

nfft = 1* fs
window = 80
step_scale = 16
zer_pad = 8
print('\n time cost:')
# measure custom spectrogram time
# for some unkown reason, the step must be scaled down to compare to scipy.
freq_min, freq_max = 0, 15
time0 = time.time()
for i in range(10):
    time_spec, beat_freq = cs.eval_beat_freq(tem, window_size=window, step_scale=step_scale, zer_pad=zer_pad)
    fmin, fmax, mask = cs.eval_mask(beat_freq, window, freq_min, freq_max, zer_pad=zer_pad)
    Sxx_custom = cs.spectrogram(sinal, window_size=window, zer_pad=zer_pad, step_scale=step_scale, freq_mask=mask)
time1 = time.time()
print('\nCUSTOM: {0} ms'.format(100 * (time1 - time0)))

plt.subplots(4, 1)
ax1 = plt.subplot(411)
plt.plot(tem, sinal)
plt.xlim(tem[0], tem[-1])
plt.ylabel('signal')
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = plt.subplot(412, sharex=ax1)
plt.pcolormesh(time_spec, beat_freq[mask], Sxx_custom)
plt.plot(tem, f, 'k', lw=2)
plt.ylim(freq_min, freq_max)
plt.ylabel('custom')
plt.setp(ax2.get_xticklabels(), visible=False)

# measure scipy spectrogram time
time0 = time.time()
for i in range(10):
    freqs, tempo, Sxx_scipy = signal.spectrogram(sinal, fs, nperseg=nfft, noverlap=99, window=signal.get_window('hann', nfft), nfft=512)
time1 = time.time()
print('\nSCIPY: {0} ms'.format(100 * (time1 - time0)))
ax3 = plt.subplot(413)
plt.pcolormesh(tempo, freqs, Sxx_scipy)
plt.ylim(freq_min, freq_max)
plt.plot(tem, f, 'k', lw=2)
plt.ylim(freq_min, freq_max)
plt.ylabel('scipy')
plt.setp(ax3.get_xticklabels(), visible=False)


# measure matplotlib spectrogram time
plt.subplot(414, sharex=ax1)
time0 = time.time()
for i in range(10):
    Sxx_malab, freqs, bins = specgram(sinal, Fs=fs, NFFT=nfft, noverlap=99, pad_to=512)
time1 = time.time()
print('\nMATPLOTLIB: {0} ms'.format(100 * (time1 - time0)))
print('\ntime resolution:\n CUSTOM: {} \t SCIPY: {} \t MATPLOTLIB: {}'.format(Sxx_custom.shape, Sxx_scipy.shape, Sxx_malab.shape))
print('\ndark line is simulated frequency\n')

plt.plot(tem, f, 'k', lw=2)
plt.pcolormesh(bins, freqs, Sxx_malab)
plt.ylim(freq_min, freq_max)
plt.ylabel('mlab')
plt.tight_layout(h_pad=0)
plt.show()

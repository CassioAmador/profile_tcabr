"""
Function: Process data for sweep frequency mode.
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: zero padding in filter should check first the size of the signal,
so it could could choose best value depending on signal size. Maybe the
same for padding in spectrogram.
Spectrogram could use a pre-sized matrix.
Some functions could be modified to work with hopping frequency.
"""

import numpy as np
from scipy import signal

from read_signal import ReadSignal
import custom_spectrogram as cs


class ProcSweep(ReadSignal):

    """Process sweep data. All data is evaluated at its current sweep."""

    def __init__(self, shot, tipo='data', save_locally=1):
        # Inherit from ReadSignal class
        ReadSignal.__init__(self, shot, tipo)
        self.save_locally = save_locally
        # evaluate size of points per sweep
        self.sweep_size = np.round(self.rate * self.sweep_dur)
        # evaluate frequency steps:
        self.find_sweep_points()
        # evaluate frequency probing window
        self.frequency_analysis()

    def find_sweep_points(self):
        """Create a list of points (called 'points') where the sweep
        started, with the 'time' channel (channel 4)."""
        self.read_channel('time', self.save_locally)
        # find min value for data, in the first 2 sweeps
        mindata = self.bindata['time'][:self.rate * 2 *
                                       (self.sweep_dur + self.interv_sweep)].min()
        # find all points with this value
        zeros = np.where(self.bindata['time'] == mindata)[0]
        # mark only the first minimum for each trigger
        singlezeros = np.where(np.diff(zeros) > 50)[0]
        # list starting from the 2nd sweep
        self.points = zeros[singlezeros + 1]
        print("Total number of sweeps: %s" % len(self.points))

    def time2sweep(self, time):
        """Convert a position in time (ms) to the correspondent
        sweep position."""
        if not hasattr(self, 'points'):
            self.find_sweep_points()
        # finds nearest sweep from a specified time
        self.sweep_cur = (abs(self.points - time * self.rate * 1e3)).argmin()
        return self.sweep_cur

    def sweep2time(self, sweep=0):
        """Convert a sweep position to its position in time (ms).
        Defaults to current sweep (sweep_cur)"""
        if not hasattr(self, 'points'):
            self.find_sweep_points()
        if sweep == 0:
            sweep = self.sweep_cur
        return (self.points[sweep] / self.rate) * 1e-3

    def frequency_analysis(self):
        """ Frequency limits for the probing frequency.
        Tests showed that 17.5 GHz and 26.5 GHz for K and Ka bands are the minimum values
        that seems to be physic relevant. But, depending on the size of the spectrogram, they are current
        out because of the window size."""
        self.probing_freq_lim = {}
        self.probing_freq_lim['K'] = (16, 100)
        self.probing_freq_lim['Ka'] = (24, 100)
        self.freqs_start = {}
        self.freqs_end = {}
        self.probing_frequency = {}
        self.mask_probe_freq = {}
        for channel in ('K', 'Ka'):
            probe_freq = np.linspace(
                self.freq_start, self.freq_end, self.sweep_size) * self.chan_factor[channel]
            self.mask_probe_freq[channel] = np.logical_and(probe_freq >= self.probing_freq_lim[
                                                           channel][0], probe_freq <= self.probing_freq_lim[channel][1])

            self.probing_frequency[channel] = probe_freq[
                self.mask_probe_freq[channel]]
            self.freqs_start[channel] = min(self.probing_frequency[channel])
            self.freqs_end[channel] = max(self.probing_frequency[channel])

    def read_single_sweep(self, channel, sweep_cur=100):
        """Read data for the current sweep (sweep_cur), for an
        specific channel"""
        self.sweep_cur = sweep_cur
        if not hasattr(self, 'points'):
            self.find_sweep_points()
        try:
            samples_ini = self.points[self.sweep_cur]
        except IndexError:
            print(len(self.points), self.sweep_cur)
            samples_ini = self.points[self.sweep_cur - self.sweep_size]
        samples_end = samples_ini + self.sweep_size
        if channel not in self.bindata.keys():
            self.read_channel(channel, self.save_locally)
        if hasattr(self, 'single_sweep_data') == False:
            self.single_sweep_data = {}
        self.single_sweep_data[channel] = self.bindata[
            channel][samples_ini:samples_end][self.mask_probe_freq[channel]]

    def single_sweep_time(self, channel, time=30):
        """same as read_single_sweep, but it reads int time (ms)"""
        self.time2sweep(time)
        self.read_single_sweep(channel, self.sweep_cur)

    def signal_filter(self, channel, freqs=(1e3, 15e3)):
        """Filter signal from specific channel. A FFT is performed in
        the signal, and the result is multiplied by a window function
        (kaiser function), which nulls the undesired beating
        frequencies, that are outside of 'freqs'. The signal is then
        recreated with an IFFT."""
        # zero padding size:
        zer_pad_filter = 4
        # alias for sweep_size
        N = int(self.sweep_size)
        # FFT with zero padding
        fft = np.fft.rfft(self.single_sweep_data[channel], zer_pad_filter * N)
        # bp=fft[:] #used for other plot functions
        fmin, fmax = freqs
        # creates the beating frequency axis, and finds the position of
        # the frequency limits in the axis
        fft_freq = np.linspace(0, self.rate * 1e3 / 2.,
                               num=zer_pad_filter * N / 2)
        cmin = (abs(fft_freq - fmin)).argmin()
        cmax = (abs(fft_freq - fmax)).argmin()
        # creates window function for filter. Second argument of kaiser
        # function must be float
        window_func = np.concatenate((np.zeros(cmin + 1),
                                      np.kaiser(cmax - cmin, 2.), np.zeros(zer_pad_filter * N / 2 - cmax)))
        # multiply the window by the signal's FFT.
        bp = np.multiply(fft, window_func)
        # Since the signal is REAL, we use IRFFT, and takes only the
        # first half, since the other half is symmetric.
        newsig = np.fft.irfft(bp)[:N]
        return newsig

    def plot_sweep(self, channel):
        """Plot binary data for an specific sweep and channel."""
        import pylab as p

        if not hasattr(self, 'sweep_freq'):
            # dict with arrays for the frequencies in each channel.
            self.sweep_freq = {}
        if channel not in self.sweep_freq.keys():
            self.sweep_freq[channel] = np.linspace(
                self.freqs_start[channel], self.freqs_end[channel], num=self.sweep_size)
        p.plot(self.sweep_freq[channel], self.single_sweep_data[
               channel], label="Channel: %s" % channel)
        p.xlabel("freq (GHz)")
        p.ylabel("beating signal")

    def spectrogram(self, channel, window=256, step_scale=16, zer_pad=8, log=0,
                    group_delay=1, figure=0, normal=0, filtered=0,
                    beating_freq_filter=(1, 15)):
        """Evaluate and plot spectrogram (SFFT) of beating signal.
        Some parameters listed (others can be found in the function):
        group_delay=1 evaluates group delay.
                    0 evaluates beating frequency
        normal=1 normalize spectrum
        log=0
            1 for log spectrogram
        """

        if not hasattr(self, 'Dt_DF'):
            self.Dt_DF = {}
            self.X = {}
            self.delta_freq = {}
            self.index_X = {}
            self.Y = {}
            self.mask_bf = {}
        if channel not in self.Dt_DF.keys():
            tem = np.linspace(0, self.sweep_dur, self.sweep_size)
            tem = tem[self.mask_probe_freq[channel]]
            self.time_spec, beat_freq = cs.eval_beat_freq(
                tem, window_size=window, step_scale=step_scale, zer_pad=zer_pad)
            fmin, fmax, mask_bf = cs.eval_mask(beat_freq, window, beating_freq_filter[
                                               0], beating_freq_filter[1], zer_pad=zer_pad)
            self.mask_bf[channel] = mask_bf
            # print(len(mask_bf),fmin,fmax,len(beat_freq))
            #print(mask_bf, beating_freq_filter, beat_freq)

            # Inverse of dF/dt sweeping rate:
            self.Dt_DF[channel] = (tem[1] - tem[0]) / \
                (self.probing_frequency[channel][1] -
                 self.probing_frequency[channel][0])

            # X is and array with the probing frequency.
            self.X[channel] = np.linspace(self.probing_frequency[channel][
                                          window / 2], self.probing_frequency[channel][-window / 2], \
                                          (len(tem) - window) * step_scale / window)

            self.delta_freq[channel] = self.X[channel][1] - self.X[channel][0]

           # Y is the beating frequency, in MHz, or the group delay, in ns
            self.Y[channel] = beat_freq[self.mask_bf[channel]]
            if group_delay == 1:
                # group delay in ns
                self.Y[channel] *= self.Dt_DF[channel]

        if filtered == 1:
            sig = self.signal_filter(channel, beating_freq_filter)
        else:
            sig = self.single_sweep_data[channel]

        mat_cs = cs.spectrogram(
            sig, window_size=window, zer_pad=zer_pad, step_scale=step_scale, normalize=normal, freq_mask=self.mask_bf[channel])

        return abs(mat_cs)

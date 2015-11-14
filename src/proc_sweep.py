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

from read_signal import ReadSignal


class ProcSweep(ReadSignal):

    """Process sweep data. All data is evaluated at its current sweep."""

    def __init__(self, shot, tipo='data', save_locally=1):
        # Inherit from ReadSignal class
        ReadSignal.__init__(self, shot, tipo)
        # multiplying factor for channel frequencies
        self.chan_factor = {"K": 2, "Ka": 3, "ref": 1, "time": 1}
        # evaluate size of points per sweep
        self.sweep_size = np.round(self.rate * self.sweep_dur)
        self.save_locally = save_locally

    def find_sweep_points(self):
        """Create a list of points (called 'points') where the sweep
        started, with the 'time' channel (channel 4)."""
        self.read_channel('time', self.save_locally)
        # find min value for data, in the first 2 sweeps
        mindata = self.bindata['time'][:self.rate * 2 * (self.sweep_dur + self.interv_sweep)].min()
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

    def sweep2time(self, sweep):
        """Convert a sweep position to its position in time (ms)."""
        if not hasattr(self, 'points'):
            self.find_sweep_points()
        return (self.points[sweep]/self.rate) * 1e-3

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
            samples_ini = self.points[self.sweep_cur-self.sweep_size]
        samples_end = samples_ini + self.sweep_size
        if channel not in self.bindata.keys():
            self.read_channel(channel, self.save_locally)
        if hasattr(self, 'single_sweep_data') == False:
            self.single_sweep_data = {}
        self.single_sweep_data[channel] = self.bindata[channel][samples_ini:samples_end]

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
        fft_freq = np.linspace(0, self.rate * 1e3 / 2., num=zer_pad_filter * N / 2)
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
            self.sweep_freq[channel] = np.linspace(self.freq_start, self.freq_end, num=self.sweep_size) * self.chan_factor[channel]
        p.plot(self.sweep_freq[channel], self.single_sweep_data[channel], label="Channel: %s" % channel)
        p.xlabel("freq (GHz)")
        p.ylabel("beating signal")

    def spectrogram(self, channel, window_scale=4, step_scale=16, log=0,
                    group_delay=1, figure=0, normal=1, filtered=1,
                    freqs=(1e3, 15e3)):
        """Evaluate and plot spectrogram (SFFT) of beating signal.
        Some parameters listed (others can be found in the function):
        group_delay=1 evaluates group delay.
                    0 evaluates beating frequency
        normal=1 normalize spectrum
        log=0
            1 for log spectrogram
        """
        # scale for zero padding
        zer_pad = 4
        # alias for sweep_size
        N = int(self.sweep_size)
        # SFFT window size
        window = N / window_scale
        # SFFT step size, by
        step = window / step_scale
        # creates the window function that will 'move' trough the signal
        # to evaluate each FFT.
        window_func = np.concatenate((np.zeros(N - window / 2), np.hanning(window), np.zeros(N - window / 2)))
        if filtered == 1:
            sig = self.signal_filter(channel, freqs)
        else:
            sig = self.single_sweep_data[channel]
        # create a matrix to receive the spectra
        matrix = np.empty(shape=(1 + (N - window) / step, 2 * N))
        # loop trough all the possible windows, and evaluates the FFT.
        for i in range(1 + int((N - window) / step)):
            t = i * step + window / 2
            new_window = window_func[N - t:N + N - t]
            new_sig = np.multiply(sig, new_window)
            # We ignore the fft's first point (DC component).
            fft_sig = np.fft.rfft(new_sig, zer_pad * N)[1:]
            if log == 1:
                fft_sig = np.log(fft_sig)
            if normal == 1:
                fft_sig *= (1. / fft_sig.max())
            matrix[i] = abs(fft_sig)

        if not hasattr(self,'Dt_DF'):
            self.Dt_DF={}
        if channel not in self.Dt_DF.keys():
            # Inverse of dF/dt sweeping rate:
            self.Dt_DF[channel] = self.sweep_dur / ((self.freq_end - self.freq_start) * self.chan_factor[channel])

        if not hasattr(self,'X'):
            self.X={}
        if channel not in self.X.keys():
            # X is and array with the probing frequency.
            self.X[channel] = np.linspace(self.freq_start, self.freq_end, num=len(matrix)) * self.chan_factor[channel]
        if not hasattr(self,'Y'):
            self.Y={}
        if channel not in self.Y.keys():
            # Y is the beating frequency, in MHz, or the group delay, in ns
            self.Y[channel] = np.linspace(0, self.rate / 2., num=zer_pad * N / 2)
            if group_delay == 1:
                    # group delay in ns
                    self.Y[channel] *= self.Dt_DF[channel]

        # transpose matrix for spectrogram.
        matrix = matrix.transpose()
        return matrix
    ARRUMAR

    def spectrogram(self, channel, window_scale=4, step_scale=16, log=0,
                    group_delay=1, ploti=1, figure=0, normal=1, filtered=1,
                    freqs=(1e3, 15e3)):
        """Evaluate and plot spectrogram (SFFT) of beating signal.
        Some parameters listed (others can be found in the function):
        group_delay=1 evaluates group delay.
                    0 evaluates beating frequency
        ploti=2 for contour plot
              1 for animation. Don't try if running trough network.
        normal=1 normalize spectrum
        log=0
            1 for log spectrogram
        """
        if ploti > 0:
            import pylab as p

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
            sig = self.singlesweep_data[channel]
        if ploti == 2:
            p.ion()
            a = input('press any key to begin plot. 0 to cancel:   ')
            if a == '0':
                ploti = 1
                p.ioff()
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
            if ploti == 2:
                p.figure(5)
                p.clf()
                p.subplot(2, 1, 1)
                p.plot(abs(fft_sig))
                p.subplot(2, 1, 2)
                p.plot(sig)
                p.plot(new_sig)
                p.twinx()
                p.plot(new_window, 'r-')
                p.draw()
                input('')
        if ploti == 2:
            p.ioff()
            pass
        # transpose matrix for spectrogram.
        matrix = matrix.transpose()
        # creates arrays with beating frequency and band frequency.
        X = np.linspace(self.freq_start, self.freq_end, num=len(matrix[0])) * self.chan_factor[channel]
        # Y is the frequency, in MHz
        Y = np.linspace(0, self.rate / 2., num=zer_pad * N / 2)
        if group_delay == 1:
            # Inverse of dF/dt sweeping rate:
            Dt_DF = self.sweep_dur / (max(X) - min(X))
            # group delay in ns
            Y *= Dt_DF
        # print(N,len(matrix),len(matrix[0]),len(X),len(Y))
        if ploti >= 1:
            if figure == 0:
                p.figure()
            else:
                p.figure(figure)
            p.contourf(X, Y, matrix)
            p.ylabel('beating freq (MHz)')
            if group_delay == 1:
                p.ylabel('group delay (ns)')
            p.xlabel('sweep freq (GHz)')
            if (filtered == 1) & (group_delay == 0):
                p.ylim(freqs[0] * 1e-3, freqs[1] * 1e-3)
            if group_delay == 1:
                p.ylim(freqs[0] * 1e-3 * Dt_DF, freqs[1] * 1e-3 * Dt_DF)
        return matrix, X, Y
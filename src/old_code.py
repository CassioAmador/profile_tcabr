# polynomial coeffs of order 15, for the abel integral inversion
abel_factor = p.array([2048 / 6435., 429 * p.pi / 4096., 1024 / 3003., 231 * p.pi / 2048.,
                       256 / 693., 63 * p.pi / 512, 128 / 315., 35 * p.pi / 256., 16 / 35.,
                       5 * p.pi / 32., 8 / 15., 3 * p.pi / 16., 2 / 3., p.pi / 4., 1, p.pi / 2.], dtype=p.float32)

class ProcProfile(ProcGroupDelay):

    """Recreate profiles from specific shot."""

    def smooth_signal(self, signal, window=7, order=3):
        """ Smooth signal. Work for scipy version >= 0.14.0."""
        from scipy import signal
        signal_smoothed = signal.savgol_filter(signal, window_length=window, polyorder=order)
        return signal_smoothed

    def profile_poly(self, order=9, all_shot=0):
        """Recreate a profile fiiting the group delay with a polynomial
        of specified order. The abel inversion is evaluated explicitly,
        since for a poly fit the integral turns out to be another
        polynomial function. More information on the manual."""
        polynomial = p.polyfit(self.freqs, self.gd, order)
        self.pol = p.poly1d(polynomial)
        if not hasattr(self, 'freq_poly'):
            self.freq_poly = p.linspace(0, self.freqs[-1], num=100)
        if not hasattr(self, 'ne_poly'):
            self.ne_poly = rf.freq2den(self.freq_poly)
        if all_shot == 0:
            # only needed if one wants to check the polynomial fit.
            self.gd_poly = self.pol(self.freq_poly)
        self.pol_coeffs = self.pol.coeffs * abel_factor[-order - 1:]
        self.r = const * p.array([p.polyval(self.pol_coeffs, freq) for freq in self.freq_poly])

    def abel_factor_2(self, f_c, f_i, f_f, order=2):
        # polynomial coeffs for any order k, for the abel integral inversion starting from non-zero frequency.
        # f_c: frequency of reflection
        # f_i: lower limit of integration
        # f_f: upper limit of integration
        self.abel = p.zeros(order+1)
        if f_c == f_f:
            self.abel[0] = p.arccos(f_i/f_c)
            for k in range(1, order+1):
                self.abel[k] = math.pow(f_c, k) * p.sqrt(p.pi) * sp_sp.gamma((1 + k)/2)/(k * sp_sp.gamma(k/2)) - math.pow(f_i, 1 + k) * sp_sp.hyp2f1(1/2, (1 + k)/2, (3 + k)/2, math.pow(f_i/f_c, 2))/(f_c*(1 + k))
        if f_c != f_f:
            if f_i == 0:
                for k in range(order+1):
                    self.abel[k] = math.pow(f_f, 1 + k) * sp_sp.hyp2f1(1, (1 + k)/2, (3 + k)/2,  math.pow(f_f/f_c, 2))/(f_c * (1 + k))
            else:
                for k in range(order+1):
                    self.abel[k] = 1/(f_c * (1 + k)) * (-math.pow(f_i, 1 + k) * sp_sp.hyp2f1(1, (1 + k)/2, (3 + k)/2, math.pow(f_i/f_c, 2)) + math.pow(f_f, 1 + k) * sp_sp.hyp2f1(1, (1 + k)/2, (3 + k)/2,  math.pow(f_f/f_c, 2)))
        # if order>=0:
        #     self.abel[0]=1.5708 -(A * sp_sp.hyp2f1(0.5,1./2,3./2,math.pow(A/B,2)))/B
        # if order>=1:
        #     self.abel[1]=math.sqrt(-math.pow(A,2)+math.pow(B,2))
        # if order>=2:
        #     self.abel[2]=0.785398 * math.pow(B,2)-(0.333333 * math.pow(A,3) * sp_sp.hyp2f1(0.5,3./2,5./2,math.pow(A/B,2)))/B
        # if order>=3:
        #     self.abel[3]=math.sqrt(-math.pow(A,2)+math.pow(B,2)) * (0.333333 * math.pow(A,2)+0.666666 * math.pow(B,2))
        # if order>=4:
        #     self.abel[4]=0.589049 * math.pow(B,4)-(0.2 * math.pow(A,5) * sp_sp.hyp2f1(0.5,5./2,7./2,math.pow(A/B,2)))/B
        # if order>=5:
        #     self.abel[5]=math.sqrt(-math.pow(A,2)+math.pow(B,2)) * (0.2 * math.pow(A,4)+0.266667 * math.pow(A,2) * math.pow(B,2)+0.533333 * math.pow(B,4))
        # if order>=6:
        #     self.abel[6]=0.490874 * math.pow(B,6)-(0.142857 * math.pow(A,7) * sp_sp.hyp2f1(0.5,7./2,9./2,math.pow(A/B,2)))/B
        # if order>=7:
        #     self.abel[7]=math.sqrt(-math.pow(A,2)+math.pow(B,2)) * (0.142857*  math.pow(A,6)+0.171429 * math.pow(A,4) * math.pow(B,2)+0.228571 * math.pow(A,2) * math.pow(B,4)+0.457143 * math.pow(B,6))
        # if order>=8:
        #     self.abel[8]=0.429515 * math.pow(B,8)-(0.111111 * math.pow(A,9) * sp_sp.hyp2f1(0.5,9./2,11./2,math.pow(A/B,2)))/B
        # if order>=9:
        #     self.abel[9]=math.sqrt(-math.pow(A,2)+math.pow(B,2)) * (0.111111 * math.pow(A,8)+0.126984 * math.pow(A,6) * math.pow(B,2)+0.152381 * math.pow(A,4) * math.pow(B,4)+0.203175 * math.pow(A,2) * math.pow(B,6)+0.406349 * math.pow(B,8))
        # if order>=10:
        #     self.abel[10]=0.386563 * math.pow(B,10)-(0.0909091 * math.pow(A,11) * sp_sp.hyp2f1(0.5,11./2,13./2,math.pow(A/B,2)))/B

    def profile_poly_6(self, order=9, all_shot=0):
        """Linear initialization, group delay dividec in two parts"""
        self.r = p.zeros(1+len(self.X))
        polynomial = p.polyfit(self.X, self.gd_m, order)
        self.pol = p.poly1d(polynomial)
        if not hasattr(self, 'freqs'):
            self.freqs = p.concatenate(([0], self.X))
        if not hasattr(self, 'ne'):
            self.ne = rf.freq2den(self.freqs*1e9)
        if all_shot == 0:
            # only needed if one wants to check the polynomial fit.
            self.gd_poly = self.pol(self.X)
            self.gd = p.concatenate(([0], self.gd_poly))
        coeff = (self.gd_m[0]/self.X[0]) * 1e-18
        self.abel_factor_2(self.X[0], 0, self.X[0], 1)
        self.r[0] = const * coeff * self.abel[1]
        for f in range(1, len(self.X)):
            self.abel_factor_2(self.X[f], 0, self.X[0], 1)
            r_0 = const * coeff * self.abel[1]
            self.abel_factor_2(self.X[f], self.X[0], self.X[f], order)
            print(r_0)
            self.r[f] = r_0+const * p.polyval(self.pol.coeffs * self.abel[::-1], self.X[f])

    def profile_poly_2(self, order=2, all_shot=0):
        if self.no_plasma:
            return
        init = 2
        f_probe = self.X
        tau = self.gd_m
        self.r = np.zeros(len(tau))
        tau2_coef = np.polyfit(f_probe, tau, order)
        if not hasattr(self, 'freq_poly'):
            self.freq_poly = self.X * 1e9
        if not hasattr(self, 'ne_poly'):
            self.ne_poly = rf.freq2den(self.freq_poly)
        # inicializao
        if init == 1:
            tau1_coef = np.array([np.polyval(tau2_coef, f_probe[0]) / f_probe[0], 0])
        elif init == 2:
            tau1_coef = np.array([np.polyval(tau2_coef, f_probe[0]) / f_probe[0] ** 2., 0, 0])
        I1 = np.zeros(len(tau1_coef))
        I2 = np.zeros(len(tau2_coef))
        I3 = np.zeros(len(tau2_coef))
        for i in range(len(f_probe)):
            for k in range(len(tau1_coef)):
                I1[k] = (np.power(f_probe[0], k + 1) / f_probe[i]) * \
                    (sp_sp.hyp2f1(1 / 2., (k + 1) / 2., (k + 3) / 2., np.power(f_probe[0] / f_probe[i], 2))) / (k + 1.)
            for k in range(len(tau2_coef)):
                I2[k] = (np.power(f_probe[0], k + 1) / f_probe[i]) * \
                    (sp_sp.hyp2f1(1 / 2., (k + 1) / 2., (k + 3) / 2., np.power(f_probe[0] / f_probe[i], 2))) / (k + 1.)
                I3[k] = np.power(f_probe[i], k) * np.sqrt(np.pi) * sp_sp.gamma(0.5 + k / 2.) / (2. * sp_sp.gamma(1 + k / 2.))
            self.r[i] = 1e-9 * (3e8 / np.pi) * (np.dot(tau1_coef[::-1], I1) + np.dot(tau2_coef[::-1], I3 - I2))

    def profile_poly_3(self, order=9, all_shot=0):
        """..."""
        # choose number of steps to start group_delay
        if not hasattr(self, 'freqs'):
            self.freqs = p.concatenate((p.zeros(1), self.X))
        if not hasattr(self, 'ne'):
            self.ne = rf.freq2den(self.freqs*1e9)
        self.gd = p.concatenate((p.zeros(1), self.gd_m))
        self.r = p.zeros(len(self.freqs))
        order = 1
        for f in range(1, len(self.freqs)-order+1):
            polynomial = p.polyfit(self.freqs[f-1:f+order+1], self.gd[f-1:f+order+1], order)
            self.pol = p.poly1d(polynomial)
            self.abel_factor_2(self.freqs[f-1], self.freqs[f-1], self.freqs[f+order-1], order)
    #       p.plot(self.freqs[f-1:f+order+1],self.gd[f-1:f+order+1],'.')
    #       p.plot(self.freqs[f-1:f+order+1],self.pol(self.freqs[f-1:f+order+1]))
    #       polynomial coefficients are in different order
            self.pol_coeffs = self.pol.coeffs * self.abel[::-1][-order-1:]
            self.r[f] = self.r[f-1] + const * p.polyval(self.pol_coeffs, self.freqs[f])
            print(self.r[f], self.r[f-1])

    def profile_poly_4(self, all_shot=0):
        """..."""
        # choose number of steps to start group_delay
        if (not hasattr(self, 'freqs')) or (self.freqs.size != (self.X.size-1)):
            self.freqs = p.concatenate((p.zeros(1), self.X))
        if (not hasattr(self, 'ne')) or (self.ne.size != (self.X.size-1)):
            self.ne = rf.freq2den(self.freqs*1e9)
        self.gd = p.concatenate((p.zeros(1), self.gd_m))
        self.r = p.zeros(len(self.freqs))
        meth = 1
        for f in range(1, len(self.freqs)):
            pos = 0
            for ff in range(f+1):
                if f != ff:
                    ss = (self.gd[ff]+self.gd[ff-1])*(self.freqs[ff]-self.freqs[ff-1])/(2*p.sqrt(math.pow(self.freqs[f], 2)-math.pow(self.freqs[ff], 2)))
                else:
                    if meth == 0:
                        ss = self.gd[ff]/1e-9
                    else:
                        ss = 1.5708 - (self.freqs[f-1] * sp_sp.hyp2f1(0.5, 0.5, 1.5, math.pow(self.freqs[f-1]/self.freqs[f], 2)))/self.freqs[f]
                        ss = self.gd[f-1]*ss + ((self.gd[f]-self.gd[f-1])/(self.freqs[f]-self.freqs[f-1])) * math.sqrt(-math.pow(self.freqs[f-1], 2) + math.pow(self.freqs[f], 2))
                pos += ss
            self.r[f] = const * pos

    def abel_transform(self, tau, f_probe, order=2, init=1):
        """Abel invertion for group delay (tau), related to the probing frequency f_probe."""
        f_probe = f_probe * 1e-9
        tau = tau * 1e9
        rc = np.zeros(len(tau))
        tau2_coef = np.polyfit(f_probe, tau, order)
        # init
        if init == 1:
            tau1_coef = np.array([np.polyval(tau2_coef, f_probe[0]) / f_probe[0], 0])
        elif init == 2:
            tau1_coef = np.array([np.polyval(tau2_coef, f_probe[0]) / f_probe[0] ** 2., 0, 0])
        I1 = np.zeros(len(tau1_coef))
        I2 = np.zeros(len(tau2_coef))
        I3 = np.zeros(len(tau2_coef))
        for i in range(len(f_probe)):
            for k in range(len(tau1_coef)):
                I1[k] = (np.power(f_probe[0], k + 1) / f_probe[i]) * (sp_sp.hyp2f1(1 / 2., (k + 1) / 2., (k + 3) / 2., np.power(f_probe[0] / f_probe[i], 2))) / (k + 1.)
            for k in range(len(tau2_coef)):
                I2[k] = (np.power(f_probe[0], k + 1) / f_probe[i]) * (sp_sp.hyp2f1(1 / 2., (k + 1) / 2., (k + 3) / 2., np.power(f_probe[0] / f_probe[i], 2))) / (k + 1.)
                I3[k] = np.power(f_probe[i], k) * np.sqrt(np.pi) * sp_sp.gamma(0.5 + k / 2.) / (2. * sp_sp.gamma(1 + k / 2.))
            rc[i] = 0.18 - 1e-9 * (3e8 / np.pi) * (np.dot(tau1_coef[::-1], I1) + np.dot(tau2_coef[::-1], I3 - I2))
        return rc


class ProcSweep(self):
    def spectrogram2(self, channel, window_scale=4, step_scale=16, log=0,
                     group_delay=1, figure=0, normal=1, filtered=1,
                     freqs=(1e3, 15e3), probing_freqs=(0, 100)):
        """Evaluate and plot spectrogram (SFFT) of beating signal.
        Some parameters listed (others can be found in the function):
        group_delay=1 evaluates group delay.
                    0 evaluates beating frequency
        normal=1 normalize spectrum
        log=0
            1 for log spectrogram
        """
        if filtered == 1:
            sig = self.signal_filter(channel, freqs)
        else:
            sig = self.single_sweep_data[channel]
        nfft = self.rate * 1.5
        f, t, Sxx = signal.spectrogram(
            sig, self.rate, nperseg=nfft, noverlap=nfft - 1, window=signal.get_window('hann', nfft), nfft=3 * nfft)
        if normal == 1:
            Sxx = Sxx / Sxx.max(axis=0)

        if not hasattr(self, 'Dt_DF'):
            self.Dt_DF = {}
        if channel not in self.Dt_DF.keys():
            # Inverse of dF/dt sweeping rate:
            self.Dt_DF[channel] = self.sweep_dur / \
                ((self.freq_end - self.freq_start) * self.chan_factor[channel])

        if not hasattr(self, 'X'):
            self.X = {}
        if not hasattr(self, 'Y'):
            self.Y = {}
        # X is and array with the probing frequency.
        self.X[channel] = (self.freq_start + (self.freq_end -
                                              self.freq_start) * t / t[-1]) * self.chan_factor[channel]
        index_X = np.logical_and(self.X[channel] > probing_freqs[
                                 0], self.X[channel] < probing_freqs[1])
        self.X[channel] = self.X[channel][index_X].copy()
        # Y is the beating frequency, in MHz, or the group delay, in ns
        self.Y[channel] = f
        if group_delay == 1:
            # group delay in ns
            self.Y[channel] *= self.Dt_DF[channel]
        return Sxx[:, index_X]
"""
Function: Evaluates phase derivative and group delay for profile reconstruction.
Reference is vessel wall, and calibration was done with a pin. After calibration,
signal from both channels are overlapped.
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: test if end of band K is more trustable than beginning of band Ka
Make a better info file with parameters used.
"""

from os import path
import numpy as np
from scipy import integrate, interpolate

import ref_functions as rf
from proc_sweep import ProcSweep

# CONSTANTS
# distance (39.5 cm) in nanoseconds from vessel wall to limiter, x2
wall_corr = 2.6351563520654012


class ProcGroupDelay(ProcSweep):

    """Evaluates phase derivative and group delay for profile reconstruction.
    Reference is vessel wall, and calibration was done with a pin. After calibration, 
    signal from both channels are overlapped."""

    def __init__(self, shot, tipo="data", save_locally=1):
        self.version = 1.0
        # inherit from ProcSweep class
        ProcSweep.__init__(self, shot, tipo, save_locally)
        self.find_sweep_points()

    def average_specgram(self, sweeps=8, sweep_ini=20, all_shot=0, window=256, step_scale=16,zer_pad=8):
        """Average spectrograms of a specified cluster of sweeps."""
        # if the class has already a mean stored, deletes it.
        if hasattr(self, "matrix_k_mean"):
            del(self.matrix_k_mean)
        if hasattr(self, "matrix_ka_mean"):
            del(self.matrix_ka_mean)
        # makes the spectrogram for each sweep.
        for sweep in range(sweep_ini, sweep_ini + sweeps):
            self.read_single_sweep("K", sweep)
            self.read_single_sweep("Ka", sweep)
            if all_shot == 0:
                print(self.sweep_cur)
            matrix_k = self.spectrogram(
                'K', figure=sweep + 1, normal=0, beating_freq_filter=(4, 15),window=window,step_scale=step_scale, zer_pad=zer_pad)
            matrix_ka = self.spectrogram(
                'Ka', figure=sweep + 1, normal=0, beating_freq_filter=(5, 17),window=window,step_scale=step_scale, zer_pad=zer_pad)
            if hasattr(self, 'matrix_k_mean'):
                self.matrix_k_mean += matrix_k
                self.matrix_ka_mean += matrix_ka
            # if there are no mean, creates it.
            else:
                self.matrix_k_mean = matrix_k.copy()
                self.matrix_ka_mean = matrix_ka.copy()

        self.matrix_k_mean /= sweeps
        self.matrix_ka_mean /= sweeps

    def reference_gd(self, sw_clustersize=20, all_shot=0):
        """Average the spectrogram and evaluates the group delay for the
        reference position. In this case, the reflection in the wall is
        stronger than in the pin, so we evaluate the time delay of the
        wall, and subtracts the position of limiter. The calibration
        is done with a pin, and stored in 'wall_cor' parameter.
        Average is done in the first and/or last sweeps of the shot,
        when there is no plasma"""
        self.average_specgram(sweeps=sw_clustersize,
                              sweep_ini=1, all_shot=all_shot)
        self.gd_k0 = self.Y['K'][self.matrix_k_mean.argmax(axis=1)] - wall_corr
        self.gd_ka0 = self.Y['Ka'][
            self.matrix_ka_mean.argmax(axis=1)] - wall_corr

    def plasma_gd(self, sweep_cur=1000, sw_clustersize=8, all_shot=0):
        """Evaluate group delay for a cluster of sweeps, and subtracts
        the reference to find the group delay from the limiter"""
        self.average_specgram(sweeps=sw_clustersize,
                              sweep_ini=sweep_cur, all_shot=all_shot)
        self.gd_k_mean = self.Y['K'][self.matrix_k_mean.argmax(axis=1)]
        self.gd_ka_mean = self.Y['Ka'][self.matrix_ka_mean.argmax(axis=1)]
        self.gd_k = self.gd_k_mean - self.gd_k0
        self.gd_ka = self.gd_ka_mean - self.gd_ka0

    def eval_freq_overlap(self):
        if not hasattr(self, 'freqs'):
            # index position of overlap in K and Ka bands
            self.cmin = None
            self.cmax = None
            # check the number of points of overlaped frquency
            self.freq_k_max = self.X['K'][-1]
            self.freq_ka_min = self.X['Ka'][0]
            delta_freq_gap = self.freq_ka_min - self.freq_k_max
            if delta_freq_gap < 0:
                # index position of overlap in K band
                self.cmax = np.where((self.X['K'] - min(self.X['Ka'])) > 0)[0][0]
                # frequency in K band that starts to overlap
                self.freq_k_max = self.X['K'][self.cmax]
                # index position of overlap in Ka band
                self.cmin = np.where(
                    (self.X['Ka'] - max(self.X['K'])) < 0)[0][-1] + 1
                # frequency in Ka band that finishs the overlap
                print(self.cmin,self.cmax)
                self.freq_ka_min = self.X['Ka'][self.cmin]
                # creates an array with all the overlapping frequencies
                self.freqs_over = np.sort(np.concatenate(
                    (self.X['K'][self.cmax:], self.X['Ka'][:self.cmin])))

            else:
                self.freqs_over = np.array([])
                # delta_freq_ka = self.X['Ka'][1] - self.X['Ka'][0]
                # self.freqs_over = np.empty()
                # # if gap is larger than a frequency step in Ka band, add points
                # if delta_freq_gap > (delta_freq_ka):
                #     self.freqs_over = np.arange(
                #         self.freq_k_max + delta_freq_ka, self.freq_ka_min, delta_freq_ka)
                # else:
                #     self.freqs_over = np.array([])
            # creates a frequency array for all bands
            self.freqs = np.concatenate(
                (self.X['K'][:self.cmax], self.freqs_over, self.X['Ka'][self.cmin:]))
            self.ne = rf.freq2den(self.freqs * 1e9)

    def eval_gd_overlap(self):
        """Overlap group delay from K and Ka bands, if they probe the
         same density layers, with interpolation of K signal into Ka band's
        and Ka signal into K band, averaging with predeterminated weights."""

        if (self.cmin, self.cmax) == (None, None):
            self.gd_over = np.array([])
        else:
            # weight for K and Ka band, for overlapping average
            weights = (2, 1)
            # interpolate the group delay in the overlap region
            func = interpolate.interp1d(self.X['K'], self.gd_k)
            gd_k_over = func(self.freqs_over)
            func = interpolate.interp1d(self.X['Ka'], self.gd_ka)
            gd_ka_over = func(self.freqs_over)
            self.gd_over = (
                gd_k_over * weights[0] + gd_ka_over * weights[1]) / sum(weights)
        self.gd = np.concatenate(
                (self.gd_k[:self.cmax], self.gd_over, self.gd_ka[self.cmin:]))

    def init_gd(self, init_type="quad", steps=6):
        """Group delay initialization:

        init_type: 'line': linear initialization
                    'quad': second order initialization. Default

        steps: number of points to extrapolate between 0 and first probed frequency. Default to 6.
        """

        if not hasattr(self, 'nfreqs'):
            self.freq_ini = np.linspace(0, self.freqs[0], num=steps)
            self.nfreqs = np.concatenate((self.freq_ini, self.freqs[1:]))
            self.freqs = self.nfreqs.copy()
            self.ne = rf.freq2den(self.nfreqs * 1e9)
        if init_type == 'line':
            self.ini_t = np.linspace(0, self.gd[0], num=steps)
        elif init_type == 'quad':
            mean_pos = steps / 2
            mean_t = self.gd[0] * 0.5 * 0.8
            polynomial = np.polyfit((0, self.freq_ini[mean_pos], self.freq_ini[-1]),
                                    (0, mean_t, self.gd[0]), 2)
            pol = np.poly1d(polynomial)
            self.ini_t = pol(self.freq_ini)
        # frequency and group delay now starts from zero
        self.gd = np.concatenate((self.ini_t[:-1], self.gd))

    def eval_phase(self):
        """ From the formula: del phi / del (2 pi* prob_freq) = group_delay
        we can find the phase """
        self.phi = integrate.cumtrapz(self.gd, x=self.freqs, initial=0)*(2*np.pi)

    def prepare_gd(self,sweep, cluster_sweeps=4,all_shot=0):
        """evaluates phase for a given sweep.
        Added just to proccess the group_delay in the right order"""
        
        self.plasma_gd(sweep, cluster_sweeps, all_shot=all_shot)
        self.eval_gd_overlap()
        self.init_gd()

#####################################################

    def save_proc_info(self, sweeps_average, initial_time, last_time):
        p.savetxt(path.join(self.prof_folder, "prof_info.dat"), [
                  self.version, sweeps_average, initial_time, last_time])

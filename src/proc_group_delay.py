"""
Function: Evaluates phase derivative and group delay for profile reconstruction.
Reference is vessel wall, and calibration was done with a pin. After calibration,
signal from both channels are overlapped.
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: give an option to overlap bands in the spectrogram.
test if end of band K is more trustable than beginning of band Ka
Make an 'ini' file.
Make a better info file with parameters used.
"""

from os import path
import pylab as p
import numpy as np

import ref_functions as rf
from proc_sweep import ProcSweep

# CONSTANTS
# 1e-9*c/pi, freq in GHz, gd in ns
const = 95426903.18473884 * 1e-9
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
        self.freq_check = 0

    def average_specgram(self, sweeps=8, sweep_ini=20, all_shot=0, freq_min=(8.4, 26), freq_max=(40, 40)):
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
            matrix_k = self.spectrogram('K', figure=sweep + 1, normal=1, freqs=(4e3, 15e3))
            matrix_ka = self.spectrogram('Ka', figure=sweep + 1, normal=1, freqs=(2e3, 17e3))
            if hasattr(self, 'matrix_k_mean'):
                self.matrix_k_mean += matrix_k[:, np.logical_and(self.X['K'] > freq_min[0], self.X['K'] < freq_max[0])]
                self.matrix_ka_mean += matrix_ka[:, np.logical_and(self.X['Ka'] > freq_min[1], self.X['Ka'] < freq_max[1])]
            # if there are no mean, creates it.
            else:
                self.matrix_k_mean = matrix_k[:, np.logical_and(self.X['K'] > freq_min[0], self.X['K'] < freq_max[0])].copy()
                self.matrix_ka_mean = matrix_ka[:, np.logical_and(self.X['Ka'] > freq_min[1], self.X['Ka'] < freq_max[1])].copy()
                if self.freq_check == 0:
                    self.freq_check = 1
                    # I hope it does not break!
                    self.X['K'] = self.X['K'][np.logical_and(self.X['K'] > freq_min[0], self.X['K'] < freq_max[0])]
                    self.X['Ka'] = self.X['Ka'][np.logical_and(self.X['Ka'] > freq_min[1], self.X['Ka'] < freq_max[1])]

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
        self.average_specgram(sweeps=sw_clustersize, sweep_ini=1, all_shot=all_shot)
        self.gd_k0 = find_max(self.matrix_k_mean, self.Y['K']) - wall_corr
        self.gd_ka0 = find_max(self.matrix_ka_mean, self.Y['Ka']) - wall_corr

    def plasma_gd(self, sweep_cur=1000, sw_clustersize=8, all_shot=0):
        """Evaluate group delay for a cluster of sweeps, and subtracts
        the reference to find the group delay from the limiter"""
        self.average_specgram(sweeps=sw_clustersize, sweep_ini=sweep_cur, all_shot=all_shot)
        self.gd_k_mean = find_max(self.matrix_k_mean, self.Y['K'])
        self.gd_ka_mean = find_max(self.matrix_ka_mean, self.Y['Ka'])
        self.gd_k = self.gd_k_mean - self.gd_k0
        self.gd_ka = self.gd_ka_mean - self.gd_ka0

    def overlap_freq(self):
        # discard first and last points (for instance, 9, and -3)
        limit_min = None
        limit_max = None
        if not hasattr(self, 'freqs_overlap'):
            if min(self.X['Ka']) <= max(self.X['K']):
                # check the number of points of overlaped frquency.
                self.cmin = (abs(self.X['K'] - min(self.X['Ka']))).argmin()
                self.cmax = (abs(self.X['Ka'] - max(self.X['K']))).argmin()
                # create an even spaced array in the overlap region
                xn = np.arange(self.X['K'][self.cmin], self.X['Ka'][self.cmax], self.X['Ka'][1] - self.X['Ka'][0])
                # frequency array for all bands.
                self.freqs_overlap = np.concatenate((self.X['K'][limit_min:self.cmin], xn, self.X['Ka'][self.cmax:limit_max]))
            else:
                # ********************* NAO FUNCIONA *************************************
                self.cmin = None
                self.cmax = None
                self.freqs_overlap = np.concatenate((xk[limit_min:self.cmin], xka[self.cmax:limit_max]))
            self.ne = rf.freq2den(self.freqs_overlap * 1e9)

    def overlap_gd(self):
        """Overlap group delay from K and Ka bands, if they probe the
         same density layers, with interpolation"""
        # discard first and last points (for instance, 9, and -3)
        limit_min = None
        limit_max = None

        # ********************* NAO FUNCIONA *************************************
        # if min(self.X['Ka']) <= max(self.X['K']):
        # interpolate the group delay in the overlap region
        #    gdkn = p.interp(xn, self.X['K'], self.gd_k)
        #    gdkan = p.interp(xn, self.X['K'], self.gd_ka)
        #    gdn = (gdkn + gdkan) / 2
        #    self.gd_m = p.concatenate((self.gd_k[limit_min:self.cmin], gdn, self.gd_ka[self.cmax:limit_max]))
        # else:
        self.gd_m = p.concatenate((self.gd_k[limit_min:self.cmin], self.gd_ka[self.cmax:limit_max]))

    def init_gd(self, tipo="line"):
        """Group delay initialization, default is first order. Second
        order is also available."""
        # choose number of steps to start group_delay
        steps = 6
        if not hasattr(self, 'freqs'):
            self.ini_f = p.linspace(0, self.freqs_overlap[0], num=steps)[:-1]
            self.freqs = p.concatenate((self.ini_f, self.freqs_overlap))
            self.ne = rf.freq2den(self.freqs * 1e9)
        if tipo == 'line':
            self.ini_t = p.linspace(0, self.gd_m[0], num=steps)[:-1]
            polynomial = p.polyfit(self.ini_f, self.ini_t, 1)
        elif tipo == 'quad':
            mean_pos = len(self.ini_f) / 2
            mean_t = self.ini_t[mean_pos] * (1.15)
            polynomial = p.polyfit((self.ini_f[0], self.ini_f[mean_pos], self.ini_f[-1]),
                                   (self.ini_t[0], mean_t, self.ini_t[-1]), 2)
        self.pol = p.poly1d(polynomial)
        self.ini_t = self.pol(self.ini_f)
        # frequency and group delay now starts from zero.
        self.gd = p.concatenate((self.ini_t, self.gd_m))

    def save_proc_info(self, sweeps_average, initial_time, last_time):
        p.savetxt(path.join(self.prof_folder, "prof_info.dat"), [self.version, sweeps_average, initial_time, last_time])

    def plot_average_spectrogram(self, colormap=0):
        """Plot the average spectrogram overlaped of both bands, at the current settings.
        Colormap set to 0 plots only contours, set to 1 plots a colormap."""
        if colormap == 0:
            cont = p.contour
        elif colormap == 1:
            cont = p.contourf
        cont(self.X['K'], self.Y['K'], self.matrix_k_mean)
        cont(self.X['Ka'], self.Y['Ka'], self.matrix_ka_mean)
        p.xlabel("freq (GHz)")
        p.ylabel("group delay (ns)")
        p.title("# %s - time: %s ms" % (self.shot, self.sweep2time(self.sweep_cur)))
        p.ylim(0.4, 13)

    def plot_groupdelay(self):
        p.plot(self.freqs, self.gd)
        p.plot(self.freq_poly * 1e-9, self.gd_poly * 1e9)
        p.xlabel("freq (GHz)")
        p.ylabel("group delay (ns)")
        p.title("# %s - time: %s ms" % (self.shot, self.sweep2time(self.sweep_cur)))


def find_max(matrix, spec_yaxis):
    """Find the curve which follow the maximum for each spectrum."""
    return spec_yaxis[matrix.argmax(axis=0)]

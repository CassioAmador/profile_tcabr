"""
Function: Recreate density profile with th Bottolier-Curtet method
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: ??
"""

import numpy as np
import pylab as p
from scipy import integrate

import ref_functions as rf
from proc_group_delay import ProcGroupDelay

factor = 1 / 2
c = 299792458
# 4*pi/c
const = 41.91690043903363
# c/4*pi
const_inv = 0.023856725796184714


class Bottollier(ProcGroupDelay):

    """Evaluates phase from group delay"""

    def __init__(self, shot, tipo="data", save_locally=1):
        self.version = 1.0
        # inherit from ProcGroupDealy class
        ProcGroupDelay.__init__(self, shot, tipo, save_locally)

    def area_prev(self, n, freqs_prob, pos):
        """Evaluates phase from known part of plasma"""

        # frequency times refraction index
        refrac_index = rf.refrac(freqs_prob[:n + 1], freqs_prob[n])
        area = freqs_prob[n] * integrate.trapz(refrac_index, pos[:n + 1])
        return area

    def pos_eval(self, n, freqs_prob, pos, phase_cur):
        """Evaluates new position"""

        area_dif = factor * freqs_prob[n] * rf.refrac(freqs_prob[n - 1], freqs_prob[n])
        if n == 1:
            phase_prev = 0
        else:
            phase_prev = self.area_prev(n, freqs_prob, pos)
        return pos[n - 1] + (phase_cur - phase_prev + np.pi / 2) / (area_dif)

    def find_pos(self, nX, phase):
        """Iterates over all frequencies to evaluate their reflection distance"""

        pos = np.zeros(len(nX))
        r0 = 0
        for n in range(1, len(nX)):
            pos[n] = self.pos_eval(n, nX, pos, phase[n])
        return pos + r0

    def profile(self, sweep, sweeps=4, all_shot=0):
        """Creates profile and the density array"""

        self.plasma_gd(sweep, sweeps, all_shot=all_shot)
        if not hasattr(self, 'nX'):
            self.nX = {}
            self.overlap_freq()
            for channel in ('K', 'Ka'):
                dif_X = (self.X[channel][1] - self.X[channel][0]) / 2
                self.nX[channel] = self.X[channel][:-1] + dif_X
                if channel == 'K':
                    self.nX[channel] = np.concatenate(([0], self.nX[channel]))
            self.nX['2bands'] = np.concatenate((self.nX['K'][:self.cmin], self.nX['Ka']))
            self.ne = rf.freq2den(self.nX['2bands'] * 1e9)

        phase_dif_k = self.gd_k / self.Dt_DF['K']
        self.phi_k = integrate.cumtrapz(phase_dif_k, initial=0)
        phase_dif_ka = self.gd_ka / self.Dt_DF['Ka']
        self.phi_ka = integrate.cumtrapz(phase_dif_ka) + self.phi_k[self.cmin]
        self.overlap_phase()
        self.pos = const_inv * self.find_pos(self.nX['2bands'], self.phase)

    def overlap_phase(self):
        """If both bands overlap, interpolates the K band phase into Ka band's
        frequency resolution, which is smaller. Then it averages both bands
        signal at the overlap."""

        if self.cmin is not None:
            phase_interp = np.interp(self.nX['Ka'][:self.cmax], self.nX['K'][self.cmin:], self.phi_k[self.cmin:])
            phase_over = (phase_interp + self.phi_ka[:self.cmax]) / 2
            self.phase = np.concatenate((self.phi_k[:self.cmin], phase_over, self.phi_ka[self.cmax:]))
        else:
            self.phase = np.concatenate((self.phi_k, self.phi_ka))

    def plot_phase(self):
        # p.plot(self.nX['Ka'][:self.cmax],npp, marker = 'o', linestyle = '-', color = 'r')
        p.plot(self.nX['K'], self.phi_k, marker='o', linestyle='-', color='b')
        p.plot(self.nX['Ka'], self.phi_ka, marker='o', linestyle='-', color='k')
        p.plot(self.nX['2bands'], self.phase, color='c')

if __name__ == "__main__":
    shot = Bottollier(32111)
    shot.reference_gd()
    shot.profile(6666, 2)

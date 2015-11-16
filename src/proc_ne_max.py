ARRUMAR
"""
Function: Evaluates the maximum density of the plasma, if it is lower than the maximum
probing frequency. Try to find a jump in the group delay
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: choose best one and erase the other
"""

from os import path, getcwd, makedirs
import pylab as p
import numpy as np
import scipy.signal

from proc_group_delay import ProcGroupDelay

    def find_ne_max(self):
        if not hasattr(self, "ne_max"):
            self.ne_max = []
        zeros_gd = p.where(self.gd_m <= 0.0)[0]
#       print(max(p.diff(self.gd_m)))
        if len(zeros_gd) != 0:
            print("pequeno")
            self.ne_max.append(rf.freq2den(1e9 * self.X[zeros_gd[0]]))
            self.gd_m[zeros_gd[0]:] = float("NaN")
        else:
            dif_gd = p.where(p.diff(self.gd_m) > 3.5)[0]
            if len(dif_gd) != 0:
                print("grande")
                self.ne_max.append(rf.freq2den(1e9 * self.X[dif_gd[-1]]))
                self.gd_m[dif_gd[0]:] = float("NaN")
            else:
                self.ne_max.append(float("NaN"))
#       zero_jumps = p.where(p.diff(zero_gd) > 1)[0]

    def find_ne_max2(self):
        """Find maximum density."""
        self.freqs_full = np.sort(np.concatenate((self.X_k, self.X_ka)))
        self.ne_full = rf.freq2den(1e9 * self.freqs_full)
        index = np.where(np.logical_or(self.gd_k < 0, self.gd_k > 0.95 * wall_corr))[0]
        if len(index) > 0:
            self.ne_max = rf.freq2den(1e9 * self.X_k[index[0]])
            # No plasma?
            if self.ne_max < 0.45e19:
                self.no_plasma = True
                self.ne_max = np.nan
                self.gd2 = []
                self.freqs2 = []
            else:
                self.no_plasma = False
                self.gd2 = self.gd_k[:index[0]]
                self.freqs2 = self.X_k[:index[0]]
        else:
            index = np.where(np.logical_or(self.gd_ka < 0, self.gd_ka > 0.95 * wall_corr))[0]
            if len(index) > 0:
                self.ne_max = rf.freq2den(1e9 * self.X_ka[index[0]])
                if index[0] == 0:
                    self.gd2 = self.gd_k[:index[0]]
                    self.freqs2 = self.X_k[:index[0]]
                    return
            else:
                # plasma density > max probing density
                index = [-1]
                self.ne_max = np.nan
            freqs = np.concatenate((self.X_k, self.X_ka[:index[0]]))
            sort_index = np.argsort(freqs)
            self.gd2 = np.concatenate((self.gd_k, self.gd_ka[:index[0]]))[sort_index]
            self.freqs2 = freqs[sort_index]
            self.no_plasma = False
        return

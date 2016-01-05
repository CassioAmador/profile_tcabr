"""
Function: Recreate density profile from group delay, with Abel Inversion.
Treats each group delay point as a first order function.
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: Incorporate other abel inversions.
"""

import numpy as np
from scipy import special as sp_sp

from proc_group_delay import ProcGroupDelay

# 1e-9*c/pi, freq in GHz, gd in ns
const = 0.09542690318473884


def eval_pos(freq_cut, freq_prev, gd_prev, incli):
    """Evaluates new position"""

    # polynomial coeffs for a first order function, for the abel integral inversion starting from non-zero frequency.
    abel0 = np.arccos(freq_prev / freq_cut)
    abel1 = np.sqrt(freq_cut**2 - freq_prev**2)
    return abel0 * gd_prev + (abel1 - abel0 * freq_prev) * incli


def find_pos(freqs_prob, gd):
    """Iterates over all frequencies to evaluate their reflection distance"""

    pos = np.zeros(len(freqs_prob))
    dif_gd_freq = np.diff(gd) / np.diff(freqs_prob)
    for n in range(1, len(freqs_prob)):
        new_pos = eval_pos(freqs_prob[n], freqs_prob[n - 1], gd[n - 1], dif_gd_freq[n - 1])
        pos[n] = pos[n - 1] + new_pos
    # group delay should be dependente of omega, not frequency, so we divide by 2pi
    return pos * const / (2 * np.pi)


def find_pos_pol(f_probe, tau, order=5, init=2):
    """Abel inversion for the polynomial fit of the group delay."""
    rc = np.zeros(len(tau))
    tau2_coef = np.polyfit(f_probe, tau, order)
    # inicialização
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
        rc[i] = 1e-9 * (3e8 / np.pi) * (np.dot(tau1_coef[::-1], I1) + np.dot(tau2_coef[::-1], I3 - I2))
    return rc


class ProcProfile(ProcGroupDelay):

    """Recreate density profile with the Abel Inversion method,
    from group delay evaluated in ProcGroupDelay"""

    def __init__(self, shot, tipo="data", save_locally=1):
        self.version = 1.0
        # inherit from ProcGroupDealy class
        ProcGroupDelay.__init__(self, shot, tipo, save_locally)

    def profile(self):
        """Creates profile array"""
        self.pos = find_pos(self.freqs, self.gd)


if __name__ == "__main__":
    shot = ProcProfile(32111)
    shot.reference_gd()
    shot.eval_freq_overlap()
    shot.prepare_gd(6666, 2)
    shot.profile()
    import matplotlib.pyplot as plt
    plt.plot(shot.pos * 1e2, shot.ne)
    print('\nplasma position (cm):')
    print(shot.pos * 1e2)

"""
Function: Recreate density profile with th Bottolier-Curtet method
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: Maybe the processing will speed up if the refraction index is evaluated only once per frequency,
      and all results stored in a matrix for recurrent use in each profile.
"""

import numpy as np
from scipy import integrate

import ref_functions as rf
from proc_group_delay import ProcGroupDelay

# 4*pi/c
const = 41.91690043903363
# c/4*pi
const_inv = 0.023856725796184714

def area_prev(n, freqs_prob, pos):
    """Evaluates phase from known part of plasma"""

    # frequency times refraction index
    refrac_index = rf.refrac(freqs_prob[:n + 1], freqs_prob[n])
    area = const * freqs_prob[n] * \
        integrate.trapz(refrac_index, pos[:n + 1])
    return area

def eval_pos(n, freqs_prob, pos, phase_cur,factor=2/3):
    """Evaluates new position"""

    area_dif = factor * freqs_prob[n] * \
        rf.refrac(freqs_prob[n - 1], freqs_prob[n])
    if n == 1:
        phase_prev = 0
    else:
        phase_prev = area_prev(n, freqs_prob, pos)
    return pos[n - 1] + const_inv * (phase_cur - phase_prev + np.pi / 2) / (area_dif)

def find_pos(freqs_prob, phase,factor=2/3):
    """Iterates over all frequencies to evaluate their reflection distance"""

    pos = np.zeros(len(freqs_prob))
    r0 = 0
    for n in range(1, len(freqs_prob)):
        pos[n] = eval_pos(n, freqs_prob, pos, phase[n],factor=factor)
    return pos + r0


class ProcProfile(ProcGroupDelay):

    """Recreate density profile with the Bottolier-Curtet method,
    from phase evaluated in ProcGroupDelay"""

    def __init__(self, shot, tipo="data", save_locally=1,factor=2/3):
        self.version = 1.0
        # inherit from ProcGroupDealy class
        self.factor=factor
        ProcGroupDelay.__init__(self, shot, tipo, save_locally)

    def profile(self):
        """Creates profile array"""

        self.eval_phase()
        self.pos = find_pos(self.freqs[1:], self.phi,factor=self.factor)

if __name__ == "__main__":
    shot = ProcProfile(32111)
    shot.reference_gd()
    shot.eval_freq_overlap()
    shot.prepare_gd(6666,2)
    shot.profile()
    print('\nplasma position (cm):')
    print(shot.pos*1e2)
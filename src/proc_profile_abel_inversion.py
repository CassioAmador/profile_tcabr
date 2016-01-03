"""
Function: Recreate density profile from group delay, with Abel Inversion.
Treats each group delay point as a first order function.
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: Incorporate other abel inversions.
"""

import numpy as np
from scipy import integrate
from scipy import special as sp_sp

import ref_functions as rf
from proc_group_delay import ProcGroupDelay

# 1e-9*c/pi, freq in GHz, gd in ns
const = 0.09542690318473884

def eval_pos(freq_cut, freq_prev,gd_prev,incli):
    """Evaluates new position"""

    # polynomial coeffs for a first order function, for the abel integral inversion starting from non-zero frequency.
    abel0 = np.arccos(freq_prev/freq_cut)
    abel1 = np.sqrt(freq_cut**2 - freq_prev**2)
    return abel0*gd_prev+(abel1-abel0*freq_prev)*incli

def find_pos(freqs_prob, gd):
    """Iterates over all frequencies to evaluate their reflection distance"""

    pos = np.zeros(len(freqs_prob))
    dif_gd_freq=np.diff(gd)/np.diff(freqs_prob)
    for n in range(1,len(freqs_prob)):
        new_pos = eval_pos(freqs_prob[n], freqs_prob[n-1], gd[n-1], dif_gd_freq[n-1])
        pos[n]=pos[n-1]+new_pos
    # group delay should be dependente of omega, not frequency, so we divide by 2pi
    return pos*const/(2*np.pi)

class ProcProfile(ProcGroupDelay):

    """Recreate density profile with the Abel Inversion method,
    from group delay evaluated in ProcGroupDelay"""

    def __init__(self, shot, tipo="data", save_locally=1):
        self.version = 1.0
        # inherit from ProcGroupDealy class
        ProcGroupDelay.__init__(self, shot, tipo, save_locally)

    def profile(self):
        """Creates profile array"""
        self.pos = find_pos(self.freqs,self.gd)




if __name__ == "__main__":
    shot = ProcProfile(32111)
    shot.reference_gd()
    shot.eval_freq_overlap()
    shot.prepare_gd(6666,2)
    shot.profile()
    import matplotlib.pyplot as plt
    plt.plot(shot.pos*1e2,shot.ne)
    print('\nplasma position (cm):')
    print(shot.pos*1e2)
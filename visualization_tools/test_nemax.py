# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, './../src/')

import proc_profile_abel_inversion as pp


def main(argv):
    """Test maximum density."""
    if len(argv) < 1:
        shot = 33772
    else:
        shot = int(argv[1])
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    shot = pp.ProcProfile(shot)
    sweeps_average = 33

    initial_sweep = shot.time2sweep(32)
    last_sweep = shot.time2sweep(140)
    shot.reference_gd(all_shot=1, sw_clustersize=sweeps_average)
    sweeps = np.arange(initial_sweep, last_sweep, sweeps_average)
    ne_max = np.zeros(len(sweeps))
    t = np.zeros(len(sweeps))

    for i in range(len(sweeps)):
        print("Sweep: %d" % sweeps[i])
        shot.plasma_gd(sweeps[i], sweeps_average, all_shot=1)
        shot.find_ne_max2()
        ne_max[i] = shot.ne_max
        t[i] = shot.sweep2time(sweeps[i])

    plt.plot(t, ne_max * 1e-19)
    plt.ylabel("ne [$10^{19}$ m$^-3$] ")
    plt.xlabel("time (ms)")
    plt.title("# %s" % (shot.shot))

    plt.show()

if __name__ == '__main__':
    sys.exit(main(sys.argv))

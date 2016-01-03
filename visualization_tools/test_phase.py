"""
Test phase and group delay for both bands.

TODO: Show overlap with different colors
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np
import scipy
import sys
sys.path.insert(0, './../src/')

import proc_profile_bottollier as ppb

if len(sys.argv) < 1:
    shot_number = int(open('shot_number.txt', 'r').read())
else:
    if (len(sys.argv) == 1) & ("py" in sys.argv[0]):
        shot_number = int(open('shot_number.txt', 'r').read())
    else:
        shot_number = int(sys.argv[1])

shot = ppb.ProcProfile(shot_number)
sweeps_average = 8

shot.reference_gd(all_shot=1, sw_clustersize=sweeps_average)
shot.eval_freq_overlap()
shot.plasma_gd(5000, sweeps_average, 1)
shot.eval_gd_overlap()
shot.init_gd()
shot.eval_phase()

fig = plt.figure()
fig.subplots_adjust(bottom=0.25)
ax1 = fig.add_subplot(111)
lns1 = ax1.plot(shot.freqs, shot.gd, marker='.',
                linestyle='-', color='k', label="gd")
lns3 = ax1.plot(shot.freqs_over, shot.gd_over, marker='o',
                linestyle='-', color='r', label="gd_overlap")
lns4 = ax1.plot(shot.freq_ini, shot.ini_t, marker='o',
                linestyle='-', color='g', label="gd_init")
plt.xlabel("probing freq (GHz)")
plt.ylabel("group delay (ns)")
plt.title("# %s - time: %s ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
ax1.set_ylim(0, 2)
ax2 = ax1.twinx()
lns2 = ax2.plot(shot.freqs, shot.phi, marker='.',
                linestyle='-', color='b', label="phase")
plt.ylabel("phase (rad)")
plt.xlim(0, 40)

# point_over=np.where(shot.freqs==shot.freqs_over[0])[0][0]
#lns4= ax1.axvline(shot.freqs[point_over],color='r',label="band overlap")
# lns3=ax1.axhline(shot.gd[point_over],color='r')


# added these three lines
lns = lns1 + lns2 + lns3+lns4
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=2)

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.13, 0.1, 0.77, 0.03], axisbg=axcolor)

sweep = Slider(axfreq, 'Sweep', 1, len(shot.points) - 1 -
               sweeps_average, valinit=5000, valfmt='%1.f')


def update(val):
    shot.plasma_gd(int(sweep.val), sweeps_average, all_shot=1)
    shot.eval_gd_overlap()
    shot.init_gd()
    shot.eval_phase()
    lns1[0].set_ydata(shot.gd)
    lns2[0].set_ydata(shot.phi)
    lns3[0].set_ydata(shot.gd_over)
    lns4[0].set_ydata(shot.ini_t)
    # lns3.set_ydata(shot.gd[point_over])
    ax2.set_ylim(min(shot.phi), max(shot.phi))
    ax1.set_title("# %s - time: %.3f ms" %
                  (shot.shot, shot.sweep2time(shot.sweep_cur)))
    fig.canvas.draw_idle()

sweep.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sweep.reset()
button.on_clicked(reset)

# plt.show()

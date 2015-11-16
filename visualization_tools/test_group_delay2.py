"""
Test group delay (raw) from each band.

Compare it with median filtered, smoothed and fitted group delay.
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np
import scipy
import sys
sys.path.insert(0, './../src/')

import proc_profile_abel_inversion as pp

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)
shot = pp.ProcProfile(33708)
sweeps_average = 33

shot.reference_gd(all_shot=1, sw_clustersize=sweeps_average)
shot.plasma_gd(5000, sweeps_average, 1)
shot.find_ne_max2()

tau_coef = np.polyfit(shot.freqs2, shot.gd2, 3)
ax1, = plt.plot(shot.freqs2, shot.gd2, 'ok', label="gd (raw)")
ax2, = plt.plot(shot.freqs2, scipy.signal.medfilt(shot.gd2, 5), 'ob', label="gd (median filt-5)")
ax3, = plt.plot(shot.freqs2, shot.smooth_signal(shot.gd2), '--c', label="gd (smothed)")
ax4, = plt.plot(shot.freqs2, np.polyval(tau_coef, shot.freqs2), '--r', label="gd (fitted)")
plt.xlabel("freq (GHz)")
plt.ylabel("group delay (ns)")
plt.title("# %s - time: %s ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
plt.legend(loc=2)
plt.xlim(15, 40)
plt.ylim(-0.1, 3)

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.13, 0.1, 0.77, 0.03], axisbg=axcolor)

sweep = Slider(axfreq, 'Sweep', 1, len(shot.points) - 1 - sweeps_average, valinit=1, valfmt='%1.f')


def update(val):
    shot.plasma_gd(int(sweep.val), sweeps_average, 1)
    shot.find_ne_max2()
    tau_coef = np.polyfit(shot.freqs2, shot.gd2, 3)
    ax1.set_xdata(shot.freqs2)
    ax1.set_ydata(shot.gd2)
    ax2.set_xdata(shot.freqs2)
    ax2.set_ydata(scipy.signal.medfilt(shot.gd2, 5))
    ax3.set_xdata(shot.freqs2)
    ax3.set_ydata(shot.smooth_signal(shot.gd2))
    ax4.set_xdata(shot.freqs2)
    ax4.set_ydata(np.polyval(tau_coef, shot.freqs2))
    ax.set_title("# %s - time: %.3f ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
    fig.canvas.draw_idle()
sweep.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sweep.reset()
button.on_clicked(reset)

plt.show()

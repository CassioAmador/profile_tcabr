import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np
import sys
sys.path.insert(0, './../src/')

import proc_profile_abel_inversion as pp
import ref_functions as rf

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)
shot = pp.ProcProfile(33708)

shot.reference_gd(all_shot=1, sw_clustersize=1)
cluster = 50
shot.plasma_gd(7000, cluster, 1)
shot.find_ne_max2()
ne = shot.ne_full*1e-19

r = np.nan*np.ones(len(ne))
if not shot.no_plasma:
    r_inver = shot.abel_transform(shot.gd2*1e-9, shot.freqs2*1e9, order=2, init=2)
    r[:len(r_inver)] = r_inver
    tau0 = np.polyval(np.polyfit(shot.freqs2, shot.gd2, 2), shot.freqs2[0])
    f_init = np.linspace(1, shot.freqs2[0], 10)
    tau_init = f_init * tau0 / (f_init[-1]**2)
    ne_init = rf.freq2den(f_init*1e9)*1e-19
    r_init = shot.abel_transform(tau_init*1e-9, f_init*1e9, order=2, init=2)

l, = plt.plot(r, ne)
m, = plt.plot(r_init, ne_init, 'r--')
plt.xlabel("r [m]")
plt.ylabel("ne [$10^{19}$ m$^-3$]")
plt.title("# %s - time: %s ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
plt.ylim(0, 2)
plt.xlim(-0.05, 0.2)

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

sweep = Slider(axfreq, 'Sweep', 1, len(shot.points) - 1 - cluster, valinit=1, valfmt='%1.f')


def update(val):
    shot.plasma_gd(int(sweep.val), cluster, 1)
    shot.find_ne_max2()
    r = np.nan*np.ones(len(ne))
    if not shot.no_plasma:
        r_inver = shot.abel_transform(shot.gd2*1e-9, shot.freqs2*1e9, order=3, init=1)
        r[:len(r_inver)] = r_inver
        tau0 = np.polyval(np.polyfit(shot.freqs2, shot.gd2, 2), shot.freqs2[0])
        f_init = np.linspace(1, shot.freqs2[0], 10)
        tau_init = f_init * tau0 / (f_init[-1]**2)
        r_init = shot.abel_transform(tau_init*1e-9, f_init*1e9, order=2, init=2)
    l.set_xdata(r)
    l.set_ydata(ne)
    m.set_xdata(r_init)
    m.set_ydata(ne_init)
    ax.set_title("# %s - time: %.3f ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
    fig.canvas.draw_idle()
sweep.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sweep.reset()
button.on_clicked(reset)

plt.show()

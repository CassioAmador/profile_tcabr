import sys
sys.path.insert(0, './../src/')

import proc_profile as pp
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)
shot = pp.ProcProfile(30407)

shot.reference_gd(all_shot=1, sw_clustersize=1)
cluster = 10
shot.plasma_gd(5000, cluster, 1)
shot.find_ne_max2()

l, = plt.plot(shot.freqs2 , shot.gd2, 'b')
plt.xlabel("freq (GHz)")
plt.ylabel("group delay (ns)")
plt.title("# %s - time: %s ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
plt.xlim(15, 40)
#plt.ylim(-1, 5)

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

sweep = Slider(axfreq, 'Sweep', 1, len(shot.points) - 1 - cluster, valinit=1, valfmt='%1.f')


def update(val):
    shot.plasma_gd(int(sweep.val), cluster, 1)
    shot.find_ne_max2()
    l.set_ydata(shot.gd2)
    l.set_xdata(shot.freqs2)
    ax.set_title("# %s - time: %.3f ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
    fig.canvas.draw_idle()
sweep.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sweep.reset()
button.on_clicked(reset)

plt.show()

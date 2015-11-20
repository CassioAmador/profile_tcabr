import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import sys
sys.path.insert(0, './../src/')

import proc_profile_bottollier as pp

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)
cluster = 25
shot = pp.Bottollier(33777, save_locally=1)
shot.reference_gd(all_shot=1, sw_clustersize=cluster)
shot.profile(3000, sweeps=cluster, all_shot=1)

l, = plt.plot(shot.pos, shot.ne * 1e2)
ax.legend(loc='upper left')
ax.set_xlabel('r (cm)')
ax.set_ylabel('density (10^19 m^-3)')
plt.title("# %s - time: %s ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
plt.axis([0, 20, shot.ne[1], shot.ne[-1] * 1.1])

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

sweep = Slider(axfreq, 'Sweep', 1, len(shot.points) - 1 - cluster, valinit=3000, valfmt='%1.f')


def update(val):
    shot.profile(int(sweep.val), sweeps=cluster, all_shot=1)
    l.set_xdata(shot.pos * 1e2)
    l.set_ydata(shot.ne)
    ax.set_title("# %s - time: %.3f ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
    fig.canvas.draw_idle()
sweep.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sweep.reset()
button.on_clicked(reset)

plt.show()

"""Test group delay (raw) from each band."""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import sys
sys.path.insert(0, './../src/')

import proc_profile_bottollier as ppb

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)
if len(sys.argv) < 1:
    shot_number=int(open('shot_number.txt','r').read())
else:
    if (len(sys.argv) == 1) & ("py" in sys.argv[0]):
        shot_number=int(open('shot_number.txt','r').read())
    else:    
        shot_number = int(sys.argv[1])
shot = ppb.ProcProfile(shot_number)

shot.reference_gd(all_shot=1, sw_clustersize=33)
cluster = 33
shot.plasma_gd(1, cluster, 1)

ax1, = plt.plot(shot.X['K'], shot.gd_k, 'ob')
ax2, = plt.plot(shot.X['Ka'], shot.gd_ka, 'or')
plt.xlabel("freq (GHz)")
plt.ylabel("group delay (ns)")
plt.title("# %s - time: %s ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
plt.ylim(-1, 5)
plt.xlim(15, 41)

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

sweep = Slider(axfreq, 'Sweep', 1, len(shot.points) - 1 - cluster, valinit=1, valfmt='%1.f')


def update(val):
    shot.plasma_gd(int(sweep.val), cluster, 1)
    ax1.set_ydata(shot.gd_k)
    ax2.set_ydata(shot.gd_ka)
    ax.set_title("# %s - time: %.3f ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
    fig.canvas.draw_idle()
sweep.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sweep.reset()
button.on_clicked(reset)

plt.show()

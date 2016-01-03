import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import sys
sys.path.insert(0, './../src/')

import proc_profile_bottollier as ppb

if len(sys.argv) < 1:
    shot_number=int(open('shot_number.txt','r').read())
else:
    if (len(sys.argv) == 1) & ("py" in sys.argv[0]):
        shot_number=int(open('shot_number.txt','r').read())
    else:    
        shot_number = int(sys.argv[1])
shot = ppb.ProcProfile(shot_number)

shot.reference_gd(all_shot=1)
cluster = 20
shot.plasma_gd(1, cluster, 1)

fig, ax = plt.subplots(1, 2)
plt.subplots_adjust(bottom=0.25)
ax[0].pcolormesh(shot.X['K'], shot.Y['K'], shot.matrix_k_mean)
ax[1].pcolormesh(shot.X['Ka'], shot.Y['Ka'], shot.matrix_ka_mean)

l, = ax[0].plot(shot.X['K'], shot.Y['K'][shot.matrix_k_mean.argmax(axis=0)], color='r', linewidth=2.0)
ax[0].set_ylim(0, 12)
ax[0].set_xlim(shot.X['K'].min(), shot.X['K'].max())
ax[0].set_ylabel("group delay (ns)")
ax[0].set_xlabel("freq (GHz)")
m, = ax[1].plot(shot.X['Ka'], shot.Y['Ka'][shot.matrix_ka_mean.argmax(axis=0)], color='r', linewidth=2.0)
ax[1].set_xlabel("freq (GHz)")
ax[1].set_ylim(0, 12)
ax[1].set_xlim(shot.X['Ka'].min(), shot.X['Ka'].max())
title = fig.suptitle("# %s - time: %s ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

sweep = Slider(axfreq, 'Sweep', 1, len(shot.points) - 1 - cluster, valinit=1, valfmt='%1.f')


def update(val):
    shot.plasma_gd(int(sweep.val), cluster, 1)
    ax[0].pcolormesh(shot.X['K'], shot.Y['K'], shot.matrix_k_mean)
    ax[1].pcolormesh(shot.X['Ka'], shot.Y['Ka'], shot.matrix_ka_mean)
    l.set_ydata(shot.Y['K'][shot.matrix_k_mean.argmax(axis=0)])
    m.set_ydata(shot.Y['Ka'][shot.matrix_ka_mean.argmax(axis=0)])
    title.set_text("# %s - time: %.3f ms" % (shot.shot, shot.sweep2time(shot.sweep_cur)))
    fig.canvas.draw_idle()
sweep.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sweep.reset()
button.on_clicked(reset)

plt.show()

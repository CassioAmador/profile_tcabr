"""Test custom spectrogram used in sweep reflectometry
Plot the average spectrogram overlaped of both bands, at the current settings.
        Colormap set to 0 plots only contours, set to 1 plots a colormap.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import sys
sys.path.insert(0, './../src/')

import proc_group_delay as pgd

if len(sys.argv) < 1:
    shot_number=int(open('shot_number.txt','r').read())
else:
    if (len(sys.argv) == 1) & ("py" in sys.argv[0]):
        shot_number=int(open('shot_number.txt','r').read())
    else:    
        shot_number = int(sys.argv[1])

shot=pgd.ProcGroupDelay(shot_number)
cluster_sweep=4
sweep=6666
zer_pad=8
window=256
step_scale=16
shot.average_specgram(sweeps=cluster_sweep, sweep_ini=sweep, all_shot=1,window=window,step_scale=step_scale, zer_pad=zer_pad)
gd_k = shot.Y['K'][shot.matrix_k_mean.argmax(axis=1)]
gd_ka = shot.Y['Ka'][shot.matrix_ka_mean.argmax(axis=1)]

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)
#ax.legend(loc='upper left')
ax.set_title("# %s - time: %.3f ms - sweeps averaged: %s \n window: %s - step: %s - zero_pad: %s" %
                 (shot.shot, shot.sweep2time(shot.sweep_cur),cluster_sweep,window,window/step_scale, zer_pad))
#plt.axis([0, 20, shot.ne[1], shot.ne[-1] * 1.1])
plt.xlabel("probing freq (GHz)")
plt.ylabel("group delay (ns)")

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
sweep = Slider(axfreq, 'Sweep', 1, len(shot.points) -
               1 - cluster_sweep, valinit=sweep, valfmt='%1.f')

colormap=1
if colormap == 0:
    cont = ax.contour
elif colormap == 1:
    cont = ax.contourf
cont(shot.X['K'], shot.Y['K'], shot.matrix_k_mean.transpose())
cont(shot.X['Ka'], shot.Y['Ka'], shot.matrix_ka_mean.transpose())
lns3=ax.plot(shot.X['K'],gd_k,'k',marker='.')
lns4=ax.plot(shot.X['Ka'],gd_ka,'k',marker='.')
plt.ylim(0.4, 12)

def update(val):
    shot.average_specgram(sweeps=cluster_sweep, sweep_ini=int(sweep.val), all_shot=1,window=window,step_scale=step_scale, zer_pad=zer_pad)
    gd_k = shot.Y['K'][shot.matrix_k_mean.argmax(axis=1)]
    gd_ka = shot.Y['Ka'][shot.matrix_ka_mean.argmax(axis=1)]
    cont(shot.X['K'], shot.Y['K'], shot.matrix_k_mean.transpose())
    cont(shot.X['Ka'], shot.Y['Ka'], shot.matrix_ka_mean.transpose())
    lns3[0].set_ydata(gd_k)
    lns4[0].set_ydata(gd_ka)
    ax.set_title("# %s - time: %.3f ms - sweeps averaged: %s \n window: %s - step: %s - zero_pad: %s" %
                 (shot.shot, shot.sweep2time(shot.sweep_cur),cluster_sweep,window,window/step_scale, zer_pad))
    fig.canvas.draw_idle()

sweep.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    sweep.reset()
button.on_clicked(reset)

plt.show()

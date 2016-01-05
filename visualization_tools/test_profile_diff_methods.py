"""Show profiles from shots recreated with the Bottollier method"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import sys
if './../src/' not in sys.path:
    sys.path.insert(0, './../src/')

import proc_profile_abel_inversion as ppa
import proc_profile_bottollier as ppb
import ref_functions as rf

if len(sys.argv) < 1:
    shot_number = int(open('shot_number.txt', 'r').read())
else:
    if (len(sys.argv) == 1) & ("py" in sys.argv[0]):
        shot_number = int(open('shot_number.txt', 'r').read())
    else:
        shot_number = int(sys.argv[1])

cluster = 8
shot = ppb.ProcProfile(shot_number, save_locally=1)
sweep_ini = shot.time2sweep(100)
shot.reference_gd(all_shot=1, sw_clustersize=cluster)
shot.eval_freq_overlap()
shot.prepare_gd(sweep_ini, cluster, all_shot=1)
shot.profile()
pos_abel = ppa.find_pos(shot.freqs, shot.gd)

fig = plt.figure()
ax1 = fig.add_subplot(212)
plt.subplots_adjust(bottom=0.25)
lineb, = ax1.plot(shot.pos * 1e2, shot.ne[1:], label="Bottollier")
linea, = ax1.plot(pos_abel * 1e2, shot.ne, label="Abel Inversion")
ax1.legend(loc='upper left')
ax1.set_xlabel('r (cm)')
ax1.set_ylabel('density (10^19 m^-3)')
plt.axis([0, 20, shot.ne[1], shot.ne[-1] * 1.1])
ax1.axhline(rf.freq2den(shot.freq_ini[-1] * 1e9), color='g', linestyle='-.')
if len(shot.freqs_over) != 0:
    ax1.axhline(rf.freq2den(shot.freqs_over[
                0] * 1e9), color='k', linestyle='-.')
    ax1.axhline(rf.freq2den(
        shot.freqs_over[-1] * 1e9), color='k', linestyle='-.')
else:
    ax1.axhline(rf.freq2den(shot.X['K'][-1] * 1e9), color='r', linestyle='-.')
    ax1.axhline(rf.freq2den(shot.X['Ka'][0] * 1e9), color='r', linestyle='-.')

ax2 = fig.add_subplot(211)
lined, = ax2.plot(abs(pos_abel[:-1] - shot.pos) * 1e2, 'k')
ax2.set_ylabel("abs diff (cm)")
ax2.set_title("# %s - time: %s ms - #sweeps average: %s " %
              (shot_number, shot.sweep2time(shot.sweep_cur), cluster))


axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

sweep = Slider(axfreq, 'Sweep', 1, len(shot.points) -
               1 - cluster, valinit=3000, valfmt='%1.f')


def update(val):
    shot.prepare_gd(int(sweep.val), cluster, all_shot=1)
    shot.profile()
    pos_abel = ppa.find_pos(shot.freqs, shot.gd)
    lineb.set_xdata(shot.pos * 1e2)
    linea.set_xdata(pos_abel * 1e2)
    lined.set_ydata(abs(pos_abel[:-1] - shot.pos) * 1e2)
ax2.set_title("# %s - time: %s ms - #sweeps average: %s " % (shot_number, shot.sweep2time(shot.sweep_cur), cluster))    fig.canvas.draw_idle()
sweep.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sweep.reset()
button.on_clicked(reset)

plt.show()

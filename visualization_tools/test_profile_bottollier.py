"""Show profiles from shots recreated with the Bottollier method"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import sys
if './../src/' not in sys.path:
    sys.path.insert(0, './../src/')

import proc_profile_bottollier as ppb
import ref_functions as rf

if len(sys.argv) < 1:
    shot_number=int(open('shot_number.txt','r').read())
else:
    if (len(sys.argv) == 1) & ("py" in sys.argv[0]):
        shot_number=int(open('shot_number.txt','r').read())
    else:    
        shot_number = int(sys.argv[1])

factors={'1':'1','2/3':'2/3','1/2':'1/2'}
pos={}
cluster = 20
shot= ppb.ProcProfile(shot_number, save_locally=1)
sweep_ini = shot.time2sweep(100)
shot.reference_gd(all_shot=1, sw_clustersize=cluster)
shot.eval_freq_overlap()
shot.prepare_gd(sweep_ini, cluster, all_shot=1)
shot.eval_phase()
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)
lines={}

for f in factors.keys():
    pos[f]=ppb.find_pos(shot.freqs[1:], shot.phi,factor=eval(f))
    lines[f], = plt.plot(pos[f] * 1e2, shot.ne[1:],label=factors[f])

ax.legend(loc='upper left',title='Factor:')
ax.set_xlabel('r (cm)')
ax.set_ylabel('density (10^19 m^-3)')
plt.title("# %s - time: %s ms" % (shot_number, shot.sweep2time(shot.sweep_cur)))
plt.axis([0, 20, shot.ne[1], shot.ne[-1] * 1.1])
plt.axhline(rf.freq2den(shot.freq_ini[-1]*1e9),color='g',linestyle='-.')
if len(shot.freqs_over)!=0:
    plt.axhline(rf.freq2den(shot.freqs_over[0]*1e9),color='k',linestyle='-.')
    plt.axhline(rf.freq2den(shot.freqs_over[-1]*1e9),color='k',linestyle='-.')
else:
    plt.axhline(rf.freq2den(shot.X['K'][-1]*1e9),color='r',linestyle='-.')
    plt.axhline(rf.freq2den(shot.X['Ka'][0]*1e9),color='r',linestyle='-.')

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

sweep = Slider(axfreq, 'Sweep', 1, len(shot.points) -
               1 - cluster, valinit=3000, valfmt='%1.f')


def update(val):
    shot.prepare_gd(int(sweep.val), cluster, all_shot=1)
    shot.eval_phase()
    for f in factors.keys():
        pos[f] = ppb.find_pos(shot.freqs[1:], shot.phi,factor=eval(f))
        lines[f].set_xdata(pos[f] * 1e2)

    # l.set_ydata(shot.ne_phi)
    ax.set_title("# %s - time: %.3f ms" %
                 (shot.shot, shot.sweep2time(shot.sweep_cur)))
    fig.canvas.draw_idle()
sweep.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sweep.reset()
button.on_clicked(reset)

plt.show()

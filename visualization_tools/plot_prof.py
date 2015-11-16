"""
Given a shot_number, it reads all files in 'prof_folder' to plot it.
With a slider, it is possible to change time of profile with a mouse.
author: Cassio Amador
TODO: read from 'ini' file.
???
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from os import chdir, listdir, path, getcwd


shot_number = 33708
prof_folder = path.join(getcwd(), "..", "PROC_FILES", "%s" % shot_number)
chdir(prof_folder)
ne = np.load('ne.npy')
info = np.loadtxt('prof_info.dat')

prof_list = listdir(prof_folder)
prof_list.sort()
position = np.empty(shape=((len(prof_list)), len(ne)))
times = np.empty(len(prof_list))

i = 0
for r_file in prof_list:
    name = r_file.strip('.npy')
    if name.isdigit():
        position[i] = np.load(r_file) * 1e2
        times[i] = name
        i += 1

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
l, = plt.plot(position[0], ne, lw=2, color='blue')
plt.axis([0, 20, ne[1], ne[-1] * 1.1])
ax.legend(loc='upper left')
ax.set_xlabel('r (cm)')
ax.set_ylabel('density (10^19 m^-3)')
ax.set_title("# %s" % shot_number)

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

stime = Slider(axfreq, 'time', info[1], info[2], valinit=info[1] + (info[2] - info[1]) * 0.2)


def update(val):
    time = stime.val
    i = (abs(time * 1e3 - times)).argmin()
    l.set_xdata(position[i])
    fig.canvas.draw_idle()
stime.on_changed(update)

plt.show()

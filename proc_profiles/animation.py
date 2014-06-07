"""
Create animation from for recreated profiles.
modifications: Cassio Amador
TODO: check automatically number of files, to correct parameter 'sw_clustersize'.
==========================
ORIGINAL:
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
===========================
"""
from os import path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

from proc_profile import ProcProfile

prof_folder=path.join(getcwd(), "..", "PROC")
shot_number=28061
sw_clustersize=4
shot=ProcProfile(shot_number)
ne=np.loadtxt(path.join(prof_folderl,"ne.dat"))

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(-0.02, 0.25), ylim=(0., 2e19))
#put legend. coordinates are in scale to window size
time_text = ax.text(0.02, 0.9, '', transform=ax.transAxes)
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

#Animation function.  This is called sequentially, defined by 'frames'
#in the 'FuncAnimation' call
def animate(i):
    time=shot.sweep2time(i*sw_clustersize)
    r = np.loadtxt(path.join(prof_folder,"%d.dat" % (time*1e3)))
    line.set_data(r, ne)
    time_text.set_text("#%s \n %.2f ms" % (shot_number,time))
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(shot.points)/sw_clustersize, interval=200, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('basic_animation.mp4', fps=10)#, extra_args=['-vcodec', 'libx264'])

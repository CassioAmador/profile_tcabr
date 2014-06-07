from os import chdir,listdir

import pylab as p

from input_explorer import inputExplorer
from proc_profile import ProcProfile

prof_folder="/home/cassio/fisica/Reflectometria/TCABR/camador_software/proc_raw/teste"
chdir(prof_folder)
ne=p.loadtxt('ne.dat')
info=p.loadtxt('prof_info.dat')
shot_number=28061
shot_number=28749


prof_list=listdir(prof_folder)
prof_list.sort()
position=p.empty(shape=((len(prof_list)-2),len(ne)))
times=p.empty(len(prof_list)-2)
print len(times)

i=0
for r_file in prof_list:
    name=r_file.strip('.dat')
    if name not in ('prof_info','ne'):
        position[i]=p.loadtxt(r_file)*1e2
        times[i]=name
        i+=1
        
fig,ax = p.subplots(1)
ax.hold(False)
 
def plotDynamics(time):
    #ax.clear()
    i=(abs(time*1e3-times)).argmin()
    #print i,time
    ax.plot(position[i],ne,label='%.2f ms' % time)
    ax.legend(loc='upper left')
    ax.set_xlabel('r (cm)')
    ax.set_ylabel('density (10^19 m^-3)')
    ax.set_xlim(0,20)
    ax.set_title("# %s" % shot_number)
    fig.canvas.draw()
 
sliders = [ { 'label' :  label,  'valmin': info[1] , 'valmax': info[2],
                'valinit': info[1]+(info[2]-info[1])*0.2 }
         for label in [ 'time' ] ]
 
inputExplorer(plotDynamics,sliders)

import time 
import pylab as p 
from os import path
from sys import argv
import proc_profile as pp

#input shot_number
if len(argv)>1:
    shot_number=int(argv[1])
else:
    shot_number = 32210
save_locally=0
time0 = time.time()
shot = pp.ProcProfile(shot_number,save_locally=save_locally)
#sweeps_average = 10
#Average of sweeps in time set in sweeps_average_time
sweeps_average_time = 1 # in ms
sweeps_average = int(sweeps_average_time*1e3/(shot.sweep_dur + shot.interv_sweep))
print sweeps_average
initial_sweep = 0
last_sweep = len(shot.points)
initial_time = 0
last_time = 170
initial_sweep = shot.time2sweep(initial_time)
last_sweep = shot.time2sweep(last_time)
#'all_shot' set to 1 avoids printing unnecessary information.
shot.reference_gd(all_shot=1)
print("time for reading files: %s s" % (time.time() - time0))
time1 = time.time()
#don't go until the end or some problems could appear for the last average.
for sweep in p.arange(initial_sweep, last_sweep-sweeps_average, sweeps_average):
    # print sweep
    shot.plasma_gd(sweep, sweeps_average, all_shot=1)
    shot.overlap_gd()
    shot.init_gd()
    shot.find_ne_max()
    # choose poly fit order
    shot.profile_poly_2(order=2, all_shot=1)
    # save profile in file with time in microsseconds
    if save_locally==1:
        p.savetxt(path.join(shot.prof_folder, "%06d.dat" % (shot.sweep2time(sweep) * 1e3)), shot.r)
shot.ne_max=p.array(shot.ne_max)
if save_locally==1:
    # separate density in other file, to save space.
    shot.save_ne()
    # save info file with parameters used to evaluate profiles.
    shot.save_proc_info(sweeps_average,initial_time,last_time)
    # save maximum estimated density.
    p.savetxt(path.join(shot.prof_folde, "ne_max.dat"), shot.ne_max)
print("time for Processing: %s s" % (time.time() - time1))
p.plot(shot.ne_max,'r')
p.plot(shot.ne_max,'b.')
p.show()
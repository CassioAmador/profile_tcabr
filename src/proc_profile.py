"""
Function: Recreate density profile for each milissecond of shot.
Stores it on MDSPlus database
First uses the reference, then evaluates overlap.
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: Make a better info file with parameters used.
Have an option to choose average by time
"""

import numpy as np
from os.path import join

from proc_profile_bottollier import ProcProfile
from read_shots_folder import check_prof_folder


def save_density(prof_folder, ne):
    np.save(join(prof_folder, "ne"), ne)


def save_profile(prof_folder, pos, time):
    np.save(join(prof_folder, "%06d" % (time)), pos)


def save_info(prof_folder, sweeps_average, initial_time, last_time):
    np.savetxt(join(prof_folder, "prof_info.dat"), [sweeps_average, initial_time, last_time])

if __name__ == "__main__":
    # change the shot number here
    shot_number = 33708
    # save data locally?
    save = 1
    time_on = 1
    shot = ProcProfile(shot_number, save_locally=save)
    shot.reference_gd()
    shot.eval_freq_overlap()
    if time_on == 1:
        import time
        time0 = time.time()
    initial_time = 29
    last_time = 150
    initial_sweep = shot.time2sweep(initial_time)
    last_sweep = shot.time2sweep(last_time)
    if last_sweep > len(shot.points):
        last_sweep = len(shot.points)
    # 'all_shot' set to 1 avoids printing unnecessary information.
    shot.reference_gd(all_shot=1)
    # sweeps average for each profile
    sweeps_average = 20
    # take out last sweep just to be sure.
    sweeps_array = np.arange(initial_sweep, last_sweep, sweeps_average)[:-1]
    # Evaluates once to check number of points per profile.
    shot.profile(1, 1, all_shot=1)
    if save == 1:
        prof_folder = check_prof_folder(shot_number)
        save_info(prof_folder, sweeps_average, initial_time, last_time)
        save_density(prof_folder, shot.ne)

    for sweep in sweeps_array:
        print("Sweep: ", sweep)
        shot.profile(sweep, sweeps_average, all_shot=1)
        if save == 1:
            save_profile(prof_folder, shot.pos, shot.sweep2time(sweep) * 1e3)

    if time_on == 1:
        delta_time = time.time() - time0
        print("it took: ", delta_time)
        print("it took per profile, in average: ", delta_time / len(sweeps_array))
        print("each profile had # of spectrograms: ", sweeps_average)

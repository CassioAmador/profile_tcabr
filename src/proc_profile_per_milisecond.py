"""
Function: Recreate density profile for each milissecond, with Bottolier-Curtet method.
Stores it on MDSPlus database.
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: Make a better info file with parameters used.
"""

import numpy as np
import MDSplus as mds

from proc_profile_bottollier import Bottollier

if __name__ == "__main__":
    # change the shot number here
    shot_number = 32211
    shot = Bottollier(shot_number, save_locally=0)
    initial_time = 0
    last_time = 170
    initial_sweep = shot.time2sweep(initial_time)
    last_sweep = shot.time2sweep(last_time)
    if last_sweep > len(shot.points):
        last_sweep = len(shot.points)
    # 'all_shot' set to 1 avoids printing unnecessary information.
    shot.reference_gd(all_shot = 1)
    # sweeps average each milissecond
    sweeps_average = int(1*1e3/(shot.sweep_dur + shot.interv_sweep))
    #take out last sweep just to be sure.
    sweeps_array=np.arange(initial_sweep, last_sweep, sweeps_average)[:-1]
    #Evaluates profile once to check number of points per profile.
    shot.profile(1,1,all_shot = 1)
    #matrix of position by time and density.
    shot.matrix = np.empty((len(sweeps_array),len(shot.pos)))

    i = 0
    for sweep in sweeps_array:
        print(sweep)
        shot.profile(sweep, sweeps_average, all_shot=1)
        shot.matrix[i] = shot.pos
        i += 1

    #Save to tree
    tree_name = "tcabr_ref"
    tree = MDSplus.Tree(tree_name, self.shot)
    # time array.
    node = tree.getNode("\\prof_time.signal")
    data = MDSplus.Float32Array(np.arange(initial_time,last_time,interval_time),dtype=np.float32)
    data.setUnits("ms")
    # density array. It will be the same for all shot.
    node = tree.getNode("\\prof_density.signal")
    data = MDSplus.Float32Array(shot.ne_poly,dtype=np.float32)
    # position matrix
    node = tree.getNode("\\prof_position.signal")
    data = MDSplus.Float32Array(shot.matrix,dtype=np.float32)
    data.setUnits("m")
import pylab as p
from proc_profile import ProcProfile

class Plasma_Behavior(ProcProfile):
    """Evaluate plasma behavior of whole shot"""
    def __init__(self, shot, tipo="data"):
        ProcProfile.__init__(self,shot,tipo)
        #colocar para testar se tem ou nao dados computados. Se nao tiver, processar todos.
        #ARRUMAR

    def proc_shot(self):
        shot.sweeps_average = 10
        initial_sweep = 0
        last_sweep = len(shot.points)
        initial_time = 58
        last_time = 95
        initial_sweep = shot.time2sweep(initial_time)
        last_sweep = shot.time2sweep(last_time)
        #'all_shot' set to 1 avoids printing unnecessary information.
        shot.reference_gd(all_shot=1)
        print("time for reading files: %s s" % (time.time() - time0))
        time1 = time.time()
        for sweep in p.arange(initial_sweep, last_sweep, sweeps_average):
            # print sweep
            shot.plasma_gd(sweep, sweeps_average, all_shot=1)
            shot.overlap_gd()
            shot.init_gd()
            # choose poly fit order
            shot.profile_poly(order=7, all_shot=1)
            # save profile in file with time in microsseconds
            p.savetxt(path.join(prof_folder, "%06d.dat" % (shot.sweep2time(sweep) * 1e3)), shot.r)
        # separate density in other file, to save space.
        p.savetxt(path.join(prof_folder, "ne.dat"), shot.ne_poly)
        # save info file with parameters used to evaluate profiles.
        p.savetxt(path.join(prof_folder, "prof_info.dat"), [sweeps_average, initial_time, last_time])
        print("time for Processing: %s s" % (time.time() - time1))

    def plot_variance(self,X,mat):
        var=p.var(mat,0, dtype=p.float64)
        p.plot(X,var)
        p.show()

    def contour(self,mat):
        pass

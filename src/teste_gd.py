# coding: utf-8
import pylab as p
import proc_profile as pp
p.ion()
shot=pp.ProcProfile(27922)
shot.reference_gd(all_shot=1)
shot.time2sweep(100)
shot.plasma_gd(6666,4,0)
shot.rate
shot.sweep_dur
shot.overlap_gd()
shot.profile_poly_4()
shot.plot_spectrogram()
p.twinx()
p.plot(shot.X,shot.gd_m,'g.-')
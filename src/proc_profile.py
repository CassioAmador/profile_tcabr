"""
Function: Recreate density profile from beating signal. Reference is
vessel wall, and calibration was done with a pin.
First uses the reference, then evaluates overlap.
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: give an option to overlap bands in the spectrogram.
Incorporate other abel inversions.
Have a better way to set folder to store profiles.
Make an 'ini' file.
"""

from os import path, getcwd

import pylab as p
from scipy import interpolate

import ref_functions as rf
from proc_sweep import ProcSweep

#c/pi
const=95426903.18473884
#distance (39.5 cm) in nanoseconds from vessel wall to limiter, x2
wall_corr=2.6351563520654012
#polynomial coeffs of order 15, for the abel integral inversion
abel_factor=p.array([2048/6435.,429*p.pi/4096.,1024/3003.,231*p.pi/2048.,
            256/693.,63*p.pi/512,128/315.,35*p.pi/256.,16/35.,
            5*p.pi/32.,8/15.,3*p.pi/16.,2/3.,p.pi/4.,1,p.pi/2.],dtype=p.float32)
    

class ProcProfile(ProcSweep):
    """recreate profiles from specific shot"""
    def __init__(self,shot,tipo="data"):
        #inherit from ProcSweep class
        ProcSweep.__init__(self,shot,tipo)
        self.mark_sweep_points()
        
    def average_specgram(self,sweeps=8,sweep_ini=20,all_shot=0):
        """Average spectrograms of a specified cluster of sweeps."""
        #if the classe has already a mean stored, deletes it.
        if hasattr(self, "gd_k"):
            del(self.matrix_k_mean)
        if hasattr(self, "gd_ka"):
            del(self.matrix_ka_mean)
        #makes the spectrogram for each sweep.
        for sweep in range(sweep_ini,sweep_ini+sweeps):
            self.read_single_sweep("K",sweep)
            self.read_single_sweep("Ka",sweep)
            if all_shot==0:
                print self.sweep_cur
            matrix_k,X_k,Y_k=self.spectrogram('K',figure=sweep+1,ploti=0)
            matrix_ka,X_ka,Y_ka=self.spectrogram('Ka',figure=sweep+1,ploti=0)
            if hasattr(self,'matrix_k'):
                self.matrix_k_mean+=matrix_k
                self.matrix_ka_mean+=matrix_ka
            #if there are no averages, creates it. 
            else:
                self.matrix_k_mean=matrix_k.copy()
                self.X_k=X_k.copy()
                self.Y_k=Y_k.copy()
                self.matrix_ka_mean=matrix_ka.copy()
                self.X_ka=X_ka.copy()
                self.Y_ka=Y_ka.copy()
            
        self.matrix_k_mean/=sweeps
        self.matrix_ka_mean/=sweeps

    def reference_gd(self,sw_clustersize=20,all_shot=0):
        """Average the spectrogram and evaluates the group delay for the
         reference position. In this case, the reflection in the wall is
          stronger than in the pin, so we evaluate the time delay of the
          wall, and subtracts the position of limiter. The calibration 
          is done with a pin, and stored in 'wall_cor' parameter.
          Average is done in the first and last sweeps of the shot, 
          when there is no plasma"""
        self.average_specgram(sweeps=sw_clustersize,sweep_ini=1,all_shot=all_shot)
        self.average_specgram(sweeps=sw_clustersize,sweep_ini=-sw_clustersize,all_shot=all_shot)
        self.gd_k0=find_max(self.matrix_k_mean,self.Y_k)-wall_corr
        self.gd_ka0=find_max(self.matrix_ka_mean,self.Y_ka)-wall_corr

    def plasma_gd(self,sweep_cur=1000,sw_clustersize=8,all_shot=0):
        """Evaluates group delay for a cluster of sweeps, and subtracts
        the reference to find the group delay from the limiter"""
        self.average_specgram(sweeps=sw_clustersize,sweep_ini=sweep_cur,all_shot=all_shot)
        self.gd_k_mean=find_max(self.matrix_k_mean,self.Y_k)
        self.gd_ka_mean=find_max(self.matrix_ka_mean,self.Y_ka)
        self.gd_k=self.gd_k_mean-self.gd_k0
        self.gd_ka=self.gd_ka_mean-self.gd_ka0

    def overlap_gd(self):
        """overlap group delay from K and Ka bands, when they probed the
         same density layers"""
        #change name for convenience. X represents the probing frequency.
        xk,xka=self.X_k,self.X_ka
        #check the number of points of overlaped frquency.
        cmin=(abs(xk-min(xka))).argmin()
        cmax=(abs(xka-max(xk))).argmin()
        #create an even spaced array in the overlap region
        xn=p.arange(self.X_k[cmin],xka[cmax],xka[1]-xka[0])
        #frequency array for all bands. 
        #first and last points are discarded because of signal noise.
        self.X=p.concatenate((xk[9:cmin],xn,xka[cmax:-3]))
        #interpolate the group delay in the overlap region
        gdkn=p.interp(xn,xk,self.gd_k)
        gdkan=p.interp(xn,xka,self.gd_ka)
        gdn=(gdkn+gdkan)/2
        self.gd_m=p.concatenate((self.gd_k[9:cmin],gdn,self.gd_ka[cmax:-3]))

    def init_gd(self,tipo="reta"):
        """Group delay initialization, default is first order. Second 
        order is also available."""
        #choose number of steps to start group_delay
        steps=6
        self.ini_f=p.linspace(0,self.X[0],num=steps)
        self.ini_t=p.linspace(0,self.gd_m[0],num=steps)
        polynomial = p.polyfit(self.ini_f,self.ini_t, 1)
        if tipo=='quad':
            mean_pos=len(self.ini_f)/2
            mean_t=self.ini_t[mean_pos]*(1.15)
            polynomial = p.polyfit((self.ini_f[0],self.ini_f[mean_pos],self.ini_f[-1]),
                                      (self.ini_t[0],mean_t,self.ini_t[-1]), 2)
        self.pol=p.poly1d(polynomial)
        self.ini_t=self.pol(self.ini_f)
        #frequency and group delay now starts from zero.
        self.freqs=p.concatenate((self.ini_f,self.X))
        self.gd=p.concatenate((self.ini_t,self.gd_m))
        if not hasattr(self,'ne'):
                self.ne=rf.freq2den(self.freqs)

    def profile_poly(self,order=9,all_shot=0):
        """Recreates a profile fiiting the group delay with a polynomial
        of specified order. The abel inversion is evaluated explicitly,
        since for a poly fit the integral turns out to be another
        polynomial function. More information on the manual."""
        polynomial = p.polyfit(self.freqs*1e9,self.gd*1e-9, order)
        self.pol=p.poly1d(polynomial)
        #self.gd_poly=pol(self.freqs*1e9)
        self.freq_poly=p.linspace(0,self.freqs[-1]*1e9,num=100)
        if not hasattr(self,'ne_poly'):
            self.ne_poly=rf.freq2den(self.freq_poly)
        if all_shot==0:
            #only needed if one wants to check the polynomial fit.
            self.gd_poly=self.pol(self.freq_poly)
        self.pol_coeffs=self.pol.coeffs*abel_factor[-order-1:]
        self.r=const*p.array([poly_eval(self.pol_coeffs,freq) for freq in self.freq_poly])

def find_max(matrix,spec_yaxis,ploti=0,spec_xaxis=0):
    """find the curve which follow the maximum for each spectrum"""
    matrix=matrix.transpose()
    gd=p.empty(len(matrix))
    for r in range(len(matrix)):
        gd[r]=spec_yaxis[p.where(matrix[r]==matrix[r].max())[0]]
    if ploti==1:
        p.figure(4)
        p.plot(spec_xaxis,gd,'b-')
        #p.ylim(freq_lim_min,freq_lim_max)
    return gd

def poly_eval(coeffs,x):
    """evaluates polynomial for a point 'x', with custom coeffs."""
    order=len(coeffs)
    s=[coeffs[i]*x**(order-i-1) for i in range(order)]
    return sum(s)

if __name__=="__main__":
    #change the shot number here
    shot_number=28061
    shot=ProcProfile(shot_number)
    #choose folder to store profiles
    lugar=path.join(getcwd(), "../post/")
    sweeps_average=4
    initial_sweep=0
    last_sweep=len(shot.points)
    #'all_shot' set to 1 avoids printing unnecessary information.
    shot.reference_gd(all_shot=1)
    #change alias for desired method.
    profile=shot.profile_poly
    for sweep in p.arange(initial_sweep,last_sweep,sweeps_average):
        print sweep
        shot.plasma_gd(sweep,sweeps_average,all_shot=1)
        shot.overlap_gd()
        shot.init_gd()
        #choose poly fit order
        profile(order=7,all_shot=1)
        #save profile in file with time in microsseconds
        p.savetxt(path.join(lugar,"%d.dat" % (shot.sweep2time(sweep)*1e3)),shot.r)
    #separate density in other file, to save space.
    p.savetxt(path.join(lugar,"ne.dat"),shot.ne_poly)

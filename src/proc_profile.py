"""
Function: Recreate density profile from beating signal. Reference is
vessel wall, and calibration was done with a pin.
First uses the reference, then evaluates overlap.
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: give an option to overlap bands in the spectrogram.
evaluate abel inversion separatedly from gd and initialization.
Incorporate other abel inversions.
Have a better way to set folder to store profiles.
Make an 'ini' file.
Make a better info file with parameters used.
"""

from os import path, getcwd, makedirs
import time
import pylab as p
import numpy as np
import math
from scipy import special as sp_sp
#os.add_path("/home/cassio/Dropbox/tcabr_ref/perfil/group_delay.py")


import ref_functions as rf
from proc_sweep import ProcSweep

# 1e-9*c/pi, freq in GHz, gd in ns
const = 95426903.18473884*1e-9
# distance (39.5 cm) in nanoseconds from vessel wall to limiter, x2
wall_corr = 2.6351563520654012
# polynomial coeffs of order 15, for the abel integral inversion
abel_factor = p.array([2048 / 6435., 429 * p.pi / 4096., 1024 / 3003., 231 * p.pi / 2048.,
                       256 / 693., 63 * p.pi / 512, 128 / 315., 35 * p.pi / 256., 16 / 35.,
                       5 * p.pi / 32., 8 / 15., 3 * p.pi / 16., 2 / 3., p.pi / 4., 1, p.pi / 2.], dtype=p.float32)


class ProcProfile(ProcSweep):

    """recreate profiles from specific shot"""

    def __init__(self, shot, tipo="data",save_locally=1):
        self.version=1.0
        # inherit from ProcSweep class
        ProcSweep.__init__(self, shot, tipo,save_locally)
        self.mark_sweep_points()
        self.check_prof_folder()
        
    def check_prof_folder(self):
        self.shotproc_folder = path.join(getcwd(), "..", "PROC","%s" % self.shot)
        self.prof_folder = path.join(self.shotproc_folder,"level_0")
        if path.exists(self.shotproc_folder):
            pass
        else:
            makedirs(self.shotproc_folder)
            makedirs(self.prof_folder)

    def average_specgram(self, sweeps=8, sweep_ini=20, all_shot=0):
        """Average spectrograms of a specified cluster of sweeps."""
        # if the class has already a mean stored, deletes it.
        if hasattr(self, "gd_k"):
            del(self.matrix_k_mean)
        if hasattr(self, "gd_ka"):
            del(self.matrix_ka_mean)
        # makes the spectrogram for each sweep.
        for sweep in range(sweep_ini, sweep_ini + sweeps):
            self.read_single_sweep("K", sweep)
            self.read_single_sweep("Ka", sweep)
            if all_shot == 0:
                print(self.sweep_cur)
            matrix_k, X_k, Y_k = self.spectrogram('K', figure=sweep + 1, normal=0, ploti=0,freqs=(3e3,15e3))
            matrix_ka, X_ka, Y_ka = self.spectrogram('Ka', figure=sweep + 1, normal=0, ploti=0, freqs=(2e3,15e3))
            if hasattr(self, 'matrix_k_mean'):
                self.matrix_k_mean += matrix_k
                self.matrix_ka_mean += matrix_ka
            # if there are no averages, creates it.
            else:
                self.matrix_k_mean = matrix_k.copy()
                self.X_k = X_k.copy()
                self.Y_k = Y_k.copy()
                self.matrix_ka_mean = matrix_ka.copy()
                self.X_ka = X_ka.copy()
                self.Y_ka = Y_ka.copy()

        self.matrix_k_mean /= sweeps
        self.matrix_ka_mean /= sweeps

    def reference_gd(self, sw_clustersize=20, all_shot=0):
        """Average the spectrogram and evaluates the group delay for the
         reference position. In this case, the reflection in the wall is
          stronger than in the pin, so we evaluate the time delay of the
          wall, and subtracts the position of limiter. The calibration
          is done with a pin, and stored in 'wall_cor' parameter.
          Average is done in the first and last sweeps of the shot,
          when there is no plasma"""
        self.average_specgram(sweeps=sw_clustersize, sweep_ini=1, all_shot=all_shot)
        self.gd_k0 = find_max(self.matrix_k_mean, self.Y_k) - wall_corr
        self.gd_ka0 = find_max(self.matrix_ka_mean, self.Y_ka) - wall_corr

    def plasma_gd(self, sweep_cur=1000, sw_clustersize=8, all_shot=0):
        """Evaluates group delay for a cluster of sweeps, and subtracts
        the reference to find the group delay from the limiter"""
        self.average_specgram(sweeps=sw_clustersize, sweep_ini=sweep_cur, all_shot=all_shot)
        self.gd_k_mean = find_max(self.matrix_k_mean, self.Y_k)
        self.gd_ka_mean = find_max(self.matrix_ka_mean, self.Y_ka)
        self.gd_k = self.gd_k_mean - self.gd_k0
        self.gd_ka = self.gd_ka_mean - self.gd_ka0

    def overlap_gd(self):
        """overlap group delay from K and Ka bands, when they probed the
         same density layers"""
        if min(self.X_ka)<=max(self.X_k):
            # change name for convenience. X represents the probing frequency.
            xk, xka = self.X_k, self.X_ka
            # check the number of points of overlaped frquency.
            cmin = (abs(xk - min(xka))).argmin()
            cmax = (abs(xka - max(xk))).argmin()
            # create an even spaced array in the overlap region
            xn = p.arange(self.X_k[cmin], xka[cmax], xka[1] - xka[0])
            # frequency array for all bands.
            # first and last points are discarded because of signal noise.
            #self.X = p.concatenate((xk[9:cmin], xn, xka[cmax:-3]))
            self.X = p.concatenate((xk[:cmin], xn, xka[cmax:]))
            # interpolate the group delay in the overlap region
            gdkn = p.interp(xn, xk, self.gd_k)
            gdkan = p.interp(xn, xka, self.gd_ka)
            gdn = (gdkn + gdkan) / 2
            #self.gd_m = p.concatenate((self.gd_k[9:cmin], gdn, self.gd_ka[cmax:-3]))
            self.gd_m = p.concatenate((self.gd_k[:cmin], gdn, self.gd_ka[cmax:]))
        else:
            self.X = p.concatenate((xk[9:cmin], xka[cmax:-3]))
            self.gd_m = p.concatenate((self.gd_k[9:cmin], self.gd_ka[cmax:-3]))
        #self.adjust_gd()

    def init_gd(self, tipo="line"):
        """Group delay initialization, default is first order. Second
        order is also available."""
        # choose number of steps to start group_delay
        steps = 6
        self.ini_f = p.linspace(0, self.X[0], num=steps)[:-1]
        self.ini_t = p.linspace(0, self.gd_m[0], num=steps)[:-1]
        polynomial = p.polyfit(self.ini_f, self.ini_t, 1)
        if tipo == 'quad':
            mean_pos = len(self.ini_f) / 2
            mean_t = self.ini_t[mean_pos] * (1.15)
            polynomial = p.polyfit((self.ini_f[0], self.ini_f[mean_pos], self.ini_f[-1]),
                                   (self.ini_t[0], mean_t, self.ini_t[-1]), 2)
        self.pol = p.poly1d(polynomial)
        self.ini_t = self.pol(self.ini_f)
        # frequency and group delay now starts from zero.
        self.freqs = p.concatenate((self.ini_f, self.X))
        self.gd = p.concatenate((self.ini_t, self.gd_m))
        if not hasattr(self, 'ne'):
            self.ne = rf.freq2den(self.freqs*1e9)

    def find_ne_max(self):
        if not hasattr(self,"ne_max"):
            self.ne_max=[]
        zeros_gd = p.where(self.gd_m <= 0.0)[0]
        #print max(p.diff(self.gd_m))
        if len(zeros_gd)!=0:
            print "pequeno"
            self.ne_max.append(rf.freq2den(1e9*self.X[zeros_gd[0]]))
            self.gd_m[zeros_gd[0]:]=float("NaN")
        else:
            dif_gd = p.where(p.diff(self.gd_m) > 3.5)[0]
            if len(dif_gd)!=0:
                print "grande"
                self.ne_max.append(rf.freq2den(1e9*self.X[dif_gd[-1]]))
                self.gd_m[dif_gd[0]:]=float("NaN")
            else:
                self.ne_max.append(float("NaN"))
        #zero_jumps = p.where(p.diff(zero_gd) > 1)[0]

    def adjust_gd(self):
        zeros = p.where(self.X < 18)[0]
        if len(zeros)!=0:
            cut=(abs(self.X - 18)).argmin()
            self.X=self.X[cut:]
            self.gd_m=self.gd_m[cut:]
        # find all points lower than zero
        zero_gd = p.where(self.gd_m <= 0)[0]
        #zero_gd = p.concatenate((zero_gd,zero_gd+10,zero_gd+15))
        while(len(zero_gd)!=0):
            zero_jumps = p.where(p.diff(zero_gd) > 1)[0]
            if len(zero_jumps)==0:
                zeros = zero_gd
            else:
                zeros = zero_gd[:zero_jumps[0]+1]
            print zero_gd,zeros
            if 0 in zeros:
                self.X=self.X[zeros[-1]+1:]
                self.gd_m=self.gd_m[zeros[-1]+1:]
            else:
                poly = p.polyfit((self.X[zeros[0]-1],self.X[zeros[-1]+1]),(self.gd_m[zeros[0]-1],self.gd_m[zeros[-1]+1]),1)
                pol = p.poly1d(poly)
                #p.plot(self.X[zeros[0]-1:zeros[-1]+2],self.gd_m[zeros[0]-1:zeros[-1]+2],'.-')
                self.gd_m[zeros] = pol(self.X[zeros])
                #p.plot(self.X[zeros[0]-1:zeros[-1]+2],self.gd_m[zeros[0]-1:zeros[-1]+2],'.-')
            zero_gd = p.where(self.gd_m <= 0)[0]

    def profile_poly(self, order=9, all_shot=0):
        """Recreates a profile fiiting the group delay with a polynomial
        of specified order. The abel inversion is evaluated explicitly,
        since for a poly fit the integral turns out to be another
        polynomial function. More information on the manual."""
        polynomial = p.polyfit(self.freqs, self.gd, order)
        self.pol = p.poly1d(polynomial)
        if not hasattr(self, 'freq_poly'):
            self.freq_poly = p.linspace(0, self.freqs[-1], num=100)
        if not hasattr(self, 'ne_poly'):
            self.ne_poly = rf.freq2den(self.freq_poly)
        if all_shot == 0:
            # only needed if one wants to check the polynomial fit.
            self.gd_poly = self.pol(self.freq_poly)
        self.pol_coeffs = self.pol.coeffs * abel_factor[-order - 1:]
        self.r = const * p.array([poly_eval(self.pol_coeffs, freq) for freq in self.freq_poly])

    def abel_factor_2(self,f_c,f_i,f_f,order=2):
        # polynomial coeffs for any order k, for the abel integral inversion starting from non-zero frequency.
        # f_c: frequency of reflection
        # f_i: lower limit of integration
        # f_f: upper limit of integration
        self.abel = p.zeros(order+1)
        if f_c == f_f:
            self.abel[0] = p.arccos(f_i/f_c)
            for k in range(1,order+1):
                self.abel[k]=math.pow(f_c,k)*p.sqrt(p.pi)*sp_sp.gamma((1 + k)/2)/(k*sp_sp.gamma(k/2)) - math.pow(f_i,1 + k) * sp_sp.hyp2f1(1/2, (1 + k)/2, (3 + k)/2, math.pow(f_i/f_c,2))/(f_c*(1 + k))
                #print 'f_c == f_f',f_c,f_i,f_f,k,self.abel[k]
        if f_c != f_f:
            if f_i==0:
                for k in range(order+1):
                    #self.abel[k] = (math.pow(B, k + 1) / A) * (sp_sp.hyp2f1(1 / 2., (k + 1) / 2., (k + 3) / 2., np.pow(B / A, 2))) / (k + 1.)
                    self.abel[k]= math.pow(f_f,1 + k) * sp_sp.hyp2f1(1, (1 + k)/2, (3 + k)/2,  math.pow(f_f/f_c,2))/(f_c * (1 + k))
                    #print 'f_c != f_f, f_i==0',f_c,f_i,f_f,k,self.abel[k]
            else:
                for k in range(order+1):
                    #self.abel[k] = (math.pow(B, k + 1) / A) * (sp_sp.hyp2f1(1 / 2., (k + 1) / 2., (k + 3) / 2., np.pow(B / A, 2))) / (k + 1.)
                    self.abel[k]= 1/(f_c * (1 + k)) * (-math.pow(f_i,1 + k) * sp_sp.hyp2f1(1, (1 + k)/2, (3 + k)/2, math.pow(f_i/f_c,2)) + math.pow(f_f,1 + k) * sp_sp.hyp2f1(1, (1 + k)/2, (3 + k)/2,  math.pow(f_f/f_c,2)))
                    #print 'f_c != f_f, f_i!=0',f_c,f_i,f_f,k,self.abel[k]
        # if order>=0:
        #     self.abel[0]=1.5708 -(A * sp_sp.hyp2f1(0.5,1./2,3./2,math.pow(A/B,2)))/B
        # if order>=1:
        #     self.abel[1]=math.sqrt(-math.pow(A,2)+math.pow(B,2))
        # if order>=2:
        #     self.abel[2]=0.785398 * math.pow(B,2)-(0.333333 * math.pow(A,3) * sp_sp.hyp2f1(0.5,3./2,5./2,math.pow(A/B,2)))/B
        # if order>=3:
        #     self.abel[3]=math.sqrt(-math.pow(A,2)+math.pow(B,2)) * (0.333333 * math.pow(A,2)+0.666666 * math.pow(B,2))
        # if order>=4:
        #     self.abel[4]=0.589049 * math.pow(B,4)-(0.2 * math.pow(A,5) * sp_sp.hyp2f1(0.5,5./2,7./2,math.pow(A/B,2)))/B
        # if order>=5:
        #     self.abel[5]=math.sqrt(-math.pow(A,2)+math.pow(B,2)) * (0.2 * math.pow(A,4)+0.266667 * math.pow(A,2) * math.pow(B,2)+0.533333 * math.pow(B,4))
        # if order>=6:
        #     self.abel[6]=0.490874 * math.pow(B,6)-(0.142857 * math.pow(A,7) * sp_sp.hyp2f1(0.5,7./2,9./2,math.pow(A/B,2)))/B
        # if order>=7:
        #     self.abel[7]=math.sqrt(-math.pow(A,2)+math.pow(B,2)) * (0.142857*  math.pow(A,6)+0.171429 * math.pow(A,4) * math.pow(B,2)+0.228571 * math.pow(A,2) * math.pow(B,4)+0.457143 * math.pow(B,6))
        # if order>=8:
        #     self.abel[8]=0.429515 * math.pow(B,8)-(0.111111 * math.pow(A,9) * sp_sp.hyp2f1(0.5,9./2,11./2,math.pow(A/B,2)))/B
        # if order>=9:
        #     self.abel[9]=math.sqrt(-math.pow(A,2)+math.pow(B,2)) * (0.111111 * math.pow(A,8)+0.126984 * math.pow(A,6) * math.pow(B,2)+0.152381 * math.pow(A,4) * math.pow(B,4)+0.203175 * math.pow(A,2) * math.pow(B,6)+0.406349 * math.pow(B,8))
        # if order>=10:
        #     self.abel[10]=0.386563 * math.pow(B,10)-(0.0909091 * math.pow(A,11) * sp_sp.hyp2f1(0.5,11./2,13./2,math.pow(A/B,2)))/B

    def profile_poly_6(self, order=9, all_shot=0):
        """linear initialization, group delay dividec in two parts"""
        self.r=p.zeros(1+len(self.X))
        polynomial = p.polyfit(self.X, self.gd_m, order)
        self.pol = p.poly1d(polynomial)
        if not hasattr(self, 'freqs'):
            self.freqs = p.concatenate(([0], self.X))
        if not hasattr(self, 'ne'):
            self.ne = rf.freq2den(self.freqs*1e9)
        if all_shot == 0:
            # only needed if one wants to check the polynomial fit.
            self.gd_poly = self.pol(self.X)
            self.gd=p.concatenate(([0],self.gd_poly))
        coeff = (self.gd_m[0]/self.X[0]) * 1e-18
        self.abel_factor_2(self.X[0],0,self.X[0],1)
        self.r[0] = const * coeff * self.abel[1]
        for f in range(1,len(self.X)):
            self.abel_factor_2(self.X[f],0,self.X[0],1)
            r_0 = const * coeff * self.abel[1]
            self.abel_factor_2(self.X[f],self.X[0],self.X[f],order)
            print r_0
            self.r[f] = r_0+const * poly_eval(self.pol.coeffs * self.abel[::-1], self.X[f])

    def profile_poly_2(self, order=2, all_shot=0):
        init = 2
        f_probe = self.X
        tau = self.gd_m
        self.r = np.zeros(len(tau))
        tau2_coef = np.polyfit(f_probe, tau, order)
        if not hasattr(self, 'freq_poly'):
            self.freq_poly = self.X * 1e9
        if not hasattr(self, 'ne_poly'):
            self.ne_poly = rf.freq2den(self.freq_poly)
        # inicializao
        if init == 1:
            tau1_coef = np.array([np.polyval(tau2_coef, f_probe[0]) / f_probe[0], 0])
        elif init == 2:
            tau1_coef = np.array([np.polyval(tau2_coef, f_probe[0]) / f_probe[0] ** 2., 0, 0])
        I1 = np.zeros(len(tau1_coef))
        I2 = np.zeros(len(tau2_coef))
        I3 = np.zeros(len(tau2_coef))
        for i in range(len(f_probe)):
            for k in range(len(tau1_coef)):
                I1[k] = (np.power(f_probe[0], k + 1) / f_probe[i]) * \
                    (sp_sp.hyp2f1(1 / 2., (k + 1) / 2., (k + 3) / 2., np.power(f_probe[0] / f_probe[i], 2))) / (k + 1.)
            for k in range(len(tau2_coef)):
                I2[k] = (np.power(f_probe[0], k + 1) / f_probe[i]) * \
                    (sp_sp.hyp2f1(1 / 2., (k + 1) / 2., (k + 3) / 2., np.power(f_probe[0] / f_probe[i], 2))) / (k + 1.)
                I3[k] = np.power(f_probe[i], k) * np.sqrt(np.pi) * sp_sp.gamma(0.5 + k / 2.) / (2. * sp_sp.gamma(1 + k / 2.))
            self.r[i] = 1e-9 * (3e8 / np.pi) * (np.dot(tau1_coef[::-1], I1) + np.dot(tau2_coef[::-1], I3 - I2))

    def profile_poly_3(self, order=9, all_shot=0):
        """..."""
        # choose number of steps to start group_delay
        if not hasattr(self, 'freqs'):
            self.freqs = p.concatenate((p.zeros(1), self.X))
        if not hasattr(self, 'ne'):
            self.ne = rf.freq2den(self.freqs*1e9)
        self.gd = p.concatenate((p.zeros(1), self.gd_m))
        self.r=p.zeros(len(self.freqs))
        order=1
        for f in range(1,len(self.freqs)-order+1):
            polynomial = p.polyfit(self.freqs[f-1:f+order+1],self.gd[f-1:f+order+1],order)
            self.pol = p.poly1d(polynomial)
            self.abel_factor_2(self.freqs[f-1],self.freqs[f-1],self.freqs[f+order-1],order)
            #p.plot(self.freqs[f-1:f+order+1],self.gd[f-1:f+order+1],'.')
            #p.plot(self.freqs[f-1:f+order+1],self.pol(self.freqs[f-1:f+order+1]))
            #polynomial coefficients are in different order
            self.pol_coeffs = self.pol.coeffs * self.abel[::-1][-order-1:]
            self.r[f] = self.r[f-1] + const * poly_eval(self.pol_coeffs, self.freqs[f])
            print self.r[f],self.r[f-1]

    def profile_poly_4(self, all_shot=0):
        """..."""
        # choose number of steps to start group_delay
        if (not hasattr(self, 'freqs')) or (self.freqs.size!=(self.X.size-1)):
            self.freqs = p.concatenate((p.zeros(1), self.X))
        if (not hasattr(self, 'ne')) or (self.ne.size!=(self.X.size-1)):
            self.ne = rf.freq2den(self.freqs*1e9)
        self.gd = p.concatenate((p.zeros(1), self.gd_m))
        self.r=p.zeros(len(self.freqs))
        meth=1
        for f in range(1,len(self.freqs)):
            pos=0
            for ff in range(f+1):
                if f!=ff:
                    ss=(self.gd[ff]+self.gd[ff-1])*(self.freqs[ff]-self.freqs[ff-1])/(2*p.sqrt(math.pow(self.freqs[f],2)-math.pow(self.freqs[ff],2)))
                else:
                    if meth==0:
                        ss=self.gd[ff]/1e-9
                    else:
                        ss=1.5708 -(self.freqs[f-1] * sp_sp.hyp2f1(0.5,0.5,1.5,math.pow(self.freqs[f-1]/self.freqs[f],2)))/self.freqs[f]
                        ss=self.gd[f-1]*ss + ((self.gd[f]-self.gd[f-1])/(self.freqs[f]-self.freqs[f-1])) * math.sqrt(-math.pow(self.freqs[f-1],2)+math.pow(self.freqs[f],2))
                pos+=ss
            self.r[f]=const * pos

    def profile_int_part(self,all_shot=0):
        # choose number of steps to start group_delay
        steps = 3
        self.ini_f = p.linspace(0, self.X[0], num=steps)[:-1]
        self.ini_t = p.linspace(0, self.gd_m[0], num=steps)[:-1]
        self.r=p.zeros(len(self.ini_f)+len(self.X))
        if not hasattr(self, 'freqs'):
            self.freqs = p.concatenate((self.ini_f * 1e9, self.X * 1e9))
            self.freqs_dif=p.diff(self.freqs)
        self.gd=p.concatenate((self.ini_t * 1e-9,self.gd_m * 1e-9))
        if not hasattr(self, 'ne_poly'):
            self.ne = rf.freq2den(self.freqs)

        f0=0
        self.gd_dif=p.diff(self.gd)/self.freqs_dif
        for f in range(1,len(self.r)):
            freq=self.freqs[f]
            integral=sum([self.gd_dif[ff]*p.sqrt(freq*freq - self.freqs[ff]*self.freqs[ff])*(self.freqs[ff]-self.freqs[ff-1]) for ff in range(f) ])
            self.r[f]=self.gd[f]*(p.sqrt(freq*freq - f0*f0))/self.freqs[f] - integral
            self.r[f]*=const

    def save_ne(self):
        p.savetxt(path.join(self.prof_folder, "ne.dat"), self.ne_poly)

    def save_proc_info(self,sweeps_average,initial_time,last_time):
        p.savetxt(path.join(self.prof_folder, "prof_info.dat"), [self.version, sweeps_average, initial_time, last_time])

    def plot_spectrogram(self, colormap=0):
        """plots the mean spectrogram overlaped of both bands, at the
        current time evaluated in the class. Colormap set to 0 plots
        only contours, set to 1 plots a colormap."""
        if colormap == 0:
            cont = p.contour
        elif colormap == 1:
            cont = p.contourf
        cont(self.X_k, self.Y_k, self.matrix_k_mean)
        cont(self.X_ka, self.Y_ka, self.matrix_ka_mean)
        p.xlabel("freq (GHz)")
        p.ylabel("group delay (ns)")
        p.title("# %s - time: %s ms" % (self.shot, self.sweep2time(self.sweep_cur)))
        p.ylim(0.4, 13)

    def plot_groupdelay(self):
        p.plot(self.freqs, self.gd)
        p.plot(self.freq_poly * 1e-9, self.gd_poly * 1e9)
        p.xlabel("freq (GHz)")
        p.ylabel("group delay (ns)")
        p.title("# %s - time: %s ms" % (self.shot, self.sweep2time(self.sweep_cur)))


def find_max(matrix, spec_yaxis, ploti=0, spec_xaxis=0):
    """find the curve which follow the maximum for each spectrum"""
    matrix = matrix.transpose()
    gd = p.empty(len(matrix))
    for r in range(len(matrix)):
        gd[r] = spec_yaxis[p.where(matrix[r] == matrix[r].max())[0]]
    if ploti == 1:
        p.figure(4)
        p.plot(spec_xaxis, gd, 'b-')
        # p.ylim(freq_lim_min,freq_lim_max)
    return gd


def poly_eval(coeffs, x):
    """evaluates polynomial for a point 'x', with custom coeffs."""
    order = len(coeffs)
    s = [coeffs[i] * (x ** (order - i - 1)) for i in range(order)]
    return sum(s)

if __name__ == "__main__":
    # change the shot number here
    shot_number = 32210
    # shot_number=28749
    time0 = time.time()
    shot = ProcProfile(shot_number,save_locally=0)
    # choose folder to store profiles
    sweeps_average = 10
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
    for sweep in p.arange(initial_sweep, last_sweep, sweeps_average):
        # print sweep
        shot.plasma_gd(sweep, sweeps_average, all_shot=1)
        shot.overlap_gd()
        shot.init_gd()
        # choose poly fit order
        shot.profile_poly_2(order=2, all_shot=1)
        # save profile in file with time in microsseconds
        p.savetxt(path.join(shot.prof_folder, "%06d.dat" % (shot.sweep2time(sweep) * 1e3)), shot.r)
    # separate density in other file, to save space.
    p.savetxt(path.join(shot.prof_folder, "ne.dat"), shot.ne_poly)
    # save info file with parameters used to evaluate profiles.
    p.savetxt(path.join(shot.prof_folder, "prof_info.dat"), [sweeps_average, initial_time, last_time])
    print("time for Processing: %s s" % (time.time() - time1))

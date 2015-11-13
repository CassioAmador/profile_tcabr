import pylab as p
from scipy import integrate

import ref_functions as rf
from proc_group_delay import ProcGroupDelay

factor=2/3
c=299792458
# 2*pi/c, frequency in GHz
const = 41.91690043903363
const_fac=const*factor

class Phase(ProcGroupDelay):
    """Evaluates phase from group delay"""

    def __init__(self, shot, tipo="data", save_locally=1):
        self.version = 1.0
        # inherit from ProcGroupDealy class
        ProcGroupDelay.__init__(self, shot, tipo, save_locally)

    def overlap_phase(self):
        pass

    def area_prev(self,n,freqs_prob,pos):
        #frequency times refraction index
        refrac_index=rf.refrac(freqs_prob[:n-1],freqs_prob[n])
        #print(refrac_index)
        area=freqs_prob[n]*integrate.trapz(refrac_index,pos[:n-1])
        return area

    def pos_eval(self,n,freqs_prob,pos,phase_cur):
        #print(rf.refrac(freqs_prob[n-1],freqs_prob[n])*freqs_prob[n])
        area_dif=const_fac*freqs_prob[n]*rf.refrac(freqs_prob[n-1],freqs_prob[n])
        if n==1:
            phase_prev=0
        else:
            phase_prev=self.area_prev(n,freqs_prob,pos)*const
        print(phase_prev,phase_cur,pos[n-1])
        return pos[n-1]+(phase_cur-phase_prev+p.pi/2)/(area_dif)

    def find_pos(self,nX,phase):
        pos=p.zeros(len(nX))
        r0=0
        for n in range(1,len(nX)):
            pos[n]=self.pos_eval(n,nX,pos,phase[n])
        return pos+r0

    def profile(self,sweep,sweeps=4):
        if not hasattr(self,'gd_k0'):
            self.reference_gd()
        self.plasma_gd(sweep,sweeps)
        if not hasattr(self,'nX'):
            self.nX={}
            self.overlap_freq()
            for channel in ('K','Ka'):
                dif_X=(self.X[channel][1]-self.X[channel][0])/2
                self.nX[channel]=self.X[channel][:-1]+dif_X
                if channel=='K':
                    self.nX[channel]=p.concatenate(([0],self.nX[channel]))
            self.nX['2bands']=p.concatenate((self.nX['K'][:self.cmin],self.nX['Ka']))
            self.ne = rf.freq2den(self.nX['2bands']*1e9)

        phase_dif=self.gd_k*1e6/self.Dt_DF['K']
        self.phi=integrate.cumtrapz(phase_dif,initial=0)
        #39.5e-2*self.nX['K']/c
        phase_difa=self.gd_ka*1e6/self.Dt_DF['Ka']
        self.phia=integrate.cumtrapz(phase_difa)+self.phi[self.cmin]

        np=p.interp(self.nX['Ka'][:self.cmax], self.nX['K'][self.cmin:], self.phi[self.cmin:])
        npp=(np+self.phia[:self.cmax])/2
        phase=p.concatenate((self.phi[:self.cmin],npp,self.phia[self.cmax:]))
        pos=self.find_pos(self.nX['2bands'],phase)

        p.plot(self.nX['2bands'],phase)
        p.plot(self.nX['Ka'][:self.cmax],npp, marker='o', linestyle='-', color='r')
        p.plot(self.nX['K'],self.phi, marker='o', linestyle='-', color='b')
        p.plot(self.nX['Ka'],self.phia, marker='o', linestyle='-', color='k')
        # p.figure()
        # p.plot(pos,self.ne)

if __name__ == "__main__":
    pass
    shot.reference_gd()
    shot.plasma_gd(6666,2)
    shot.overlap_gd()
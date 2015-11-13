import pylab as p
from scipy import integrate

phi_cal=0
factor=2/3

def refrac(w0,w1):
    #zero=wp2[abs(wp2-ww).argmin()]
    return p.sqrt(1-pow(w0/w1,2))

def area_trap(f0,r0,f1,r1):
    return (r1+r0)*p.pi*f1*refrac(f0,f1)

def phase_eval(phase_dif):
    phase=integrate.cumtrapz(phase_dif)-phi_cal
    return phase

def phi_0(n,freqs_prob,pos):
    area=0
    method=1
    if method==1:    
        area=integrate.simps(freqs_prob[n]*refrac(freqs_prob[n],freqs_prob[:n-1]))
    elif method==2:
        for f in range(1,n):
            area+=area_trap(freqs_prob[f-1],pos[f-1],freqs_prob[f],pos[f])
    return area

def pos_new(r0,n,freqs_prob,pos,factor,phase):
    alt=2*p.pi*freqs_prob[n]*refrac(freqs_prob[n-1],freqs_prob[n])
    area_prev=phi_0(n,freqs_prob,pos)
    return r0+factor*(phase[n]-area_prev+p.pi/2)/alt
"""
Function: Process data for sweep frequency mode.
Author: Cassio Amador (cassioamador at yahoo.com.br)
Thanks for G. Ronchi for debugging
TODO: zero padding in filter should check first the size of the signal,
so it could could choose best value depending on signal size. Maybe the
same for padding in spectrogram.
Spectrogram could use a pre-sized matrix.
Some functions could be modified to work with hopping frequency.
"""

import numpy as np

from read_signal import ReadSignal

class ProcSweep(ReadSignal):
    """Process sweep data. All data is evaluated at its current sweep."""
    def __init__(self,shot,tipo='data'):
        #Inherit from ReadSignal class
        ReadSignal.__init__(self,shot,tipo)
        #multiplying factor for channel frequencies
        self.chan_factor={"K":2,"Ka":3,"ref":1,"time":1}
        #evaluate size of points per sweep
        self.sweep_size=self.rate*self.sweep_dur

    def mark_sweep_points(self):
        """Creates a list of points (called 'points') where the sweep
        started, with the 'time' channel (channel 4)."""
        self.read_channel('time')
        #find min value for data, in the first 2 sweeps
        mindata=self.bindata['time'][:self.rate*2*(self.sweep_dur+self.interv_sweep)].min()
        #find all points with this value
        zeros=np.where(self.bindata['time']==mindata)[0]
        #mark only the first minimum for each trigger
        singlezeros=np.where(np.diff(zeros)>50)[0]
        self.points=np.append(zeros[0],zeros[singlezeros+1])
        print "Total number of sweeps: %s" % len(self.points)

    def time2sweep(self,time):
        """Converts a position in time (ms) to the correspondent
        sweep position."""
        if not hasattr(self,'points'):
            self.mark_sweep_points()
        #finds nearest sweep from a specified time
        self.sweep_cur=(abs(self.points-time*self.rate*1e3)).argmin()
        return self.sweep_cur
    
    def sweep2time(self,sweep):
        """Converts a sweep position to its position in time (ms)."""
        return sweep*(self.sweep_dur+self.interv_sweep)*1e-3

    def read_single_sweep(self,channel,sweep_cur=100):
        """Reads data for the current sweep (sweep_cur), for an
        specific channel"""
        self.sweep_cur=sweep_cur
        if not hasattr(self,'points'):
            self.mark_sweep_points()
        samples_ini=self.points[self.sweep_cur]
        samples_end=samples_ini+self.sweep_size
        if channel not in self.bindata.keys():
            self.read_channel(channel)
        if hasattr(self,'singlesweep_data')==False:
            self.singlesweep_data={}
        self.singlesweep_data[channel]=self.bindata[channel][samples_ini:samples_end]

    def single_sweep_time(self,channel,time=30):
        """same as read_single_sweep, but it reads int time (ms)"""
        self.time2sweep(time)
        self.read_single_sweep(channel,self.sweep_cur)

    def signal_filter(self,channel,freqs=(1e3,15e3)):
        """Filter signal from specific channel. A FFT is performed in
        the signal, and the result is multiplied by a window function
        (kaiser function), which nulls the undesired beating
        frequencies, that are outside of 'freqs'. The signal is then
        recreated with an IFFT."""
        #zero padding size:
        zer_pad_filter=4
        #alias for sweep_size
        N=int(self.sweep_size)
        #FFT with zero padding
        fft=np.fft.rfft(self.singlesweep_data[channel],zer_pad_filter*N)
        #bp=fft[:] #used for other plot functions
        fmin,fmax=freqs
        #creates the beating frequency axis, and finds the position of
        #the frequency limits in the axis
        fft_freq=np.linspace(0,self.rate*1e3/2.,num=zer_pad_filter*N/2)
        cmin=(abs(fft_freq-fmin)).argmin()
        cmax=(abs(fft_freq-fmax)).argmin()
        #creates window function for filter. Second argument of kaiser
        #function must be float
        window_func=np.concatenate((np.zeros(cmin+1),
                       np.kaiser(cmax-cmin,2.),np.zeros(zer_pad_filter*N/2-cmax)))
        #multiply the window by the signal's FFT.
        bp=np.multiply(fft,window_func)
        #Since the signal is REAL, we use IRFFT, and takes only the
        #first half, since the other half is symmetric.
        newsig=np.fft.irfft(bp)[:N]
        return newsig

    def plot_sweep(self,channel):
        """Plot binary data for an especific channel."""
        import pylab as p
        
        if hasattr(self,'sweepfreq'):
            #dict with arrays for the frequencies in each channel.
            self.sweep_freq={}
        if channel not in self.sweep_freq.keys():
            self.sweep_freq[channel]=np.linspace(self.freq_start,self.freq_end,num=self.sweep_size)*self.chan_factor[channel]
        p.plot(self.sweep_freq[channel],self.sweep_bindata[channel],label="Channel: %s" % channel)
        p.xlabel("freq (GHz)")
        p.ylabel("beating signal")

    def spectrogram(self,channel,window_scale=4,step_scale=16,log=0,
                    group_delay=1,ploti=1,figure=0,normal=1,filtered=1,
                    freqs=(1e3,15e3)):
        """Evaluates and plot spectrogram (SFFT) of beating signal.
        Some parameters listed (others can be found in the function):
        group_delay=1 evaluates group delay.
                    0 evaluates beating frequency
        ploti=2 for contour plot
              1 for animation. Don't try if running trough network.
        normal=1 normalize spectrum
        log=0 
            1 for log spectrogram
        """
        if ploti>0:
            import pylab as p

        #scale for zero padding
        zer_pad=4
        #alias for sweep_size
        N=int(self.sweep_size)
        #SFFT window size
        window=N/window_scale
        #SFFT step size, by 
        step=window/step_scale
        #creates the window function that will 'move' trough the signal
        #to evaluate each FFT.
        window_func=np.concatenate((np.zeros(N-window/2),np.hanning(window),np.zeros(N-window/2)))
        if filtered==1:
            sig=self.signal_filter(channel,freqs)
        else:
            sig=self.singlesweep_data[channel]
        if ploti==2:
            p.ion()
            a=raw_input('press any key to begin plot. 0 to cancel:   ')
            if a=='0':
                ploti=1
                p.ioff()
        #creates a matrix to receive the spectra
        matrix=np.empty(shape=(1+(N-window)/step,2*N))
        #loops trough all the possible windows, and evaluates the FFT.
        for i in range(1+(N-window)/step):
            t=i*step+window/2
            new_window=window_func[N-t:N+N-t]
            new_sig=np.multiply(sig,new_window)
            #We ignore the fft's first point (DC component).
            fft_sig=np.fft.rfft(new_sig,zer_pad*N)[1:]
            if log==1:
                fft_sig=np.log(fft_sig)
            if normal==1:
                fft_sig *= (1./fft_sig.max())
            #Tries to join fft arrays in a matrix.
            #If it is the first one, it creates the matrix.
            matrix[i]=abs(fft_sig)
            if ploti==2:
                p.figure(5)
                p.clf()
                p.subplot(2,1,1)
                p.plot(abs(fft_sig))
                p.subplot(2,1,2)
                p.plot(sig)
                p.plot(new_sig)
                p.twinx()
                p.plot(new_window,'r-')
                p.draw()
                #raw_input('')
        if ploti==2:
            p.ioff()
        #transpose matrix for spectrogram.
        matrix=matrix.transpose()
        #creates arrays with beating frequency and band frequency.
        X=np.linspace(self.freq_start,self.freq_end,num=len(matrix[0]))*self.chan_factor[channel]
        #Y is the frequency, in MHz
        Y=np.linspace(0,self.rate/2.,num=zer_pad*N/2)
        if group_delay==1:
            #group delay in ns
            Y*=self.sweep_dur/(max(X)-min(X))
        #print N,len(matrix),len(matrix[0]),len(X),len(Y)
        if ploti>=1:
            if figure==0:
                p.figure()
            else:
                p.figure(figure)
            ctr=p.contourf(X,Y,matrix)
            #creates colorbar
            #colorbar=p.colorbar(ctr, shrink=0.8, extend='both',format='%.3g')
            p.ylabel('beating freq (MHz)')
            if group_delay==1:
                p.ylabel('group delay (ns)')
            p.xlabel('sweep freq (GHz)')
            if (filtered==1) & (group_delay==0):
                p.ylim(freqs[0]*1e-3,freqs[1]*1e-3)
            #fig.subplots_adjust(left=0.15)
            #fig.subplots_adjust(bottom=0.15)
            #p.title(channel)
        return matrix,X,Y

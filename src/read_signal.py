"""
Function: Read info file for a shot and binary data from specific 
channel from Reflectometry. Check data in the following order: local 
folder defined in shots_loc.py, MDSPlus database from TCABR internal
network, than MDSPlus from public machine "tcabrcl".
Can plot channel data, downsampled.
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at gmail.com)
TODO: clean and fix names.
Make it possible to store locally data from MDSPlus, if wanted.
Write specific unities for info read.
"""

from os.path import isfile,join
from sys import exc_info
import numpy as np
import sqlite3

try:
    import MDSplus as mds
except ImportError:    
    print("No python module for MDSplus.")

import shots_folder

class ReadSignal():
    def __init__(self,shot='0',acq_mode='data'):
        if shot=='0':
            print("Select a shot number")
            return
        self.shot_folder=shots_folder.find_folder(acq_mode)
        self.acq_mode=acq_mode
        self.shot=shot
        self.__readMDSplus__=False
        self.__LocalMDSplusSever__=False
        self.bindata={}
        self.read_info()

    def read_info(self):
        ''' Reads reflectometer configuration.
            Try to read local file, or from MDSplus server (local or
            external).'''
        #dict of frequency modes abbreviations
        mode_names={'sf':'sweep freq.',
                    'ff':'fixed freq.','hf':'hopping freq.',}
        #checks if 'info' file is stored locally. If not, activates
        # MDSPlus mode.
        if isfile(join(self.shot_folder,'%s_info.dat' % self.shot)):
            infofile=open(join(self.shot_folder,'%s_info.dat' % self.shot))
            infolines=infofile.readlines()
            infofile.close()
            #frequency mode abbreviation
            self.mode=infolines[0].rstrip('\n')
            #full frequency mode name
            self.mode_name=mode_names[self.mode]
            #date of acquisition
            self.date=infolines[1].rstrip('\n')
            #create a shot parameter (class attribute) for each line in
            # 'info' file.
            for line in infolines[2:]:
                line=line.rstrip('\n')
                #if there is a '[', it is a list of frequencies. 
                if '[' in line:
                    info1=line.split(':')[0]
                    info2=line.split('[')[1].split(']')[0].split(',')
                    setattr(self,info1,str(info2))
                else:
                    info1=line.split(':')[0]
                    if info1=='sweep':
                        info1='sweep_dur'
                    info2=line.split(':')[1]
                    setattr(self,info1,float(info2))
            print("Read from local file")
        else:
            self.__readMDSplus__= True
            try:
                # Local MDSplus Server
                tree = mds.Tree('tcabr_ref', self.shot)
                mds_get=tree.getNode
                self.__LocalMDSplusSever__= True
            except:
                # External MDSplus Server
                conn = mds.Connection('tcabrcl.if.usp.br')
                try:
                    conn.openTree('tcabr_ref', self.shot)
                except mds._mdsshr.MdsException:
                    print("No reflectometry data available for this shot.")
                    exit()
                mds_get=conn.get
            #Common parameters
            self.mode=mds_get('\\REFPARAMETER.REFMODE').data().strip()
            self.mode_name=mode_names[self.mode]
            self.angle=mds_get('\\REFPARAMETER.ANGLE').data() # degrees
            self.rate=np.int(mds_get('\\REFPARAMETER.RATE').data())/1e6   # MHz
            self.time_dur=1e-3*mds_get('\\REFPARAMETER.SAMPLES').data()/self.rate
            if self.mode=='ff':
                #print("Fixed Freq. enable")
                self.freq=mds_get('\\FIXEDFREQ.FREQUENCY').data()
            elif self.mode=='hf':
                #print("Hopping Freq. enable")
                self.freq=mds_get('\\HOPPINGFREQ.FREQ_TABLE').data()           # GHz
                self.restart_time=mds_get('\\HOPPINGFREQ.RESTART_TIME').data() # us
                self.time_step=mds_get('\\HOPPINGFREQ.TIME_STEP').data()       # us
            elif self.mode=='sf':    
                #print('Sweep Freq. enable')
                self.freq_start=np.float(mds_get('\\SWEEPFREQ.FREQ_START').data())  # GHz
                self.freq_end=np.float(mds_get('\\SWEEPFREQ.FREQ_END').data())       # GHz
                self.sweep_dur=np.int(mds_get('\\SWEEPFREQ.SWEEP_TIME').data())          # us
                self.interv_sweep=np.int(mds_get('\\SWEEPFREQ.INTERV_SWEEP').data()) # us
            if self.__LocalMDSplusSever__:
                print("Read from local MDSplus server")
            else:
                conn.closeAllTrees()
                print("Read from external MDSplus server")
        if hasattr(self,'mode')==False:
            print("No reflectometry data available for this shot.")
            return

    def read_channel(self,chan):
        '''Read bin data from the channel specified.'''
        #dict to exchange between number and name of channel.
        channels={1:'K',2:'Ka',3:'ref',4:'time'}
        bands={'K':1,'Ka':2,'ref':3,'time':4}
        try:
            channel=channels[chan]
        except KeyError:
            channel=chan
            chan=bands[channel]
        if self.__readMDSplus__:
            signals={1:'\\KBAND.SIGNAL',2:'\\KABAND.SIGNAL',
                        3:'\\MIRNOV.SIGNAL',4:'\\TRIGGER.SIGNAL'}
            try:
                # Local MDSplus Server
                tree = mds.Tree('tcabr_ref', self.shot)
                self.bindata[channel] = tree.getNode(signals[chan]).data()
            except:
                # External MDSplus Server
                conn = mds.Connection('tcabrcl.if.usp.br')
                conn.openTree('tcabr_ref', self.shot)
                self.bindata[channel]= np.array(conn.get(signals[chan]))
                conn.closeAllTrees()
        else:
            self.arq=join(self.shot_folder,'%s_%s.bin' % (self.shot,chan))
            self.bindata[channel]=np.fromfile(self.arq,np.int16)
            print("binary file: %s\n#Samples: %s" % (self.arq,len(self.bindata[channel])))
        self.datasize=len(self.bindata[channel])

#    def save_channel(self,chan):

    def plot_channel(self,chan,factor=200):
        """Plots channel data, but downsampled by a 'factor' because of
        memory issues."""
        import pylab as p
        self.times=np.arange(0,self.datasize)/(factor*self.rate*1e3)
        p.plot(self.times,self.bindata[chan][::factor])

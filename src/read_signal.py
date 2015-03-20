"""
Function: Read info file for a shot and binary data from specific
channel from Reflectometry. Check data in the following order: local
folder defined in shots_folder.py, MDSPlus database from TCABR internal
network, than MDSPlus from public machine "tcabrcl".
Can plot channel data, downsampled.
Authors:
    Cassio Amador (cassioamador at yahoo.com.br)
    Gilson Ronchi (gronchi at if.usp.br)
TODO: clean and fix names.
Write specific unities for info read.
"""

from os import path
from sys import exit
import numpy as np

try:
    import MDSplus as mds
except ImportError:
    print("No python module for MDSplus.")

import shots_folder

class ReadSignal():

    def __init__(self, shot='0', acq_mode='data'):
        if shot == '0':
            print("Select a shot number")
            return
        self.shot_folder = shots_folder.find_folder(acq_mode)
        self.acq_mode = acq_mode
        self.shot = shot
        self.__readMDSplus__ = False
        self.__LocalMDSplusSever__ = False
        self.bindata = {}
        self.read_info()

    def read_info(self):
        """ Reads reflectometer configuration.
            Try to read local file, or from MDSplus server (local or
            external)."""
        # dict of frequency modes abbreviations
        mode_names = {'sf': 'sweep freq.',
                      'ff': 'fixed freq.', 'hf': 'hopping freq.', }
        # checks if 'info' file is stored locally. If not, activates
        # MDSPlus mode.
        if path.isfile(path.join(self.shot_folder, '%s_info.dat' % self.shot)):
            infofile = open(path.join(self.shot_folder, '%s_info.dat' % self.shot))
            infolines = infofile.readlines()
            infofile.close()
            # frequency mode abbreviation
            self.mode = infolines[0].rstrip('\n')
            # full frequency mode name
            self.mode_name = mode_names[self.mode]
            # date of acquisition
            self.date = infolines[1].rstrip('\n')
            # create a shot parameter (class attribute) for each line in
            # 'info' file.
            for line in infolines[2:]:
                line = line.rstrip('\n')
                # if there is a '[', it is a list of frequencies.
                if '[' in line:
                    info1 = line.split(':')[0]
                    info2 = line.split('[')[1].split(']')[0].split(',')
                    setattr(self, info1, str(info2))
                else:
                    info1 = line.split(':')[0]
                    if info1 == 'sweep':
                        info1 = 'sweep_dur'
                    info2 = line.split(':')[1]
                    setattr(self, info1, float(info2))
            print("Read from local file")
        else:
            self.__readMDSplus__ = True
            try:
                # Local MDSplus Server
                tree = mds.Tree('tcabr_ref', self.shot)
                mds_get = tree.getNode
                self.__LocalMDSplusSever__ = True
            except:
                # External MDSplus Server
                conn = mds.Connection('tcabrcl.if.usp.br')
                try:
                    conn.openTree('tcabr_ref', self.shot)
                except mds._mdsshr.MdsException:
                    print("No reflectometry data available for this shot.")
                    exit()
                mds_get = conn.get
            # Common parameters
            self.mode = mds_get('\\REFPARAMETER.REFMODE').data().decode('utf-8')
            self.mode_name = mode_names[self.mode]
            self.angle = mds_get('\\REFPARAMETER.ANGLE').data()  # degrees
            self.rate = np.int(mds_get('\\REFPARAMETER.RATE').data()) / 1e6  # MHz
            self.time_dur = 1e-3 * mds_get('\\REFPARAMETER.SAMPLES').data() / self.rate
            if self.mode == 'ff':
                #print("Fixed Freq. enable")
                self.freq = mds_get('\\FIXEDFREQ.FREQUENCY').data()
            elif self.mode == 'hf':
                #print("Hopping Freq. enable")
                self.freq = mds_get('\\HOPPINGFREQ.FREQ_TABLE').data()            # GHz
                self.restart_time = mds_get('\\HOPPINGFREQ.RESTART_TIME').data()  # us
                self.time_step = mds_get('\\HOPPINGFREQ.TIME_STEP').data()        # us
            elif self.mode == 'sf':
                #print('Sweep Freq. enable')
                self.freq_start = np.float(mds_get('\\SWEEPFREQ.FREQ_START').data())    # GHz
                self.freq_end = np.float(mds_get('\\SWEEPFREQ.FREQ_END').data())        # GHz
                self.sweep_dur = np.int(mds_get('\\SWEEPFREQ.SWEEP_TIME').data())       # us
                self.interv_sweep = np.int(mds_get('\\SWEEPFREQ.INTERV_SWEEP').data())  # us
            if self.__LocalMDSplusSever__:
                print("Read from local MDSplus server")
            else:
                conn.closeAllTrees()
                print("Read from external MDSplus server")
        if hasattr(self, 'mode') is False:
            print("No reflectometry data available for this shot.")
            return

    def read_channel(self, chan,save_locally=1):
        """Try to read bin data from the channel specified, first locally, then from MDSPlus"""
        # dict to exchange between number and name of channel.
        channels = {1: 'K', 2: 'Ka', 3: 'ref', 4: 'time'}
        bands = {'K': 1, 'Ka': 2, 'ref': 3, 'time': 4}
        try:
            channel = channels[chan]
        except KeyError:
            channel = chan
            chan = bands[channel]
        if (chan==4) & (self.mode=="ff"):
            print "no channel 4, data for fixed frequency shots"
            return
        try:
            self.arq = path.join(self.shot_folder, '%s_%s.bin' % (self.shot, chan))
            self.bindata[channel] = np.fromfile(self.arq, np.int16)
            print("binary file: %s\n#Samples: %s" % (self.arq, len(self.bindata[channel])))
        except IOError:
            print "no local data, connecting to MDSPlus"
            signals = {1: '\\KBAND.SIGNAL', 2: '\\KABAND.SIGNAL',
                       3: '\\MIRNOV.SIGNAL', 4: '\\TRIGGER.SIGNAL'}
            try:
                # Local MDSplus Server
                tree = mds.Tree('tcabr_ref', self.shot)
                self.bindata[channel] = tree.getNode(signals[chan]).data()
            except:
                # External MDSplus Server
                conn = mds.Connection('tcabrcl.if.usp.br')
                conn.openTree('tcabr_ref', self.shot)
                self.bindata[channel] = np.array(conn.get(signals[chan]))
                conn.closeAllTrees()
            if save_locally==1:
                self.save_channel(chan)
                self.save_info_file()
        self.datasize = len(self.bindata[channel])

    def save_channel(self, chan):
        channels = {1: 'K', 2: 'Ka', 3: 'ref', 4: 'time'}
        arq = path.join(self.shot_folder, '%s_%s.bin' % (self.shot, chan))
        if not path.isfile(arq):
            self.bindata[channels[chan]].tofile(arq)

    def save_info_file(self):
        arq = path.join(self.shot_folder, '%s_info.dat' % self.shot)
        if not path.isfile(arq):
            arq = open(arq, 'w')
            arq.write('%s\n' % (self.mode))
            arq.write('00:00:00 January 01 1970\n')  # fake timestamp
            if self.mode == 'ff':
                arq.write('freq: %s\n' % (self.freq))
            elif self.mode == 'hf':
                arq.write('freq: %s\n' % (self.freq))
                arq.write('restart_time: %s\n' % (self.restart_time))
                arq.write('time_step: %s\n' % (self.time_step))
            elif self.mode == 'sf':
                arq.write('sweep: %s\n' % (self.sweep_dur))
                arq.write('interv_sweep: %s\n' % (self.interv_sweep))
                arq.write('freq_start: %s\n' % (self.freq_start))
                arq.write('freq_end: %s\n' % (self.freq_end))
            arq.write('angle: %s\n' % (self.angle))
            arq.write('rate: %s\n' % (self.rate))
            arq.write('time_dur: %s\n' % (1e3 * self.time_dur))
            arq.close()

    def plot_channel(self, chan, factor=200):
        """Plots channel data, but downsampled by a 'factor' because of
        memory issues."""
        import pylab as p
        self.times = np.arange(0, 1 + self.datasize, factor) / (self.rate * 1e3)
        #print(len(self.times), len(self.bindata[chan][::factor]))
        p.plot(self.times, self.bindata[chan][::factor])

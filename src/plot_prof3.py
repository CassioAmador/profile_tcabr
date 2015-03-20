import Pmw
import Tkinter as tk
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from os import chdir,listdir,path,getcwd

import read_signal as rs

class PlotWindow():
    def __init__(self,shot_number=0,text_box=0):
        self.shot_number=shot_number
        self.st=text_box
        self.proc_folder=path.join(getcwd(), "..", "PROC")
        chdir(self.proc_folder)
        self.read_data()
        self.plot_profiles()
    
    def changeshot(self,shot_number):
        if shot_number!=self.shot_number:
            self.shot_number=shot_number
            self.read_data()
            self.plot_update(0)

    def read_data(self):
        shotfile=rs.ReadSignal(self.shot_number)
        self.st.clear()
        info={}
        self.st.appendtext("\nAcq. Rate: %s" % shotfile.rate)
        self.st.appendtext("\nK band: %s" % shotfile.freq_end)
        self.st.appendtext("\nSweep_dur: %s " % shotfile.sweep_dur)
        self.st.appendtext("\nInterv sweep: %s" % shotfile.interv_sweep)
        #self.info=shotfile.info
        #print shotfile.info
        prof_folder=path.join(self.proc_folder,"%s" % self.shot_number,"level_0")
        chdir(prof_folder)
        self.ne=np.loadtxt('ne.dat')
        self.prof_info=np.loadtxt('prof_info.dat')
        print self.prof_info

        prof_list=listdir(prof_folder)
        prof_list.sort()
        self.position=np.empty(shape=((len(prof_list)-2),len(self.ne)))
        self.times=np.empty(len(prof_list)-2)

        i=0
        for r_file in prof_list:
            name=r_file.strip('.dat')
            if name not in ('prof_info','ne'):
                self.position[i]=np.loadtxt(r_file)*1e2
                self.times[i]=name
                i+=1        
    
    def plot_profiles(self):
        self.fig, self.ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)
        self.l, = plt.plot(self.position[0],self.ne, lw=2, color='blue')
        plt.axis([0, 20, self.ne[0],self.ne[-1]])
        self.ax.legend(loc='upper left')
        self.ax.set_xlabel('r (cm)')
        self.ax.set_ylabel('density (10^19 m^-3)')
        self.ax.set_title("# %s" % self.shot_number)

        axcolor = 'lightgoldenrodyellow'
        self.axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

        self.stime = Slider(self.axfreq, 'time', self.prof_info[1], self.prof_info[2], valinit=self.prof_info[1]+(self.prof_info[2]-self.prof_info[1])*0.2)

        self.stime.on_changed(self.plot_update)

    def plot_update(self,val):
        time = self.stime.val
        i=(abs(time*1e3-self.times)).argmin()
        self.l.set_xdata(self.position[i])
        self.fig.canvas.draw_idle()



class PlotProf():
    def __init__(self,root):
        self.root=root
        self.draw()
    def draw(self):
        self.infobox(70)
        self.st.appendtext('Input shot')
        self.entries_caption=['Shot:']
        self.entries_default=['32214']
        self.keys_entries=['shot']

        self.number_entries = len(self.keys_entries)
        self.entries = {}
        for e in range(self.number_entries):
            self.makeentry(self.keys_entries[e], self.entries_caption[e], self.entries_default[e])
        Pmw.alignlabels(self.entries.values())

        self.buttons()
        self.entries['shot'].configure(command=self.choose_shot)
        self.entries['shot'].focus_set()
        self.plot_frame=tk.Frame()
        self.plot_frame.pack()

    def infobox(self,height=150):
        # Create the ScrolledText with headers.
        #'Helvetica', 'Times', 'Fixed', 'Courier' or 'Typewriter'
#        fixedFont = Pmw.logicalfont('Fixed',size=12)
        self.st = Pmw.ScrolledText(self.root,
                borderframe = 1,
                usehullsize = 1,
                hull_height = height,
                text_wrap='word',
#                text_font = fixedFont,
                text_padx = 4,
                text_pady = 4,
        )
        self.st.pack(side='bottom', padx = 5, pady = 5, fill = 'both', expand = 1)
        # Prevent users' modifying text
        self.st.configure(text_state = 'disabled')

    def buttons(self):
        # Create the button box
        self.buttonBox = Pmw.ButtonBox(self.root,labelpos = 'nw')
        self.buttonBox.pack(fill = 'x', expand = 1, padx = 3, pady = 3)

        # Add buttons to the ButtonBox.
        self.buttonBox.add('About', command=self.about_gui)
        self.buttonBox.add('Close', command = self.close)
        self.buttonBox.add('Set', command = self.choose_shot)

        # Make all the buttons the same width.
        self.buttonBox.alignbuttons()

    def makeentry(self,key,caption, default):
        self.entries[key]=Pmw.EntryField(self.root,
            labelpos = 'w',
            label_text = caption,
            value= default)
        self.entries[key].pack(side= 'top', fill= 'x', expand=1, padx=10, pady=5)

    def choose_shot(self):
        self.shot_number=int(self.entries['shot'].get())
        self.st.clear()
        self.st.appendtext('Reflectometry parameters:')
        if hasattr(self,'plotwindow'):
            self.plotwindow.changeshot(self.shot_number)
        else:
            self.plotwindow=PlotWindow(self.shot_number,self.st)

    def close(self):
        self.root.destroy()

    def about_gui(self):
        Pmw.aboutversion('1.0\n Mar 15 2015')
        Pmw.aboutcopyright('Author: Cassio H. S. Amador')
        Pmw.aboutcontact(
            'For more informations/bug reporting:\n' +
            '  email: cassioamador@yahoo.com.br'
        )
        self.about = Pmw.AboutDialog(self.root, applicationname = 'Ref Setup')
        self.about.withdraw()
        self.about.show()

if __name__== '__main__':
    root = tk.Tk()
    root.title('Ref Profiles')
    Pmw.initialise(fontScheme='pmw2')
    plotprof=PlotProf(root)
    plt.ion()
    root.mainloop()
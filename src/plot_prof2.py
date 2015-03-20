import Pmw
import Tkinter as tk
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from os import chdir,listdir,path,getcwd

class PlotWindow():
    def __init__(self,shot=0,plot_frame=0):
        self.shot=shot
        self.frame=plot_frame
        self.plot_profiles()
    def changeshot(self,shot):
        self.shot=shot
    def plot_profiles(self):
        prof_folder=path.join(getcwd(), "..", "PROC","%s" % self.shot,"level_0")
        chdir(prof_folder)
        ne=np.loadtxt('ne.dat')
        info=np.loadtxt('prof_info.dat')

        prof_list=listdir(prof_folder)
        prof_list.sort()
        position=np.empty(shape=((len(prof_list)-2),len(ne)))
        times=np.empty(len(prof_list)-2)

        i=0
        for r_file in prof_list:
            name=r_file.strip('.dat')
            if name not in ('prof_info','ne'):
                position[i]=np.loadtxt(r_file)*1e2
                times[i]=name
                i+=1

        self.fig=plt.figure()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.get_tk_widget().grid(column=0,row=1)
        ax = self.fig.add_subplot(111)
        plt.subplots_adjust(left=0.25, bottom=0.25)
        l, = ax.plot(position[0],ne, lw=2, color='blue')
        ax.axis([0, 20, ne[0],ne[-1]])
        ax.legend(loc='upper left')
        ax.set_xlabel('r (cm)')
        ax.set_ylabel('density (10^19 m^-3)')
        ax.set_title("# %s" % self.shot)

        axcolor = 'lightgoldenrodyellow'
        axfreq = mpl.axes.Axes(self.fig,[0.25, 0.4, 0.65, 0.03], axisbg=axcolor)

        stime = Slider(axfreq, 'time', info[1], info[2], valinit=info[1]+(info[2]-info[1])*0.2)

        def update(val):
            time = stime.val
            i=(abs(time*1e3-times)).argmin()
            l.set_xdata(position[i])
            self.fig.canvas.draw_idle()
        stime.on_changed(update)

        plt.show()


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
        self.shot=int(self.entries['shot'].get())
        self.st.clear()
        self.st.appendtext('Reflectometry parameters:')
        if hasattr(self,'plotwindow'):
            self.plotwindow.changeshot(self.shot)
        else:
            self.plotwindow=PlotWindow(self.shot,self.plot_frame)

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
    root.mainloop()
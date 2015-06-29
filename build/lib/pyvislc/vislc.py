from Tkinter import *
import Tkconstants, tkFileDialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.colors import Normalize, BoundaryNorm, LogNorm
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import sys, lsp, argparse, os, re, defaults
import numpy as np
#from scipy.signal import lombscargle   # < -- can't get this to work properly...
import pandas as pd 
from math import *
from textwrap import *
import readhatlc as rhlc
from processing import *
from utils import *

cmaps = [ 'jet', 'binary', 'hot' ]

default_custom_plot_params = {
    'plot-type' : 'scatter',
    'cmap' : 'jet',
    '3rd-dim' : None,
    'xlabel' : 'X',
    'ylabel' : 'Y',
    'zlabel' : 'Z',
    'figtitle' : 'title',
    'figsize' : (4, 4),
}
SCATTER_MARKER = '.'

#categorical_types = [ '' ]

class NavigationToolbar(NavigationToolbar2TkAgg):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar2TkAgg.toolitems if
                 t[0] in ('Home', 'Pan', 'Zoom', 'Save')]

def ax_lomb_scargle(ax, periods, powers, opts):
    """ Adds a lomb-scargle plot to a matplotlib.Axis instance 

        inputs:
            ax          matplotlib.Axis instance to which a LS plot will be added
            periods     A set of periods (1-d array)
            powers      A set of corresponding Lomb-Scargle powers
            opts        Options

    """

    ax.set_title("Lomb-Scargle")
    ax.set_ylabel("Power")
    ax.set_xlabel("Period (days)")
    #ax.set_yscale('log')
    ax.set_xscale('log')

    data1 = ax.plot(periods, powers, color='k')
    if 'window_lsp' in opts: window_data = ax.plot(opts['window_lsp']['periods'], opts['window_lsp']['powers'],color='r',alpha=0.5) 
    else: window_data = None
    '''
    ls_opts = {}
    for o in opts:
        ls_opts[o] = opts[o]
    for o in default_ls_opts:
        if o not in ls_opts:
            ls_opts[o] = default_ls_opts[o]
    '''
    ax.set_xlim(min(periods), max(periods))
    #ax.scatter(opts['periods_to_label'], opts['powers_to_label'], 
    #    color='r', alpha=defaults.settings['alpha'], marker=SCATTER_MARKER)
    labellist = []
    labelvals = []
    for i in range(0,len(opts['periods_to_label'])):
        ax.plot(opts['periods_to_label'][i], opts['powers_to_label'][i], 'o', color='r', alpha = defaults.settings['alpha'])
        labellist.append("%d"%(i+1))
        labelvals.append("%.2e days"%(opts['periods_to_label'][i]))
        period = opts['periods_to_label'][i]
        power = opts['powers_to_label'][i]
        ax.text(period, power,"%d"%(i+1), fontsize=12, ha='center', va='bottom')
    if len(labellist) > 0:
        txt = ""
        for i in range(0,len(opts['periods_to_label'])):
            if not txt == "": txt = "%s\n"%(txt)
            txt = "%s%s: %s"%(txt, labellist[i], labelvals[i])
        txt = "Best periods\n%s"%(txt)
        ax.text(0.95, 0.95, txt,fontsize=10, transform=ax.transAxes, ha='right', va='top', bbox=dict(facecolor='red', alpha=0.3))
        #ax.legend(labellist, labelvals, loc='best', fontsize=10)
    ax.format_coord = lambda x, y : ''
def ax_phase_folded(ax, t, y, opts):
    """ Adds a phase-folded plot to a matplotlib.Axis instance

        inputs:
            ax          matplotlib.Axis instance to which phase-folded plot will be added
            t           Central values of phase bins
            y           Average magnitude in each phase bin
            opts        Options (required, should contain 'ylabel', 'period', 'yerr' )

    """
    invert_y = True
    if 'invert' in opts: invert_y = opts['invert']
    if defaults.settings['bin-phase']:
        ax.plot(t, y, label="P=%.4f days"%(opts['period']),color='b',lw=2)
        ax.fill_between(t, y-opts['yerr'],y+opts['yerr'],facecolor='b',alpha=defaults.settings['alpha'])
    else:
        ax.scatter(t, y, label="P=%.4f days"%(opts['period']),c='b'
            ,alpha=defaults.settings['alpha'], marker=SCATTER_MARKER)
    
    ax.set_ylabel(opts['ylabel'])
    ax.set_xlabel("Phase")
    ax.set_title("Phase-folded")
    ax.set_xlim(0,1)
    ax.legend(loc='upper right', fontsize=12)
    if invert_y: ax.invert_yaxis()

    ax.format_coord = lambda x, y : ''
def ax_raw_lc(ax, t, y, opts):
    """ Adds a lomb-scargle plot to a matplotlib.Axis instance 

        inputs:
            ax          matplotlib.Axis instance to which a LS plot will be added
            t           Observation times
            y           Magnitude measurement at each time
            opts        Options (required, should contain 'xlabel', 'ylabel')

    """
    invert_y = True
    if 'invert' in opts: invert_y = opts['invert']
    ax.scatter(t, y, alpha=defaults.settings['alpha'],marker=SCATTER_MARKER)
    ax.set_ylabel(opts['ylabel'])
    ax.set_xlabel(opts['xlabel'])
    ax.set_xlim(min(t),max(t))
    ax.set_title("Raw lightcurve")
    if invert_y: ax.invert_yaxis()
def ax_default(ax, t, y, opts):
    ax.scatter(t, y, alpha=defaults.settings['alpha'],marker=SCATTER_MARKER)
def ax_custom(ax, t, y, opts={}):
    # Are we dealing with numbers or categories?
    names = ['t', 'y']
    arrs = {
        't' : t,
        'y' : y
    }
    types = {}
    cats = {}

    for p in default_custom_plot_params:
        if not p in opts:
            opts[p] = default_custom_plot_params[p]

    for n in arrs:
        types[n] = data_kind(arrs[n])
        if types[n] == "categorical":
            cats[n], arrs[n] = encode(arrs[n])

    
    if opts['plot-type'] == 'scatter':
        if not opts['3rd-dim'] is None:
            z = opts['3rd-dim']
            z_type = data_kind(opts['3rd-dim'])
            if z_type == "categorical":
                z_cat, z = encode(z)
            norm = Normalize(min(z), max(z))
            cbar_obj = ax.scatter(arrs['t'],arrs['y'],c=z,
                cmap=opts['cmap'], norm=norm, alpha=defaults.settings['alpha'],
                marker = SCATTER_MARKER)
            plt.colorbar(cbar_obj)
        else:
            ax.scatter(arrs['t'],arrs['y'],alpha=defaults.settings['alpha'],
                marker=SCATTER_MARKER)
    else:
        ax.plot(arrs['t'],arrs['y'])
    if types['t'] == "categorical":
        ax.set_xticks(np.arange(0,len(cats)))
        ax.set_xticklabels(cats['t'])
    if types['y'] == "categorical":
        ax.set_yticks(np.arange(0,len(cats)))
        ax.set_yticklabels(cats['y'])

    ax.set_ylabel(opts['ylabel'])
    ax.set_xlabel(opts['xlabel'])
def open_help_window():
    window = Toplevel()
    window.wm_title("Abbreviations")

    frame = Frame(window)
    scrollbar = Scrollbar(frame, orient='vertical')
    lbox = Listbox(frame, width=50,height=100,yscrollcommand=scrollbar.set)
    scrollbar.config(command=lbox.yview)
    abbrs = []
    descs = []
    for col in sorted(rhlc.TEXTLC_OUTPUT_COLUMNS):
        abbrs.append(col)
        descs.append(rhlc.TEXTLC_OUTPUT_COLUMNS[col][0])
    

    for i in range(0,len(abbrs)):
        lbox.insert("end","%-10s = %-100s"%(abbrs[i], descs[i]))
        #abbrs[i].grid(row=i+1,column=0, sticky = "W")
        #Label(frame,text=" = ", font=defaults.settings['value-font']).grid(row=i+1,column = 1,sticky="NEWS")
        #descs[i].grid(row=i+1,column=2, sticky = "W")
    scrollbar.pack(side="left", fill="y")
    lbox.pack(side="right", fill="both", expand=1)
    frame.pack()
def open_custom_plot_window(lc):
    window = Toplevel()
    window.wm_title("Custom plot for %s"%(lc['hatid']))

    figure_frame = Frame(window)
    button_frame = Frame(window)
    title = Label(button_frame, text="Plotting options", font=defaults.settings['label-font'], 
            bg=defaults.settings['label-bg'])
    plot_types = [ 'scatter', 'line' ]
    var_types = [ v for v in rhlc.TEXTLC_OUTPUT_COLUMNS if v in lc ]
    #print defaults.settings['label-font']

    variables = {
        'x_type' : StringVar(),
        'y_type' : StringVar(),
        'z_type' : StringVar(),
        'cmap' : StringVar(),
        'plot_type' : StringVar()
    }
    variables['x_type'].set(defaults.settings['time-type'])
    variables['y_type'].set('%s%s'%(defaults.settings['mag-type'],defaults.settings['aperture'] ))
    variables['z_type'].set(defaults.settings['na-str'])
    variables['cmap'].set('jet')
    variables['plot_type'].set('scatter')
    allowed_values = {
        'x_type' : var_types,
        'y_type' : var_types,
        'z_type' : var_types + [defaults.settings['na-str']],
        'cmap' : cmaps,
        'plot_type' : plot_types
    }
    labels = {
        'plot_type' : Label(button_frame,text="Plot:", font = defaults.settings['label-font']),
        'x_type' : Label(button_frame,text="x=", font = defaults.settings['label-font']),
        'y_type' : Label(button_frame,text="y=", font = defaults.settings['label-font']),
        'z_type' : Label(button_frame,text="z=", font = defaults.settings['label-font']),
        'cmap' : Label(button_frame,text="Colormap:", font = defaults.settings['label-font'])
    }
    def set_and_replot():
        values = {}
        for v in variables:
            values[v] = str(variables[v].get())
        popts = {
            'cmap' : values['cmap'],
            'xlabel' : values['x_type'],
            'ylabel' : values['y_type'],
            'zlabel' : values['z_type'],
            'plot-type' : values['plot_type']
        }
        if values['z_type'] not in var_types:
            popts['3rd-dim'] = None
        else:
            popts['3rd-dim'] = lc[values['z_type']]
        t, x = lc[values['x_type']], lc[values['y_type']]
        popts['plot'] = ax_custom
        contents['figure_canvas'].replot(t,x,options=popts)

    contents = {
        'figure_canvas' : PlotFrame(figure_frame, [], []),
        'plot_button' : Button(button_frame,text="Plot", command=set_and_replot)
    }
    for v in variables:
        contents[v] = OptionMenu(button_frame,variables[v],*(allowed_values[v]))
        
    rcs = {
        'x_type' : (1, 1, 1, 1, "W"),
        'y_type' : (2, 1, 1, 1, "W"),
        'z_type' : (3, 1, 1, 1, "W"),
        'cmap' : (4, 1, 1, 1, "W"),
        'plot_type' : (5, 1, 1, 1, "W"),
        'plot_button' : (6, 1, 1, 1, "NEWS")
    }

    
    for c in contents:
        if c == 'figure_canvas' : continue
        if c in labels:
            labels[c].grid(row=rcs[c][0], column=(rcs[c][1] - 1), rowspan=rcs[c][2], 
                            columnspan=rcs[c][3], sticky="E")
        contents[c].grid(row=rcs[c][0], column=rcs[c][1], rowspan=rcs[c][2], 
                            columnspan=rcs[c][3], sticky=rcs[c][4])
        
    
    title.grid(row=0, column=0, rowspan=1, columnspan=2, sticky="NEWS")
    figure_frame.grid(row=0,column=0, sticky="NEWS")
    button_frame.grid(row=0,column=1, sticky="N")
XPOS = "250"
YPOS = "80"
class FileList:
    def __init__(self, parent):
        self.parent = parent
        self.window = Toplevel()
        self.window.geometry("+0+%s"%(YPOS))
        self.window.wm_title("File list")
        self.frame = Frame(self.window)
        
        
        self.buttonframe = Frame(self.window)
        self.title = Label(self.buttonframe,text="File list",font=defaults.settings['label-font'],bg=defaults.settings['label-bg'])
        self.scrollbar = Scrollbar(self.frame, orient='vertical')
        self.lbox = Listbox(self.frame, width=30,height=40,yscrollcommand=self.scrollbar.set, font=defaults.settings['file-list-font'])
        self.scrollbar.config(command=self.lbox.yview)
        self.lbox.bind("<Double-Button-1>", self.select_hid)
        self.flags_to_show=None
        self.set_lbox()
        self.set_buttons()
        
        
        self.select_all()
        self.title.grid(row=0,column=0,columnspan=3, sticky="NEWS")
        self.scrollbar.grid(row=1, column=0)
        self.lbox.grid(row=1,column=1)
        self.buttonframe.grid(row=0, column=0)
        self.frame.grid(row=1,column=0)
    def set_buttons(self):
        self.selection_buttons = {   
            'Select All' : Button(self.buttonframe, text="Select All", command=self.select_all),
            'Flagged' : Button(self.buttonframe, text="Flagged", command=self.select_flagged),
            'Unselect All' : Button(self.buttonframe, text="Unselect All", command=self.unselect_all),
        }
        self.ckboxes = {}
        self.ckboxes_vars = {}
        for flag in self.all_flags:
            self.ckboxes_vars[flag] = IntVar()
            self.ckboxes[flag] = Checkbutton(self.buttonframe, text=flag, variable=self.ckboxes_vars[flag])

        self.apply = Button(self.buttonframe, text='Apply', command=self.apply_selection)
        self.num_items_var = StringVar()
        self.num_items_var.set('%d files listed'%(len(self.hids_shown)))
        self.num_items = Label(self.buttonframe, textvariable=self.num_items_var)
        dim = 2
        rownum = 0
        for i,sb in enumerate(self.selection_buttons):
            self.selection_buttons[sb].grid(row=i/dim + 1, column=i%dim, sticky="W")
            rownum = i/dim + 1

        for i,cb in enumerate(self.ckboxes):
            self.ckboxes[cb].grid(row=rownum + 1 + i/dim, column=i%dim, sticky="W")
            rownum2 = rownum + 1 + i/dim
        rownum = rownum2 
        self.apply.grid(row=rownum + 1, column=0, columnspan = 1)
        self.num_items.grid(row=rownum + 1, column = 1, columnspan = 1)

    def apply_selection(self):
        flags_to_show = []
        for flag in self.ckboxes:
            is_selected = (self.ckboxes_vars[flag].get() == 1)
            if is_selected: flags_to_show.append(flag)
        if len(flags_to_show) == 0:
            print "Select at least one category"
        else:
            self.flags_to_show = flags_to_show
            self.lbox.delete("0","end")
            self.set_lbox()

    def select_all(self):
        for flag in self.ckboxes:
            self.ckboxes[flag].select()

    def select_flagged(self):
        for flag in self.ckboxes:
            if flag == 'Unflagged': self.ckboxes[flag].deselect()
            else: self.ckboxes[flag].select()

    def unselect_all(self):
        for flag in self.ckboxes:
            self.ckboxes[flag].deselect()
    def get_prev_lc_index(self):
        inds = []
        for hid in self.hids_shown:
            if self.indices_dict[hid] < self.parent.lc_index: 
                inds.append(self.indices_dict[hid])
        if len(inds) == 0:
            inds = []
            for hid in self.hids_shown:
                if self.indices_dict[hid] > self.parent.lc_index: 
                    inds.append(self.indices_dict[hid])
            if len(inds) == 0: return None
            else: return max(inds)
        return max(inds)
        
    def get_next_lc_index(self):
        inds = []
        for hid in self.hids_shown:
            if self.indices_dict[hid] > self.parent.lc_index: 
                inds.append(self.indices_dict[hid])
        if len(inds) == 0:
            inds = []
            for hid in self.hids_shown:
                if self.indices_dict[hid] < self.parent.lc_index: 
                    inds.append(self.indices_dict[hid])
            if len(inds) == 0: return None
            else: return min(inds)
        return min(inds)
    def set_lbox(self):
        lc_files = self.parent.lc_files

        self.hids = []
        self.indices_dict = {}
        self.flag_dict = {}
        self.hidre = re.compile("(HAT-[0-9]{3}-[0-9]+)")
        for i,lcfile in enumerate(lc_files):
            result = self.hidre.search(lcfile)
            hid = result.groups(1)[0]
            self.hids.append(hid)
            self.indices_dict[hid] = i
            if hid in self.parent.log_info:
                self.flag_dict[hid] = self.parent.log_info[hid]
            else:
                self.flag_dict[hid] = "Unflagged"
        self.all_flags = np.unique([ self.flag_dict[hid] for hid in self.flag_dict  ])
        if self.flags_to_show is None:
            self.flags_to_show = self.all_flags
        

        self.hids_shown = []
        for hid in self.hids:
            if self.flag_dict[hid] in self.flags_to_show:
                self.hids_shown.append(hid)
                self.lbox.insert("end","%-15s %-15s"%(hid, self.flag_dict[hid]))
        for i,hid in enumerate(self.hids_shown):
            if self.flag_dict[hid] is not "Unflagged":
                self.lbox.itemconfig(i, {'fg' : 'forest green'})
        for hid in self.indices_dict:
            if self.indices_dict[hid] == self.parent.lc_index:
                if hid not in self.hids_shown:
                    self.parent.lc_index = self.indices_dict[self.hids_shown[0]]
                    self.selection_index = 0
                    self.parent.save()
                    self.parent.set()

                else:
                    self.selection_index = self.hids_shown.index(hid)
        self.lbox.selection_anchor(self.selection_index)
        self.lbox.selection_set(self.selection_index)


    def update(self):
        self.lbox.delete("0","end")
        self.set_lbox()
        self.num_items_var.set('%d files listed'%(len(self.hids_shown)))
        #self.set_buttons()
        
    def select_hid(self, event):
        selection = self.lbox.curselection()
        selection = int(selection[0])
        listentry = self.lbox.get(selection)
        hid = self.hidre.search(listentry).groups(1)[0]
        #print hid, " selected."
        self.parent.lc_index = self.indices_dict[hid]
        self.parent.save()
        self.parent.set()
class PlotFrame:
    def __init__(self, frame,t,y, options=None, figsize=(4,4)):
        #frame = Frame(master)
        self.frame = frame
        self.options = options
        self.fig = Figure(figsize=figsize,tight_layout={'pad' : 1.5})
        self.ax = self.fig.add_subplot(111)
        if self.options is None or not 'plot' in self.options: 
            ax_default(self.ax, t, y, self.options)
        else:
            self.fig_data = self.options['plot'](self.ax, t, y, self.options)
        
        self.canvas = FigureCanvasTkAgg(self.fig,master=self.frame)
        
        #self.toolbar_canvas = Frame(self.frame)
        NavigationToolbar2TkAgg(self.canvas, self.frame)
        self.canvas._tkcanvas.config(background="#FFFFFF", borderwidth=0, highlightthickness=0) 
        #self.canvas.get_tk_widget().grid(row=0,column=0,sticky="NEWS")
        self.canvas.get_tk_widget().pack()
        #self.toolbar_canvas.grid(row=1,column=0,sticky="W")
        self.canvas.show()
    def replot(self,t,y, options=None):
       
        self.ax.clear()
        if not (options is None): self.options = options

        if self.options is None or not ( 'plot' in self.options):
            ax_default(self.ax, t, y, self.options)
        else:
            if 'invert' not in self.options: self.options['invert'] = False
            self.options['plot'](self.ax, t, y, self.options)
        self.canvas.draw()  
class InfoFrame:
    def __init__(self,frame,parent):
        self.labelfg = 'black'
        self.parent = parent
        self.frame = frame
        self.title_var = StringVar()
        self.title_var.set(self.parent.lc['hatid'])
        self.title = Label(self.frame, textvariable=self.title_var, font=defaults.settings['label-font'],
            bg=defaults.settings['label-bg'])
        # Custom foreground colors
        self.fgs = {
            'hatid' : 'red',
            'flag' : 'forest green'
        }
        #self.labelbg = 'white'
        self.full_names = {
            'twomassid' : '2MASS ID', 'flag' : 'Note',
            'radec' : 'Ra, Dec',
            'vmag' : 'V','rmag' : 'R','jmag' : 'J',
            'imag' : 'I','hmag' : 'H','kmag' : 'K'
        }
        self.mag_index = {
            'vmag' : 0,'rmag' : 1,'jmag' : 3,
            'imag' : 2,'hmag' : 4,'kmag' : 5
        }
        self.set_custom_data()

        # Set uncustomized foreground colors to the default foreground color
        for n in self.full_names:
            if not n in self.fgs: self.fgs[n] = self.labelfg

        # Initialize non-magnitude strings
        self.strlist = [ fn for fn in self.full_names if not "mag" in fn ]
        self.strs = {}
        for s in self.strlist:
            self.strs[s] = StringVar()

        # Initialize magnitude strings
        self.maglist = [ fn for fn in self.full_names if "mag" in fn]

        self.mag_strs = {}
        for m in self.maglist:
            self.mag_strs[m] = StringVar()

        rowdim = 4
        last_row = (len(self.strlist) + 1)/rowdim

        rows = [ r%rowdim + 1 for r in range(len(self.strlist)) ]
        columns = [ 2*(r/rowdim) for r in range(len(self.strlist)) ]
        colspans = [ 1 + 5*((r/rowdim)/last_row) for r in range(len(self.strlist)) ]

        mag_rows = [ len(rows)%rowdim + 1 for i in range(len(self.maglist)) ]
        mag_columns = [ 2*(len(rows)/rowdim) + r for r in range(len(self.maglist)) ]
        mag_colspans = [ 1 for r in range(len(self.maglist)) ]

        total_width = max(mag_columns) + 2

        self.objects = []
        self.placement = {
            'row' : [],
            'column' : [],
            'columnspan' : [],
            'sticky' : []
        }

        # Write labels for non-magnitude data
        columnspan = 1
        for i,s in enumerate(self.strlist):
            rowno = rows[i]
            colno = columns[i]
            
            self.strs[s].set(defaults.settings['na-str'])
            self.objects.append(Label(self.frame,text=self.full_names[s], fg=self.labelfg, font=defaults.settings['label-font']))
            self.placement['row'].append(rowno)
            self.placement['column'].append(colno)
            self.placement['columnspan'].append(columnspan)
            self.placement['sticky'].append("W")

            
            
        # Write labels for magnitudes
        rowno = mag_rows[0]
        colno = mag_columns[0]

        self.objects.append(Label(self.frame,text="Mags", fg = self.labelfg, font=defaults.settings['label-font']))
        
        self.placement['row'].append(rowno)
        self.placement['column'].append(colno)
        self.placement['columnspan'].append(columnspan)
        self.placement['sticky'].append("W")

        for i,m in enumerate([ 'vmag', 'rmag', 'imag', 'jmag', 'hmag', 'kmag']):
            rowno = mag_rows[i]
            colno = mag_columns[i] + 1
            
            self.mag_strs[m].set(defaults.settings['na-str'])
            self.objects.append(Label(self.frame,text = self.full_names[m],fg=self.labelfg, font=defaults.settings['label-font']))

            self.placement['row'].append(rowno)
            self.placement['column'].append(colno)
            self.placement['columnspan'].append(columnspan)
            self.placement['sticky'].append("NEWS")

        self.objects.append(self.title)
        self.placement['row'].append(0)
        self.placement['column'].append(0)
        self.placement['columnspan'].append(total_width)
        self.placement['sticky'].append("NEWS")
        
        # Write values for non-magnitude data
        self.labels = {}

        for i,s in enumerate(self.strlist):
            rowno = rows[i]
            colno = columns[i] + 1
            columnspan = colspans[i]

            self.labels[s] = Label(self.frame,textvariable=self.strs[s],fg=self.fgs[s], font=defaults.settings['value-font'])
            self.objects.append(self.labels[s])
            self.placement['row'].append(rowno)
            self.placement['column'].append(colno)
            self.placement['columnspan'].append(columnspan)
            self.placement['sticky'].append("E")


        # Write values for magnitude data
        columnspan = 1
        for i,m in enumerate([ 'vmag', 'rmag', 'imag', 'jmag', 'hmag', 'kmag']):
            rowno = mag_rows[i] + 1
            colno = mag_columns[i] + 1

            self.labels[m] = Label(self.frame,textvariable=self.mag_strs[m],fg=self.fgs[m], font=defaults.settings['value-font'])
            self.objects.append(self.labels[m])
            self.placement['row'].append(rowno)
            self.placement['column'].append(colno)
            self.placement['columnspan'].append(columnspan)
            self.placement['sticky'].append("W")
        self.set_grid()
    def set_custom_data(self):
        if 'ra' in self.parent.lc:
            self.parent.lc['radec'] = "%.5f, %.5f"%(self.parent.lc['ra'], self.parent.lc['dec'])
        if 'FLT' in self.parent.lc:
            self.full_names['filters'] = "Filters"
            self.parent.lc['filters'] = ""
            for filt in np.unique(self.parent.lc['FLT']):
                nobs = len([ 1 for fil in self.parent.lc['FLT'] if fil == filt ])
                self.parent.lc['filters'] = "%s %s (%d)"%(self.parent.lc['filters'], filt, nobs)
        if 'NET' in self.parent.lc:
            self.full_names['nets'] = "Network"
            self.parent.lc['nets'] = ""
            for entry in np.unique(self.parent.lc['NET']):
                nobs = len([ 1 for this_entry in self.parent.lc['NET'] if str(this_entry) == str(entry) ])
                self.parent.lc['nets'] = "%s %s (%d)"%(self.parent.lc['nets'], str(entry), nobs)
        if 'STF' in self.parent.lc:
            self.full_names['stfs'] = "Station"
            self.parent.lc['stfs'] = ""
            for entry in np.unique(self.parent.lc['STF']):
                nobs = len([ 1 for this_entry in self.parent.lc['STF'] if str(this_entry) == str(entry) ])
                self.parent.lc['stfs'] = "%s %s (%d)"%(self.parent.lc['stfs'], str(entry), nobs)
    def unset_grid(self):
        for o in self.objects:
            o.grid_forget()

    def set_grid(self):
        for i in range(len(self.objects)):
            row, column = self.placement['row'][i], self.placement['column'][i]
            columnspan, sticky =  self.placement['columnspan'][i], self.placement['sticky'][i]
            self.objects[i].grid(row=row, column=column, columnspan=columnspan, sticky=sticky)

    def update(self):
        lc = self.parent.lc
        if lc is None: return
        if self.parent.lc_index is None: return
        self.title_var.set(self.parent.lc['hatid'])
        self.unset_grid()
        self.set_custom_data()

        # Load magnitude values from lightcurve
        self.mag_values = {}
        for i,m in enumerate(self.maglist):
        
            self.mag_values[m] = lc['mags'][self.mag_index[m]]
        
        # Set values for non-magnitude data
        for i,s in enumerate(self.strlist):
            if s == 'flag':
                if lc['hatid'] in self.parent.log_info:
                    self.strs[s].set(self.parent.log_info[lc['hatid']])
                else:
                    self.strs[s].set(defaults.settings['na-str'])
            elif not lc[s] is None:
                self.strs[s].set(str(lc[s]))
        
        # Set values for magnitude data
        for i,m in enumerate(self.maglist):
            if not self.mag_values[m] is None or np.isnan(self.mag_values[m]):
                self.mag_strs[m].set(str(self.mag_values[m]))

        self.set_grid()
class CommentBox:
    def __init__(self,  visualizer_instance, frame=None):
        self.visualizer_instance = visualizer_instance
        if frame is None:
            self.window = Toplevel()
            self.window.wm_title("Comments for %s"%(self.visualizer_instance.lc['hatid']))
            self.frame = Frame(self.window)
        else:
            self.frame = frame
            self.window = None
        self.title = Label(self.frame,text="Comments",font=defaults.settings['label-font'],bg=defaults.settings['label-bg'])
        self.scrollbar = Scrollbar(self.frame, orient='vertical')
        self.textbox = Text(self.frame, wrap="word",width=40, height=5, yscrollcommand=self.scrollbar.set,state="normal")
        
        self.scrollbar.config(command=self.textbox.yview)
        self.textbox.insert('1.0',self.visualizer_instance.comments)
        self.submit = Button(self.frame, text="Submit", command=self.submit)
        
        self.title.grid(row=0,column=0, columnspan=3, sticky="NEWS")
        self.scrollbar.grid(row=1, column=0)
        self.textbox.grid(row=1, column=1)
        self.submit.grid(row=1, column=2,columnspan=1,sticky="NW")
        self.frame.grid(row=0,column=0,sticky="NW")
    def update(self):
        #self.visualizer_instance = visualizer_instance
        if self.window is not None:
            self.window.wm_title("Comments for %s"%(self.visualizer_instance.lc['hatid']))
        self.textbox.delete('1.0','end')
        self.textbox.insert('1.0', self.visualizer_instance.comments)
    def submit(self):
        comments = self.textbox.get('1.0','end')
        print comments
        self.visualizer_instance.comments = comments
        self.visualizer_instance.save()
class Toolbox:
    def __init__(self,frame,parent, more_buttons=None):
        
        self.parent = parent
        self.dphase = defaults.settings['dphase']

        # Set variables

        # P: string value of the period
        self.P = StringVar()
        
        # phase: phase offset; 0 to 1
        self.phase = StringVar()

        # dperiod: amount by which you can finely tune the period
        self.dperiod = StringVar()
        self.dperiod_value = defaults.settings['dperiod']
        
        # mag_type: type of reduced magnitude to use
        self.mag_type = StringVar()
        self.mag_choices = [ "TF", "IM", "RM", "EP" ]

        # aperture_num: HAT aperture number to use
        self.aperture_num = StringVar()
        self.aperture_choices = [ "1", "2", "3" ]
        
        # nbins: number of bins to use for the phase-folded lightcurve
        self.nbins = StringVar()
        self.nbins.set(str(defaults.settings['nbins']))

        # Set up a menu
        self.top_menu = Menu(self.parent.master)
        self.top_menu.add_command(label="Open", command=self.openlc)
        self.top_menu.add_command(label="Close", command=self.closelc)
        self.top_menu.add_command(label="Exit", command=sys.exit)

        # Define the contents of this toolbox
        self.contents = {
            'options_title' : Label(frame,text="Options",font=defaults.settings['label-font'],bg=defaults.settings['label-bg']),
            'period_entry' : Entry(frame,textvariable=self.P,width=defaults.settings['entry-width']),
            'enter_period' : Button(frame,text="Set period", command=self.set_period),
            'double_period' : Button(frame,text="2 x P", command=self.double_period),
            'halve_period' : Button(frame,text="1/2 x P", command=self.halve_period),
            'inc_period' : Button(frame,text="> P", command=self.increase_period),
            'dec_period' : Button(frame,text="< P", command=self.decrease_period),
            'entry_dt' : Entry(frame, textvariable=self.dperiod,width=defaults.settings['entry-width']),
            'enter_dt' : Button(frame, text="Set dP",command=self.set_dperiod),
            'find_period' : Button(frame,text="Find period", command=self.find_period),
            'next_lightcurve' : Button(frame,text="> LC", command=self.up_lc),
            'last_lightcurve' : Button(frame,text="< LC", command=self.down_lc),
            'mag_type_label'    : Label(frame, text="Mag type:"),
            'mag_type'         : OptionMenu(frame,self.mag_type,*(self.mag_choices)),
            'aperture_label'    : Label(frame, text="Aperture:"),
            'aperture_num'    : OptionMenu(frame,self.aperture_num,*(self.aperture_choices)),
            'entry_nbins' : Entry(frame,textvariable=self.nbins, width=defaults.settings['entry-width']),
            'enter_nbins' : Button(frame,text="Set nbins",command=self.set_nbins),
            'entry_phase' : Entry(frame,textvariable=self.phase, width=defaults.settings['entry-width']),
            'enter_phase' : Button(frame,text="Set phase",command=self.set_phase),
            'inc_phase' : Button(frame,text="> Phase", command=self.inc_phase),
            'dec_phase' : Button(frame,text="< Phase", command=self.dec_phase),
            'custom_plot' : Button(frame,text="Custom plot", 
                command= lambda : open_custom_plot_window(self.parent.lc)),
            'file_list' : Button(frame,text="Files", command=self.open_file_list),
            'abbrs' : Button(frame,text="Abbrs", command=open_help_window),
            'save_plot' :  Button(frame,text="Save plot", command=self.parent.saveplot),
            #'comments' : Button(frame,text="Comments", command=self.open_comment_box)
            #'open_button': Button(frame,text="Open", command=self.openlc)
        }
        # Set defaults
        self.mag_type.set(defaults.settings['mag-type'])
        self.aperture_num.set(defaults.settings['aperture'])
        self.dperiod.set(str(self.dperiod_value))
        self.P.set(str(parent.period))
        self.phase.set("0")

        # Functions that handle changes to the mag type or aperture
        def change_aperture(*args):
            ap_num = self.aperture_num.get()
            self.parent.options['aperture'] = int(ap_num)
            self.parent.set()
        def change_mag_type(*args):
            m_type = self.mag_type.get()
            self.parent.options['magtype'] = m_type
            self.parent.set()

        # Link aperture/magnitude changes to function calls
        self.aperture_num.trace('w', change_aperture)
        self.mag_type.trace('w', change_mag_type)

        # (row, column, rowspan, columnspan, sticky) values for each element
        self.rcs = {
            'options_title'     : (0, 0, 1, 4, "NEWS"),
            'period_entry'      : (1, 0, 1, 1, "E"),
            'enter_period'      : (1, 1, 1, 1, "W"),
            'entry_nbins'       : (1, 2, 1, 1, "E"),
            'enter_nbins'       : (1, 3, 1, 1, "W"),
            'dec_period'        : (2, 0, 1, 1, "E"),
            'inc_period'        : (2, 1, 1, 1, "W"),
            'entry_dt'          : (2, 2, 1, 1, "E"),
            'enter_dt'          : (2, 3, 1, 1, "W"),
            'halve_period'      : (3, 0, 1, 1, "NEWS"),
            'double_period'     : (3, 1, 1, 1, "NEWS"),
            'entry_phase'       : (3, 2, 1, 1, "E"),
            'enter_phase'       : (3, 3, 1, 1, "W"),
            'find_period'       : (4, 0, 1, 2, "NEWS"),
            'dec_phase'         : (4, 2, 1, 1, "E"),
            'inc_phase'         : (4, 3, 1, 1, "W"),
            'mag_type_label'    : (5, 0, 1, 1, "E"),
            'mag_type'          : (5, 1, 1, 1, "W"),
            'aperture_label'    : (6, 0, 1, 1, "E"),
            'aperture_num'      : (6, 1, 1, 1, "W"),
            'last_lightcurve'   : (5, 2, 2, 1, "E"),
            'next_lightcurve'   : (5, 3, 2, 1, "W"),
            'custom_plot'       : (7, 2, 1, 1, "NEWS"),
            #'open_button'       : (7, 0, 1, 1, "NEWS"),
            'abbrs'             : (7, 1, 1, 1, "NEWS"),
            #'comments'          : (8, 0, 1, 1, "NEWS"),
            'save_plot'         : (7, 0, 1, 1, "NEWS"),
            'file_list'         : (7, 3, 1, 1, "NEWS")
            #''
        }
        last_row = max( [ self.rcs[n][0] for n in self.rcs ])
        tot_cols = max( [ self.rcs[n][1] for n in self.rcs]) + 1
        def keypress(event):
            key = event.char
            if key in self.parent.flag_shortcuts:
                self.flag_lc(self.parent.flag_shortcuts[key])

        self.parent.master.bind('<Shift-Right>', self.up_lc)
        self.parent.master.bind('<Shift-Left>', self.down_lc)
        self.parent.master.bind('<Key>', keypress)

        if len(self.parent.flags) > 0:
            self.contents['flag_title'] = Label(frame,text="Flags",font=defaults.settings['label-font'],bg=defaults.settings['label-bg'])
            self.rcs['flag_title'] = (last_row + 1, 0, 1, tot_cols,"NEWS")
            last_row += 1
            for i,name in enumerate(self.parent.flags):
                self.contents[name] = Button(frame, text=name, command=lambda : self.flag_lc(name))
                self.rcs[name] = (last_row+2+i/tot_cols, i%tot_cols, 1, 1, "NEWS")
        # Place contents in their correct location
        for i,c in enumerate(self.contents):
            self.contents[c].grid(row=self.rcs[c][0], 
                column=self.rcs[c][1], rowspan=self.rcs[c][2], columnspan=self.rcs[c][3], sticky=self.rcs[c][4])

        self.parent.master.config(menu=self.top_menu)

        # Options for the file-loading menu
        self.file_opt = {
            'title' : "Load HAT lightcurve"
        }
    def flag_lc(self, flag):
        #print "Added %s to parent"%(flag)
        hid = self.parent.lc['hatid']
        self.parent.log_info[hid] = flag
        self.parent.info_box.update()
    def open_file_list(self):
        if isclosed( self.parent.file_list_window.window ): 
            self.parent.file_list_window = FileList(self.parent)

    def open_comment_box(self):
        if isclosed( self.parent.commentbox.window ): 
            self.parent.commentbox = CommentBox(self.parent)

    def inc_phase(self, event=None):
        self.parent.phase_offset += self.dphase
        self.phase.set(self.parent.phase_offset)
        self.parent.set_phase_folded()

    # Decrease phase
    def dec_phase(self, event=None):
        self.parent.phase_offset -= self.dphase
        self.phase.set(self.parent.phase_offset)
        self.parent.set_phase_folded()

    # Set custom phase
    def set_phase(self, event=None):
        _phase = nums(self.phase.get())
        if not _phase is None: 
            self.parent.phase_offset = _phase

            self.parent.set_phase_folded()

    # Set custom number of phase-folding bins
    def set_nbins(self, event=None):
        _nbins = nums(self.nbins.get())
        if not _nbins is None: 
            self.parent.nbins=int(_nbins)
            self.parent.set_phase_folded()

    # Double/halve/increase/decrease period
    def double_period(self, event=None):
        period = 2*self.parent.period
        self.P.set(str(period))
        self.parent.set_phase_folded(period = period)
    def halve_period(self, event=None):
        period = 0.5*self.parent.period
        self.P.set(str(period))
        self.parent.set_phase_folded(period = period)
    def increase_period(self, event=None):
        self.parent.period += self.dperiod_value
        self.P.set(str(self.parent.period))
        self.parent.set_phase_folded()
    def decrease_period(self, event=None):
        period = self.parent.period - self.dperiod_value
        if not (period < 0): 
            self.parent.period = period
            self.P.set(str(self.parent.period))
            self.parent.set_phase_folded()

    def update_button_values(self):
        self.P.set(str(self.parent.period))
        self.phase.set(str(self.parent.phase_offset))
        self.dperiod.set(str(self.dperiod_value))
        self.nbins.set(str(self.parent.nbins))


    # Set custom period
    def set_period(self, event=None):
        period = nums(self.P.get())
        if not period is None:
            self.P.set(str(period))
            self.parent.period = period
            self.parent.set_phase_folded()
        else:
            print "'",period,"' is not a valid period"

    # Set dperiod
    def set_dperiod(self, event=None):
        _dperiod = nums(self.dperiod.get())
        if not _dperiod is None:
            self.dperiod_value = _dperiod

    # Toggle next lightcurve
    def up_lc(self, event=None):
        
        self.parent.next_lightcurve()
        self.find_period()

    # Toggle previous lightcurve
    def down_lc(self, event=None):
        
        self.parent.last_lightcurve()
        self.find_period()
    # Pick period with largest Lomb-Scargle power
    def find_period(self, event=None):
        self.parent.period = self.parent.best_period
        self.parent.set_phase_folded()
        self.P.set(str(self.parent.period))

    # Open a lightcurve
    def openlc(self, event=None):
        fname = tkFileDialog.askopenfilename( **self.file_opt)
        if not fname is None:
            self.parent.add_lc(fname)

    # Close a lightcurve (does nothing at the moment)
    def closelc(self):
        print "closed lightcurve"


def isclosed(window):
    try:
        state = window.state()
    except:
        return True
    return False


class Visualizer:
    def __init__(self, master, lcs=None, lc_index=0, logfile=None, flag_shortcuts=None, flags=None):
        self.comments = ""

        if logfile is None:
            self.logfile = defaults.settings['logfile']
        else:
            self.logfile = logfile
        if lcs is None:
            lcs = []
        #self.filters = {
        #    'IQ1' : 
        #}
        frame = Frame(master)
        self.initialized=False

        self.options = {
            'aperture' : int(defaults.settings['aperture']),
            'magtype' : defaults.settings['mag-type'],
            'time' : defaults.settings['time-type']
        }

        self.nbins = defaults.settings['nbins']
        self.phase_offset = 0.0
        self.master = master
        self.flag_shortcuts = flag_shortcuts
        self.flags = flags

        if self.flag_shortcuts is None: self.flag_shortcuts = {}
        if self.flags is None: self.flags = []

        self.lc_files = lcs
        self.lc_index = lc_index
        self.log_info = {}
        if not self.logfile is None:
            if os.path.exists(self.logfile):
                if os.path.getsize(self.logfile):
                    try:
                        log_info = np.loadtxt(self.logfile,dtype=dt_log)
                        if hasattr(log_info,'__iter__'):
                            for i,hid in enumerate(log_info['hatid']):
                                self.log_info[hid] = log_info[i]['flag']
                        else:
                            self.log_info[log_info['hatid']] = log_info['flag']
                    except:
                        print "Can't read contents of %s"%(self.logfile)

        lcs_flagged = []
        lcs_unflagged = []
        for lc in self.lc_files:
            if lc in self.log_info:
                lcs_flagged.append(lc)
            else:
                lcs_unflagged.append(lc)
        self.lc_files = lcs_unflagged 
        self.lc_files.extend(lcs_flagged)
        self.file_list_window = FileList(self)
        
        # Set up frames
        self.raw_plot_frame = Frame(master,bd=1)#,relief='sunken')
        self.phase_plot_frame = Frame(master,bd=1)#,relief='sunken')
        self.lsp_plot_frame = Frame(master,bd=1)#,relief='sunken')
        self.info_box_frame = Frame(master,bd=1)#,relief='sunken')
        self.toolbox_frame = Frame(master,bd=1)#,relief='raised')
        self.comment_frame = Frame(master, bd=1)    
        self.raw_plot_figsize = (12.5,3)
        if len(self.lc_files) > 0: 
            # Initialize everything, load + process first lightcurve
            self.set()
             # Set frame contents
        if len(self.lc_files) > 0:
            self.raw_plot = PlotFrame( self.raw_plot_frame, self.t, self.y, options = self.opts_raw, figsize=self.raw_plot_figsize)
            self.phase_plot = PlotFrame( self.phase_plot_frame, self.phases, self.phase_mags, options = self.opts_phase_folded)
            self.lsp_plot = PlotFrame( self.lsp_plot_frame, self.periods, self.lsp_powers, options = self.opts_ls)
            self.info_box = InfoFrame(self.info_box_frame,self)
            self.toolbox = Toolbox(self.toolbox_frame,self)
            self.commentbox = CommentBox(self, frame=self.comment_frame)
            self.info_box.update()
        else:
            self.period = 0
            self.lc = None
            self.initialized = True
            self.raw_plot = PlotFrame( self.raw_plot_frame, [], [], figsize=self.raw_plot_figsize)
            self.phase_plot = PlotFrame( self.phase_plot_frame, [], [])
            self.lsp_plot = PlotFrame( self.lsp_plot_frame, [],[])
            self.info_box = InfoFrame(self.info_box_frame,self)
            self.toolbox = Toolbox(self.toolbox_frame,self)
            self.info_box.update()
        
        # Place frames

        self.info_box_frame.grid(row=0,column=0, columnspan=2, sticky="NEWS")
        self.comment_frame.grid(row=0,column=2, columnspan=1, sticky="NEWS")
        self.raw_plot_frame.grid(row=1,column=0, columnspan=3, sticky="NEWS")
        
        self.lsp_plot_frame.grid(row=2,column=0, sticky="NEWS")
        self.phase_plot_frame.grid(row=2,column=1, sticky="NEWS")
        self.toolbox_frame.grid(row=2,column=2, rowspan=2, sticky='NEWS')
        

        #self.plot_frames = [  self.phase_plot_frame, self.lsp_plot_frame, self.raw_plot_frame]

    def set_phase_folded(self, period=None):
        if not period is None:
            self.period = period
        elif not self.initialized:
            self.period = self.best_period
        self.phases, self.phase_mags, errs = phase_fold(self.t,self.y,self.period,nbins=self.nbins, phase_offset=self.phase_offset)
        self.opts_phase_folded = { 'plot' : ax_phase_folded, 'yerr' : errs, 'period' : self.period, 'ylabel' : self.ytype}
        if self.initialized:
            self.phase_plot.replot( self.phases, self.phase_mags, options = self.opts_phase_folded)

    def set_raw(self):
        self.ttype = self.options['time']
        self.ytype = "%s%d"%(self.options['magtype'],self.options['aperture'])
        self.y = self.lc[self.ytype]
        self.t = self.lc[self.ttype]
        self.filts = self.lc['FLT']
        ok_inds = [ i for i in range(len(self.y)) if not isnan(self.y[i]) ]
        self.t, self.y = filter_normalization(self.t[ok_inds], self.y[ok_inds],self.filts[ok_inds])
        
        self.opts_raw = { 'ylabel' : self.ytype, 'xlabel' : self.ttype, 'plot' : ax_raw_lc }
        if self.initialized:
            self.raw_plot.replot( self.t, self.y, options = self.opts_raw)
    def set_lsp(self):
        self.periods, self.lsp_powers, self.best_period = get_lombscarg(self.t, self.y)
        #window = np.random.normal(1., 1E-5, len(self.t))
        #self.window_lsp_periods, self.window_lsp_powers, junk = get_lombscarg(self.t, window)
        
        #self.frequencies, self.lsp_powers = scaled_lsp(self.t, self.y)
        #self.periods = np.power(self.frequencies, -1)
        #self.best_period = self.periods[0]
        #self.highest_lsp_power = self.lsp_powers[0]
        #for i in range(len(self.periods)):
        #    if self.lsp_powers[i] > self.highest_lsp_power:
        #        self.best_period = self.periods[i]
        #        self.highest_lsp_power = self.lsp_powers[i]

        self.peak_periods, self.peak_powers = find_n_peaks(self.periods, self.lsp_powers, defaults.settings['n-peaks'] )

        self.opts_ls = { 'plot' : ax_lomb_scargle, 'periods_to_label' : self.peak_periods, 
                        'powers_to_label' : self.peak_powers }#, 'window_lsp' : 
                            #{ 'periods' : self.window_lsp_periods, 'powers' : self.window_lsp_powers }}
        if self.initialized:
            self.lsp_plot.replot(self.periods, self.lsp_powers, options = self.opts_ls)
    def set(self,period=None):

        self.t, self.y = [], []
        loaded_first_lc = False
        while (len(self.t) < 2 or not loaded_first_lc) and len(self.lc_files) > 0:
            if not os.path.exists(self.lc_files[self.lc_index]):
                file_exists = False
            else:
                file_exists = True
            unreadible = False
            cant_set_raw = False
            lc_is_none = False
            if file_exists:
                print self.lc_files[self.lc_index]
                rhlc.read_hatlc(self.lc_files[self.lc_index])
                try:

                    self.lc = rhlc.read_hatlc(self.lc_files[self.lc_index])
                except:
                    unreadible = True
                if not unreadible:
                    if self.lc is None: lc_is_none = True
                    if not lc_is_none:
                        try:
                            self.set_raw()
                        except:
                            cant_set_raw = True
            

            loaded_first_lc = True
            if len(self.t) < 2 or (unreadible or not file_exists or cant_set_raw or lc_is_none):
                if unreadible:
                    print "%s is unreadible..."%(self.lc_files[self.lc_index])
                elif not file_exists:
                    print "%s does not exist!"%(self.lc_files[self.lc_index])
                elif lc_is_none:
                    print "read_hatlc returned None"
                elif cant_set_raw:
                    print "Can't do set_raw()"
                else:
                    print "%s has no data..."%(self.lc_files[self.lc_index])
                if len(self.lc_files) == 1:
                    print "ALL lightcurves are bad!"
                    self.lc_files = []
                    self.lc_index = 0
                    self.period = 0
                    self.lc = None
                    return
                    #sys.exit()
                else:
                    dind = self.lc_index
                    self.lc_index = (self.lc_index + 1)%(len(self.lc_files) - 1)
                    self.lc_files.remove(self.lc_files[dind])
        self.commentfilename = "%s/%s.comments"%(defaults.settings['comment-file-dir'],self.lc['hatid'])
        if os.path.exists(self.commentfilename):
            with open(self.commentfilename, 'r') as commentfile:
                self.comments = commentfile.read()
        else:
            self.comments = ""
        self.set_lsp()
        self.set_phase_folded(period=self.best_period)
        
        if not isclosed( self.file_list_window.window ):
            self.file_list_window.update()
        if self.initialized:    
            self.info_box.update()
            self.toolbox.update_button_values()
            self.commentbox.update()
            #if not isclosed(self.commentbox.window):
            
                
                
        else:
            self.initialized = True

    def next_lightcurve(self):
        self.save()
        if len(self.lc_files) == 1: return

        if not isclosed(self.file_list_window.window):
            self.lc_index = self.file_list_window.get_next_lc_index()

        else: 
            self.lc_index += 1
            self.lc_index = self.lc_index%len(self.lc_files)

        self.set()
    def last_lightcurve(self):
        self.save()
        if len(self.lc_files) == 1: return
        if not isclosed(self.file_list_window.window):
            self.lc_index = self.file_list_window.get_prev_lc_index()
        else:
            self.lc_index -= 1
            self.lc_index = self.lc_index%len(self.lc_files)
        self.set()
    def update_options(self, new_opts):
        assert (isinstance(new_opts, dict))
        for n in new_opts:
            self.options[n] = new_opts[n]
        self.set()
    def add_lc(self, fname):
        # If no filenames loaded, initialize everything
        if len(self.lc_files) == 0: 
            self.lc_files.append(fname)
            self.set()
            return

        # If you already loaded the filename, just move there
        if fname in self.lc_files:
            self.lc_index = self.lc_files.index(fname)
            self.set()
            return

        # Not a special case - add the new filename
        # to the self.lcfnames list at i = self.lc_index + 1
        # then move to the next index (i.e. the file the user just opened)
        new_lc_files = []
        for i in range(0,self.lc_index+1):
            new_lc_files.append(self.lc_files[i])
        new_lc_files.append(fname)
        if self.lc_index + 2 < len(self.lc_files):
            for i in range(self.lc_index+2,len(self.lc_files)):
                new_lc_files.append(self.lc_files[i])
        self.lc_files = new_lc_files[:]
        self.next_lightcurve()
    def save(self, plot=True):
        # Save the flagged hat-ids
        if self.comments != "":
            with open(self.commentfilename,'w') as commentfile:
                commentfile.write(self.comments)
        if self.lc is None: return

        if self.lc['hatid'] not in self.log_info: return
        data = "%s %s"%(self.lc['hatid'], self.log_info[self.lc['hatid']])
        if self.logfile is None: 
            print "Logfile is none"
            print data
        else:
            with open(self.logfile, "w") as lfile:
                for hid in self.log_info:
                    lfile.write("%s %s\n"%(hid, self.log_info[hid]))
        

    def saveplot(self, fname=None):
        fig = plt.figure(figsize=(16, 8))
        ax_raw =  plt.subplot2grid((2,3), (0,0), colspan=2)
        ax_lsp = plt.subplot2grid((2,3), (1,0), colspan=1)
        ax_pf = plt.subplot2grid((2,3), (1,1), colspan=1)

        #ax_info = fig.add_subplot(236)

        mags = []
        for i in range(6):
            mags.append( "%.3f"%(float(self.lc['mags'][i]) ))

        hid = self.lc['hatid']
        if hid in self.log_info:
            linfo = self.log_info[hid]
        else:
            linfo = "---"
        #txt = '{0:15} {1}\n'.format('HATID:', hid) + \
        txt=  '{0:15} {1}\n'.format('2MASS ID:', self.lc['twomassid']) + \
              '{0:15} {1}\n'.format('FLAG:', linfo) + \
              '{0:15} {1:.6f}\n'.format('RA (J2000):', float(self.lc['ra'])) + \
              '{0:15} {1:.6f}\n\nMagnitude\n'.format('DEC (J2000):', float(self.lc['dec'])) + \
              '{0:3} {1} '.format('V:', mags[0]) + \
              '{0:3} {1} '.format('R:', mags[1]) + \
              '{0:3} {1}\n'.format('I:',mags[2]) + \
              '{0:3} {1} '.format('J:', mags[3]) + \
              '{0:3} {1} '.format('H:', mags[4]) + \
              '{0:3} {1}\n'.format('K:',mags[5])
        commenttxt = '\n'.join(wrap(self.comments,30))

        #fig.suptitle(hid)
        fig.text(0.7, 0.62, "Comments", fontproperties=FontProperties(family='sans-serif', weight='bold', size=14))
        fig.text(0.7, 0.6, commenttxt, family = 'monospace', ha='left', va='top')
        fig.text(0.7, 0.92, hid, fontproperties=FontProperties(family='sans-serif', weight='bold', size=14))
        ax_lomb_scargle(ax_lsp, self.periods, self.lsp_powers, self.opts_ls)
        ax_raw_lc(ax_raw, self.t, self.y, self.opts_raw )
        ax_phase_folded(ax_pf, self.phases, self.phase_mags, self.opts_phase_folded)
        fig.text(0.7, 0.9, txt,family = 'monospace', ha='left', va='top')
        #ax_info.set_axis_off()
        fig.set_tight_layout(True)
        #fig.subplots_adjust(top=0.8)
        if fname is None:
            fname = self.lc['hatid'] + ".png"
        fig.savefig(fname)
        plt.close(fig)


if __name__ == "__main__":


   

    # Parse command-line options
    parser = argparse.ArgumentParser(description='Visualize HAT lightcurves')
    parser.add_argument('--file', metavar='lcfile', type=str, nargs='+',
            help='A .gz file path containing a HAT lightcurve')
    parser.add_argument('--list', action='store', default=None, type=str,
            help='A file path containing a list of HAT ids or lightcurve file paths' ) 
    parser.add_argument('--flags', type=str, nargs='+',default=None,
            help='Do --flags flagname/keyboardshortcut ... -- this will allow you to flag LCs and store the results in a log file')
    parser.add_argument('--logfile', action='store', default=None, type=str,
            help='Place to store resulting LCs that you flag. Will digest the current logfile and append new results.')

    args = parser.parse_args()

    command_line_list = args.file
    list_file = args.list
    logfile = args.logfile
    shortcuts = {}
    flags = []
    if not args.flags is None:
        for flag_arg in args.flags:
            if "/" in flag_arg:
                flag,sc = flag_arg.split("/")
                shortcuts[sc] = flag
                flags.append(flag)
            else:
                flags.append(flag_arg)
    else:
        shortcuts = None
        flags = None


    # Combine the LC files specified on the command line, 
    # along with those in the list-file if one was specified
    # into one big list (all_files)
    if not command_line_list is None:
        command_line_files = command_line_list
    else:
        command_line_files = []
    if not list_file is None:
        list_file_contents = np.loadtxt(list_file,dtype='S15')
        if len(list_file_contents.shape) > 1: 
            list_file_contents = list_file_contents[:,0]
            if args.logfile is None:
                logfile = list_file

        list_files = list_file_contents

    else:
        list_files = []

    all_files = []
    for c in command_line_files: all_files.append(c)
    for l in list_files: all_files.append(l)

    # Visualize these!

    # Start Tk() instance
    root = Tk()
    root.geometry('+%s+%s'%(XPOS,YPOS)) 

    # Set window title
    root.wm_title("Lightcurve visualizer")
    app = Visualizer(root,all_files, logfile=logfile, flag_shortcuts=shortcuts, flags=flags)
    root.mainloop()


"""
Main file for Display.

Author: Wangkun Xu

The current version includes:
1. Run OPF.
2. Display the normal measurement.
3. Display both random and FDI attacks.

Limits:
1. I only write the display function with full pf and pi. 
2. The residual plot pop up automatically when running the code.

TODO:
1. In the future, we can display the estimated measurement.
2. Add button to control the pop up residual.
"""


"""
Import
"""
from src.load.case14.coordinate import *
from pypower.api import case14
from pypower.idx_brch import RATE_A, PF
from pypower.idx_bus import PD
from pypower.idx_gen import PG, QG
import tkinter as tk
from tkinter import Toplevel, ttk
from power_sys import power_env
import numpy as np
from matplotlib import cm, colors
from mea_idx import define_mea_idx_noise
import random
import matplotlib.pyplot as plt
from platform import system
from PIL import ImageTk, Image

"""
Call matplotlib-canvas backend
"""
# To connect matplotlib to tkinter canvas and also for dynamically update the plot
import matplotlib
matplotlib.use('TkAgg')    # Use Tkinter backend in matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class case_display:
    def __init__(self, display_window, case_env):
        
        """
        ATTR
        """
        self.display_window = display_window                   # The main window
        self.case_env = case_env                               # The instance of grid
        self.opf_idx = 0                                       # The number of run of opf
        self.flow_limit = case_env.case['branch'][:,RATE_A]    # The limit of active flow on each line
        self.is_pause = False
        self.att_choice = ['No Attack', 'Random Attack', 'FDI Attack']
        self.att_select = tk.StringVar()
        self.att_select.set(self.att_choice[0])   # Set the defaut as 'No Attack'

        self.residual_summary = []   # Record the residual of each run
        self.att_record = []         # Record the attack types at each opf run
        
        """
        Specify the ATTACK
        """
        # RANDOM ATT
        self.att_ratio_max = 0.15     # Maximum measurement change
        # FDI ATT
        att_spec = {}
        ang_no = np.random.randint(5,self.case_env.no_bus-1)
        att_spec['ang_posi'] = random.sample(self.case_env.non_ref_index, ang_no)
        att_spec['ang_str'] = -0.5+0.5*2*np.random.rand(ang_no)        
        mag_no = np.random.randint(self.case_env.no_bus-1)
        att_spec['mag_posi'] = random.sample(self.case_env.non_ref_index, mag_no)
        att_spec['mag_str'] = -0.001+0.002*np.random.rand()    
        self.att_spec = att_spec

        """
        Display settings
        """
        r = 15                                  # Parameter to control the ON-GRID DISPLAY, unit: px
        self.wait_time = 1000                   # Interval between two opf (in ms)
        self.color_no = 1000                    # The number of colors used in the branch loading rate
        self.color_name = 'plasma'
        self.font = 'arial'
        self.font_grid = (self.font, 16, 'bold')
        self.font_disp_legend = (self.font, 20, 'bold')
        self.font_disp_value = (self.font, 20)
        matplotlib.rcParams['lines.linewidth'] = 2    # Matplotlib settings
        matplotlib.rc('xtick', labelsize=20) 
        matplotlib.rc('ytick', labelsize=20)

        self.residual_plot_len = 20             # The No. of residuals show at the same time in the residual plot
        
        """
        ON-GRID DISPLAY
        """
        self.sys_canvas = tk.Canvas(self.display_window, width=550, height=500, bg = 'white', relief='ridge')
        
        # Record the coordinate of each bus
        for i in range(len(case_env.case['branch'])):
            f_bus = case_env.f_bus[i]
            t_bus = case_env.t_bus[i]
            x1 = coordinate[f'{f_bus+1}'][0]
            y1 = coordinate[f'{f_bus+1}'][1]
            x2 = coordinate[f'{t_bus+1}'][0]
            y2 = coordinate[f'{t_bus+1}'][1]

            # Draw branch
            self.sys_canvas.create_line(x1,y1,x2,y2,fill='black', width=5, tags=f'line_{i}')
            
            # Branch flow display
            x = (x1+x2)/2
            y = (y1+y2)/2
            self.sys_canvas.create_text(x,y, tags = f'line_text_{i}')  # The tag is used to refer and change the on-grid text
            
        for i in range(len(case_env.case['bus'])):
            x1 = coordinate[f'{i+1}'][0]
            y1 = coordinate[f'{i+1}'][1]
            x2 = coordinate[f'{i+1}'][0]
            y2 = coordinate[f'{i+1}'][1]
            self.sys_canvas.create_oval(x1-r,y1-r,x2+r,y2+r,fill='red')
            self.sys_canvas.create_text(x1,y1, text=f'{i+1}', fill='black', font=self.font_grid)

        """
        OFF-GRID DISPLAY
        """

        # Measurement display frame
        self.frame_mea = tk.Frame(self.display_window, width=430, height=500, bd = 0)  # Construct a frame to hold all the off-grid text 
        self.frame_mea.grid_propagate(0)

        # Label
        self.label_att = tk.Label(self.frame_mea, font = self.font_disp_legend, text = 'Please Choose the Attack Type:', anchor='w', padx=20)
        
        self.label_power_active = tk.Label(self.frame_mea, font = self.font_disp_legend, text = 'Total Active Load(MW):', width = '22', anchor='w', padx=20)
        self.label_power_reactive = tk.Label(self.frame_mea, font = self.font_disp_legend, text = 'Total Reactive Load(MVAr):',width = '22',anchor='w', padx=20)
        self.label_power_loss = tk.Label(self.frame_mea, font = self.font_disp_legend, text = 'Power Loss(MW):',width = '22', anchor='w', padx=20)
        self.label_residual = tk.Label(self.frame_mea, font = self.font_disp_legend, text = 'BDD Residual:', width = '22', anchor = 'w', padx=20)
        
        # Measure
        self.label_power_loss_ = tk.Label(self.frame_mea, font = self.font_disp_value, width = 7,anchor='w')
        self.label_power_active_ = tk.Label(self.frame_mea, font = self.font_disp_value, width = 7, anchor='w')
        self.label_power_reactive_ = tk.Label(self.frame_mea, font = self.font_disp_value, width = 7, anchor='w')
        self.label_residual_ = tk.Label(self.frame_mea, font = self.font_disp_value, width = 7, anchor = 'w')
        
        """
        Pause Button
        """
        self.button_pause = tk.Button(self.frame_mea, text = 'Pause', font = self.font_disp_legend, command=self.pause, relief=tk.RAISED)

        """
        Attack Combo Button
        """
        self.combo_att = ttk.Combobox(self.frame_mea, values = self.att_choice, textvariable=self.att_select, state = 'readonly', font=self.font_disp_value)

        """
        Authorship
        """
        cap_logo = Image.open('src/cap_logo.jpeg').resize((70,70), Image.ANTIALIAS)
        cap_logo = ImageTk.PhotoImage(cap_logo)
        print(cap_logo)
        label_cap_logo = tk.Label(master=self.frame_mea,image=cap_logo)
        label_cap_logo.image = cap_logo

        """
        Layout
        """
        # Main
        self.sys_canvas.grid(row=0,column=0)
        self.frame_mea.grid(row = 0,column=1)

        # Frame
        #label_cap_logo.place(anchor='nw')
        label_cap_logo.grid(row = 0, column = 1, padx=(0,0), pady = (0,40))
        self.label_att.grid(row = 1, column = 0, columnspan=1)
        self.combo_att.grid(row=2,column=0, columnspan=2, pady = (0,70))

        self.label_power_active.grid(row = 3, column = 0)
        self.label_power_reactive.grid(row = 4, column = 0)
        self.label_power_loss.grid(row = 5, column = 0)
        self.label_residual.grid(row = 6, column = 0)

        # self.button_residual.grid(row=6, column=0, pady = (70,0))
        self.button_pause.grid(row=7,column=0,columnspan=1, padx = (200,0), pady = (60,0))

        self.label_power_active_.grid(row = 3, column = 1)
        self.label_power_reactive_.grid(row = 4, column = 1)
        self.label_power_loss_.grid(row = 5, column = 1)
        self.label_residual_.grid(row = 6, column = 1)
        
        """
        TopLevel: Plot Residual
        """
        self.resid_window = Toplevel(master=self.display_window)
        self.resid_window.geometry("1200x800")
        self.resid_window.title('Residual Plot')
        figure = Figure()                            # Construct a figure
        self.sub = figure.add_subplot(111)
        self.sub.plot(self.residual_summary)
        self.figure_canvas = FigureCanvasTkAgg(figure, self.resid_window)    # Connect this figure to the Tkinter canvas
        self.figure_canvas.draw()
        self.figure_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)    # Place the figure

        # Call color map: Display different color based on the power flow
        self.color_map()
        # Call opf: Initial run the OPF 
        self.opf_loop()
        # Call display measurement
        self.show_mea()    

    def opf_loop(self):
        """
        For each display circle, call the opf_loop func to calculate the opf and display the result
        If there is an attack, then we display the attacked measurement
        """
        # Run OPF and receive the results
        result = self.case_env.run_opf(opf_idx = self.opf_idx)
        print(f'Is OPF success: {result["success"]}')             # Display in the terminal
        # Move the opf running idx by 1
        self.opf_idx+=1
        # Record the ACTUAL(NORMAL) measurement
        self.pf_actual = result['branch'][:,PF]
        self.load_active_actual = result['bus'][:,PD]

        # Generate the measurement
        load_active, load_reactive, pf, pl, residual = self.gen_mea(result)
            
        # Display
        self.load_active = load_active         # Active Load Power
        self.load_reactive = load_reactive     # Reactive Load Power
        self.pf = pf                           # Active from side power flow
        self.pl = pl                           # Active power loss
        self.residual = residual               # BDD residual
        self.residual_summary.append(residual) # Record the residual
        # print(self.residual_summary)

        # Draw and set the residual plot
        self.plot_residual()

        print('*'*50)
    
    def gen_mea(self, result):
        """
        Extract the measurement when there is an attack
        """
        # Generate the measurement
        z, z_noise, vang_ref, vmag_ref = self.case_env.construct_mea(result)
        if self.att_select.get() == self.att_choice[0]:
            """
            NO ATTACK
            """
            pass

        if self.att_select.get() == self.att_choice[1]:
            """
            RANDOM ATTACK
            """
            
            z_att_noise = self.case_env.gen_ran_att(z_noise, self.att_ratio_max)
            # (z_att_noise.shape)
            # Find the detection residual
            z_noise = z_att_noise     # pass as z_noise for consistency
            
        
        elif self.att_select.get() == self.att_choice[2]:
            """
            FDI ATTACK
            """
            # Do state estimation
            v_est = self.case_env.ac_se_pypower(z_noise, vang_ref, vmag_ref)
            z_est = self.case_env.h_x_pypower(v_est)
            # Generate FDI attack
            v_att = self.case_env.gen_fdi_att(v_est, self.att_spec)
            # Calculate the attacked measurement
            z_att_est = self.case_env.h_x_pypower(v_att)
            z_att_noise = z_noise + z_att_est - z_est
            z_noise = z_att_noise    # pass as z_noise for consistency

        # Calculate the residual
        residual = self.case_env.bdd_residual(z_noise, vang_ref, vmag_ref)

        # Define the attacked measurement (not real only shown to the system operator)
        # z = [pf, pt, pi, vang, qf, qt, qi, vmag] in the current setting
        # Measurement no
        pf_no = len(self.case_env.idx['pf'])
        pt_no = len(self.case_env.idx['pt'])
        pi_no = len(self.case_env.idx['pi'])
        vang_no = len(self.case_env.idx['vang'])
        qf_no = len(self.case_env.idx['qf'])
        qt_no = len(self.case_env.idx['qt'])
        qi_no = len(self.case_env.idx['qi'])
        vmag_no = len(self.case_env.idx['vmag'])

        z_noise = z_noise.squeeze(1)       # Reduce the last dimension
        pf = z_noise[:pf_no] 
        pt = z_noise[pf_no:pf_no+pt_no] 
        pi = z_noise[pf_no+pt_no:pf_no+pt_no+pi_no] 
        qi = z_noise[pf_no+pt_no+pi_no+vang_no+qf_no+qt_no:pf_no+pt_no+pi_no+vang_no+qf_no+qt_no+qi_no] 
        pg = result['gen'][:,PG]
        qg = result['gen'][:,QG]
        
        load_active = self.case_env.Cg@pg - pi*self.case_env.case['baseMVA']
        load_reactive = self.case_env.Cg@qg - qi*self.case_env.case['baseMVA']
        pl = np.sum(np.abs(pf+pt))
        
        return load_active, load_reactive, pf*self.case_env.case['baseMVA'], pl*self.case_env.case['baseMVA'], residual

    def plot_residual(self):
        """
        Plot the estimation residual
        """ 
        
        # Draw and set the residual plot
        self.att_record.append(self.att_select.get())
        # print(self.att_record)
        self.sub.clear()
        ticks_start = 0
        ticks_end = len(self.residual_summary)
        if len(self.residual_summary) >= self.residual_plot_len:
            residual_summary = self.residual_summary[-self.residual_plot_len:]
            ticks_start, ticks_end = len(self.residual_summary) - self.residual_plot_len, len(self.residual_summary)
        else:
            residual_summary = self.residual_summary

        self.sub.plot(residual_summary, label = 'Residual', color = 'blue')
        self.sub.hlines(xmin = 0, xmax = len(residual_summary)-1, y = self.case_env.bdd_threshold, label = 'Threshold', color = 'red')
        #plt.text(1, 50, 'Some Text', ha='center', va='center',rotation='vertical')
        plt.setp(self.sub, xticks = np.arange(len(residual_summary)), xticklabels = np.arange(ticks_start, ticks_end))
        self.sub.legend(prop = {'size': 15}, loc = 'center left')
        self.sub.set_xlabel('OPF Instance', fontsize=20)
        self.sub.set_ylabel('Residual', fontsize = 20)
        self.sub.set_ylim(0,np.max([np.max(residual_summary), self.case_env.bdd_threshold+10]))
        self.figure_canvas.draw()

    def show_mea(self):
        """
        Display the measurements
        """
        
        # Update off-grid Measurement
        self.label_power_active_['text'] = f'{np.round(np.sum(self.load_active),2)}'
        self.label_power_reactive_['text'] = f'{np.round(np.sum(self.load_reactive),2)}'
        self.label_power_loss_['text'] = f'{np.round(self.pl,2)}'
        #print(self.pl)
        self.label_residual_['text'] = f'{np.round(self.residual,2)}'
        
        # Change the residual text color
        if self.residual >= self.case_env.bdd_threshold:
            self.label_residual_['fg'] = 'red'
        else:
            self.label_residual_['fg'] = 'black'
        
        # Update the on-grid Measurement
        for i in range(case_env.no_brh):
            # Set the branch color according to the active power flow ratio
            brh_ratio = np.abs(self.pf[i]/self.flow_limit[i])           # The load rate of branch
            brh_ratio_actual = np.abs(self.pf_actual[i]/self.flow_limit[i])
            # Set the branch coloring according to the value of loading rate
            if brh_ratio >= 1:
                color_idx = -1                              # Pick the last color in the colormap
            else:
                color_idx = int(self.color_no*brh_ratio)    # Pick the color proportional to the power load ratio   
            
            self.sys_canvas.itemconfig(tagOrId = f'line_{i}', fill = self.color_map_list[color_idx])
            
            # Change the on-grid active flow load value
            if self.att_select.get() == self.att_choice[0]:
                # No att
                self.sys_canvas.itemconfig(tagOrId = f'line_text_{i}', fill = 'black', text = f'{int(brh_ratio*100)}' + '%', font = self.font_grid)
            else:
                # With attack
                # Also display the actual measurement 
                self.sys_canvas.itemconfig(tagOrId = f'line_text_{i}', fill = 'black', text = f'{int(brh_ratio*100)}%' + '\n' + f'{int(brh_ratio_actual*100)}%', font = self.font_grid)

        # Call OPF
        if self.is_pause == False:
            # Continuely call the opf
            self.opf_loop()
            self.display_window.after(1000, self.show_mea)  # Call the display function every .ms
        else:
            # Hold the current display
            pass
    
    """
    DISPLAY FUNCS
    """

    def color_map(self, color_name = 'magma'):
        """
        Create the color map used for different line active power flow loading
        """
        cmap = cm.get_cmap(self.color_name, self.color_no)
        color_map = []
        for i in range(cmap.N):
            rgba = cmap(i)
            color_map.append(colors.rgb2hex(rgba))
        self.color_map_list = color_map
        
    def pause(self):
        if self.is_pause == False:
            # Run -> Pause
            self.is_pause = True
            self.button_pause['text'] = 'Continue'
        else:
            # Pause -> Run
            self.is_pause = False
            self.button_pause['text'] = 'Pause'
            # Recall the function to continue
            self.show_mea()
    
    # def visual_resid(self):
    #     """
    #     If the visualizing residual button is clicked, call matplotlib to draw the residual plot.
    #     """

    #     def draw_residual(i):
    #         self.sub.clear()
    #         self.sub.plot(self.residual_summary)
        
    #     top = tk.Toplevel()   # pop up a new window
    #     top.title('Residual Plot')
    #     figure = Figure()                            # Construct a figure
    #     self.sub = figure.add_subplot(111)
    #     self.sub.plot(self.residual_summary)
    #     animation.FuncAnimation(figure, draw_residual, interval = 1000)
    #     figure_canvas = FigureCanvasTkAgg(figure, top)    # Connect this figure to the Tkinter canvas
    #     figure_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)    # Place the figure
        
    # def att(self):
    #     if self.is_att == False:
    #         # no att -> att
    #         self.is_att = True
    #     else:
    #         # att -> no att
    #         self.is_att = False
            
def os_config():
    platformD = system()
    print(platformD)
    if platformD == 'Darwin':
        logo_image = 'src/power_plant.icns'
    elif platformD == 'Windows':
        logo_image = 'src/power_plant.ico'
    
    return logo_image

if __name__ == "__main__":
    # Instance power env
    case_name = 'case14'
    case = case14()
    
    # Define measurement index
    mea_idx, no_mea, noise_sigma = define_mea_idx_noise(case, 'RTU_POWER')
    
    # Instance the class
    case_env = power_env(case = case, case_name = case_name, noise_sigma = noise_sigma, idx = mea_idx, fpr = 0.05)
    
    # Define the Tkinter mainloop
    root = tk.Tk()
    root.title('IEEE ' + case_name + ' system')

    #logo_image = os_config()
    #root.iconbitmap(logo_image)

    os = system()
    if os == 'Darwin':
        # mac os
        img = tk.Image("photo", file = 'src/bulb_logo.png')
        root.iconphoto(True, img)
    elif os == 'Windows':
        root.iconbitmap('src/bulb.ico')
    
    case_display(root, case_env)
    root.mainloop()

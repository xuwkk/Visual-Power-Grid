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
from gen_data import gen_case, gen_load
from src.case14.coordinate import *
from pypower.api import case14
from pypower.idx_brch import RATE_A, PF, BR_X
from pypower.idx_bus import PD
from pypower.idx_gen import PG, QG
import tkinter as tk
from tkinter import DISABLED, Toplevel, ttk
from power_sys import power_env
import numpy as np
from matplotlib import cm, colors
from config_mea_idx import define_mea_idx_noise
import random
import matplotlib.pyplot as plt
from platform import system
from PIL import ImageTk, Image
from copy import deepcopy

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
        
        self.delay_time = 1500       # The time interval between two OPFs

        # Copy the case_env whose reactance never changes for the attacker
        self.case_env_ori = deepcopy(case_env)

        """
        Measurement No.
        """
        # z = [pf, pt, pi, vang, qf, qt, qi, vmag] in the current setting
        # Measurement no
        self.pf_no = len(self.case_env.idx['pf'])
        self.pt_no = len(self.case_env.idx['pt'])
        self.pi_no = len(self.case_env.idx['pi'])
        self.vang_no = len(self.case_env.idx['vang'])
        self.qf_no = len(self.case_env.idx['qf'])
        self.qt_no = len(self.case_env.idx['qt'])
        self.qi_no = len(self.case_env.idx['qi'])
        self.vmag_no = len(self.case_env.idx['vmag'])

        """
        Attack Settings
        """

        # RANDOM ATT
        self.att_ratio_max = 0.01     # Maximum measurement change

        # FDI ATT
        self.fdi_min_posi_no = 4      # Maximum and minimum number of attacked buses
        self.fdi_max_posi_no = 8

        """
        MTD settings
        """
        self.max_reac_ratio = 0.3
        self.min_reac_ratio = 0.1

        """
        Display settings
        """
        r = 15                                  # Parameter to control the ON-GRID DISPLAY, unit: px
        self.r_brch = 10                        # shifting parameters for load flow text, unit: px
        self.wait_time = 1000                   # Interval between two opf (in ms)
        self.color_no = 1000                    # The number of colors used in the branch loading rate
        self.color_name = 'plasma'
        self.font = 'arial'
        self.font_grid = (self.font, 16, 'bold')
        self.font_disp_legend = (self.font, 20, 'bold')
        self.font_disp_value = (self.font, 20)
        matplotlib.rcParams['lines.linewidth'] = 2    # Matplotlib settings
        matplotlib.rc('xtick', labelsize=16) 
        matplotlib.rc('ytick', labelsize=16)

        self.residual_plot_len = 15             # The No. of residuals show at the same time in the residual plot
        
        """
        ON-GRID DISPLAY
        """
        self.sys_canvas = tk.Canvas(self.display_window, width=570, height=500, bg = 'white', relief='ridge')
        
        # Coordinate of each bus
        for i in range(len(case_env.case['branch'])):
            f_bus = case_env.f_bus[i]
            t_bus = case_env.t_bus[i]
            x1 = coordinate[f'{f_bus+1}'][0]
            y1 = coordinate[f'{f_bus+1}'][1]
            x2 = coordinate[f'{t_bus+1}'][0]
            y2 = coordinate[f'{t_bus+1}'][1]

            # Draw BRANCHES
            self.sys_canvas.create_line(x1,y1,x2,y2,fill='black', width=5, tags=f'line_{i}')
            
            # Branch flow display
            x = (x1+x2)/2
            y = (y1+y2)/2

            # The tag is used to change the color of different line loading
            # Design the load flow text:
            # Below the line: normal measurements
            # Above the line: Attacked measurements
            self.sys_canvas.create_text(x,y - self.r_brch, tags = f'line_text_above_{i}')    # Below: normal 
            self.sys_canvas.create_text(x,y + self.r_brch, tags = f'line_text_below_{i}')    # Above: attack             
        
        # Draw BUSES
        for i in range(len(case_env.case['bus'])):
            x1 = coordinate[f'{i+1}'][0]
            y1 = coordinate[f'{i+1}'][1]
            x2 = coordinate[f'{i+1}'][0]
            y2 = coordinate[f'{i+1}'][1]
            # The tag is used to refer and change the color of the attacked bus
            self.sys_canvas.create_oval(x1-r,y1-r,x2+r,y2+r,fill='green', tags=f'bus_{i}')        
            self.sys_canvas.create_text(x1,y1, text=f'{i+1}', fill='black', font=self.font_grid)

        """
        OFF-GRID DISPLAY
        """

        # Measurement display frame
        self.frame_mea = tk.Frame(self.display_window, width=430, height=500, bd = 0)  # Construct a frame to hold all the off-grid text 
        self.frame_mea.grid_propagate(0)

        # Label
        self.label_att = tk.Label(self.frame_mea, font = self.font_disp_legend, text = 'Choose the Attack Type:', width = '25', anchor='w')
        
        self.label_power_active = tk.Label(self.frame_mea, font = self.font_disp_legend, text = 'Total Active Load(MW):', width = '25', anchor='w')
        self.label_power_reactive = tk.Label(self.frame_mea, font = self.font_disp_legend, text = 'Total Reactive Load(MVAr):',width = '25',anchor='w')
        self.label_power_loss = tk.Label(self.frame_mea, font = self.font_disp_legend, text = 'Power Loss(MW):',width = '25', anchor='w')
        self.label_residual = tk.Label(self.frame_mea, font = self.font_disp_legend, text = 'BDD Residual:', width = '25', anchor = 'w')
        
        # Measure
        self.label_power_loss_ = tk.Label(self.frame_mea, font = self.font_disp_value, width = 7,anchor='w')
        self.label_power_active_ = tk.Label(self.frame_mea, font = self.font_disp_value, width = 7, anchor='w')
        self.label_power_reactive_ = tk.Label(self.frame_mea, font = self.font_disp_value, width = 7, anchor='w')
        self.label_residual_ = tk.Label(self.frame_mea, font = self.font_disp_value, width = 7, anchor = 'w')
        
        """
        Pause Button
        """
        self.button_pause = tk.Button(self.frame_mea, text = 'Pause', font = self.font_disp_legend, command=self.pause, relief=tk.RAISED, width=10)

        """
        Attack Combo Button
        """
        self.combo_att = ttk.Combobox(self.frame_mea, values = self.att_choice, textvariable=self.att_select, \
            state = 'readonly', font=self.font_disp_value, width ='15', justify=tk.LEFT)

        """
        MTD Check Button
        """
        self.mtd_on = tk.IntVar()
        self.check_mtd = tk.Checkbutton(self.frame_mea, text = 'MTD Trigger', variable = self.mtd_on, \
            onvalue=1, offvalue=0, state = DISABLED,font = self.font_disp_legend, width='25', anchor='w')

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
        self.sys_canvas.grid(row=0,column=0)   # The canvas to show the grid
        self.frame_mea.grid(row = 0,column=1)  # The frame for user choice and diaplay measurement

        # Frame
        #label_cap_logo.place(anchor='nw')
        label_cap_logo.grid(row = 0, column = 1, padx=(0,0), pady = (0,30))
        self.label_att.grid(row = 1, column = 0, padx = 20)
        self.combo_att.grid(row = 2, column = 0, padx = 20, sticky='w')
        self.check_mtd.grid(row = 3, column = 0, pady = (0,70), padx = 20)

        self.label_power_active.grid(row = 4, column = 0, padx = 20)
        self.label_power_reactive.grid(row = 5, column = 0, padx = 20)
        self.label_power_loss.grid(row = 6, column = 0, padx = 20)
        self.label_residual.grid(row = 7, column = 0, padx = 20)

        # self.button_residual.grid(row=6, column=0, pady = (70,0))
        self.button_pause.grid(row = 8, column = 0, columnspan=2, pady = (50,0))

        self.label_power_active_.grid(row = 4, column = 1)
        self.label_power_reactive_.grid(row = 5, column = 1)
        self.label_power_loss_.grid(row = 6, column = 1)
        self.label_residual_.grid(row = 7, column = 1)
        
        """
        TopLevel: Plot Residual
        """
        self.resid_window = Toplevel(master=self.display_window)
        #self.resid_window.geometry("600x500")
        self.resid_window.title('Residual Plot')
        figure = Figure(figsize = (9,6), dpi = 100)                            # Construct a figure
        self.sub = figure.add_subplot(111)
        self.sub.plot(self.residual_summary)
        
        self.figure_canvas = FigureCanvasTkAgg(figure, master = self.resid_window)    # Connect this figure to the Tkinter canvas
        self.figure_canvas.draw()
        self.figure_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)    # Place the figure
        # side=tk.TOP, fill=tk.BOTH, expand=1
        # Call color map: Display different color based on the power flow
        self.color_map()
        # Call opf: Initial run the OPF 
        self.opf_loop()
        # Call display measurement
        self.show_mea()
        # Call MTD trigger
        self.trigger_mtd()    

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
        load_active, load_reactive, pf_actual, pf, pl, residual = self.gen_mea(result)
            
        # Display
        self.load_active = load_active          # Active Load Power
        self.load_reactive = load_reactive      # Reactive Load Power
        self.pf_actual = pf_actual              # Actual (Nromal) power flow
        self.pf = pf                            # Active from side power flow
        self.pl = pl                            # Active power loss
        self.residual = residual                # BDD residual
        self.residual_summary.append(residual)  # Record the residual
        # print(self.residual_summary)

        # Draw and set the residual plot
        self.plot_residual()

        print('*'*60)
    
    def gen_mea(self, result):
        """
        Extract the measurement when there is an attack
        """
        # Generate the measurement
        z, z_noise, vang_ref, vmag_ref = self.case_env.construct_mea(result)
        
        # Record the ACTUAL(NORMAL) measurement
        pf_actual = result['branch'][:,PF]
        load_active_actual = result['bus'][:,PD]

        if self.att_select.get() == self.att_choice[0]:
            """
            NO ATTACK
            """
            # Deselect the MTD check button and disabled the selection
            self.check_mtd.deselect()
            self.check_mtd.config(state = tk.DISABLED)

        elif self.att_select.get() == self.att_choice[1]:
            """
            RANDOM ATTACK
            """
            z_att_noise = self.case_env.gen_ran_att(z_noise, self.att_ratio_max)
            
            z_noise = z_att_noise     # pass as z_noise for consistency

            # Deselect the MTD check button and disabled the selection
            self.check_mtd.deselect()
            self.check_mtd.config(state = tk.DISABLED)
        
        elif self.att_select.get() == self.att_choice[2]:
            """
            FDI ATTACK
            """
            # Set the MTD trigger check box to NORMAL
            self.check_mtd.config(state = tk.NORMAL)

            # FDI ATT
            att_spec = {}
            # Voltage angle
            ang_no = np.random.randint(self.fdi_min_posi_no, self.fdi_max_posi_no+1)     # Attack bus number
            att_spec['ang_posi'] = random.sample(self.case_env.non_ref_index, ang_no)    # Attack bus position
            att_spec['ang_str'] = -0.5+0.5*2*np.random.rand(ang_no)                      # Attack value
            # Voltage magnitude
            mag_no = np.random.randint(self.case_env.no_bus-1)
            att_spec['mag_posi'] = random.sample(self.case_env.non_ref_index, mag_no)
            att_spec['mag_str'] = -0.001+0.002*np.random.rand()                        # The voltage attack on manitude is very small
            # self.att_spec = att_spec

            # Do state estimation
            # NOTE: case_env_ori is used
            v_est, _ = self.case_env_ori.ac_se_pypower(z_noise, vang_ref, vmag_ref)
            residual = self.case_env_ori.bdd_residual(z_noise, v_est)
            print(f'Attacker residual: {residual}')
            z_est = self.case_env_ori.h_x_pypower(v_est)
            
            # Generate FDI attack
            v_att, self.att_posi = self.case_env_ori.gen_fdi_att(v_est, att_spec)
            
            # Calculate the attacked measurement
            z_att_est = self.case_env_ori.h_x_pypower(v_att)
            z_att_noise = z_noise + z_att_est - z_est
            
            z_noise = z_att_noise    # pass as z_noise for consistency

            # Print the attack position
            print(f'The FDI attack position is: {self.att_posi}')

        # Calculate the residual
        # NOTE: based on the changed reatance
        v_est, _ = self.case_env.ac_se_pypower(z_noise, vang_ref, vmag_ref)
        residual = self.case_env.bdd_residual(z_noise, v_est)

        z_noise = z_noise.squeeze(1)       
        pf = z_noise[:self.pf_no] 
        pt = z_noise[self.pf_no:self.pf_no+self.pt_no] 
        pi = z_noise[self.pf_no+self.pt_no:self.pf_no+self.pt_no+self.pi_no] 
        qi = z_noise[self.pf_no+self.pt_no+self.pi_no+self.vang_no+self.qf_no+self.qt_no:self.pf_no+self.pt_no+self.pi_no+self.vang_no+self.qf_no+self.qt_no+self.qi_no] 
        pg = result['gen'][:,PG]
        qg = result['gen'][:,QG]
        
        load_active = self.case_env.Cg@pg - pi*self.case_env.case['baseMVA']
        load_reactive = self.case_env.Cg@qg - qi*self.case_env.case['baseMVA']
        pl = np.sum(np.abs(pf+pt))
        
        # print(self.mtd_on.get())

        return load_active, load_reactive, pf_actual, pf*self.case_env.case['baseMVA'], pl*self.case_env.case['baseMVA'], residual

    def trigger_mtd(self):
        """
        Update/Set to original reactance in the case file
        """

        if self.mtd_on.get() == 1:
            # MTD is triggered, update the reactance
            increase_decrease = np.random.randint(0,2,size = self.case_env.no_brh)*2-1   # -1 or 1
            ratio = self.min_reac_ratio + (self.max_reac_ratio - self.min_reac_ratio)*np.random.rand(self.case_env.no_brh,)
            reac_new = self.case_env.reactance_ori * (1+ratio*increase_decrease)
            self.case_env.update_reactance(reac_new)
        
        else:
            # print('here!')
            # MTD is turned off, back to the default situation
            self.case_env.update_reactance(self.case_env.reactance_ori)

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
        self.sub.set_ylim(0,np.max([np.max(residual_summary), self.case_env.bdd_threshold+10]))
        self.sub.set_xlabel('OPF Instance', fontsize=15)
        self.sub.set_ylabel('Residual', fontsize = 15)
        self.figure_canvas.draw()

    def show_mea(self):
        if self.is_pause == False:
            # Continue
            self.opf_loop()       # Run OPF
            self.trigger_mtd()    # Determine whether the MTD should be triggered or not
            print(f'MTD statue: {self.mtd_on.get()}')
            # print(self.case_env.case['branch'][:,BR_X])

            self.display_update() # Update the display
            self.display_window.after(self.delay_time, self.show_mea)  # Refresh the whole window
        else:
            # Pause
            pass

    def display_update(self):
        """
        OFF-grid measurement update
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
        
        """
        ON-grid measurement update
        """
        
        # Branch
        for i in range(self.case_env.no_brh):
            # Set the branch color according to the active power flow ratio
            brh_ratio = np.abs(self.pf[i]/self.flow_limit[i])           # The load rate of branch
            brh_ratio_actual = np.abs(self.pf_actual[i]/self.flow_limit[i])
            
            # Set the branch coloring according to the value of loading rate
            if brh_ratio_actual >= 1:
                color_idx_actual = -1
            else:
                color_idx_actual = int(self.color_no*brh_ratio_actual)
            
            if brh_ratio >= 1:
                color_idx = -1                              # Pick the last color in the colormap
            else:
                color_idx = int(self.color_no*brh_ratio)    # Pick the color proportional to the power load ratio   
            
            self.sys_canvas.itemconfig(tagOrId = f'line_{i}', fill = self.color_map_list[color_idx])
            
            # Change the on-grid active flow load value
            if self.att_select.get() == self.att_choice[0]:
                # No att
                # Normal operation
                self.sys_canvas.itemconfig(tagOrId = f'line_text_below_{i}', fill = self.color_map_list[color_idx_actual], text = f'{int(brh_ratio_actual*100)}' + '%', font = self.font_grid)

                # Clear the Attack text: above
                self.sys_canvas.itemconfig(tagOrId = f'line_text_above_{i}', fill = self.color_map_list[color_idx_actual], text = f'', font = self.font_grid)

            else:
                # With attacks
                # Also display the actual measurement
                # Normal 
                self.sys_canvas.itemconfig(tagOrId = f'line_text_below_{i}', fill = self.color_map_list[color_idx_actual], text = f'{int(brh_ratio_actual*100)}%', font = self.font_grid)
                # Attack
                self.sys_canvas.itemconfig(tagOrId = f'line_text_above_{i}', fill = self.color_map_list[color_idx], text = f'{int(brh_ratio*100)}%', font = self.font_grid)

        # Bus
        if self.att_select.get() == self.att_choice[-1]:
            # Change the bus color if this bus is attacked by FDI attacks
            for i in range(self.case_env.no_bus):
                if i in set(self.att_posi):
                    self.sys_canvas.itemconfig(tagOrId = f'bus_{i}', fill = 'red')
                else:
                    self.sys_canvas.itemconfig(tagOrId = f'bus_{i}', fill = 'green')

        # # Call OPF
        # if self.is_pause == False:
        #     # Continuely call the opf
        #     self.opf_loop()
        #     self.display_window.after(self.delay_time, self.show_mea)  # Call the display function every .ms
        # else:
        #     # Hold the current display
        #     pass
    
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
    case = gen_case(case, 'case14')  # Modify the case
    
    # Define measurement index
    mea_idx, no_mea, noise_sigma = define_mea_idx_noise(case, 'RTU')
    
    # Instance the class
    case_env = power_env(case = case, case_name = case_name, noise_sigma = noise_sigma, idx = mea_idx, fpr = 0.05)
    
    # Generate load if it does not exist
    _, _ = gen_load(case, 'case14')

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

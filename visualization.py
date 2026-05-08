"""
Main file for Display.

Author: Wangkun Xu
Copyright (c) 2026 Wangkun Xu. All rights reserved.

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
from tkinter import ttk
from power_sys import power_env
import numpy as np
from matplotlib import cm, colors
from config_mea_idx import define_mea_idx_noise
import random
from platform import system
from PIL import ImageTk, Image
from copy import deepcopy

try:
    RESAMPLE_LANCZOS = Image.Resampling.LANCZOS
except AttributeError:
    RESAMPLE_LANCZOS = Image.ANTIALIAS

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
        self.event_markers = []      # Record attack/MTD start indices for the residual plot
        self.prev_attack_mode = self.att_select.get()
        self.prev_mtd_state = 0
        
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
        r = 14                                  # Parameter to control the ON-GRID DISPLAY, unit: px
        self.r_brch = 9                         # shifting parameters for load flow text, unit: px
        self.wait_time = 1000                   # Interval between two opf (in ms)
        self.color_no = 256                     # The number of colors used in the branch loading rate
        self.color_name = 'viridis'
        self.font = 'Arial'
        self.font_title = (self.font, 20, 'bold')
        self.font_section = (self.font, 13, 'bold')
        self.font_label = (self.font, 11)
        self.font_value = (self.font, 17, 'bold')
        self.font_badge = (self.font, 11, 'bold')
        self.font_panel_section = (self.font, 15, 'bold')
        self.font_panel_label = (self.font, 13)
        self.font_panel_value = (self.font, 20, 'bold')
        self.font_panel_button = (self.font, 13, 'bold')
        self.font_panel_footer = (self.font, 10)
        self.font_grid = (self.font, 11, 'bold')
        self.font_grid_small = (self.font, 10, 'bold')
        self.palette = {
            'bg': '#eef1f4',
            'panel': '#ffffff',
            'panel_edge': '#d5dce3',
            'canvas': '#f9fbfc',
            'text': '#1f2933',
            'muted': '#627080',
            'normal': '#256f86',
            'warning': '#c77900',
            'alarm': '#c92a2a',
            'attack': '#7c3aed',
            'ok_bg': '#dff3ea',
            'ok_fg': '#17603a',
            'warn_bg': '#fff0d1',
            'warn_fg': '#8a5200',
            'alarm_bg': '#ffe1df',
            'alarm_fg': '#a21d1d',
            'idle_bg': '#e6ebf0',
            'idle_fg': '#465461',
        }
        matplotlib.rcParams['lines.linewidth'] = 2    # Matplotlib settings
        matplotlib.rc('xtick', labelsize=10)
        matplotlib.rc('ytick', labelsize=10)

        self.residual_plot_len = 15             # The No. of residuals show at the same time in the residual plot

        self.display_window.configure(bg=self.palette['bg'])
        self.display_window.minsize(1060, 720)
        self._configure_style()

        self.app_frame = tk.Frame(self.display_window, bg=self.palette['bg'], padx=10, pady=10)
        self.app_frame.grid(row=0, column=0, sticky='nsew')
        self.display_window.rowconfigure(0, weight=1)
        self.display_window.columnconfigure(0, weight=1)
        self.app_frame.rowconfigure(1, weight=1)
        self.app_frame.columnconfigure(0, weight=1)

        self.header_frame = tk.Frame(self.app_frame, bg=self.palette['bg'])
        self.header_frame.grid(row=0, column=0, sticky='ew', pady=(0, 8))
        self.header_frame.columnconfigure(0, weight=1)

        self.title_label = tk.Label(
            self.header_frame,
            text='IEEE 14-Bus System',
            font=self.font_title,
            bg=self.palette['bg'],
            fg=self.palette['text'],
            anchor='w',
        )
        self.title_label.grid(row=0, column=0, sticky='w')

        self.run_status_label = tk.Label(self.header_frame, font=self.font_badge, padx=8, pady=4)
        self.run_status_label.grid(row=0, column=1, padx=(6, 0))
        self.bdd_status_label = tk.Label(self.header_frame, font=self.font_badge, padx=8, pady=4)
        self.bdd_status_label.grid(row=0, column=2, padx=(6, 0))
        self.opf_label = tk.Label(
            self.header_frame,
            text='OPF 0',
            font=self.font_label,
            bg=self.palette['bg'],
            fg=self.palette['muted'],
        )
        self.opf_label.grid(row=0, column=3, padx=(10, 0))

        self.content_frame = tk.Frame(self.app_frame, bg=self.palette['bg'])
        self.content_frame.grid(row=1, column=0, sticky='nsew')
        self.content_frame.rowconfigure(0, weight=1)
        self.content_frame.columnconfigure(0, weight=1)

        self.grid_panel = tk.Frame(
            self.content_frame,
            bg=self.palette['panel'],
            highlightbackground=self.palette['panel_edge'],
            highlightthickness=1,
            padx=8,
            pady=8,
        )
        self.grid_panel.grid(row=0, column=0, sticky='nsew', padx=(0, 10))
        self.grid_panel.rowconfigure(1, weight=1)
        self.grid_panel.columnconfigure(0, weight=1)

        self.grid_header = tk.Frame(self.grid_panel, bg=self.palette['panel'])
        self.grid_header.grid(row=0, column=0, sticky='ew', pady=(0, 6))
        self.grid_header.columnconfigure(0, weight=1)
        tk.Label(
            self.grid_header,
            text='Grid State',
            font=self.font_section,
            bg=self.palette['panel'],
            fg=self.palette['text'],
            anchor='w',
        ).grid(row=0, column=0, sticky='w')
        self.attack_status_label = tk.Label(self.grid_header, font=self.font_badge, padx=7, pady=3)
        self.attack_status_label.grid(row=0, column=1, padx=(6, 0))
        self.mtd_status_label = tk.Label(self.grid_header, font=self.font_badge, padx=7, pady=3)
        self.mtd_status_label.grid(row=0, column=2, padx=(6, 0))
        
        """
        ON-GRID DISPLAY
        """
        self.sys_canvas = tk.Canvas(
            self.grid_panel,
            width=570,
            height=468,
            bg=self.palette['canvas'],
            highlightthickness=1,
            highlightbackground=self.palette['panel_edge'],
        )
        self.sys_canvas.grid(row=1, column=0, sticky='nsew')
        self.grid_diagram_offset = (0, 0)
        self.grid_legend_bottom = 42
        
        # Coordinate of each bus
        for i in range(len(case_env.case['branch'])):
            f_bus = case_env.f_bus[i]
            t_bus = case_env.t_bus[i]
            x1 = coordinate[f'{f_bus+1}'][0]
            y1 = coordinate[f'{f_bus+1}'][1]
            x2 = coordinate[f'{t_bus+1}'][0]
            y2 = coordinate[f'{t_bus+1}'][1]

            # Draw BRANCHES
            self.sys_canvas.create_line(
                x1,
                y1,
                x2,
                y2,
                fill='#9ba8b3',
                width=4,
                capstyle=tk.ROUND,
                tags=(f'line_{i}', 'grid_diagram'),
            )
            
            # Branch flow display
            x = (x1+x2)/2
            y = (y1+y2)/2

            # The tag is used to change the color of different line loading.
            self.sys_canvas.create_text(x, y - self.r_brch, tags=(f'line_text_above_{i}', 'grid_diagram'))
            self.sys_canvas.create_text(x, y + self.r_brch, tags=(f'line_text_below_{i}', 'grid_diagram'))
        
        # Draw BUSES
        for i in range(len(case_env.case['bus'])):
            x1 = coordinate[f'{i+1}'][0]
            y1 = coordinate[f'{i+1}'][1]
            x2 = coordinate[f'{i+1}'][0]
            y2 = coordinate[f'{i+1}'][1]
            # The tag is used to refer and change the color of the attacked bus
            self.sys_canvas.create_oval(
                x1-r,
                y1-r,
                x2+r,
                y2+r,
                fill='#ffffff',
                outline=self.palette['normal'],
                width=2,
                tags=(f'bus_{i}', 'grid_diagram'),
            )
            self.sys_canvas.create_text(
                x1,
                y1,
                text=f'{i+1}',
                fill=self.palette['text'],
                font=self.font_grid,
                tags=('grid_diagram', 'bus_label'),
            )

        self._draw_grid_legend()
        self.grid_diagram_bbox = self.sys_canvas.bbox('grid_diagram')
        self.sys_canvas.bind('<Configure>', self._center_grid_diagram)
        self.display_window.after_idle(self._center_grid_diagram)

        """
        OFF-GRID DISPLAY
        """

        # Measurement display frame
        self.frame_mea = tk.Frame(
            self.content_frame,
            width=322,
            bg=self.palette['panel'],
            highlightbackground=self.palette['panel_edge'],
            highlightthickness=1,
            padx=14,
            pady=12,
        )  # Construct a frame to hold all the off-grid text
        self.frame_mea.grid(row=0, column=1, sticky='ns')
        self.frame_mea.grid_propagate(False)
        self.frame_mea.columnconfigure(0, weight=1)
        self.frame_mea.rowconfigure(3, weight=1)

        logo_header = tk.Frame(self.frame_mea, bg=self.palette['panel'])
        logo_header.grid(row=0, column=0, sticky='ew', pady=(0, 12))
        logo_header.columnconfigure(0, weight=1)

        tk.Label(
            logo_header,
            text='Operations',
            font=self.font_panel_section,
            bg=self.palette['panel'],
            fg=self.palette['text'],
            anchor='w',
        ).grid(row=0, column=0, sticky='w')

        """
        Attack Combo Button
        """
        control_frame = tk.Frame(self.frame_mea, bg=self.palette['panel'])
        control_frame.grid(row=1, column=0, sticky='ew')
        control_frame.columnconfigure(0, weight=1)

        self.label_att = tk.Label(
            control_frame,
            font=self.font_panel_label,
            text='Attack Type',
            bg=self.palette['panel'],
            fg=self.palette['muted'],
            anchor='w',
        )
        self.label_att.grid(row=0, column=0, sticky='ew')

        self.combo_att = ttk.Combobox(
            control_frame,
            values=self.att_choice,
            textvariable=self.att_select,
            state='readonly',
            font=self.font_panel_label,
            width=18,
            justify=tk.LEFT,
        )
        self.combo_att.grid(row=1, column=0, sticky='ew', pady=(3, 8))

        """
        MTD Check Button
        """
        self.mtd_on = tk.IntVar()
        self.check_mtd = ttk.Checkbutton(
            control_frame,
            text='MTD Trigger',
            variable=self.mtd_on,
            onvalue=1,
            offvalue=0,
            state=tk.DISABLED,
        )
        self.check_mtd.grid(row=2, column=0, sticky='w', pady=(0, 8))

        self.button_pause = ttk.Button(
            control_frame,
            text='Pause',
            command=self.pause,
            style='Primary.TButton',
        )
        self.button_pause.grid(row=3, column=0, sticky='ew', pady=(0, 12))

        metric_frame = tk.Frame(self.frame_mea, bg=self.palette['panel'])
        metric_frame.grid(row=2, column=0, sticky='ew')
        metric_frame.columnconfigure(0, weight=1)
        tk.Label(
            metric_frame,
            text='Live Metrics',
            font=self.font_panel_section,
            bg=self.palette['panel'],
            fg=self.palette['text'],
            anchor='w',
        ).grid(row=0, column=0, sticky='ew', pady=(0, 6))

        self.label_power_active_ = self._metric_row(metric_frame, 1, 'Active Load', 'MW')
        self.label_power_reactive_ = self._metric_row(metric_frame, 2, 'Reactive Load', 'MVAr')
        self.label_power_loss_ = self._metric_row(metric_frame, 3, 'Power Loss', 'MW')
        self.label_residual_ = self._metric_row(metric_frame, 4, 'BDD Residual', '')
        self.label_threshold_ = self._metric_row(metric_frame, 5, 'BDD Threshold', '')

        """
        Authorship
        """
        cap_logo = Image.open('src/cap_logo.jpeg').resize((70,70), RESAMPLE_LANCZOS)
        cap_logo = ImageTk.PhotoImage(cap_logo)
        label_cap_logo = tk.Label(master=logo_header, image=cap_logo, bg=self.palette['panel'])
        label_cap_logo.image = cap_logo
        label_cap_logo.grid(row=0, column=1, sticky='e')

        copyright_label = tk.Label(
            self.frame_mea,
            text='Power System Visualization Demo\nCopyright (c) 2026 Wangkun Xu\nAll rights reserved.',
            font=self.font_panel_footer,
            bg=self.palette['panel'],
            fg=self.palette['muted'],
            justify=tk.LEFT,
            anchor='sw',
        )
        copyright_label.grid(row=4, column=0, sticky='sew')

        """
        Residual Plot
        """
        self.chart_panel = tk.Frame(
            self.app_frame,
            bg=self.palette['panel'],
            highlightbackground=self.palette['panel_edge'],
            highlightthickness=1,
            padx=10,
            pady=8,
        )
        self.chart_panel.grid(row=2, column=0, sticky='ew', pady=(8, 0))
        self.chart_panel.columnconfigure(0, weight=1)
        chart_header = tk.Frame(self.chart_panel, bg=self.palette['panel'])
        chart_header.grid(row=0, column=0, sticky='ew', pady=(0, 4))
        chart_header.columnconfigure(0, weight=1)
        tk.Label(
            chart_header,
            text='BDD Residual History',
            font=self.font_section,
            bg=self.palette['panel'],
            fg=self.palette['text'],
            anchor='w',
        ).grid(row=0, column=0, sticky='w')
        self._draw_residual_key(chart_header)
        figure = Figure(figsize=(8.8, 1.1), dpi=100, facecolor=self.palette['panel'])
        self.sub = figure.add_subplot(111)
        figure.subplots_adjust(left=0.055, right=0.985, top=0.92, bottom=0.20)
        self.sub.plot(self.residual_summary)
        
        self.figure_canvas = FigureCanvasTkAgg(figure, master=self.chart_panel)
        self.figure_canvas.draw()
        self.figure_canvas.get_tk_widget().grid(row=1, column=0, sticky='ew')
        # Call color map: Display different color based on the power flow
        self.color_map()
        self._set_badge(self.run_status_label, 'Starting', 'idle')
        self._set_badge(self.bdd_status_label, 'BDD Clear', 'ok')
        self._set_badge(self.attack_status_label, 'No Attack', 'idle')
        self._set_badge(self.mtd_status_label, 'MTD Off', 'idle')
        self.label_threshold_['text'] = f'{np.round(self.case_env.bdd_threshold, 2)}'
        # Call display measurement
        self.show_mea()
        # Call MTD trigger
        self.trigger_mtd()

    def _configure_style(self):
        style = ttk.Style(self.display_window)
        try:
            style.theme_use('clam')
        except tk.TclError:
            pass
        style.configure('TCombobox', padding=4)
        style.configure('TCheckbutton', background=self.palette['panel'], font=self.font_panel_label)
        style.configure('Primary.TButton', font=self.font_panel_button, padding=(9, 7))

    def _metric_row(self, parent, row, name, unit):
        row_frame = tk.Frame(parent, bg=self.palette['panel'])
        row_frame.grid(row=row, column=0, sticky='ew', pady=(0, 6))
        row_frame.columnconfigure(0, weight=1)

        tk.Label(
            row_frame,
            text=name,
            font=self.font_panel_label,
            bg=self.palette['panel'],
            fg=self.palette['muted'],
            anchor='w',
        ).grid(row=0, column=0, sticky='w')

        value = tk.Label(
            row_frame,
            text='--',
            font=self.font_panel_value,
            bg=self.palette['panel'],
            fg=self.palette['text'],
            anchor='e',
        )
        value.grid(row=0, column=1, sticky='e')

        if unit:
            tk.Label(
                row_frame,
                text=unit,
                font=self.font_panel_label,
                bg=self.palette['panel'],
                fg=self.palette['muted'],
                anchor='e',
            ).grid(row=0, column=2, sticky='e', padx=(5, 0))

        return value

    def _set_badge(self, label, text, state):
        colors = {
            'ok': (self.palette['ok_bg'], self.palette['ok_fg']),
            'warn': (self.palette['warn_bg'], self.palette['warn_fg']),
            'alarm': (self.palette['alarm_bg'], self.palette['alarm_fg']),
            'idle': (self.palette['idle_bg'], self.palette['idle_fg']),
            'attack': ('#eadcff', self.palette['attack']),
        }
        bg, fg = colors[state]
        label.config(text=text, bg=bg, fg=fg)

    def _center_grid_diagram(self, event=None):
        if not hasattr(self, 'grid_diagram_bbox') or self.grid_diagram_bbox is None:
            return

        canvas_width = event.width if event is not None else self.sys_canvas.winfo_width()
        canvas_height = event.height if event is not None else self.sys_canvas.winfo_height()
        if canvas_width <= 1 or canvas_height <= 1:
            return

        x1, y1, x2, y2 = self.grid_diagram_bbox
        diagram_width = x2 - x1
        diagram_height = y2 - y1
        x_margin = 12
        y_margin = 4
        available_top = self.grid_legend_bottom + y_margin
        available_bottom = canvas_height - y_margin
        available_height = max(0, available_bottom - available_top)

        target_x = max(x_margin, (canvas_width - diagram_width) / 2)
        target_y = available_top
        if available_height > diagram_height:
            target_y = available_top + (available_height - diagram_height) / 2

        target_offset = (target_x - x1, target_y - y1)
        dx = target_offset[0] - self.grid_diagram_offset[0]
        dy = target_offset[1] - self.grid_diagram_offset[1]
        if abs(dx) < 0.1 and abs(dy) < 0.1:
            return

        self.sys_canvas.move('grid_diagram', dx, dy)
        self.grid_diagram_offset = target_offset

    def _draw_grid_legend(self):
        x1, y1, x2, y2 = 16, 8, 554, 42
        self.sys_canvas.create_rectangle(
            x1,
            y1,
            x2,
            y2,
            fill=self.palette['canvas'],
            outline=self.palette['panel_edge'],
            width=1,
            tags='grid_legend',
        )
        self.sys_canvas.create_text(
            x1 + 12,
            y1 + 17,
            text='Legend',
            anchor='w',
            fill=self.palette['text'],
            font=self.font_badge,
            tags='grid_legend',
        )

        legend_items = [
            (94, 'Normal', self.palette['normal'], 'line'),
            (208, 'High', self.palette['warning'], 'line'),
            (306, 'Alarm', self.palette['alarm'], 'line'),
            (416, 'FDI bus', self.palette['alarm'], 'bus'),
        ]
        y = y1 + 17
        for x, text, color, kind in legend_items:
            if kind == 'line':
                self.sys_canvas.create_line(
                    x,
                    y,
                    x + 34,
                    y,
                    fill=color,
                    width=4,
                    capstyle=tk.ROUND,
                    tags='grid_legend',
                )
            else:
                self.sys_canvas.create_oval(
                    x + 9,
                    y - 8,
                    x + 23,
                    y + 6,
                    fill='#ffffff',
                    outline=color,
                    width=3,
                    tags='grid_legend',
                )
            self.sys_canvas.create_text(
                x + 44,
                y,
                text=text,
                anchor='w',
                fill=self.palette['muted'],
                font=self.font_label,
                tags='grid_legend',
            )

    def _draw_residual_key(self, parent):
        key_canvas = tk.Canvas(
            parent,
            width=492,
            height=24,
            bg=self.palette['panel'],
            highlightthickness=0,
        )
        key_canvas.grid(row=0, column=1, sticky='e')
        key_canvas.create_line(4, 12, 42, 12, fill=self.palette['normal'], width=3)
        key_canvas.create_text(
            50,
            12,
            text='Residual',
            anchor='w',
            fill=self.palette['muted'],
            font=self.font_label,
        )
        key_canvas.create_line(
            138,
            12,
            176,
            12,
            fill=self.palette['alarm'],
            width=3,
            dash=(7, 4),
        )
        key_canvas.create_text(
            184,
            12,
            text='Threshold',
            anchor='w',
            fill=self.palette['muted'],
            font=self.font_label,
        )
        event_items = [
            (296, self.att_choice[1]),
            (360, self.att_choice[2]),
            (418, 'MTD'),
        ]
        for x, event_type in event_items:
            label, color, dash = self._event_style(event_type)
            key_canvas.create_line(
                x,
                5,
                x,
                19,
                fill=color,
                width=2,
                dash=dash,
            )
            key_canvas.create_text(
                x + 10,
                12,
                text=label,
                anchor='w',
                fill=self.palette['muted'],
                font=self.font_label,
            )

    def _event_style(self, event_type):
        styles = {
            self.att_choice[1]: ('Random', '#9333ea', (2, 2)),
            self.att_choice[2]: ('FDI', '#4f46e5', (5, 2)),
            'MTD': ('MTD', self.palette['warning'], (1, 2)),
        }
        return styles.get(event_type, (str(event_type), self.palette['muted'], (3, 2)))

    def _record_plot_events(self):
        current_idx = len(self.residual_summary) - 1
        attack_mode = self.att_select.get()

        if current_idx < 0:
            self.prev_attack_mode = attack_mode
            return

        if attack_mode in self.att_choice[1:] and attack_mode != self.prev_attack_mode:
            self.event_markers.append({'index': current_idx, 'type': attack_mode})

        self.prev_attack_mode = attack_mode

    def _draw_event_markers(self, ticks_start, ticks_end):
        for event in self.event_markers:
            if not ticks_start <= event['index'] < ticks_end:
                continue
            _, color, dash = self._event_style(event['type'])
            self.sub.axvline(
                x=event['index'] - ticks_start,
                color=color,
                linewidth=1.5,
                linestyle=(0, dash),
                alpha=0.9,
                zorder=2,
            )

    def _loading_color(self, ratio):
        if ratio >= 1:
            return self.palette['alarm']
        if ratio >= 0.85:
            return self.palette['warning']
        color_idx = min(self.color_no - 1, max(0, int(ratio * (self.color_no - 1))))
        return self.color_map_list[color_idx]

    def opf_loop(self):
        """
        For each display circle, call the opf_loop func to calculate the opf and display the result
        If there is an attack, then we display the attacked measurement
        """
        # Run OPF and receive the results
        result = self.case_env.run_opf(opf_idx = self.opf_idx)
        self.opf_success = bool(result["success"])
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

        # Draw and set the residual plot
        self.plot_residual()
    
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
            # Deselect the MTD check button and disable the selection.
            self.mtd_on.set(0)
            self.check_mtd.config(state=tk.DISABLED)

        elif self.att_select.get() == self.att_choice[1]:
            """
            RANDOM ATTACK
            """
            z_att_noise = self.case_env.gen_ran_att(z_noise, self.att_ratio_max)
            
            z_noise = z_att_noise     # pass as z_noise for consistency

            # Deselect the MTD check button and disable the selection.
            self.mtd_on.set(0)
            self.check_mtd.config(state=tk.DISABLED)
        
        elif self.att_select.get() == self.att_choice[2]:
            """
            FDI ATTACK
            """
            # Set the MTD trigger check box to NORMAL
            self.check_mtd.config(state=tk.NORMAL)

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
            self.attacker_residual = residual
            z_est = self.case_env_ori.h_x_pypower(v_est)
            
            # Generate FDI attack
            v_att, self.att_posi = self.case_env_ori.gen_fdi_att(v_est, att_spec)
            
            # Calculate the attacked measurement
            z_att_est = self.case_env_ori.h_x_pypower(v_att)
            z_att_noise = z_noise + z_att_est - z_est
            
            z_noise = z_att_noise    # pass as z_noise for consistency

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

        return load_active, load_reactive, pf_actual, pf*self.case_env.case['baseMVA'], pl*self.case_env.case['baseMVA'], residual

    def trigger_mtd(self):
        """
        Update/Set to original reactance in the case file
        """

        mtd_state = self.mtd_on.get()
        if mtd_state == 1 and self.prev_mtd_state != 1:
            self.event_markers.append({'index': len(self.residual_summary), 'type': 'MTD'})

        if mtd_state == 1:
            # MTD is triggered, update the reactance
            increase_decrease = np.random.randint(0,2,size = self.case_env.no_brh)*2-1   # -1 or 1
            ratio = self.min_reac_ratio + (self.max_reac_ratio - self.min_reac_ratio)*np.random.rand(self.case_env.no_brh,)
            reac_new = self.case_env.reactance_ori * (1+ratio*increase_decrease)
            self.case_env.update_reactance(reac_new)
        
        else:
            # MTD is turned off, back to the default situation
            self.case_env.update_reactance(self.case_env.reactance_ori)

        self.prev_mtd_state = mtd_state

    def plot_residual(self):
        """
        Plot the estimation residual
        """ 
        
        # Draw and set the residual plot
        self.att_record.append(self.att_select.get())
        self._record_plot_events()
        self.sub.clear()
        ticks_start = 0
        ticks_end = len(self.residual_summary)
        if len(self.residual_summary) >= self.residual_plot_len:
            residual_summary = self.residual_summary[-self.residual_plot_len:]
            ticks_start, ticks_end = len(self.residual_summary) - self.residual_plot_len, len(self.residual_summary)
        else:
            residual_summary = self.residual_summary

        self.sub.set_facecolor(self.palette['panel'])
        self.sub.set_axisbelow(True)
        self._draw_event_markers(ticks_start, ticks_end)
        self.sub.plot(residual_summary, label='Residual', color=self.palette['normal'], zorder=3)
        self.sub.hlines(
            xmin=0,
            xmax=max(len(residual_summary)-1, 0),
            y=self.case_env.bdd_threshold,
            label='Threshold',
            color=self.palette['alarm'],
            linestyles='--',
            zorder=2,
        )
        self.sub.set_xticks(np.arange(len(residual_summary)))
        self.sub.set_xticklabels(np.arange(ticks_start, ticks_end))
        self.sub.set_ylim(0, np.max([np.max(residual_summary), self.case_env.bdd_threshold+10]))
        self.sub.set_xlabel('')
        self.sub.set_ylabel('')
        self.sub.tick_params(colors=self.palette['muted'])
        self.sub.grid(True, color='#e4e9ee', linewidth=0.8)
        self.sub.spines['top'].set_visible(False)
        self.sub.spines['right'].set_visible(False)
        self.sub.spines['left'].set_color(self.palette['panel_edge'])
        self.sub.spines['bottom'].set_color(self.palette['panel_edge'])
        self.figure_canvas.draw()

    def show_mea(self):
        if self.is_pause == False:
            # Continue
            self.opf_loop()       # Run OPF
            self.trigger_mtd()    # Determine whether the MTD should be triggered or not

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
        self.label_power_active_['text'] = f'{np.round(np.sum(self.load_active), 2)}'
        self.label_power_reactive_['text'] = f'{np.round(np.sum(self.load_reactive), 2)}'
        self.label_power_loss_['text'] = f'{np.round(self.pl, 2)}'
        self.label_residual_['text'] = f'{np.round(self.residual, 2)}'

        is_alarm = self.residual >= self.case_env.bdd_threshold
        self.label_residual_.config(fg=self.palette['alarm'] if is_alarm else self.palette['text'])
        self._set_badge(self.run_status_label, 'Running' if self.opf_success else 'OPF Failed', 'ok' if self.opf_success else 'alarm')
        self._set_badge(self.bdd_status_label, 'BDD Alarm' if is_alarm else 'BDD Clear', 'alarm' if is_alarm else 'ok')

        attack_mode = self.att_select.get()
        if attack_mode == self.att_choice[0]:
            self._set_badge(self.attack_status_label, 'No Attack', 'idle')
        elif attack_mode == self.att_choice[1]:
            self._set_badge(self.attack_status_label, 'Random Attack', 'attack')
        else:
            self._set_badge(self.attack_status_label, 'FDI Attack', 'attack')

        if self.mtd_on.get() == 1:
            self._set_badge(self.mtd_status_label, 'MTD Active', 'warn')
        else:
            self._set_badge(self.mtd_status_label, 'MTD Off', 'idle')

        load_index = getattr(self.case_env, 'current_load_index', '--')
        self.opf_label.config(text=f'OPF {self.opf_idx} | Load {load_index}')
        
        """
        ON-grid measurement update
        """
        
        # Branch
        for i in range(self.case_env.no_brh):
            # Set the branch color according to the active power flow ratio
            brh_ratio = np.abs(self.pf[i]/self.flow_limit[i])           # The load rate of branch
            brh_ratio_actual = np.abs(self.pf_actual[i]/self.flow_limit[i])

            actual_color = self._loading_color(brh_ratio_actual)
            measured_color = self._loading_color(brh_ratio)
            line_ratio = brh_ratio if attack_mode != self.att_choice[0] else brh_ratio_actual
            line_color = measured_color if attack_mode != self.att_choice[0] else actual_color
            line_width = 6 if line_ratio >= 1 else 5 if line_ratio >= 0.85 else 4
            self.sys_canvas.itemconfig(tagOrId=f'line_{i}', fill=line_color, width=line_width)
            
            # Change the on-grid active flow load value
            if attack_mode == self.att_choice[0]:
                # No att
                # Normal operation
                self.sys_canvas.itemconfig(
                    tagOrId=f'line_text_below_{i}',
                    fill=actual_color,
                    text=f'{int(brh_ratio_actual*100)}%',
                    font=self.font_grid_small,
                )

                # Clear the Attack text: above
                self.sys_canvas.itemconfig(
                    tagOrId=f'line_text_above_{i}',
                    fill=actual_color,
                    text='',
                    font=self.font_grid_small,
                )

            else:
                # With attacks
                # Also display the actual measurement
                # Normal 
                self.sys_canvas.itemconfig(
                    tagOrId=f'line_text_below_{i}',
                    fill=actual_color,
                    text=f'{int(brh_ratio_actual*100)}%',
                    font=self.font_grid_small,
                )
                # Attack
                self.sys_canvas.itemconfig(
                    tagOrId=f'line_text_above_{i}',
                    fill=self.palette['attack'],
                    text=f'{int(brh_ratio*100)}%',
                    font=self.font_grid_small,
                )

        # Bus
        attacked_buses = set(self.att_posi) if attack_mode == self.att_choice[-1] else set()
        for i in range(self.case_env.no_bus):
            if i in attacked_buses:
                self.sys_canvas.itemconfig(
                    tagOrId=f'bus_{i}',
                    fill=self.palette['alarm_bg'],
                    outline=self.palette['alarm'],
                    width=3,
                )
            else:
                self.sys_canvas.itemconfig(
                    tagOrId=f'bus_{i}',
                    fill='#ffffff',
                    outline=self.palette['normal'],
                    width=2,
                )
    
    """
    DISPLAY FUNCS
    """

    def color_map(self, color_name = 'magma'):
        """
        Create the color map used for different line active power flow loading
        """
        try:
            cmap = matplotlib.colormaps[self.color_name].resampled(self.color_no)
        except AttributeError:
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
            self.button_pause.config(text='Continue')
            self._set_badge(self.run_status_label, 'Paused', 'idle')
        else:
            # Pause -> Run
            self.is_pause = False
            self.button_pause.config(text='Pause')
            # Recall the function to continue
            self.show_mea()
            
def os_config():
    platformD = system()
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
    
    # Generate load if it does not exist
    _, _ = gen_load(case, 'case14')
    
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

import numpy as np 
import tkinter as tk
import astropy.units as u
import matplotlib.pyplot as plt
import astropy.constants as const
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from binary_model import PulsarBinaryModel

# run:
# binary = InteractiveOrbit(period=5, pulsar_mass=1.4, companion_mass=2, eccentricity=0.5, inclination=90, LoAN=0, AoP=0)
# binary.plot_orbit()

class InteractiveOrbit(PulsarBinaryModel):
    def __init__(self, period=5, pulsar_mass=1.4, companion_mass=2, eccentricity=0.5, inclination=90, LoAN=0, AoP=0):
        super().__init__(period=period/3600, pulsar_mass=pulsar_mass, companion_mass=companion_mass, eccentricity=eccentricity,
                         inclination=inclination, LoAN=LoAN, AoP=AoP)
        self.phase_i = 0
        self.sliders = []
        self.sim_phase = False
        self.sim_scale = 1#u.s.to(u.h)
        self.traj = self.true_anomaly(np.linspace(0, self.period, 1000), time_scale=self.sim_scale)

        self.root = tk.Tk()
        fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
        canvas = FigureCanvasTkAgg(fig, master=self.root)
        canvas.get_tk_widget().grid(row=0, column=0, rowspan=1, columnspan=2, sticky='nsew')

        fig1, ax1 = plt.subplots(figsize=(6,2))
        canvas1 = FigureCanvasTkAgg(fig1, master=self.root)
        canvas1.get_tk_widget().grid(row=1, column=2, rowspan=1, columnspan=6, sticky='sew')
        self.fig = [fig, fig1]
        self.ax = [ax, ax1]
        self.canvas = [canvas, canvas1]
        
    def update_vals(self, values):
        self.mass_p = values[0]
        self.mass_c = values[1]
        self.e = values[2]
        self.I = np.deg2rad(values[3])
        self.LoAN = np.deg2rad(values[4])
        self.AoP = np.deg2rad(values[5])
        
    def plot_data(self, values):
        self.update_vals(values)
        self.ax[0].clear()
        theta_range = np.linspace(0, 2*np.pi, 100)
        self.ax[0].plot(*self.star_coord(theta_range, star='pulsar')/self.a, 'b-')
        self.ax[0].plot(*self.star_coord(theta_range, star='companion')/self.a, 'r-')
        self.ax[0].plot(*self.star_coord(self.traj[self.phase_i], star='pulsar')/self.a, 'b.', markersize=8)
        self.ax[0].plot(*self.star_coord(self.traj[self.phase_i], star='companion')/self.a, 'r.', markersize=8)
        view_range = [-1.2, 1.2]
        self.ax[0].set_xlim(view_range)
        self.ax[0].set_ylim(view_range)
        self.ax[0].set_zlim(view_range)
        self.ax[0].set_zlabel(r'to observer $\longrightarrow$')
        self.ax[0].set_xlabel('RA (J2000)')
        self.ax[0].set_ylabel('DEC (J2000)')
        self.ax[0].set_xticklabels([])
        self.ax[0].set_yticklabels([])
        self.ax[0].set_zticklabels([])
        self.canvas[0].draw()
    
    def plot_delay(self, values):
        self.update_vals(values)
        self.ax[1].clear()
        t_range = np.linspace(0, self.period, 1000)
        delay = -self.star_coord(self.true_anomaly(t_range, time_scale=self.sim_scale))[2]/const.c.value
        self.ax[1].plot(t_range/self.period, delay, 'b-')
        self.ax[1].plot(self.phase_i/1000, -self.star_coord(self.true_anomaly(self.phase_i/1000*self.period, time_scale=self.sim_scale))[2]/const.c.value,
                    'b.', markersize=8)
        self.ax[1].set_ylabel('Römer delay (s)')
        self.ax[1].set_xlabel(f'Phase ({self.period} hr)')
        self.fig[1].subplots_adjust(bottom=0.3)
        self.canvas[1].draw()

    def animate_phase(self):
        if self.sim_phase:
            self.phase_i += round(10 * 20/self.period)
            if self.phase_i == 1000:
                self.phase_i = 0
            self.plot_data([slider.get() for slider in self.sliders])
            self.plot_delay([slider.get() for slider in self.sliders])
        self.root.after(10, self.animate_phase)
    
    def create_slider(self, variable, label, low, high, res, col_i, color='black'):
        slider = tk.Scale(self.root, from_=low, to=high, length=200, orient='vertical', 
                            command=self.update_plot, resolution=res, label=label, foreground=color)
        slider.grid(row=0, column=col_i,rowspan=1, padx=0, pady=0, sticky='ns')
        slider.set(variable)
        self.sliders.append(slider)

    def update_plot(self, val):
        if self.sliders:
            self.traj = self.true_anomaly(np.linspace(0, self.period, 1000), time_scale=self.sim_scale)
            self.plot_data([slider.get()for slider in self.sliders])
            self.plot_delay([slider.get()for slider in self.sliders])
    
    def toggle(self, val):
        self.sim_phase = not self.sim_phase
        
    def plot_orbit(self):

        self.create_slider(self.mass_p, r'Mp', 2.2, 0.8, 0.1, 2, color='blue')
        self.create_slider(self.mass_c, r'Mc', 5, 0.2, 0.1, 3, color='red')
        self.create_slider(self.e, 'ecc', 0.95, 0, 0.01, 4)
        self.create_slider(np.rad2deg(self.I), 'Inc', 90, 0, 1, 5)
        self.create_slider(np.rad2deg(self.LoAN), 'LoAN', 360, 0, 1, 6)
        self.create_slider(np.rad2deg(self.AoP), 'AoP', 360, 0, 1, 7)
        
        self.plot_data([slider.get()for slider in self.sliders])
        self.plot_delay([slider.get()for slider in self.sliders])
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)
        self.root.bind("<KeyPress-space>", self.toggle)
        self.animate_phase()

        self.root.mainloop()

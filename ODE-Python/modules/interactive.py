import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as lin
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, TextBox

import modules.PATHDEF as PATHDEF
import modules.CRSCDEF as CRSCDEF
from modules.spring import Spring
from modules.utils import deg2rad
import modules.materials as materials

from typing import Optional
from datetime import datetime

class Interactive_Spring(Spring):
    def __init__(self, path: PATHDEF.RadiallyEndedPolynomial, 
                       crsc: CRSCDEF.Piecewise_Ic_Control, 
                       matl: materials.Material,
                       torqueCapacity=3000,
                       name=None):
        # print("YOU SHOULD SEE THIS")
        self.version = 0
        self.baseName = name
        self.stacker=0.03
        self.center = (1-0.25)/2

        # self.path     = path
        # self.crsc     = crsc
        # self.material = matl

        super().__init__(path, crsc, matl, torqueCapacity=torqueCapacity, name=name)

        self.fig, self.ax = plt.subplots(figsize=(15,15))
        
        self.adjusters = []

        self.x = self.path.pts[:,0]
        self.y = self.path.pts[:,1]
        
        # pts = np.empty([len(self.x),2])
        for i in range(len(self.x)):
            self.path.pts[i,0] = self.x[i]
            self.path.pts[i,1] = self.y[i]

        # self.distance1 = np.sqrt((self.x[1]-self.x[0])**2+(self.y[1]-self.y[0])**2)
        # self.distance2 = np.sqrt((self.x[2]-self.x[1])**2+(self.y[2]-self.y[1])**2)
        self.initialize_plot()

        self.sliderHeigths = []
        # Create the necessary adjustors
        for index, pathPt in np.ndenumerate(path.XYFactors):
            name = "arcLenStep"+str(index[0]+1)
            a, b = self.add_adjuster(name, [0, 1], pathPt, columns=[3, 0])
            setattr(self, name+"Slider", a)
            setattr(self, name+"Text", b)
        for index, IcPt in np.ndenumerate(self.crsc.IcFactors):
            name = "IcArcLen"+str(index[0]+1)
            a, b = self.add_adjuster(name, [0, 1], IcPt, columns=[3, 0])
            setattr(self, name+"Slider", a)
            setattr(self, name+"Text", b)
        self.IcSliders = []
        self.sliderHeigths.append(self.stacker)
        self.stacker = 0.03
        for index, IcVal in np.ndenumerate(self.crsc.IcPts):
            name = "IcValue"+str(index[0]+1)
            max = 3*IcVal
            a, b = self.add_adjuster(name, [0, max], IcVal, scale='linear', columns=[3, 1])
            setattr(self, name+"Slider", a)
            setattr(self, name+"Text", b)
            self.IcSliders.append(self.adjusters[-1])
            self.sliderHeigths.append(self.stacker)
        self.sliderHeigths.append(self.stacker)
        self.stacker = 0.03
        for index, alphaVal in np.ndenumerate(self.path.alphaAngles):
            name = "AlphaAngle"+str(index[0]+1)
            max = 90*deg2rad
            a, b = self.add_adjuster(name, [0, max], alphaVal, scale='linear', columns=[3,2])
            setattr(self, name+"Slider", a)
            setattr(self, name+"Text", b)

        # ax_checkbox=plt.axes([self.center-0.1, bottomIc, 0.1, self.stacker-bottomIc])
        # ax_checkbox.patch.set_alpha(0)
        # self.IcCheckBoxes=CheckButtons(ax_checkbox, ["","",""])
        # self.IcCheckBoxes.on_clicked(self.link_sliders)

        # Make room for however many adjustors you made
        plt.subplots_adjust(bottom=np.max(self.sliderHeigths)+0.03)

        self.free_pts = [i for i in range(1, len(self.x) - 1)]
        self.constr_pts = [len(self.x)-1]
        self.dragging = False
        self.index = None

        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)

        # Create the print button
        ax_button = plt.axes([0.05, 0.9, 0.1, 0.05])  # Position of the button
        self.print_button = Button(ax_button, 'Print Points')
        self.print_button.on_clicked(self.print_points)  # Connect to the print method

        # Create the surface export button
        ax_button1 = plt.axes([0.05, 0.84, 0.1, 0.05])  # Position of the button
        self.export_button = Button(ax_button1, 'Export Surfaces')
        self.export_button.on_clicked(self.export_surfaces)  # Connect to the export method

        # Create the save spring button
        ax_buttonS = plt.axes([0.05, 0.78, 0.1, 0.05])
        self.save_button = Button(ax_buttonS, 'Save Spring')
        self.save_button.on_clicked(self.export_parameters)  # Connect to the export method

        # Create the deform button
        ax_button2 = plt.axes([0.85, 0.9, 0.1,0.05])
        self.deform_button = Button(ax_button2, 'Deform Spring')
        self.deform_button.on_clicked(self.deform_interactive)

        # Create deform output text boxes
        ax_deform_text = plt.axes([0.8,0.8,0.1,0.03])
        ax_stress_text = plt.axes([0.8,0.8-0.06,0.1,0.03])
        self.deflection_textbox = TextBox(ax_deform_text, "deformation angle", initial=0)
        self.max_stress_textbox = TextBox(ax_stress_text, "maximum stress   ", initial=0)
        self.fig.text(0.8, 0.8-0.1, "design stress:"+str(self.designStress))

    def export_surfaces(self, event):
        # print("YOU SHOULD ALSO SEE THIS")
        date = datetime.now().strftime("%Y%m%d")
        self.name = date+self.baseName+str(self.version)
        self.version+=1
        # print(self.name)
        self.A, self.B = self.crsc.get_outer_geometry(self.resl)
        return super().export_surfaces()

    def export_parameters(self):
        date = datetime.now().strftime("%Y%m%d")
        self.name = date+self.baseName+str(self.version)
        return super().export_parameters()

    def deform_interactive(self, event):
        print("attempting to deform current spring")
        print(self.torqueCapacity)
        self.deform_by_torque_predict_forces(self.torqueCapacity, self.deform_ODE)
        
        print("Deflection:", self.dBeta/deg2rad)
        self.update_text_box(self.dBeta/deg2rad, self.deflection_textbox)
        # self.update_text_box(self.dBeta, self.deform_textbox)
        if hasattr(self, "ax_deform_plot"):
            self.ax_deform_plot.clear()
            self.ax_deform_plot.remove()
        self.ax_deform_plot = plt.axes([0.75, np.max(self.sliderHeigths)+0.03, 0.25-.03, 0.7-(np.max(self.sliderHeigths)+0.06)])
        self.plot_deform(True, targetAxes=self.ax_deform_plot)
        print("Max Stress:", self.maxStress)
        self.update_text_box(self.maxStress, self.max_stress_textbox)

    def initialize_plot(self):
        curveMesh = np.linspace(0,self.path.arcLen,self.path.measureResolution)

        self.path.xy_poly()
        self.curve, = self.ax.plot(self.path.get_xy_n(curveMesh, 'x'), self.path.get_xy_n(curveMesh, 'y'))
        self.a_surface, self.b_surface = self.crsc.get_outer_geometry(len(curveMesh))
        self.inner, = self.ax.plot(self.a_surface[:,0], self.a_surface[:,1])
        self.outer, = self.ax.plot(self.b_surface[:,0], self.b_surface[:,1])

        self.ax.plot([0,self.outerRadius*np.cos(2*np.pi/self.n)],
                     [0,self.outerRadius*np.sin(2*np.pi/self.n)],'k')
        self.ax.plot([0,self.outerRadius],[0,0],'k')

        self.line, = self.ax.plot(self.x, self.y, 'o-')
        outerCircle = plt.Circle([0,0],self.outerRadius,color ="k",fill=False)
        innerCircle = plt.Circle([0,0],self.innerRadius,color ="k",fill=False)
        self.ax.add_patch(outerCircle)
        self.ax.add_patch(innerCircle)
        
        self.outerCircle = self.semicircle
        self.curveMesh = curveMesh
        
    def add_adjuster(self, label, range, value, scale='linear', columns = None):
        
        if columns is None:
            center = self.center
        else:
            center = (1-0.25*columns[0])/2+0.25*columns[1]

        ax_slider = plt.axes([center, self.stacker, 0.10, 0.03])
        if scale=='log':
            slider = Sliderlog(ax_slider, label, range[0], range[1], valinit=value)
        else:
            slider = Slider(ax_slider, label, range[0], range[1], valinit=value)
        slider.valtext.set_visible(False)
        ax_text = plt.axes([center+0.10+.02, self.stacker, .05, .03])
        textbox = TextBox(ax_text, label, initial=str(slider.valinit))
        textbox.label.set_visible(False)

        slider.on_changed(self.update_curve)  # Connect slider to update method
        slider.on_changed(lambda val: self.update_text_box(val, textbox))
        textbox.on_submit(lambda text: self.update_slider(text, slider, textbox))
        
        self.stacker+=.05

        self.adjusters.append(slider)

        return slider, textbox

    def update_slider(self, text, slider, textbox):
        try: 
            val = float(text)
            if slider.valmin <= val <= slider.valmax:
                slider.set_val(val)
            else:
                textbox.set_val(str(slider.val))
        except ValueError:
            textbox.set_val(str(slider.val))
        
    def update_text_box(self, val, textbox):
        try:
            form = "{set: .7f}"
            textbox.set_val(form.format(set = val))
        except:
            pass

    def link_sliders(self):
        active_buttons = self.IcCheckBoxes.get_status()
        checked_indexes = [i for i, checked in enumerate(active_buttons) if checked]
        for j in checked_indexes:
            self.IcSliders[j].on_changed(self.IcSliders[j].set_val(self.IcSliders[0]))

    def print_points(self, event):
        # for i in range(len(self.x)):
        #     print(f'Point {i}: ({self.x[i]}, {self.y[i]})')
        # return self.x, self.y
        self.calc_angles_radii()
        for i in range(len(self.x)):
            print(f'Point {i}: (radius {self.radii[i]}, at angle {self.angles[i]/deg2rad})')
    
    def calc_angles_radii(self):
        self.radii  = []
        self.angles = []
        for i in range(len(self.x)):
            self.radii.append(lin.norm([self.x[i],self.y[i]]))
            self.angles.append(np.atan2(self.y[i],self.x[i]))
            # print(f'Point {i}: (radius {self.radii[i]}, at angle {self.angles[i]/deg2rad})')
        return self.radii, self.angles
        
    def update_curve(self, _):

        # print("YOU SHOULD SEE THIS")
        vals = [slider.val for slider in self.adjusters]
        assert len(vals)==len(self.path.XYFactors)+len(self.crsc.IcFactors)+len(self.crsc.IcPts)+len(self.path.alphaAngles)

        # update control points based on drag
        for i in range(len(self.x)):
            self.path.pts[i,0] = self.x[i]
            self.path.pts[i,1] = self.y[i]

        # update distortion effects based on sliders
        for index in range(len(self.path.XYFactors)):
            # print(self.path.XYParamLens)
            self.path.XYParamLens[index+1] = vals[index]*self.path.arcLen
            # print(self.path.XYParamLens)
        # update Ic positions based on sliders
        for index in range(len(self.crsc.IcFactors)):
            self.crsc.IcParamLens[index+1] = vals[index+len(self.path.XYFactors)]*self.crsc.arcLen
        # update Ic locations based on sliders
        for index in range(len(self.crsc.IcPts)):
            self.crsc.IcPts[index] = vals[index+len(self.path.XYFactors)+len(self.crsc.IcFactors)]
        # update alpha angles based on sliders
        for index in range(len(self.path.alphaAngles)):
            self.path.alphaAngles[index] = vals[index+len(self.path.XYFactors)+len(self.crsc.IcFactors)+len(self.crsc.IcPts)]

        # recalcualte necessary stuff for plotting
        self.path.XCoeffs, self.path.YCoeffs  = self.path.xy_poly()
        self.crsc.IcCoeffs, self.crsc.domains = self.crsc.Ic_multiPoly()

        self.curveMesh = np.linspace(0,self.path.arcLen,500)       
        self.a_surface, self.b_surface = self.crsc.get_outer_geometry(len(self.curveMesh))
        
        
        self.curve.set_data(self.path.get_xy_n(self.curveMesh, 'x'), self.path.get_xy_n(self.curveMesh, 'y'))   
        self.inner.set_data(self.a_surface[:,0], self.a_surface[:,1])
        self.outer.set_data(self.b_surface[:,0], self.b_surface[:,1])

        self.fig.canvas.draw_idle()

    def on_press(self, event):
        if event.inaxes != self.ax:
            return
        contains, ind = self.line.contains(event)
        if contains:
            self.index = ind['ind'][0]
            if self.index in self.free_pts or self.index in self.constr_pts:
                self.dragging = True

    def on_release(self, event):
        self.dragging = False
        self.index = None

    def on_motion(self, event):
        if self.dragging and self.index is not None:
            if self.index in self.constr_pts:
                # print("you touched a constr point")
                self.x[self.index] = event.xdata
                self.y[self.index] = self.semicircle(self.x[self.index])
            else:
                self.x[self.index] = event.xdata
                self.y[self.index] = event.ydata

            self.line.set_data(self.x, self.y)
            self.update_curve("garbage")
            # self.ax.relim()  # Recompute limits
            # self.ax.autoscale_view()  # Autoscale the view
            self.fig.canvas.draw_idle()

    def semicircle(self, x):
        # print("called")
        return np.sqrt(self.path.outerRadius**2-x**2)

class Sliderlog(Slider):

    """Logarithmic slider.

    Takes in every method and function of the matplotlib's slider.

    Set slider to *val* visually so the slider still is lineat but display 10**val next to the slider.

    Return 10**val to the update function (func)"""

    def set_val(self, val):

        if self.orientation == 'vertical':
            self.poly.set_height(val - self.poly.get_y())
            self._handle.set_ydata([val])
        else:
            self.poly.set_width(val - self.poly.get_x())
            self._handle.set_xdata([val])
        self.valtext.set_text(self._format(val) % 10**val)
        if self.drawon:
            self.ax.figure.canvas.draw_idle()
        self.val = val
        if self.eventson:
            self._observers.process('changed', 10**val)
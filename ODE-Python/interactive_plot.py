import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as lin
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, TextBox

import modules.PATHDEF as PATHDEF
import modules.CRSCDEF as CRSCDEF
from modules.utils import deg2rad

class DraggableSpring:
    def __init__(self, path: PATHDEF.RadiallyEndedPolynomial, crsc: CRSCDEF.Piecewise_Ic_Control):
        
        self.stacker=0.03
        self.center = (1-0.25)/2

        self.path = path
        self.crsc = crsc

        self.fig, self.ax = plt.subplots(figsize=(8,10))
        
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

        # Create the necessary adjustors
        for index, pathPt in np.ndenumerate(path.XYFactors):
            name = "arcLenStep"+str(index[0]+1)
            a, b = self.add_adjuster(name, [0, 1], pathPt)
            setattr(self, name+"Slider", a)
            setattr(self, name+"Text", b)
        for index, IcPt in np.ndenumerate(self.crsc.IcFactors):
            name = "IcArcLen"+str(index[0]+1)
            a, b = self.add_adjuster(name, [0, 1], IcPt)
            setattr(self, name+"Slider", a)
            setattr(self, name+"Text", b)
        bottomIc = self.stacker
        self.IcSliders = []
        for index, IcVal in np.ndenumerate(self.crsc.IcPts):
            name = "IcValue"+str(index[0]+1)
            max = 3*IcVal
            a, b = self.add_adjuster(name, [0, max], IcVal, scale='linear')
            setattr(self, name+"Slider", a)
            setattr(self, name+"Text", b)
            self.IcSliders.append(self.adjusters[-1])

        # ax_checkbox=plt.axes([self.center-0.1, bottomIc, 0.1, self.stacker-bottomIc])
        # ax_checkbox.patch.set_alpha(0)
        # self.IcCheckBoxes=CheckButtons(ax_checkbox, ["","",""])
        # self.IcCheckBoxes.on_clicked(self.link_sliders)

        # Make room for however many adjustors you made
        plt.subplots_adjust(bottom=self.stacker+0.03)

        self.free_pts = [i for i in range(1, len(self.x) - 1)]
        self.constr_pts = [len(self.x-1)]
        self.dragging = False
        self.index = None

        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)

        # Create the print button
        ax_button = plt.axes([0.05, 0.9, 0.1, 0.05])  # Position of the button
        self.print_button = Button(ax_button, 'Print Points')
        self.print_button.on_clicked(self.print_points)  # Connect to the print method

    def initialize_plot(self):
        curveMesh = np.linspace(0,self.path.arcLen,self.path.measureResolution)

        self.path.xy_poly()
        self.curve, = self.ax.plot(self.path.get_xy_n(curveMesh, 'x'), self.path.get_xy_n(curveMesh, 'y'))
        self.a_surface, self.b_surface = self.crsc.get_outer_geometry(len(curveMesh))
        self.inner, = self.ax.plot(self.a_surface[:,0], self.a_surface[:,1])
        self.outer, = self.ax.plot(self.b_surface[:,0], self.b_surface[:,1])

        self.line, = self.ax.plot(self.x, self.y, 'o-')
        
        self.outerCircle = self.semicircle
        self.curveMesh = curveMesh
        
    def add_adjuster(self, label, range, value, scale='linear'):
        
        ax_slider = plt.axes([self.center, self.stacker, 0.25, 0.03])
        if scale=='log':
            slider = Sliderlog(ax_slider, label, range[0], range[1], valinit=value)
        else:
            slider = Slider(ax_slider, label, range[0], range[1], valinit=value)
        slider.valtext.set_visible(False)
        ax_text = plt.axes([self.center+0.25+.02, self.stacker, .1, .03])
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

        vals = [slider.val for slider in self.adjusters]
        assert len(vals)==len(self.path.XYFactors)+len(self.crsc.IcFactors)+len(self.crsc.IcPts)

        # update control points based on drag
        for i in range(len(self.x)):
            self.path.pts[i,0] = self.x[i]
            self.path.pts[i,1] = self.y[i]

        # update distortion effects based on sliders
        for index in range(len(self.path.XYFactors)):
            self.path.XYParamLens[index+1] = vals[index]*self.path.arcLen
        # update Ic positions based on sliders
        for index in range(len(self.crsc.IcFactors)):
            self.crsc.IcParamLens[index+1] = vals[index+len(self.path.XYFactors)]*self.crsc.arcLen
        # update Ic locations based on sliders
        for index in range(len(self.crsc.IcPts)):
            self.crsc.IcPts[index] = vals[index+len(self.path.XYFactors)+len(self.crsc.IcFactors)]

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
                self.x[self.index] = event.xdata
                self.y[self.index] = self.outerCircle(self.x[self.index])
            else:
                self.x[self.index] = event.xdata
                self.y[self.index] = event.ydata

            self.line.set_data(self.x, self.y)
            self.update_curve("piss")
            # self.ax.relim()  # Recompute limits
            # self.ax.autoscale_view()  # Autoscale the view
            self.fig.canvas.draw_idle()

    def semicircle(self, x):
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


path = PATHDEF.RadiallyEndedPolynomial(2,1)
crsc = CRSCDEF.Piecewise_Ic_Control(path, 0.375, IcPts=np.array([0.008,0.0005,0.008]), IcParamLens=np.array([0.6]))

x=[]
y=[]
for i in range(len(path.pts)):
    x.append(path.pts[i,0])
    y.append(path.pts[i,1])

# x = [1,0.5,.16,-1.5]
# y = [0,0.5,1.52,2]

# def dumbCurve(x):
#     return abs(x)

radius = path.outerRadius

# plot = DraggableSpring(x, y, np.linspace(0,path.arcLen,500), [1, 2], [3], semicircle)
plot = DraggableSpring(path, crsc)
plot.ax.set_xlim(-radius*1.1,radius*1.1)
plot.ax.set_ylim(-radius*1.1,radius*1.1)
plot.ax.set_aspect('equal')
plt.show()
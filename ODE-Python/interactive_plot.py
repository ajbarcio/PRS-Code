import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as lin
from matplotlib.widgets import Slider, Button, TextBox

import modules.PATHDEF as PATHDEF
import modules.CRSCDEF as CRSCDEF
from modules.utils import deg2rad

class DraggableSpring:
    def __init__(self, x, y, curveMesh, free_pts, constr_pts, outerCircle):
        
        self.stacker=0

        self.fig, self.ax = plt.subplots()
        plt.subplots_adjust(bottom=0.3)
        
        self.adjusters = []

        self.x = x
        self.y = y
        
        # pts = np.empty([len(self.x),2])
        for i in range(len(self.x)):
            path.pts[i,0] = self.x[i]
            path.pts[i,1] = self.y[i]

        # self.distance1 = np.sqrt((self.x[1]-self.x[0])**2+(self.y[1]-self.y[0])**2)
        # self.distance2 = np.sqrt((self.x[2]-self.x[1])**2+(self.y[2]-self.y[1])**2)
        self.initialize_plot(curveMesh, free_pts, constr_pts, outerCircle)

        # Create the necessary adjustors
        for index, pathPt in np.ndenumerate(path.XYFactors):
            name = "arcLenStep"+str(index[0]+1)
            a, b = self.add_adjuster(name, [0, 1], pathPt)
            setattr(self, name+"Slider", a)
            setattr(self, name+"Text", b)
        for index, IcPt in np.ndenumerate(crsc.IcFactors):
            name = "IcArcLen"+str(index[0]+1)
            a, b = self.add_adjuster(name, [0, 1], IcPt)
            setattr(self, name+"Slider", a)
            setattr(self, name+"Text", b)
        for index, IcVal in np.ndenumerate(crsc.IcPts):
            name = "IcValue"+str(index[0]+1)
            max = 10*IcVal
            a, b = self.add_adjuster(name, [0, max], IcVal)
            setattr(self, name+"Slider", a)
            setattr(self, name+"Text", b)
        
        self.free_pts = free_pts
        self.constr_pts = constr_pts
        self.dragging = False
        self.index = None

        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)

        # Create the button
        ax_button = plt.axes([0.05, 0.9, 0.1, 0.05])  # Position of the button
        self.print_button = Button(ax_button, 'Print Points')
        self.print_button.on_clicked(self.print_points)  # Connect to the print method

    def initialize_plot(self, curveMesh, free_pts, constr_pts, outerCircle):
        path.xy_poly()
        self.curve, = self.ax.plot(path.get_xy_n(curveMesh, 'x'), path.get_xy_n(curveMesh, 'y'))
        self.a_surface, self.b_surface = crsc.get_outer_geometry(len(curveMesh))
        self.inner, = self.ax.plot(self.a_surface[:,0], self.a_surface[:,1])
        self.outer, = self.ax.plot(self.b_surface[:,0], self.b_surface[:,1])

        self.line, = self.ax.plot(self.x, self.y, 'o-')
        
        self.outerCircle = outerCircle
        self.curveMesh = curveMesh
        
    def add_adjuster(self, label, range, value):
        center = (1-0.25)/2

        ax_slider = plt.axes([center, self.stacker, 0.25, 0.03])
        slider = Slider(ax_slider, label, range[0], range[1], valinit=value)
        slider.valtext.set_visible(False)
        ax_text = plt.axes([center+0.25+.02, self.stacker, .1, .03])
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

    def print_points(self, event):
        # for i in range(len(self.x)):
        #     print(f'Point {i}: ({self.x[i]}, {self.y[i]})')
        # return self.x, self.y
        self.calc_angles_radii()
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
        assert len(vals)==len(path.XYFactors)+len(crsc.IcFactors)+len(crsc.IcPts)

        # update control points based on drag
        for i in range(len(self.x)):
            path.pts[i,0] = self.x[i]
            path.pts[i,1] = self.y[i]

        # update distortion effects based on sliders
        for index in range(len(path.XYFactors)):
            path.XYParamLens[index+1] = vals[index]*path.arcLen
        # update Ic positions based on sliders
        for index in range(len(crsc.IcFactors)):
            crsc.IcParamLens[index+1] = vals[index+len(path.XYFactors)]*crsc.arcLen
        # update Ic locations based on sliders
        for index in range(len(crsc.IcPts)):
            crsc.IcPts[index] = vals[index+len(path.XYFactors)+len(crsc.IcFactors)]

        # recalcualte necessary stuff for plotting
        path.XCoeffs, path.YCoeffs  = path.xy_poly()
        crsc.IcCoeffs, crsc.domains = crsc.Ic_multiPoly()

        self.curveMesh = np.linspace(0,path.arcLen,500)       
        self.a_surface, self.b_surface = crsc.get_outer_geometry(len(self.curveMesh))
        
        
        self.curve.set_data(path.get_xy_n(self.curveMesh, 'x'), path.get_xy_n(self.curveMesh, 'y'))   
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


path = PATHDEF.RadiallyEndedPolynomial(2,1)
crsc = CRSCDEF.Piecewise_Ic_Control(path, 0.375)

x=[]
y=[]
for i in range(len(path.pts)):
    x.append(path.pts[i,0])
    y.append(path.pts[i,1])

# x = [1,0.5,.16,-1.5]
# y = [0,0.5,1.52,2]

# def dumbCurve(x):
#     return abs(x)

def semicircle(x):
    return np.sqrt(path.outerRadius**2-x**2)

radius = path.outerRadius

plot = DraggableSpring(x, y, np.linspace(0,path.arcLen,500), [1, 2], [3], semicircle)
plot.ax.set_xlim(-radius*1.1,radius*1.1)
plot.ax.set_ylim(-radius*1.1,radius*1.1)
plot.ax.set_aspect('equal')
plt.show()
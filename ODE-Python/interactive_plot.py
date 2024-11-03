import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as lin
from matplotlib.widgets import Slider  # Import the Slider
from matplotlib.widgets import Button

import modules.PATHDEF as PATHDEF
import modules.CRSCDEF as CRSCDEF
from modules.utils import deg2rad

class DraggablePlot:
    def __init__(self, x, y, curveMesh, free_pts, constr_pts, outerCircle):
        
        self.fig, self.ax = plt.subplots()
        
        self.x = x
        self.y = y
        
        # pts = np.empty([len(self.x),2])
        for i in range(len(self.x)):
            path.pts[i,0] = self.x[i]
            path.pts[i,1] = self.y[i]

        # self.distance1 = np.sqrt((self.x[1]-self.x[0])**2+(self.y[1]-self.y[0])**2)
        # self.distance2 = np.sqrt((self.x[2]-self.x[1])**2+(self.y[2]-self.y[1])**2)

        path.xy_poly()
        self.curve, = self.ax.plot(path.get_xy_n(curveMesh, 'x'), path.get_xy_n(curveMesh, 'y'))
        self.a_surface, self.b_surface = crsc.get_outer_geometry(len(curveMesh))
        self.inner, = self.ax.plot(self.a_surface[:,0], self.a_surface[:,1])
        self.outer, = self.ax.plot(self.b_surface[:,0], self.b_surface[:,1])

        self.line, = self.ax.plot(self.x, self.y, 'o-')
        
        self.outerCircle = outerCircle
        self.curveMesh = curveMesh
        self.free_pts = free_pts
        self.constr_pts = constr_pts
        self.dragging = False
        self.index = None

        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)

        # Create the slider
        ax_slider = plt.axes([0.2, 0.1, 0.65, 0.03])  # Position of the slider
        self.slider = Slider(ax_slider, 'Second Param', 0.0, 1.0, valinit=0.5)  # Set range and initial value
        self.slider.on_changed(self.update_curve)  # Connect slider to update method

        # Create the button
        ax_button = plt.axes([0.8, 0.01, 0.1, 0.05])  # Position of the button
        self.print_button = Button(ax_button, 'Print Points')
        self.print_button.on_clicked(self.print_points)  # Connect to the print method
    
    def print_points(self, event):
        # for i in range(len(self.x)):
        #     print(f'Point {i}: ({self.x[i]}, {self.y[i]})')
        # return self.x, self.y
        self.print_angles_radii()
    
    def print_angles_radii(self):
        radii  = []
        angles = []
        for i in range(len(self.x)):
            radii.append(lin.norm([self.x[i],self.y[i]]))
            angles.append(np.atan2(self.y[i],self.x[i]))
            print(f'Point {i}: (radius {radii[i]}, at angle {angles[i]/deg2rad})')
        return radii, angles
        
    def update_curve(self, val):

        for i in range(len(self.x)):
            path.pts[i,0] = self.x[i]
            path.pts[i,1] = self.y[i]
        path.XYParamLens = np.array([0,val-0.1,val+0.1,1])*path.arcLen
        path.XCoeffs, path.YCoeffs = path.xy_poly()
        
        self.curveMesh = np.linspace(0,path.arcLen,500)
        self.curve.set_data(path.get_xy_n(self.curveMesh, 'x'), path.get_xy_n(self.curveMesh, 'y'))
        self.a_surface, self.b_surface = crsc.get_outer_geometry(len(self.curveMesh))
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
            self.update_curve(self.slider.val)
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

plot = DraggablePlot(x, y, np.linspace(0,path.arcLen,500), [1, 2], [3], semicircle)
plot.ax.set_xlim(-radius*1.1,radius*1.1)
plot.ax.set_ylim(-radius*1.1,radius*1.1)
plot.ax.set_aspect('equal')
plt.show()
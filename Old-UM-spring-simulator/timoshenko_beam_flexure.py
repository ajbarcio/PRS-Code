import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from materials import *


prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
R = lambda gamma: np.array([[np.cos(gamma), np.sin(gamma)], [-np.sin(gamma), np.cos(gamma)]])
def forward_sum(array):
    ret = np.array(array)
    for i in range(len(ret))[1:]:
        ret[i]+=ret[i-1]
    return ret

IGNORE_TIMOSHENKO_SHEAR = True

class TimoshenkoBeamFlexure(object):
    """ A beam which can be
        geometrically specified
        flexed by a force/torque in local tip coordinates
        measured, post-deflection, including global force (rotated and displaced)
        plotted
        """
    def __init__(self, material):
        self.N=20
        self.Ny=9 # for von-mises graphic
        self.material = material
        self.insertion_depth = 6.080
        self.xgrid = np.zeros((self.N, self.Ny))
        self.ygrid = np.zeros((self.N, self.Ny))
        self.zgrid = np.zeros((self.N, self.Ny))

        self.set_geometry_mm([25.92,20,10,5,0], [1.65,1.4494,1.0249,0.7247,0.81/2], 6.35) # 1/4 in thichness, not 5 mm
        self.flex_by_local_tip_force(-20,90)


    def set_geometry_mm(self, xs, ys, z):
        # Beam definition
        self.z_depth = z # mm thickness
        self.tip_r = ys[-1] #0.81/2
        self.length = xs[0] # 25.92 # mm
        self.xs = np.linspace(0,self.length, self.N)
        self.ys = np.interp(self.xs, self.length-np.array(xs),np.array(ys))

        self.Iyys = self.z_depth*1e-3 *1/12.* (2*self.ys*1e-3)**3 # second moments of area for each segment in m^4
        self.As = self.z_depth*1e-3 * (self.ys*1e-3) # area in m^2

    def flex_by_local_tip_force(self, fx, fy, m=0):
        # construct these by translating the force/moment backwards from the tip
        self.local_xy_force = np.zeros((self.N,2))
        self.local_moment = np.zeros((self.N,))
        self.dphis = np.array(self.xs*0+.00)
        self.dzs = np.array(self.xs*0+0.00)
        self.dxs = np.array(self.xs*0+0.00)
        self.local_xy_force[-1,:]=[fx,fy] # IC (for reverse integration)
        self.local_moment[-1]=m # IC (for reverse integration)
        self.safety_margins = np.zeros((self.N,))

        kappa = 5./6. # Timoshenko shear coefficient
        for i in list(reversed(range(self.N)[1:])):
            dx = (self.xs[i]-self.xs[i-1]) *1e-3 # meters
            self.dphis[i] = self.local_moment[i]/(self.material.E*self.Iyys[i])*dx
            self.dzs[i] = self.local_xy_force[i,1] / (self.material.G*kappa*self.As[i])*dx*1e3 # in units of mm
            if IGNORE_TIMOSHENKO_SHEAR: self.dzs[i]*=0
            self.dxs[i] = self.local_xy_force[i,0] / (self.material.E*self.As[i])*dx*1e3 # in units of mm
            if IGNORE_TIMOSHENKO_SHEAR: self.dxs[i]*=0
            self.local_xy_force[i-1,:] = self.local_xy_force[i,:].dot(R(self.dphis[i]).T)
            self.local_moment[i-1] = self.local_moment[i] + dx*self.local_xy_force[i,1]
            # calculate von mises yield
            sigma_xx_bending = self.local_moment[i]/self.Iyys[i]*(self.ys[i]*1e-3)
            sigma_xx_force = self.local_xy_force[i,0]/self.As[i]
            sigma_xy_force = self.local_xy_force[i,1]/self.As[i]
            sigma_vm = np.sqrt((abs(sigma_xx_bending)+abs(sigma_xx_force))**2+3*sigma_xy_force**2)
            for j, y in enumerate(np.linspace(-self.ys[i], self.ys[i], self.Ny)):
                self.zgrid[i,j] = np.sqrt((sigma_xx_bending*y/self.ys[i]+sigma_xx_force)**2+3*sigma_xy_force**2)
            self.safety_margins[i]=self.material.yield_stress/sigma_vm-1
        self.total_safety_margin = min(self.safety_margins[1:])

        self.phis = forward_sum(self.dphis)
        self.xy_mid = np.zeros((self.N,2))
        self.xy_mid[0,:] = [0,0] # IC for forward integration
        for i in range(self.N)[1:]:
            self.xy_mid[i,:] = self.xy_mid[i-1,:]+np.array([(self.xs[i]-self.xs[i-1])+self.dxs[i],self.dzs[i]]).dot(R(self.phis[i]))
            for j, y in enumerate(np.linspace(-self.ys[i], self.ys[i], self.Ny)):
                self.xgrid[i,j]=self.xy_mid[i,0]-np.sin(self.phis[i])*y
                self.ygrid[i,j]=self.xy_mid[i,1]+np.cos(self.phis[i])*y
        # convert the mid-tip point into polar gear coordinates.
        # print(self.xy_mid[-1,0]-(self.insertion_depth + self.length) , self.xy_mid[-1,1])
        self.global_flexure_tip_xy = self.xy_mid[-1,:]
        flexure_tip_xy = self.xy_mid[-1,0]-(self.insertion_depth + self.length) , self.xy_mid[-1,1]
        flexure_tip_r = np.sqrt((flexure_tip_xy[0])**2+(flexure_tip_xy[1])**2)
        self.global_xy_tip_force = self.local_xy_force[-1,:].dot(R(self.phis[-1] ))
        global_tip_force_angle = np.arctan2(self.global_xy_tip_force[1], self.global_xy_tip_force[0])
        flexure_tip_angle = np.arctan2(flexure_tip_xy[1], flexure_tip_xy[0])
        flexure_tip_force_pressure_angle = np.pi/2-(global_tip_force_angle-flexure_tip_angle)
        #print(self.phis)
        self.final_phi = self.phis[-1]
        flexure_tip_rotation = self.final_phi
        return flexure_tip_r, flexure_tip_force_pressure_angle, flexure_tip_angle, self.global_xy_tip_force, flexure_tip_rotation

    def global_to_local_force(self, global_xy_tip_force):
        local_xy_tip_force = global_xy_tip_force.dot(R(self.phis[-1]).T)
        return local_xy_tip_force

    def draw_tip_force(self, scale=0.01, **kwargs):
        plt.arrow(self.global_flexure_tip_xy[0], self.global_flexure_tip_xy[1], scale*self.global_xy_tip_force[0], scale*self.global_xy_tip_force[1], **kwargs)


    def draw_beam_graphics(self, transform=lambda x,y: (x, y), legend=False, draw_tip=False, draw_unflexed=False):
        # integrate forward to find the position of every point on the beam

        xy_topside = np.zeros((self.N,2))
        xy_botside = np.zeros((self.N,2))
        for i in range(len(self.xs)):
            xy_topside[i,:]=self.xy_mid[i,:]+np.array([0,self.ys[i]]).dot(R(self.phis[i]))
            xy_botside[i,:]=self.xy_mid[i,:]+np.array([0,-self.ys[i]]).dot(R(self.phis[i]))

        # OG flexure shape
        if draw_unflexed:
            plt.plot(*transform(self.xs, self.ys), ":",lw=.5, color=colors[0])
            tip_thetas = np.linspace(-np.pi/2, np.pi/2, 100)
            if draw_tip: plt.plot(*transform(self.length+self.tip_r*np.cos(tip_thetas), self.tip_r*np.sin(tip_thetas)), ":",lw=.5, color=colors[0])
            plt.plot(*transform(self.xs, -self.ys), ":",lw=.5, color=colors[0])

        # Flexed flexure shape
        plt.plot(*transform(xy_topside[:,0], xy_topside[:,1]), zorder=-5, color=colors[0])
        tip_thetas = np.linspace(-np.pi/2+self.phis[-1], np.pi/2+self.phis[-1], 100)
        if draw_tip: plt.plot(*transform(self.xy_mid[-1,0]+self.tip_r*np.cos(tip_thetas),self.xy_mid[-1,1]+self.tip_r*np.sin(tip_thetas)), zorder=-5, color=colors[0])
        plt.plot(*transform(xy_botside[:,0], xy_botside[:,1]), zorder=-5, color=colors[0])

        gear_x_center=self.length+self.insertion_depth
        xy_topside[:,0]-= gear_x_center
        tip_xy = np.array([[self.xy_mid[-1,0]-gear_x_center+self.tip_r*np.cos(theta),self.xy_mid[-1,1]+self.tip_r*np.sin(theta)] for theta in tip_thetas])
        xy_botside[:,0]-= gear_x_center

        # for gamma in np.linspace(0,2*np.pi,24+1)[1:-1]:
        # #     involute_xy.dot(R(gamma))
        #     plt.plot(*transform(gear_x_center + xy_topside.dot(R(gamma))[:,0], xy_topside.dot(R(gamma))[:,1]), zorder=-6, color=colors[0])
        #     plt.plot(*transform(gear_x_center + tip_xy.dot(R(gamma))[:,0], tip_xy.dot(R(gamma))[:,1]), zorder=-6, color=colors[0])
        #     plt.plot(*transform(gear_x_center + xy_botside.dot(R(gamma))[:,0], xy_botside.dot(R(gamma))[:,1]), zorder=-6, color=colors[0])

        if True: # fancy vonmises
            levels = MaxNLocator(nbins=15).tick_values(self.zgrid.min(), self.zgrid.max())


            # pick the desired colormap, sensible levels, and define a normalization
            # instance which takes data values and translates those into levels.
            cmap = plt.get_cmap('viridis')#PiYG
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

            # fig, (ax0, ax1) = plt.subplots(nrows=2)

            im = plt.pcolormesh(*transform(self.xgrid, self.ygrid), self.zgrid, cmap=cmap, norm=norm)
            if legend: plt.gcf().colorbar(im, ax=plt.gca())


def plot_beamwheel(beam):
    x0=-beam.length-beam.insertion_depth
    for gamma in np.linspace(0,2*np.pi,24+1)[1:-1]:
        (Rxx, Rxy), (Ryx, Ryy) = R(gamma)
        beam.draw_beam_graphics(transform=lambda X, Y: (Rxx*(X+x0) + Rxy*Y-x0, Ryx*(X+x0)+Ryy*Y))

class GearShaft(object):
    def __init__(self, pressure_angle_deg=24):
        self.little_r = 10.6/2
        self.big_r = 13/2
        self.key_r = 12.2/2
        self.pressure_angle=pressure_angle_deg*np.pi/180.
        self.addendum_r = self.key_r*np.cos(self.pressure_angle)
        # q = lambda theta: key_r*np.tan(pressure_angle+theta)-key_r*np.tan(pressure_angle)-key_r*theta
        # involute = lambda theta: (key_r/np.cos(pressure_angle+theta)-q(theta)*np.cos(np.pi-pressure_angle-theta) - key_r/np.cos(pressure_angle), -q(theta)*np.sin(np.pi-pressure_angle-theta))
        self.involute = InvoluteProfile(self.addendum_r, self.key_r, self.pressure_angle)


    def setup_topfillet_profile(self, tip_r, plot=True):

        topfillet_OG_thetas = np.linspace(0,np.pi/2, 1000)
        self.topfillet_r = .25
        topfillet_xy = np.array([[self.big_r-self.topfillet_r+self.topfillet_r*np.sin(theta), self.topfillet_r*np.cos(theta)] for theta in topfillet_OG_thetas])
        self.topfillet_rs = np.array([np.sqrt(x**2+y**2) for x,y in topfillet_xy])
        topfillet_thetas = np.arcsin(((self.big_r-self.topfillet_r)*np.sin(topfillet_OG_thetas)+self.topfillet_r)/self.topfillet_rs)
        print(topfillet_thetas.shape)
        print(topfillet_OG_thetas.shape)
        assert(topfillet_thetas.shape==topfillet_OG_thetas.shape)


        involute_OG_thetas = np.linspace(-self.pressure_angle,.6,1000)
        self.involute_xy = np.array([self.involute(theta) for theta in involute_OG_thetas])
        mz, imz, mq, imq = -9999, -1, 99999, -1
        for i,(x,y) in enumerate(self.involute_xy):
            z = np.sqrt((x+self.key_r)**2+y**2)
            if z < self.addendum_r:
                mz=z
                imz=i
            if z > self.big_r:
                if imq==-1:
                    imq=i
                    mq = z
        print(imz, imq, involute_OG_thetas[0], involute_OG_thetas[imq])
        involute_OG_thetas = np.linspace(involute_OG_thetas[0], involute_OG_thetas[imq],1000)
        self.involute_xy = np.array([self.involute(theta) for theta in involute_OG_thetas])
        self.reverse_involute_xy = np.array(self.involute_xy)
        self.reverse_involute_xy[:,1]*=-1
        self.involute_xy+=np.hstack([np.ones((len(involute_OG_thetas),1))*(self.key_r), 0*np.ones((len(involute_OG_thetas),1))])
        self.reverse_involute_xy+=np.hstack([np.ones((len(involute_OG_thetas),1))*(self.key_r), 0*np.ones((len(involute_OG_thetas),1))])
        involute_r = np.array([np.sqrt(x**2+y**2) for x,y in self.involute_xy])
        involute_thetas = np.arccos(self.addendum_r/involute_r)

        test1 = involute_r-np.interp(involute_thetas, topfillet_thetas, self.topfillet_rs)
        involute_ndx = test1.searchsorted(0)-1
        test2 = np.interp(topfillet_thetas, involute_thetas, involute_r) -self.topfillet_rs
        topfillet_ndx = test2.searchsorted(0)

        self.combo_thetas = np.hstack([involute_thetas[:involute_ndx],topfillet_thetas[topfillet_ndx:]])
        self.combo_r = np.hstack([involute_r[:involute_ndx], self.topfillet_rs[topfillet_ndx:]])
        topfillet_correction_angle = np.arctan2(self.involute_xy[involute_ndx,1], self.involute_xy[involute_ndx,0])-np.arctan2(topfillet_xy[topfillet_ndx,1], topfillet_xy[topfillet_ndx,0])
        self.combo_xy = np.vstack([self.involute_xy[:involute_ndx,:], topfillet_xy[topfillet_ndx:,:].dot(R(topfillet_correction_angle))])
        combo_polar_angle = np.arctan2(self.combo_xy[:,1],self.combo_xy[:,0])


        # center-tip: this is the lookup we use to find what angle the gear is at to exactly touch the deflected flexure
        self.center_tip_r = np.sqrt(self.combo_r**2+tip_r**2-2*self.combo_r*tip_r*np.cos(self.combo_thetas+np.pi/2))
        self.center_tip_theta = np.arccos(self.combo_r*np.cos(self.combo_thetas)/self.center_tip_r)
        self.center_tip_xy=self.combo_xy+np.array([[tip_r*np.cos(angle), tip_r*np.sin(angle)] for angle in (np.pi/2-self.combo_thetas)+combo_polar_angle])
        self.center_tip_polar_angle = np.arctan2(self.center_tip_xy[:,1],self.center_tip_xy[:,0])


        # topfillet_ndx = 
        if plot:
            old_ax = plt.gca()
            plt.figure()
            plt.plot(involute_thetas, involute_r)
            plt.plot(topfillet_thetas, self.topfillet_rs)
            plt.plot(np.linspace(0,np.pi/2,2), np.array([self.big_r, self.big_r]))
            plt.plot(involute_thetas,test1 )
            plt.plot(involute_thetas[involute_ndx], test1[involute_ndx],"o")
            plt.plot(topfillet_thetas, test2)
            plt.plot(topfillet_thetas[topfillet_ndx], self.topfillet_rs[topfillet_ndx],"o")
            plt.plot(self.center_tip_theta, self.center_tip_r, lw=1, zorder=-3)
            plt.plot(self.combo_thetas, self.combo_r, lw=8, zorder=-3)
            plt.sca(old_ax)

        return self.center_tip_r, self.center_tip_theta, self.center_tip_polar_angle



    def draw(self, gear_x_center, gear_angle, tooth_offset=np.pi/12*.295):
        # plt.plot(gear_x_center + self.involute_xy[:,0], self.involute_xy[:,1], zorder=20, color=colors[2])
        # plt.plot(gear_x_center + self.reverse_involute_xy[:,0], self.reverse_involute_xy[:,1], zorder=20, color=colors[2])
        # plt.plot(gear_x_center + self.topfillet_xy[:,0], self.topfillet_xy[:,1], zorder=20, color=colors[2])
        # plt.plot(gear_x_center + self.combo_xy[:,0], self.combo_xy[:,1], zorder=20, color=colors[3])

        rev_center_tip_xy = np.array(self.center_tip_xy)
        rev_center_tip_xy[:,1]*=-1
        rev_combo_xy = np.array(self.combo_xy)
        rev_combo_xy[:,1]*=-1

        rev_center_tip_xy=rev_center_tip_xy.dot(R(gear_angle))
        rev_combo_xy=rev_combo_xy.dot(R(gear_angle))

        center_tip_xy = self.center_tip_xy.dot(R(gear_angle+np.pi))
        plt.plot(gear_x_center + rev_combo_xy[:,0], rev_combo_xy[:,1], zorder=20, color=colors[3])
        plt.plot(gear_x_center + rev_center_tip_xy[:,0], rev_center_tip_xy[:,1], ":", zorder=20, color=colors[3])



        # normal visualization
        if False:
            skip=100
            for (x,y), theta, r, polar_angle in list(zip(rev_center_tip_xy, self.center_tip_theta, self.center_tip_r, self.center_tip_polar_angle))[::skip]:
                psi = gear_angle-polar_angle-np.pi/2+theta
                plt.arrow(gear_x_center+x, y, .25*np.cos(psi), .25*np.sin(psi), width=.0025)

        combo_xy = self.combo_xy.dot(R(gear_angle))
        gamma_1 = np.arctan2(combo_xy[-1,1], combo_xy[-1,0])
        gamma_2 = np.arctan2(rev_combo_xy[-1,1], rev_combo_xy[-1,0])
        print (gamma_2-gamma_1)
        tooth_offset = gamma_2-gamma_1
        # involute_xy=self.involute_xy.dot(R(np.pi))
        # reverse_involute_xy=self.reverse_involute_xy.dot(R(np.pi))
        # involute_xy=involute_xy.dot(R(-.05))
        # reverse_involute_xy=reverse_involute_xy.dot(R(-.05-.160))
        # plt.plot(gear_x_center  + involute_xy[:,0], involute_xy[:,1], zorder=20, color=colors[2])
        # plt.plot(gear_x_center  + reverse_involute_xy[:,0], reverse_involute_xy[:,1], zorder=20, color=colors[2])
        for gamma in np.linspace(0,2*np.pi,24+1)[:-1]:
        #     involute_xy.dot(R(gamma))
            plt.plot(gear_x_center + rev_combo_xy.dot(R(gamma))[:,0], rev_combo_xy.dot(R(gamma))[:,1], color=colors[1])
            plt.plot(gear_x_center + self.combo_xy.dot(R(gamma+gear_angle+tooth_offset))[:,0], self.combo_xy.dot(R(gamma+gear_angle+tooth_offset))[:,1], color=colors[1])
        #     plt.plot(gear_x_center + reverse_involute_xy.dot(R(gamma))[:,0], reverse_involute_xy.dot(R(gamma))[:,1], color=colors[1])
        full_thetas = np.linspace(0, 2*np.pi, 300)
        plt.plot(gear_x_center + self.big_r*np.cos(full_thetas), self.big_r*np.sin(full_thetas), "--", lw=.75, color=colors[1])
        plt.plot(gear_x_center + self.little_r*np.cos(full_thetas), self.little_r*np.sin(full_thetas), "--", lw=.75, color=colors[1])
        plt.plot(gear_x_center + self.key_r*np.cos(full_thetas), self.key_r*np.sin(full_thetas), ":", lw=.5, color=colors[1])
        plt.plot(gear_x_center + self.addendum_r*np.cos(full_thetas), self.addendum_r*np.sin(full_thetas), ":", lw=.5, color=colors[1])

class InvoluteProfile(object):
    def __init__(self, addendum_r, key_r, pressure_angle):    
        Rad = addendum_r
        Rkey = key_r
        phi = pressure_angle
        x0 = Rkey*np.sin(phi)
        sin = np.sin
        cos = np.cos
        self.involute = lambda theta: ( -Rkey + Rad*cos(phi+theta) + (theta*Rad+x0)*sin(phi+theta),-Rad*sin(phi+theta)+ (theta*Rad+x0)*cos(phi+theta) )

    def __call__(self, theta):
        return self.involute(theta)

def main():

    material = MaragingSteelC300()
    print("Maraging Steel yield strain", MaragingSteelC300().yield_stress/MaragingSteelC300().E)
    # print("Aluminum7075 yield strain", Aluminum7075().yield_stress/Aluminum7075().E)
    # print("steel yield strain", SS410().yield_stress/SS410().E)
    # print("titanium yield strain", Titanium5().yield_stress/Titanium5().E)
    gear = GearShaft()
    beam = TimoshenkoBeamFlexure(material) #-12.32542434, 41.70018387
    # beam.insertion_depth = 6.35#
    # beam.set_geometry_mm(
    #     [ 23,  10,    1,      0], 
    #     [0.8, .45,  .25, 1.01/2], 5)
    beam.flex_by_local_tip_force(-89.14754078, 70.52246586)
    beam.local_moment

    center_tip_r, center_tip_theta, center_tip_polar_angle = gear.setup_topfillet_profile(beam.tip_r, plot=False)

    deflections = []
    tip_angles = []
    torques = []
    vm_safety_margins=[]
    hertz_safety_margins=[]
    local_force = np.array([0,1]) # unit vector in the right direction
    j = 1
    for force in np.logspace(math.log2(0.5),math.log2(500),101,base=2.0): # Are these in pounds? # original: np.linspace(0.5,240,51)
        local_force=force*local_force/np.linalg.norm(local_force)
        for i in range(5000):
            flexure_tip_r, flexure_tip_force_pressure_angle, flexure_tip_angle, flexure_tip_force, flexure_tip_rotation = beam.flex_by_local_tip_force(local_force[0], local_force[1])
            ndx = center_tip_r.searchsorted(flexure_tip_r)
            if ndx==len(center_tip_polar_angle):
                # this is the case where the spring is so deflected that it does not contact the gear.
                # beam.draw_beam_graphics()
                # local_force*=np.exp(.02*(center_tip_r[-1]-flexure_tip_r))
                ndx-=1
                # plt.axis("equal")
                # plt.show()
            # print("ndx", ndx, center_tip_theta[ndx], center_tip_polar_angle[ndx])
            # print(flexure_tip_force)
            # print(flexure_tip_force.dot(R(-flexure_tip_angle)))

            standard_force = flexure_tip_force.dot(R(-flexure_tip_angle))
            torque = 24*standard_force[1]*flexure_tip_r*1e-3   #Is there a unit conversion here? What do these numbers represent??

            # original strategy for converging
            # new_x = standard_force[1]*np.tan(-center_tip_theta[ndx])
            # next_force = np.array(standard_force)
            # next_force[0]=new_x

            # new strategy
            # unit_vector = np.array([np.cos(-center_tip_theta[ndx]), np.sin(-center_tip_theta[ndx])])
            # next_force = np.linalg.norm(standard_force)*unit_vector*.5+.5*standard_force
            # global_next_force = next_force.dot(R(flexure_tip_angle))

            # even newer strategy
            gear_angle = flexure_tip_angle+center_tip_polar_angle[ndx]
            theta = center_tip_theta[ndx]
            polar_angle = center_tip_polar_angle[ndx]
            psi = gear_angle-polar_angle-np.pi/2+theta
            global_unit_vector = np.array([np.cos(psi), np.sin(psi)])

            alpha=.125 * (10./(i+10)) 
            global_next_force = (1-alpha)*flexure_tip_force+alpha*global_unit_vector*np.linalg.norm(flexure_tip_force)
            if ndx==(len(center_tip_polar_angle)-1):
                alpha=.125
                psi = center_tip_polar_angle[ndx]
                global_unit_vector = np.array([np.cos(psi), np.sin(psi)])
                global_next_force = alpha*global_unit_vector*global_unit_vector.dot(flexure_tip_force)*(1-alpha)+(1-alpha)*flexure_tip_force
                print(np.linalg.norm(flexure_tip_force))

            new_local_force = beam.global_to_local_force(global_next_force)
            if np.linalg.norm(new_local_force-local_force)<1e-4:
                print("done iterating!", new_local_force)
                break
            local_force=new_local_force
        vm_safety_margins.append(beam.total_safety_margin)
        # just some stuff to hide data beyond beam failure
        if (beam.total_safety_margin < 1)&(j == 1):
            xlim = -(flexure_tip_angle+center_tip_polar_angle[ndx]-3.2157024285455895)*180/np.pi
            j = 2
            # define stiffness as end of quasi-linear behavior (just before yield)
            stiffness = torques[-1]/deflections[-1]
        hertz_safety_margins.append(beam.total_safety_margin)
        deflections.append(flexure_tip_angle+center_tip_polar_angle[ndx]-3.2157024285455895) # Where does this number come from??
        torques.append(torque)
        tip_angles.append(flexure_tip_rotation)
        print("gear rotation angle", flexure_tip_angle+center_tip_polar_angle[ndx]-3.2157024285455895, "gear_torque", torque)
            # plt.arrow(beam.insertion_depth + beam.length+flexure_tip_r, 0, 0.05*standard_force[0], 0.05*standard_force[1], width=.25)
            # plt.arrow(beam.insertion_depth + beam.length+flexure_tip_r, 0, 0.05*new_x, 0.05*standard_force[1], width=.25)
            # plt.arrow(beam.insertion_depth + beam.length+flexure_tip_r*np.cos(flexure_tip_angle), flexure_tip_r*np.sin(flexure_tip_angle),
            #     0.05*global_next_force[0], 0.05*global_next_force[1], width=.25)
        print("tip rotation angle: ", flexure_tip_rotation)
    beam.draw_tip_force(scale=0.05, width=.25)

    print("stiffness: (Nm/deg)",stiffness)

    beam.draw_beam_graphics()
    plot_beamwheel(beam)
    gear.draw(beam.insertion_depth + beam.length, flexure_tip_angle+center_tip_polar_angle[ndx])
    plt.axis("equal")
    fig,axs = plt.subplots(2,1,sharex=True)
    axs[0].plot(-np.array(deflections)*180/np.pi, -np.array(torques))
    axs[0].plot(np.array(tip_angles)*180/np.pi, -np.array(torques))
    axs[-1].set_xlabel("degrees")
    axs[0].set_ylabel("Nm")
    axs[1].set_ylabel("von mises yield safety margin")
    axs[1].plot(-np.array(deflections)*180/np.pi, vm_safety_margins)
    axs[1].plot(-np.array(deflections)*180/np.pi, hertz_safety_margins)
    axs[1].plot(-np.array(deflections)*180/np.pi, [1]*len(vm_safety_margins))
    axs[1].set_ylim([0,5])
    axs[1].set_xlim([0,xlim])

    plt.show()

if __name__ == '__main__':
    main()
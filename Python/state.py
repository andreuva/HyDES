import numpy as np
from calc import deriv2D, midval

class state:
    def __init__(self, params):
        # initialize some of the variables
        self.gamm = 5./3.
        self.p00 = 1.15
        self.Um00 = 0.8
        self.cs00 = np.sqrt(self.gamm*self.p00/self.Um00)

        # set the equilibrium velocity
        if str(params["type"]).lower() == "shock":
            self.v00 = params["magNum"]*self.cs00
        else:
            self.v00 = params["v0"]*self.cs00

    def initilice_conditions(self, pert, pert_vx, pert_vy, ampl=1):

        # compute the initial states as eq. + perturbation  
        self.Uminit   = self.Um00  + self.Um00          *ampl*pert
        self.presinit = self.p00   + self.gamm*self.p00 *ampl*pert
        self.vxinit   = self.v00   + self.cs00          *ampl*pert_vx
        self.vyinit   = self.v00   + self.cs00          *ampl*pert_vy

        # this calculates the initial values of the 'densities' for momentum and energy:
        self.Upxinit  = self.Uminit*self.vxinit
        self.Upyinit  = self.Uminit*self.vyinit 
        self.Ueinit = self.presinit/(self.gamm-1) +\
                            self.Uminit*(self.vxinit**2+self.vyinit**2)/2.

        # this calculates the sound speed array for the initial condition:
        self.csinit = np.sqrt(self.gamm*self.presinit/self.Uminit)

        # Set the state of the simulation to the initial conditions
        # the 'densities', in the sense of the conservation laws, 
        # are 'Um', 'Upx', 'Upy' and  'Ue'.
        self.vx   = self.vxinit
        self.vy   = self.vyinit
        self.cs   = self.csinit
        self.pres = self.presinit  
        self.Um   = self.Uminit
        self.Ue   = self.Ueinit
        self.Upx  = self.Upxinit
        self.Upy  = self.Upyinit
        self.Umn  = self.Um
        self.Uen  = self.Ue
        self.Upxn = self.Upx
        self.Upyn = self.Upy
    
    def advance_step(self, domain, dt):
        # initialaice the arrays wich will contain the derivates of the fluxes  
        self.fmxdx  = self.fmx;      self.fmydy  = self.fmy
        self.fexdx  = self.fex;      self.feydy  = self.fey
        self.fpxxdx = self.fpxx;     self.fpyydy = self.fpyy
        self.fpyxdx = self.fpyx;     self.fpxydy = self.fpxy
  
        # calculate the derivates of the mas, energy and momentum fluxes in x and y
        self.fmxdx = deriv2D(self.fmx, dx=domain.dx);     self.fmydy = deriv2D(self.fmy, dy=domain.dy)
        self.fexdx = deriv2D(self.fex, dx=domain.dx);     self.feydy = deriv2D(self.fey, dy=domain.dy)
        self.fpxxdx = deriv2D(self.fpxx, dx=domain.dx);   self.fpyydy = deriv2D(self.fpyy, dy=domain.dy)
        self.fpyxdx = deriv2D(self.fpyx, dx=domain.dx);   self.fpxydy = deriv2D(self.fpxy, dy=domain.dy)
  
        # compute the densities in the time t + dt/2
        self.Umn[:-1,:-1]  = midval(self.Um)  - (dt/2)*(self.fmxdx  + self.fmydy)
        self.Uen[:-1,:-1]  = midval(self.Ue)  - (dt/2)*(self.fexdx  + self.feydy)
        self.Upxn[:-1,:-1] = midval(self.Upx) - (dt/2)*(self.fpxxdx + self.fpxydy)
        self.Upyn[:-1,:-1] = midval(self.Upy) - (dt/2)*(self.fpyydy + self.fpyxdx)
  
        # calculate the fluxes in t + dt/2
        new_fluxes = self.compute_fluxes(self.Umn[:-1,:-1],self.Upxn[:-1,:-1],self.Upyn[:-1,:-1],self.Uen[:-1,:-1])
        self.fmxn,self.fmyn,   self.fpxxn,self.fpxyn,self.fpyyn,self.fpyxn,  self.fexn,self.feyn = new_fluxes
  
        # calculate the derivates of the mas, energy and momentum in x and y in t+dt/2
        self.fmxdx = deriv2D(self.fmxn, dx=domain.dx);     self.fmydy = deriv2D(self.fmyn, dy=domain.dy)
        self.fexdx = deriv2D(self.fexn, dx=domain.dx);     self.feydy = deriv2D(self.feyn, dy=domain.dy)
        self.fpxxdx = deriv2D(self.fpxxn, dx=domain.dx);   self.fpyydy = deriv2D(self.fpyyn, dy=domain.dy)
        self.fpyxdx = deriv2D(self.fpyxn, dx=domain.dx);   self.fpxydy = deriv2D(self.fpxyn, dy=domain.dy)
  
        # compute the densities in the time t + dt 
        self.Umn[1:-1,1:-1] =  self.Um[1:-1,1:-1] - dt*(self.fmxdx  + self.fmydy)
        self.Uen[1:-1,1:-1] =  self.Ue[1:-1,1:-1] - dt*(self.fexdx  + self.feydy)
        self.Upxn[1:-1,1:-1] = self.Upx[1:-1,1:-1] - dt*(self.fpxxdx + self.fpxydy)
        self.Upyn[1:-1,1:-1] = self.Upy[1:-1,1:-1] - dt*(self.fpyydy + self.fpyxdx)


    def compute_fluxes(self, Um, Upx, Upy, Ue):
        vx, vy, pres = self.compute_primitive_quantities(Um, Upx, Upy, Ue)

        fmx = Upx
        fmy = Upy

        fpxx = Upx*Upx/Um + pres
        fpxy = Upx*Upy/Um 
        fpyy = Upy*Upy/Um + pres
        fpyx = Upy*Upx/Um 
        
        fex = (Ue + pres) * vx
        fey = (Ue + pres) * vy

        return fmx, fmy, fpxx, fpxy, fpyy, fpyx, fex, fey

    def compute_time_step(self, clf, dx, dy):
        # compute the characteristic velocities
        vcharac1 = np.abs(np.sqrt(self.Upx*self.Upx + self.Upy*self.Upy)/self.Um + self.cs)
        vcharac2 = np.abs(np.sqrt(self.Upx*self.Upx + self.Upy*self.Upy)/self.Um - self.cs)

        # compute the proper delta t with the cflparam
        dt = clf * np.min([dx,dy]) / np.max([vcharac1, vcharac2])

        return dt


    def compute_primitive_quantities(self, Um, Upx, Upy, Ue):
        vx   = Upx/Um
        vy   = Upy/Um
        pres = (self.gamm-1.)*(Ue - (Upx*Upx + Upy*Upy)/(2*Um))

        return vx, vy, pres

    def apply_boundary_conditions(self, type, border):
        # print("Applying boundary conditions")
        # print("Type: ", type)
        # print("Border: ", border)

        if border == 'top':
            self.Umn[0,:] = self.Umn[-2,:]
            self.Uen[0,:] = self.Uen[-2,:]
            self.Upxn[0,:] = self.Upxn[-2,:]
            self.Upyn[0,:] = self.Upyn[-2,:]
        elif border == 'bottom':
            self.Umn[-1,:] = self.Umn[1,:]
            self.Uen[-1,:] = self.Uen[1,:]
            self.Upxn[-1,:] = self.Upxn[1,:]
            self.Upyn[-1,:] = self.Upyn[1,:]
        elif border == 'left':
            self.Umn[:,0] = self.Umn[:,-2]
            self.Uen[:,0] = self.Uen[:,-2]
            self.Upxn[:,0] = self.Upxn[:,-2]
            self.Upyn[:,0] = self.Upyn[:,-2]
        elif border == 'right':
            self.Umn[:,-1] = self.Umn[:,1]
            self.Uen[:,-1] = self.Uen[:,1]
            self.Upxn[:,-1] = self.Upxn[:,1]
            self.Upyn[:,-1] = self.Upyn[:,1]
        else:
            raise ValueError('Border must be "left", "right", "top" or "bottom"')

    def update(self):
        self.Um = self.Umn
        self.Upx = self.Upxn
        self.Upy = self.Upyn
        self.Ue = self.Uen
        self.pres = self.presn
        self.vx = self.vxn
        self.vy = self.vyn 

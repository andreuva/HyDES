import numpy as np
from .calc import deriv2D, midval

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
        self.vx   = self.vxinit.copy()
        self.vy   = self.vyinit.copy()
        self.cs   = self.csinit.copy()
        self.pres = self.presinit.copy()
        self.Um   = self.Uminit.copy()
        self.Ue   = self.Ueinit.copy()
        self.Upx  = self.Upxinit.copy()
        self.Upy  = self.Upyinit.copy()
        self.Umn  = self.Um.copy()
        self.Uen  = self.Ue.copy()
        self.Upxn = self.Upx.copy()
        self.Upyn = self.Upy.copy()
    
    def advance_step(self, domain, dt):
        # initialaice the arrays wich will contain the derivates of the fluxes  
        self.fmxdx  = self.fmx.copy();      self.fmydy  = self.fmy.copy()
        self.fexdx  = self.fex.copy();      self.feydy  = self.fey.copy()
        self.fpxxdx = self.fpxx.copy();     self.fpyydy = self.fpyy.copy()
        self.fpyxdx = self.fpyx.copy();     self.fpxydy = self.fpxy.copy()
  
        # calculate the derivates of the mas, energy and momentum fluxes in x and y
        self.fmxdx = deriv2D(self.fmx, dx=domain.dx, axis=0);     self.fmydy = deriv2D(self.fmy, dy=domain.dy, axis=1)
        self.fexdx = deriv2D(self.fex, dx=domain.dx, axis=0);     self.feydy = deriv2D(self.fey, dy=domain.dy, axis=1)
        self.fpxxdx = deriv2D(self.fpxx, dx=domain.dx, axis=0);   self.fpyydy = deriv2D(self.fpyy, dy=domain.dy, axis=1)
        self.fpyxdx = deriv2D(self.fpyx, dx=domain.dx, axis=0);   self.fpxydy = deriv2D(self.fpxy, dy=domain.dy, axis=1)
  
        # compute the densities in the time t + dt/2
        self.Umn[:-1,:-1]  = midval(self.Um)  - (dt/2)*(self.fmxdx  + self.fmydy)
        self.Uen[:-1,:-1]  = midval(self.Ue)  - (dt/2)*(self.fexdx  + self.feydy)
        self.Upxn[:-1,:-1] = midval(self.Upx) - (dt/2)*(self.fpxxdx + self.fpxydy)
        self.Upyn[:-1,:-1] = midval(self.Upy) - (dt/2)*(self.fpyydy + self.fpyxdx)
  
        # calculate the fluxes in t + dt/2
        new_fluxes = self.compute_fluxes(self.Umn[:-1,:-1],self.Upxn[:-1,:-1],self.Upyn[:-1,:-1],self.Uen[:-1,:-1])
        self.fmxn,self.fmyn,   self.fpxxn,self.fpxyn,self.fpyyn,self.fpyxn,  self.fexn,self.feyn = new_fluxes
  
        # calculate the derivates of the mas, energy and momentum in x and y in t+dt/2
        self.fmxdx = deriv2D(self.fmxn, dx=domain.dx, axis=0);     self.fmydy = deriv2D(self.fmyn, dy=domain.dy, axis=1)
        self.fexdx = deriv2D(self.fexn, dx=domain.dx, axis=0);     self.feydy = deriv2D(self.feyn, dy=domain.dy, axis=1)
        self.fpxxdx = deriv2D(self.fpxxn, dx=domain.dx, axis=0);   self.fpyydy = deriv2D(self.fpyyn, dy=domain.dy, axis=1)
        self.fpyxdx = deriv2D(self.fpyxn, dx=domain.dx, axis=0);   self.fpxydy = deriv2D(self.fpxyn, dy=domain.dy, axis=1)
  
        # compute the densities in the time t + dt 
        self.Umn[1:-1,1:-1] =  self.Um[1:-1,1:-1] - dt*(self.fmxdx  + self.fmydy)
        self.Uen[1:-1,1:-1] =  self.Ue[1:-1,1:-1] - dt*(self.fexdx  + self.feydy)
        self.Upxn[1:-1,1:-1] = self.Upx[1:-1,1:-1] - dt*(self.fpxxdx + self.fpxydy)
        self.Upyn[1:-1,1:-1] = self.Upy[1:-1,1:-1] - dt*(self.fpyydy + self.fpyxdx)

        print()


    def compute_fluxes(self, Um, Upx, Upy, Ue):
        vx, vy, pres = self.compute_primitive_quantities(Um, Upx, Upy, Ue)

        fmx = Upx.copy()
        fmy = Upy.copy()

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
        if type == "periodic":
            if border == 'top':
                self.Umn[0,:] = self.Umn[-2,:].copy()
                self.Uen[0,:] = self.Uen[-2,:].copy()
                self.Upxn[0,:] = self.Upxn[-2,:].copy()
                self.Upyn[0,:] = self.Upyn[-2,:].copy()
            elif border == 'bottom':
                self.Umn[-1,:] = self.Umn[1,:].copy()
                self.Uen[-1,:] = self.Uen[1,:].copy()
                self.Upxn[-1,:] = self.Upxn[1,:].copy()
                self.Upyn[-1,:] = self.Upyn[1,:].copy()
            elif border == 'left':
                self.Umn[:,0] = self.Umn[:,-2].copy()
                self.Uen[:,0] = self.Uen[:,-2].copy()
                self.Upxn[:,0] = self.Upxn[:,-2].copy()
                self.Upyn[:,0] = self.Upyn[:,-2].copy()
            elif border == 'right':
                self.Umn[:,-1] = self.Umn[:,1].copy()
                self.Uen[:,-1] = self.Uen[:,1].copy()
                self.Upxn[:,-1] = self.Upxn[:,1].copy()
                self.Upyn[:,-1] = self.Upyn[:,1].copy()
            else:
                raise ValueError('Border must be "left", "right", "top" or "bottom"')
        elif type == "0deriv":
            if border == 'top':
                self.Umn[0,:] = self.Umn[1,:].copy()
                self.Uen[0,:] = self.Uen[1,:].copy()
                self.Upxn[0,:] = self.Upxn[1,:].copy()
                self.Upyn[0,:] = self.Upyn[1,:].copy()
            elif border == 'bottom':
                self.Umn[-1,:] = self.Umn[-2,:].copy()
                self.Uen[-1,:] = self.Uen[-2,:].copy()
                self.Upxn[-1,:] = self.Upxn[-2,:].copy()
                self.Upyn[-1,:] = self.Upyn[-2,:].copy()
            elif border == 'left':
                self.Umn[:,0] = self.Umn[:,1].copy()
                self.Uen[:,0] = self.Uen[:,1].copy()
                self.Upxn[:,0] = self.Upxn[:,1].copy()
                self.Upyn[:,0] = self.Upyn[:,1].copy()
            elif border == 'right':
                self.Umn[:,-1] = self.Umn[:,-2].copy()
                self.Uen[:,-1] = self.Uen[:,-2].copy()
                self.Upxn[:,-1] = self.Upxn[:,-2].copy()
                self.Upyn[:,-1] = self.Upyn[:,-2].copy()
            else:
                raise ValueError('Border must be "left", "right", "top" or "bottom"')
        else:
            raise ValueError('Type must be "periodic" or "0deriv"')

    def update(self):
        self.Um = self.Umn.copy()
        self.Upx = self.Upxn.copy()
        self.Upy = self.Upyn.copy()
        self.Ue = self.Uen.copy()
        self.pres = self.presn.copy()
        self.vx = self.vxn.copy()
        self.vy = self.vyn.copy()

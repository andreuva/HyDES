import numpy as np

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
        # are 'Um', 'Upx', 'Upz' and  'Ue'.
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
    
    def update(self):
        pass
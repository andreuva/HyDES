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

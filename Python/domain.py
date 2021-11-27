import numpy as np

# Class that represents the domain of the simulation
class domain:
    def __init__(self, params):

        # Store the domain parameters
        self.xmin = params["xmin"]
        self.xmax = params["xmax"]
        self.ymin = params["ymin"]
        self.ymax = params["ymax"]
        self.xres = params["xres"]
        self.yres = params["yres"]

        # Initialize the domain
        self.x = np.linspace(self.xmin, self.xmax, self.xres)
        self.y = np.linspace(self.ymin, self.ymax, self.yres)

        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]

        self.xmesh, self.ymesh = np.meshgrid(self.x, self.y)

        # self.xvec = np.reshape(self.xmesh, (self.xres * self.yres, 1))
        # self.yvec = np.reshape(self.ymesh, (self.xres * self.yres, 1))

        self.N = self.xres * self.yres

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

        self.xnum = self.xres + 2
        self.ynum = self.yres + 2

        self.dx = (self.xmax - self.xmin) / self.xres
        self.dy = (self.ymax - self.ymin) / self.yres

        # Initialize the domain
        self.x = np.linspace(self.xmin - self.dx/2, self.xmax + self.dx/2, self.xnum)
        self.y = np.linspace(self.ymin - self.dy/2, self.ymax + self.dy/2, self.ynum)

        self.xmesh, self.ymesh = np.meshgrid(self.x, self.y)
        self.N = self.xnum * self.ynum

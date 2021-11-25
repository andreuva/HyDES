import numpy as np

# Class that represents the domain of the simulation
class domain:
    def __init__(self, xmin, xmax, ymin, ymax, xres, yres):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.xres = xres
        self.yres = yres

        self.x = np.linspace(self.xmin, self.xmax, self.xres)
        self.y = np.linspace(self.ymin, self.ymax, self.yres)

        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]

        self.xmesh, self.ymesh = np.meshgrid(self.x, self.y)

        # self.xvec = np.reshape(self.xmesh, (self.xres * self.yres, 1))
        # self.yvec = np.reshape(self.ymesh, (self.xres * self.yres, 1))

        self.N = self.xres * self.yres

import numpy as np
import matplotlib.pyplot as plt
import os

# Class that represents the domain of the simulation
class display:
    def __init__(self, domain, state, params):

        # store the type of plot to be used
        self.plot_type = params["plotType"]
        self.savplotcad = params["savePlotCad"]

        # store the path to save the plots (if needed)
        if params["savePath"]:
            self.savepath = params["savePathPlots"]
        else:
            self.savepath = None

        if params["plotType"] == "quiver":
            self.nrows = 1
            self.ncols = 1
            self.skip = (slice(None,None,10),slice(None,None,10))
        else:
            self.nrows = 2
            self.ncols = 2


        self.simfig, self.sim_ax_list = plt.subplots(self.nrows, self.ncols, figsize=(6,6))
        plt.rcParams['image.cmap'] = 'plasma'
        self.simfig.suptitle('2D Hydrodinamic Simulation')

        self.umamp = np.max(np.abs((state.Uminit  - state.Um00)/state.Um00))
        self.pamp = np.max(np.abs((state.presinit - state.p00)/state.p00))
        
        vvinit = np.sqrt(state.vxinit*state.vxinit + state.vyinit*state.vyinit)
        v_ampinit = vvinit/state.cs00
        med   = (np.max(v_ampinit) + np.min(v_ampinit))/2
        range = np.max(v_ampinit) - np.min(v_ampinit)
        self.vmin = med - range
        self.vmax = med + range

        self.numplots = 0
        self.update(domain, state)

    def update(self, domain, state, itt=0, saveplot=True, makeplot=False):
        # Plot the state
        if self.plot_type == "contourf":
            ax_list = self.sim_ax_list.ravel()
            ax_list[0].set_title('Density')
            ax_list[0].contourf(domain.x,  domain.y, (state.Um-state.Um00)/state.Um00, 20, vmin=-self.umamp , vmax=self.umamp)
            ax_list[1].set_title('Pressure')
            ax_list[1].contourf(domain.x,  domain.y, (state.pres-state.p00)/state.p00, 20, vmin=-self.pamp , vmax=self.pamp)
            ax_list[2].set_title('X velocity')
            ax_list[2].contourf(domain.x,  domain.y, state.vx/state.cs00, 20, vmin=self.vmin, vmax=self.vmax)
            ax_list[3].set_title('Y velocity')
            ax_list[3].contourf(domain.x,  domain.y, state.vy/state.cs00, 20, vmin=self.vmin, vmax=self.vmax)

        elif self.plot_type == "color":
            ax_list = self.sim_ax_list.ravel()
            ax_list[0].set_title('Density')
            ax_list[0].imshow((state.Um-state.Um00)/state.Um00,
                            extent=[domain.xmin,domain.xmax, domain.ymin,domain.ymax],
                            vmin=-self.umamp , vmax=self.umamp)
            ax_list[1].set_title('Pressure')
            ax_list[1].imshow((state.pres-state.p00)/state.p00,
                            extent=[domain.xmin,domain.xmax, domain.ymin,domain.ymax],
                            vmin=-self.pamp , vmax=self.pamp)
            ax_list[2].set_title('X velocity')
            ax_list[2].imshow(state.vx/state.cs00,
                            extent=[domain.xmin,domain.xmax, domain.ymin,domain.ymax],
                            vmin=self.vmin, vmax=self.vmax)
            ax_list[3].set_title('Y velocity')
            ax_list[3].imshow(state.vy/state.cs00,
                          extent=[domain.xmin,domain.xmax, domain.ymin,domain.ymax],
                          vmin=self.vmin, vmax=self.vmax)
        elif self.plot_type == "quiver":
            self.sim_ax_list.clear()
            self.sim_ax_list.set_title('Density with velocity')
            self.sim_ax_list.quiver(domain.xmesh[self.skip], domain.ymesh[self.skip],
                                    (state.vx/state.cs00)[self.skip], (state.vy/state.cs00)[self.skip],
                                    ((state.pres-state.p00)/state.p00)[self.skip], angles='xy', scale_units='xy', scale=1.)
            self.sim_ax_list.set_xlabel('x')
            self.sim_ax_list.set_ylabel('y')
            self.sim_ax_list.set_aspect('equal')
        else:
            raise ValueError("Plot type not recognized")

        if self.plot_type != "quiver":
            for ax in ax_list:
                ax.label_outer()
        
        if makeplot:
            plt.draw()
            plt.show(block=False)
            plt.pause(0.01)

        if saveplot and self.savepath is not None:
            self.simfig.savefig(os.path.join(self.savepath, f'plot_{1000000 + itt}.png'))

        self.numplots += 1

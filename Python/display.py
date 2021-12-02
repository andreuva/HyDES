import numpy as np
import matplotlib.pyplot as plt
import os
import time

# Class that represents the domain of the simulation
class display:
    def __init__(self, domain, state, params):

        # store the type of plot to be used
        self.plot_type = params["plotType"]
        self.savplotcad = params["savePlotCad"]

        # store the path to save the plots
        if not os.path.exists(params["savePath"]):
            os.mkdir(params["savePath"])
        self.savepath = os.path.join(params["savePath"], f'plots_sim_{params["simName"]}_{time.strftime("%Y%m%d-%H%M")}')
        if not os.path.exists(self.savepath):
            os.mkdir(self.savepath)

        self.simfig, self.sim_ax_list = plt.subplots(2, 2, figsize=(6,6))
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

    def update(self, domain, state, time=0):
        # Plot the state
        ax_list = self.sim_ax_list.ravel()
        
        if self.plot_type == "contourf":
            ax_list[0].set_title('Density')
            ax_list[0].contourf(domain.x,  domain.y, (state.Um-state.Um00)/state.Um00, 20, vmin=-self.umamp , vmax=self.umamp)
            ax_list[1].set_title('Pressure')
            ax_list[1].contourf(domain.x,  domain.y, (state.pres-state.p00)/state.p00, 20, vmin=-self.pamp , vmax=self.pamp)
            ax_list[2].set_title('X velocity')
            ax_list[2].contourf(domain.x,  domain.y, state.vx/state.cs00, 20, vmin=self.vmin, vmax=self.vmax)
            ax_list[3].set_title('Y velocity')
            ax_list[3].contourf(domain.x,  domain.y, state.vy/state.cs00, 20, vmin=self.vmin, vmax=self.vmax)

        elif self.plot_type == "color":
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
        else:
            raise ValueError("Plot type not recognized")

        for ax in ax_list:
            ax.label_outer()
        
        plt.draw()
        plt.pause(0.00001)

        if self.numplots % self.savplotcad == 0:
            self.simfig.savefig(os.path.join(self.savepath, f'plot_{time}.png'))

        self.numplots += 1
        

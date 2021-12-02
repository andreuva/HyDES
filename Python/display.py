import numpy as np
import matplotlib.pyplot as plt

# Class that represents the domain of the simulation
class display:
    def __init__(self, domain, state, params):

        # store the type of plot to be used
        self.plot_type = params["plotType"]
        self.simfig, self.sim_ax_list = plt.subplots(2, 2)
        self.simfig.suptitle('State of the simulation')
        self.update(domain, state)

    def update(self, domain, state, time=0):
        # Plot the state
        ax_list = self.sim_ax_list.ravel()
        
        if self.plot_type == "contourf":
            ax_list[0].set_title('Density')
            ax_list[0].contourf(domain.x,  domain.y, (state.Um-state.Um00)/state.Um00, 20)
            ax_list[1].set_title('Pressure')
            ax_list[1].contourf(domain.x,  domain.y, (state.pres-state.p00)/state.p00, 20)
            ax_list[2].set_title('X velocity')
            ax_list[2].contourf(domain.x,  domain.y, state.vx/state.cs00, 20)
            ax_list[3].set_title('Y velocity')
            ax_list[3].contourf(domain.x,  domain.y, state.vy/state.cs00, 20)
        elif self.plot_type == "color":
            ax_list[0].set_title('Density')
            ax_list[0].imshow((state.Um-state.Um00)/state.Um00,
                            extent=[domain.xmin,domain.xmax, domain.ymin,domain.ymax])
            ax_list[1].set_title('Pressure')
            ax_list[1].imshow((state.pres-state.p00)/state.p00,
                            extent=[domain.xmin,domain.xmax, domain.ymin,domain.ymax])
            ax_list[2].set_title('X velocity')
            ax_list[2].imshow(state.vx/state.cs00,
                            extent=[domain.xmin,domain.xmax, domain.ymin,domain.ymax])
            ax_list[3].set_title('Y velocity')
            ax_list[3].imshow(state.vy/state.cs00,
                          extent=[domain.xmin,domain.xmax, domain.ymin,domain.ymax])
        else:
            raise ValueError("Plot type not recognized")

        for ax in ax_list:
            ax.label_outer()
        
        plt.draw()
        plt.pause(0.00001)
        

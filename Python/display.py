import numpy as np
import matplotlib.pyplot as plt

# Class that represents the domain of the simulation
class display:
    def __init__(self, domain, state, params):
        self.simfig, self.sim_ax_list = plt.subplots(2, 2)
        self.simfig.suptitle('State of the simulation')
        self.update(domain, state)

    def update(self, domain, state):
        # Plot the state
        ax_list = self.sim_ax_list.ravel()
        
        ax_list[0].set_title('Density')
        ax_list[0].contourf(domain.x,  domain.y, (state.Um-state.Um00)/state.Um00, 20)
        ax_list[1].set_title('Pressure')
        ax_list[1].contourf(domain.x,  domain.y, (state.pres-state.p00)/state.p00, 20)
        ax_list[2].set_title('X velocity')
        ax_list[2].contourf(domain.x,  domain.y, state.vx/state.cs00, 20)
        ax_list[3].set_title('Y velocity')
        ax_list[3].contourf(domain.x,  domain.y, state.vy/state.cs00, 20)

        for ax in ax_list:
            ax.label_outer()
        
        plt.draw()
        plt.pause(0.00001)
        

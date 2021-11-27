import numpy as np
import matplotlib.pyplot as plt

# Class that represents the domain of the simulation
class display:
    def __init__(self, domain, state, params):
        self.plot_state(domain, state)

    def plot_state(self, domain, state):
        # Plot the state
        plt.clf()
        plt.contourf(domain.x, domain.y, (state.presinit - state.p00)/state.p00, 100)
        plt.colorbar()
        plt.title("Pressure")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.draw()
        plt.pause(0.0001)
        

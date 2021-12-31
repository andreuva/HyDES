import json
# from importlib import resources
import os #, io
from domain import domain as dmn
from display import display as dsp
from initCond import check_inputs, compute_initial_conditions


def load_sample_params():

    # load the parameters from the input file
    with open(os.path.join(os.path.dirname(__file__),'params.json')) as file:
        params = json.load(file)

    # with resources.open_binary('parameter_files', 'params.json') as file:
    #     param_file = file.read()
    # params = json.load(io.BytesIO(param_file))

    return params



def run_sim(params=load_sample_params()):

    check_inputs(params)

    # create and compute the domain and grid
    domain = dmn(params["domain"])

    # create the state variables with the initial conditions
    state = compute_initial_conditions(domain, params["initCond"])

    # Initialize the display
    display = dsp(domain, state, params)

    itteration = 0
    time = 0.0
    while itteration < params["maxIter"] and time < params["maxTime"]:

        # compute the fluxes
        fluxes = state.compute_fluxes(state.Um, state.Upx, state.Upy, state.Ue)
        state.fmx, state.fmy, state.fpxx, state.fpxy, state.fpyy, state.fpyx, state.fex, state.fey = fluxes

        # compute the time step
        dt = state.compute_time_step(params["clf"], domain.dx, domain.dy)

        if time+dt > params["maxTime"]:
            dt = params["maxTime"]-time

        # update the state quantities
        state.advance_step(domain, dt)

        # aply the boundary conditions
        state.apply_boundary_conditions(params["boundCond"], "left")
        state.apply_boundary_conditions(params["boundCond"], "right")
        state.apply_boundary_conditions(params["boundCond"], "top")
        state.apply_boundary_conditions(params["boundCond"], "bottom")

        # calculate the primitive variables from the conserved ones
        state.vxn, state.vyn, state.presn = state.compute_primitive_quantities(state.Um, state.Upx, state.Upy, state.Ue)

        # exchange the old variables and the new ones
        state.update()

        # update the display if necessary
        if itteration % params["pltCad"] == 0:
            display.update(domain, state, time)

        # update the time
        time += dt
        itteration += 1

        print("----------------------------------")
        print("Time: ", time)
        print("Iteration: ", itteration)
        print("----------------------------------")

        # saving the state of the simulation if necessary
        if str(params["savePath"]).lower() not in ['none','no','n','false','0']:
            if itteration % params["saveSnapCad"] == 0:
                state.save_state(params["savePath"], itteration)

if __name__ == "__main__":
    run_sim()

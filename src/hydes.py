import json
import os
import numpy as np
from copy import deepcopy
from domain import domain as dmn
from display import display as dsp
from initCond import check_inputs
import state as stt

def load_sample_params():

    # load the parameters from the input file
    with open(os.path.join(os.path.dirname(__file__),'params.json')) as file:
        params = json.load(file)

    # with resources.open_binary('parameter_files', 'params.json') as file:
    #     param_file = file.read()
    # params = json.load(io.BytesIO(param_file))

    return params


def sample_perturbation(params=load_sample_params()):
    
    domain = dmn(params["domain"])
    # compute the sigma and central point in relation with the fisical domain choosed
    # sigma [sigmax, sigmay]
    sigmax, sigmay = np.mean([[domain.xmax-domain.xmin], [domain.ymax-domain.ymin]], axis=1)/8
    xp, yp = (domain.xmax + domain.xmin)/2, (domain.ymax + domain.xmin)/2

    # compute the perturbation
    # pert = np.exp(-((domain.xmesh-xp)**2/(2*sigmax**2) + (domain.ymesh-yp)**2/(2*sigmay**2)))
    pert = np.exp(-(((domain.xmesh-xp)/sigmax)**2 + ((domain.ymesh-yp)/sigmay)**2))
    pert_vx = pert*2*(domain.xmesh-xp)/sigmax
    pert_vy = pert*2*(domain.ymesh-yp)/sigmay

    return pert, pert_vx, pert_vy


def run_sim(params=load_sample_params(),
            pert=sample_perturbation()[0],
            pert_vx=sample_perturbation()[1],
            pert_vy=sample_perturbation()[2],
            restore_path=False):

    check_inputs(params)

    # create and compute the domain and grid
    domain = dmn(params["domain"])

    # create the state variables with the initial conditions
    # initialize the state class
    init_state = stt.state(params["initCond"])

    # add the perturbation to the initial conditions
    init_state.initilice_conditions(pert, pert_vx, pert_vy, params["initCond"]["ampl"])

    if restore_path:
        # restore the state from the previous run
        init_state.load_snapshot(restore_path)

    state = deepcopy(init_state)
    # state = compute_initial_conditions(domain, params["initCond"])

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


        # update the time
        time += dt
        itteration += 1

        print("----------------------------------")
        print("Time: ", time)
        print("Iteration: ", itteration)
        print("----------------------------------")

        # update the display if necessary
        if itteration % params["pltCad"] == 0:
            makeplot=True
        else:
            makeplot=False

        if itteration % params["savePlotCad"] == 0:
            saveplot=True
        else:
            saveplot=False

        if saveplot or makeplot:
            display.update(domain, state, itteration, saveplot, makeplot)

        # saving the state of the simulation if necessary
        if str(params["savePath"]).lower() not in ['none','no','n','false','0']:
            if itteration % params["saveSnapCad"] == 0:
                state.save_state(params["savePath"], itteration)

    return params["savePath"]

if __name__ == "__main__":
    run_sim()

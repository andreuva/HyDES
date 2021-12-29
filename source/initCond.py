import numpy as np
import warnings
from state import state

def check_inputs(params):
    if params["clf"] >= 0.99:
        warnings.warn(f'WARNING: Cfl param is {params["clf"]} stability may be compromised')

    # warn if the rain is activated outside a gaussian mode
    if str(params["initCond"]["rain"]["active"]).lower() in ['true','y','t','yes', '1'] and \
    str(params["initCond"]["type"]).lower() != 'gaussian':
        params["initCond"]["rain"]["active"] = True
        warnings.warn('WARNING: RAIN IS ACTIVATED IN A NON GAUSSIAN MODE')

    # Force 0 derivative at the boundaries if we are simulating a packet
    if str(params["initCond"]["type"]).lower() == 'packet':
        params["boundCond"] = "0deriv"

    # Check that the domain is valid
    if params["domain"]["xmin"] >= params["domain"]["xmax"]:
        raise ValueError('xmin must be smaller than xmax\n'+\
                        f'xmin = {params["domain"]["xmin"]}\t xmax = {params["domain"]["xmax"]}')
    if params["domain"]["ymin"] >= params["domain"]["ymax"]:
        raise ValueError('ymin must be smaller than ymax'+\
                        f'ymin = {params["domain"]["ymin"]}\t ymax = {params["domain"]["ymax"]}')
    if params["domain"]["xres"] <= 2:
        raise ValueError(f'xres must be greater than 2 and is {params["domain"]["xres"]}')
    if params["domain"]["yres"] <= 2:
        raise ValueError(f'yres must be greater than 2 and is {params["domain"]["yres"]}')

    # make sure that the wavenumbers are lees than the lenght
    if params["initCond"]["ondx"] > params["domain"]["xres"] or\
       params["initCond"]["ondy"] > params["domain"]["yres"]:
        raise ValueError('The wavenumber of the initial condition is greater than the number of intervals')

    # check the boundary conditions and puts the 'periodic' as default
    if str(params["boundCond"]).lower() not in ['periodic', '0deriv']:
        raise ValueError('boundcond must be "periodic" or "0deriv"')

    # warn if no wave is selected in the wave modes
    if params["initCond"]["ondx"] == 0 and params["initCond"]["ondy"] == 0 and \
       str(params["initCond"]["type"]).lower() in ['sound wave', 'chladni', 'packet']:
        warnings.warn(f'WARNING: There are no modes selected kx=0,ky=0 : ky=1 has been set')
        params["initCond"]["ondy"] = 1


# Calculate the arrays of density, velocity and pressure at time t=0.
def compute_initial_conditions(domain, params):
    # initialize the state class
    init_state = state(params)

    # For now, I will use only the gaussian mode

    # compute the sigma and central point in relation with the fisical domain choosed
    # sigma [sigmax, sigmay]
    sigmax, sigmay = np.mean([[domain.xmax-domain.xmin], [domain.ymax-domain.ymin]], axis=1)/8
    xp, yp = (domain.xmax + domain.xmin)/2, (domain.ymax + domain.xmin)/2

    # compute the perturbation
    # pert = np.exp(-((domain.xmesh-xp)**2/(2*sigmax**2) + (domain.ymesh-yp)**2/(2*sigmay**2)))
    pert = np.exp(-(((domain.xmesh-xp)/sigmax)**2 + ((domain.ymesh-yp)/sigmay)**2))
    pert_vx = pert*2*(domain.xmesh-xp)/sigmax
    pert_vy = pert*2*(domain.ymesh-yp)/sigmay

    # add the perturbation to the initial conditions
    init_state.initilice_conditions(pert, pert_vx, pert_vy, params["ampl"])

    # if v00 is not set to 0 warn that the gaussian will move te center over time
    if init_state.v00 != 0:
        warnings.warn('The gaussian will be moving in the direction of k because of the inital velociti v0 in the inputs')

    return init_state

import warnings, os, time

def check_inputs(params):
    if params["clf"] >= 0.99:
        warnings.warn(f'WARNING: Cfl param is {params["clf"]} stability may be compromised')

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

    # check the boundary conditions and puts the 'periodic' as default
    if str(params["boundCond"]).lower() not in ['periodic', '0deriv']:
        raise ValueError('boundcond must be "periodic" or "0deriv"')

    # check the path to save the outputs
    # store the path to save the plots (if needed)
    if str(params["savePath"]).lower() in ['none','no','n','false','0']:
        params["savePath"] = None
    else:
        # complete the save path and print it to the screen
        params["savePath"] = os.path.join(os.getcwd(), params["savePath"])

        if not os.path.exists(params["savePath"]):
            os.mkdir(params["savePath"])
        # add simulation name to the save path
        params["savePath"] = os.path.join(params["savePath"], f'sim_{params["simName"]}_{time.strftime("%Y%m%d-%H%M")}')

        if not os.path.exists(params["savePath"]):
            os.mkdir(params["savePath"])
        params["savePathPlots"] = os.path.join(params["savePath"], 'plots')
        if not os.path.exists(params["savePathPlots"]):
            os.mkdir(params["savePathPlots"])


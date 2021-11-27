import json
from domain import domain as dmn
from display import display as dsp
from initCond import check_inputs, compute_initial_conditions


# read the input file
with open('parameter_files/params.json') as file:
    params = json.load(file)

check_inputs(params)

# create and compute the domain and grid
domain = dmn(params["domain"])

# create the state variables with the initial conditions
state = compute_initial_conditions(domain, params["initCond"])

# Initialize the display
display = dsp(domain, state, params)

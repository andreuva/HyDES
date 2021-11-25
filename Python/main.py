import json
from domain import domain as dmn

# ------------------
# READ IN PARAMETERS
# ------------------
with open('parameter_files/params.json') as file:
    params = json.load(file)

# --------------------------------------
# CREATE THE DOMAIN AND COMPUTE THE GRID
# --------------------------------------
domain = dmn(params["Domain"]["x0"], params["Domain"]["xf"],
             params["Domain"]["y0"], params["Domain"]["yf"],
             params["Domain"]["Nx"], params["Domain"]["Ny"])

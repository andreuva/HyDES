import hydes as hd
from animation import animation
import os
import numpy as np


params = hd.load_sample_params()
pert, pert_vx, pert_vy = hd.sample_perturbation()

pert_vx = np.ones(pert_vx.shape)
pert_vy = np.zeros(pert_vy.shape)

x = np.linspace(params["domain"]["xmin"], params["domain"]["xmax"], params["domain"]["xres"]+2)
x,_ = np.meshgrid(x, x)
pert = np.cos(x)*1e-3

reults_path = hd.run_sim(params, pert, pert_vx, pert_vy)
frames = os.path.join(os.path.join(reults_path, 'plots'), '*')
animation(frames=frames, path_out=reults_path)

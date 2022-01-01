# HyDES

## Concept

HydroDinamical Equation Solver Code to numerically solve the hydrodynamical equations in 2 dimensions in a regular rectangular mesh.
The code performs a first order staggered Lax-Friedrich scheme in time and a central first order derivative in space.
There are several initial conditions to test and a some boundary conditions, like 0 derivative in the boundaries or periodic ones.
The code is impleted originaly in IDL in may 2019 and it was transcripted to python in december 2021.

## Installation

To install the library just pip install it with:

    pip install hydes

## Getting started

To start simulating your models just load the library, load the parameters (template loaded with load_sample_params())
compute the desired initial conditions with a perturbation and simulate the evolution of the sistem like:


    import hydes as hd

    params = hd.load_sample_params()
    pert, pert_vx, pert_vy = hd.sample_perturbation()

    hd.run_sim(params, pert, pert_vx, pert_vy)


if you also want to have a movie of the saved plots you can do:


    import hydes as hd
    from animation import animation

    params = hd.load_sample_params()
    pert, pert_vx, pert_vy = hd.sample_perturbation()

    reults_path = hd.run_sim(params, pert, pert_vx, pert_vy)
    frames = os.path.join(os.path.join(reults_path, 'plots'), '*')
    animation(frames=frames, path_out=reults_path)

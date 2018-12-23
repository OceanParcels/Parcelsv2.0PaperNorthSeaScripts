
This repository contains the files necessary to run and postprocess the simulations studying the distribution of
floating microplastic in the North West European continental shelf and its sensitivity to different physical processes, which is discussed in *The Parcels v2.0 Lagrangian framework: new field interpolation schemes*, by Philippe Delandmeter and Erik van Sebille, submitted to *Geoscientific Model Development*.

The master file is `northsea_mp.py` that calls, using the different parameters, the simulation code `run_northsea_mp.py`. The kernels defining the particle dynamics are implemented in `northsea_mp_kernels.py`.

All postprocessing scripts and useful functions are stored in the `postprocess` directory.

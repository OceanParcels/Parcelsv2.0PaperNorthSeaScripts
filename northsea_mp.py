from run_northsea_mp import run_northsea_mp

outfile = '/scratch/shared/delandmeter/northSea_plastic/test_pr486.nc'
run_northsea_mp(outfile, nemo_res='0083', cmems=True, stokes=True, diffusion=1)

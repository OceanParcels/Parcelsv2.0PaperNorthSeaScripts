#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter

Simulation master file
"""

from run_northsea_mp import run_northsea_mp

outfile = '/scratch/shared/delandmeter/northSea_plastic/test.nc'
run_northsea_mp(outfile, nemo_res='0083', cmems=True, stokes=True, diffusion=1)

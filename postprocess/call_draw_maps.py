#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter

Master postprocessing file
"""


from draw_maps import draw


path = '../ncfiles/'

filename = 'nemo0083.nc'
draw(path+filename)

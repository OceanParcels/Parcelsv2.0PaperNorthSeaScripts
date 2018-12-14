#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter

Function postprocessing floating MP particles released in Thames and Rhiver estuaries
"""


import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import postprocess as postpro
from postprocess import NWcontinental_shelf_zones
from parcels import RectilinearZGrid, Field
import matplotlib as mpl
from os import path



def draw(filename): 
    p = postpro.ParticleData(filename)
    basename = path.splitext(path.basename(filename))[0]
    print('Total number of particles: %d' % p.npart)
    
    pland = np.logical_and(abs(p.lon[:,0]-p.lon[:,-1]) < 1e-4, abs(p.lat[:,0]-p.lat[:,-1]) < 1e-4)
    pind = np.where(np.logical_not(pland))
    p.keepParticles(pind)
    wetParticles = p.npart
    print('Number of initially wet particles: %d' % wetParticles)
    
    pland = np.logical_and(abs(p.lon[:,580]-p.lon[:,500]) < 1e-4, abs(p.lat[:,580]-p.lat[:,500]) < 1e-4)
    pind = np.where(np.logical_not(pland))
    beachingParticles = wetParticles - pind[0].shape[0]
    print('Number of beaching particles: %d' % beachingParticles)
    
    p.origin = np.where(p.lon[:,0]>2, 0, 1)
    
    pind = np.logical_and(p.age > 0*365, p.age <= 3*365)
    p.keepParticles(pind, 'time')
    assert abs(np.max(p.age)-3*365) < 10
    
    glon = np.arange(-120, 120,1)
    glat = np.arange(35, 90,.5)
    glon_m, glat_m = np.meshgrid(glon, glat)
    
    utils = postpro.Utils()
    
    pxi, pyi = utils.locate(p.lon, p.lat, glon, glat)
    map_count = utils.histogram(pxi, pyi, (len(glat), len(glon)))
    map_mean_age = utils.mean_age(pxi, pyi, p.age, (len(glat), len(glon)), map_count)
    map_touched = utils.touched(pxi, pyi, (len(glat), len(glon)))

    glon_f = np.arange(-120, 120,1/4.)
    glat_f = np.arange(35, 90,.5/4.)
    glon_fm, glat_fm = np.meshgrid(glon_f, glat_f)
    zones = NWcontinental_shelf_zones(glon_fm, glat_fm)
    pxi_f, pyi_f = utils.locate(p.lon, p.lat, glon_f, glat_f)
    zone_concentration = utils.zone_concentration(zones, p.age, 2,  pxi_f, pyi_f)
    np.save('%s_zone_concentration' % basename, zone_concentration)

    grid = RectilinearZGrid(glon, glat, mesh='spherical')
    cellSizeField = Field('cellsize', np.zeros((len(glat), len(glon))), grid=grid)
    cellSizeField.data[0,:] =  cellSizeField.cell_areas()
    pind = np.logical_and(p.age > 2*365, p.age <= 3*365)
    p.keepParticles(pind, 'time')
    pxi, pyi = utils.locate(p.lon, p.lat, glon, glat)
    map_density_3rd_year = utils.histogram(pxi, pyi, (len(glat), len(glon)))
    map_density_3rd_year = map_density_3rd_year.astype(np.float32) / cellSizeField.data[0,:] #/ p.npart

    
    def plot(data, title, cmap_name='plasma_r', log=False, vmin=None, vmax=None, under=None, over=None, show=True, fname=''):
        print '\n\n PLOTTING %s \ninto %s' % (title, fname)

        if fname:
            fig = plt.figure(figsize=(14, 8.5), dpi=150, facecolor='w', edgecolor='k')
        else:
            fig = plt.figure(figsize=(14, 8.5), dpi=60, facecolor='w', edgecolor='k')
        ax = fig.add_axes([.05, .05, .9, .9])
        #ttl = ax.set_title(title, fontsize=20)
        
        m = Basemap(width=8e6,height=5e6,
                    resolution='l',projection='stere',\
                    lat_ts=70,lat_0=70,lon_0=5.)
        
        m.drawcoastlines(zorder=3)
        m.drawparallels(np.arange(-90, 91, 10), labels=[True, False, False, False], fontsize=15, zorder=3)
        m.drawmeridians(np.arange(-180, 181, 30), labels=[False, False, False, True], fontsize=15, zorder=3)
        m.fillcontinents(color='blanchedalmond')

        legend_text = title
        xshift = 0 if len(legend_text) > 6 else .095
        ax.text(0.92+xshift, 1.05, legend_text,
                transform=ax.transAxes,
                verticalalignment='top',
                fontsize=24)
    
        if not under and not over:
            extend = 'neither'
            under = 'white'
            over = 'white'
        elif not under:
            extend = 'max'
            under = 'white'
        elif not over:
            extend = 'min'
            over = 'white'
        else:
            extend = 'both'
        cmap = plt.get_cmap(cmap_name)
        cmap.set_under(under)
        cmap.set_over(over)
        
        xs, ys = m(glon_m, glat_m)
        #m.pcolormesh(xs, ys, data, vmin=10, vmax=100, cmap=cmap);
        #m.pcolormesh(xs, ys, (data), norm=colors.LogNorm(vmin=.1, vmax=100), cmap=cmap);
    
        if vmin is None:
            vmin = np.min(data)
        if vmax is None:
            vmax = np.max(data)
        if under == 'white':
            data = np.ma.masked_where(data < vmin+1e-15 , data)
    
        if log:
            m.pcolormesh(xs, ys, data, cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax), zorder=2)
        else:
            m.pcolormesh(xs, ys, data, cmap=cmap, vmin=vmin, vmax=vmax, zorder=2)
    
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
    
        cbar = plt.colorbar(extend=extend, cax=cax)
        cbar.ax.tick_params(labelsize=15)
    
        if show:
            plt.show()
        if fname:
            plt.savefig(fname)
        plt.close(fig)

    
    #title = run_name + ': Particle mean age (days)'
    #plot(map_mean_age, title, cmap_name='plasma_r', log=False, vmin=0, vmax=3*365, show=False, fname='nemo_c'+runid+'_mean_age.png')
    #title = run_name + ': Particle mean age (days, saturated at one year)'
    #plot(map_mean_age, title, cmap_name='plasma_r', log=False, vmin=0, vmax=365, over='green', show=False, fname='nemo_c'+runid+'_mean_age_short.png')
    title = '# Part / km$^2$'
    plot(map_density_3rd_year, title, cmap_name='hot_r', log=True, vmin=1e-10, vmax=1e-5, under='white', over='blue', show=False, fname=basename+'_3rd_yr_density.png')
    title = '%'
    plot(map_touched, title, cmap_name='Spectral_r', log=True, vmin=.1, vmax=100, under='white', show=False, fname=basename+'_touched_log.png')
    #title = run_name + ': Affected cell proportion'
    #plot(map_touched, title, cmap_name='Spectral_r', log=False, vmin=.1, vmax=100, under='white', show=False, fname='nemo_c'+runid+'_touched.png')

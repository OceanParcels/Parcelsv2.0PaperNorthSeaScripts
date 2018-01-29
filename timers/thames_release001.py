from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, AdvectionRK4, ParticleFile
from parcels import compute_curvilinearGrid_rotationAngles, ErrorCode
from argparse import ArgumentParser
import numpy as np
import pytest
from os import path

ptype = {'scipy': ScipyParticle, 'jit': JITParticle}


def run_nemo_curvilinear(mode, outfile):
    data_path = '/data2/imau/oceanparcels/hydrodynamic_data/NEMO-MEDUSA/ORCA025-N006/'

    mesh_filename = data_path + 'domain/coordinates.nc'
    rotation_angles_filename = './rotation_angles.nc'
    variables = {'cosU': 'cosU',
                 'sinU': 'sinU',
                 'cosV': 'cosV',
                 'sinV': 'sinV'}
    dimensions = {'U': {'lon': 'glamu', 'lat': 'gphiu'},
                  'V': {'lon': 'glamv', 'lat': 'gphiv'},
                  'F': {'lon': 'glamf', 'lat': 'gphif'}}
    compute_curvilinearGrid_rotationAngles(mesh_filename, rotation_angles_filename, variables, dimensions)

    filenames = {'U': data_path + 'means/ORCA025-N06_2000????d05U.nc',
                 'V': data_path + 'means/ORCA025-N06_2000????d05V.nc',
                 'cosU': rotation_angles_filename,
                 'sinU': rotation_angles_filename,
                 'cosV': rotation_angles_filename,
                 'sinV': rotation_angles_filename}
    variables = {'U': 'uos',
                 'V': 'vos',
                 'cosU': 'cosU',
                 'sinU': 'sinU',
                 'cosV': 'cosV',
                 'sinV': 'sinV'}

    dimensions = {'U': {'lon': 'nav_lon', 'lat': 'nav_lat', 'time': 'time_counter'},
                  'V': {'lon': 'nav_lon', 'lat': 'nav_lat', 'time': 'time_counter'},
                  'cosU': {'lon': 'glamu', 'lat': 'gphiu'},
                  'sinU': {'lon': 'glamu', 'lat': 'gphiu'},
                  'cosV': {'lon': 'glamv', 'lat': 'gphiv'},
                  'sinV': {'lon': 'glamv', 'lat': 'gphiv'}}

    field_set = FieldSet.from_netcdf(filenames, variables, dimensions, mesh='spherical', allow_time_extrapolation=False)
    field_set.U.grid.lon[1:,1:] = field_set.cosU.grid.lon
    field_set.U.grid.lat[1:,1:] = field_set.cosU.grid.lat
    field_set.V.grid.lon[1:,1:] = field_set.cosV.grid.lon
    field_set.V.grid.lat[1:,1:] = field_set.cosV.grid.lat


    def DeleteParticle(particle, fieldset, time, dt):
       particle.delete()
 

    nTime = 100
    nStep = 86400*2

    time0 = field_set.U.grid.time[0]
    #lon_rotterdam = 3.7
    #lat_rotterdam = 51.9
    #lon_r = [lon_rotterdam for i in range(nTime)]
    #lat_r = [lat_rotterdam for i in range(nTime)]
    time_r = [time0 + nStep*i for i in range(nTime)]
    lon_thames = 1.35
    lat_thames = 51.59
    lon_t = [lon_thames for i in range(nTime)]
    lat_t = [lat_thames for i in range(nTime)]
    time_t = [time0 + nStep*i for i in range(nTime)]
    pset = ParticleSet.from_list(field_set, ptype[mode], lon=lon_t, lat=lat_t, time=time_t)
    kernel = AdvectionRK4
    pfile = ParticleFile(outfile, pset, type="indexed")
    pfile.write(pset, pset[0].time)
    for _ in range(350):
        pset.execute(kernel, runtime=86400, dt=900, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
        pfile.write(pset, pset[0].time)


if __name__ == "__main__":
    p = ArgumentParser(description="""Chose the mode using mode option""")
    p.add_argument('--mode', choices=('scipy', 'jit'), nargs='?', default='jit',
                   help='Execution mode for performing computation')
    args = p.parse_args()

    outfile = __file__[:-3]

    run_nemo_curvilinear(args.mode, outfile)

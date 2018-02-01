from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, AdvectionRK4, ParticleFile
from parcels import compute_curvilinearGrid_rotationAngles, ErrorCode
from parcels import timer
from argparse import ArgumentParser
import numpy as np
import pytest
from os import path
from glob import glob

ptype = {'scipy': ScipyParticle, 'jit': JITParticle}

def set_nemo_fieldset(ufiles, vfiles, mesh_mask=None):
    filenames = {'U': ufiles,
                 'V': vfiles}
    if mesh_mask:
        filenames['mesh_mask'] = mesh_mask

    variables = {'U': 'uo',
                 'V': 'vo'}
    dimensions = {'U': {'lon': 'nav_lon', 'lat': 'nav_lat', 'depth': 'depthu', 'time': 'time_counter'},
                  'V': {'lon': 'nav_lon', 'lat': 'nav_lat', 'depth': 'depthv', 'time': 'time_counter'}}
    if mesh_mask:
        return FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=False)
    else:
        return FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=False)

def run_nemo_curvilinear(mode, outfile):
    timer.fieldset = timer.Timer('FieldSet', parent=timer.nemo)
    data_path = '/data2/imau/oceanparcels/hydrodynamic_data/NEMO-MEDUSA/ORCA025-N006/'
    ufiles = sorted(glob(data_path+'means/ORCA025-N06_20000???d05U.nc'))
    vfiles = sorted(glob(data_path+'means/ORCA025-N06_20000???d05V.nc'))
    timer.fset_load = timer.Timer('fset loading', parent=timer.fieldset)
    field_set = set_nemo_fieldset(ufiles[0:3], vfiles[0:3], data_path + 'domain/coordinates.nc')
    timer.fset_load.stop()

    timer.fset_clean = timer.Timer('fset cleaning', parent=timer.fieldset)
    field_set.U.grid.lon[1:,1:] = field_set.cosU.grid.lon
    field_set.U.grid.lat[1:,1:] = field_set.cosU.grid.lat
    field_set.V.grid.lon[1:,1:] = field_set.cosV.grid.lon
    field_set.V.grid.lat[1:,1:] = field_set.cosV.grid.lat
    timer.fset_clean.stop()
    timer.fieldset.stop()
    timer.pset = timer.Timer('ParticleSet creation', parent=timer.nemo)


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
    timer.pset.stop()
    timer.exe = timer.Timer('ParticleSet execution', parent=timer.nemo)
    pfile = ParticleFile(outfile, pset, type="indexed")
    pfile.write(pset, pset[0].time)
    timer.run = timer.Timer('computation', parent=timer.exe, start=False)
    timer.advanceTime = timer.Timer('fset advancetime', parent=timer.exe, start=False)
    timer.write = timer.Timer('writing', parent=timer.exe, start=False)
    for d in range(100):
        timer.run.start()
        pset.execute(kernel, runtime=86400, dt=900, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
        timer.run.stop()
        timer.advanceTime.start()
        if d%5 == 1 and d/5 > 0:
            field_set.advancetime(set_nemo_fieldset(ufiles[d/5+2], vfiles[d/5+2]))
        timer.advanceTime.stop()
        timer.write.start()
        pfile.write(pset, pset[0].time)
        timer.write.stop()
    timer.exe.stop()


if __name__ == "__main__":
    timer.root = timer.Timer('Main')
    p = ArgumentParser(description="""Chose the mode using mode option""")
    p.add_argument('--mode', choices=('scipy', 'jit'), nargs='?', default='jit',
                   help='Execution mode for performing computation')
    args = p.parse_args()

    outfile = __file__[:-3]

    timer.nemo = timer.Timer('Nemo run', parent=timer.root)
    run_nemo_curvilinear(args.mode, outfile)
    timer.nemo.stop()
    timer.root.stop()
    timer.root.print_tree()

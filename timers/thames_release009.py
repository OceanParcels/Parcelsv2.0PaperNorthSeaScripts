from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, AdvectionRK4, ParticleFile
from parcels import ErrorCode
from parcels import timer
from argparse import ArgumentParser
import numpy as np
import pytest
from os import path
from glob import glob
import xarray as xr
import time as timelib

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
    field_set = FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=False)

    dataset = xr.open_dataset(mesh_mask, decode_times=False)
    lonU = np.array(dataset.glamu).squeeze()
    latU = np.array(dataset.gphiu).squeeze()
    lonV = np.array(dataset.glamv).squeeze()
    latV = np.array(dataset.gphiv).squeeze()
    lonT = np.array(dataset.glamt).squeeze()
    latT = np.array(dataset.gphit).squeeze()

    field_set.U.grid.lon = lonU
    field_set.U.grid.lat = latU
    field_set.V.grid.lon = lonV
    field_set.V.grid.lat = latV
    if not field_set.U.grid.lon.dtype == np.float32:
        field_set.U.grid.lon = field_set.U.grid.lon.astype(np.float32)
        field_set.U.grid.lat = field_set.U.grid.lat.astype(np.float32)
        field_set.V.grid.lon = field_set.V.grid.lon.astype(np.float32)
        field_set.V.grid.lat = field_set.V.grid.lat.astype(np.float32)
    dataset.close()
    return field_set


def run_nemo_curvilinear(mode, outfile):
    timer.fieldset = timer.Timer('FieldSet', parent=timer.nemo)
    data_path = '/data2/imau/oceanparcels/hydrodynamic_data/NEMO-MEDUSA/ORCA025-N006/'
    ufiles = sorted(glob(data_path+'means/ORCA025-N06_2000????d05U.nc'))
    vfiles = sorted(glob(data_path+'means/ORCA025-N06_2000????d05V.nc'))
    field_set = set_nemo_fieldset(ufiles, vfiles, data_path + 'domain/coordinates.nc')
    timer.fieldset.stop()

    timer.pset = timer.Timer('ParticleSet creation', parent=timer.nemo)
    def DeleteParticle(particle, fieldset, time, dt):
       particle.delete()


    time0 = field_set.U.grid.time[0]
    lonv = np.arange(-10,10.1,.05)
    latv = np.arange(45,66.1,.05)
    lon_t, lat_t = np.meshgrid(lonv, latv)
    lon_t = lon_t.flatten()
    lat_t = lat_t.flatten()
    time_t = time0 * np.ones(lon_t.shape)
    pset = ParticleSet.from_list(field_set, ptype[mode], lon=lon_t, lat=lat_t, time=time_t)
    kernel = AdvectionRK4
    timer.pset.stop()
    timer.pfile = timer.Timer('ParticleSet writing', parent=timer.nemo)
    pfile = ParticleFile(outfile, pset)
    pfile.write(pset, pset[0].time)
    timer.pfile.stop()
    tic = timelib.time()
    ndays = 30
    timer.run = timer.Timer('computation', parent=timer.nemo, start=False)
    for d in range(ndays):
        print('running %d / %d: time %g' % (d, ndays, timelib.time()-tic))
        timer.run.start()
        pset.execute(kernel, runtime=86400, dt=900)#, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
        timer.run.stop()
        timer.pfile.start()
        pfile.write(pset, pset[0].time)
        timer.pfile.stop()


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

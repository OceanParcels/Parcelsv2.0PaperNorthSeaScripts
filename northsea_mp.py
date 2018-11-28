from parcels import FieldSet, Field, VectorField, ParticleFile, ParticleSet, JITParticle, Variable
from parcels import ErrorCode
import numpy as np
from glob import glob
import time as timelib
from datetime import timedelta as delta
from northsea_mp_kernels import *

def get_nemo_fieldset(res='0083'):
    data_dir = '/projects/0/topios/hydrodynamic_data/NEMO-MEDUSA/ORCA%s-N006/' % res
    ufiles = sorted(glob(data_dir+'means/ORCA%s-N06_200?????d05U.nc' % res))
    vfiles = sorted(glob(data_dir+'means/ORCA%s-N06_200?????d05V.nc' % res))
    mesh_mask = data_dir + 'domain/coordinates.nc'

    filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
                 'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles}}


    variables = {'U': 'uo',
                 'V': 'vo'}
    dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'},
                  'V': {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}}

    fieldset = FieldSet.from_nemo(filenames, variables, dimensions)
    fieldset.U.vmax = 5
    fieldset.V.vmax = 5
    fieldset.nemo_res = res
    return fieldset

def set_nemo_unbeaching(fieldset):
    files = '/home/philippe/data/ORCA%s-N006_unbeaching_vel.nc' % fieldset.nemo_res
    filenames = {'unBeachU': files,
                 'unBeachV': files,
                 'mesh_mask': files}

    variables = {'unBeachU': 'unBeachU',
                 'unBeachV': 'unBeachV'}
    dimensions = {'lon': 'glamf', 'lat': 'gphif'}
    fieldsetUnBeach = FieldSet.from_nemo(filenames, variables, dimensions, tracer_interp_method='cgrid_linear')
    UVunbeach = VectorField('UV_unbeach', fieldsetUnBeach.unBeachU, fieldsetUnBeach.unBeachV)
    fieldset.add_field(fieldsetUnBeach.unBeachU)
    fieldset.add_field(fieldsetUnBeach.unBeachV)
    fieldset.add_vector_field(UVunbeach)

def get_particle_set():

    class PlasticParticle(JITParticle):
        age = Variable('age', dtype=np.float32, initial=0.)
        # beached : 0 sea, 1 beached, 2 after non-beach dyn, 3 after beach dyn, 4 please unbeach
        beached = Variable('beached', dtype=np.int32, initial=0.)

    # meshgrid containing 11x11 points uniformly distributed in a [0,1]x[0,1] quad
    vec = np.linspace(0,1,11)
    xsi, eta = np.meshgrid(vec, vec)
    
    # get particles in Rhine estuary
    lonCorners = [2.96824026, 3.22713804, 3.26175451, 3.002671]
    latCorners = [51.60693741, 51.58454132, 51.73711395, 51.759758] 
    lon_r = (1-xsi)*(1-eta) * lonCorners[0] + xsi*(1-eta) * lonCorners[1] + \
            xsi*eta * lonCorners[2] + (1-xsi)*eta * lonCorners[3]
    lat_r = (1-xsi)*(1-eta) * latCorners[0] + xsi*(1-eta) * latCorners[1] + \
            xsi*eta * latCorners[2] + (1-xsi)*eta * latCorners[3]
    
    # get particles in Thames estuary
    lonCorners = [1.37941658, 1.63887346, 1.67183721, 1.41217935]
    latCorners = [51.58309555, 51.56196213, 51.71636581, 51.73773575]
    lon_t = (1-xsi)*(1-eta) * lonCorners[0] + xsi*(1-eta) * lonCorners[1] + \
            xsi*eta * lonCorners[2] + (1-xsi)*eta * lonCorners[3]
    lat_t = (1-xsi)*(1-eta) * latCorners[0] + xsi*(1-eta) * latCorners[1] + \
            xsi*eta * latCorners[2] + (1-xsi)*eta * latCorners[3]
    
    # gather all particles, released every day for one year
    lons = np.concatenate((lon_r.flatten(), lon_t.flatten()))
    lats = np.concatenate((lat_r.flatten(), lat_t.flatten()))
    times = np.arange(np.datetime64('2000-01-05'), np.datetime64('2001-01-05'))
    
    return ParticleSet.from_list(fieldset, PlasticParticle,
                                 lon=np.tile(lons, [len(times)]),
                                 lat=np.tile(lats, [len(times)]),
                                 time=np.repeat(times, len(lons)))


fieldset = get_nemo_fieldset()
set_nemo_unbeaching(fieldset)
pset = get_particle_set()

kernel = AdvectionRK4 + pset.Kernel(BeachTesting) + pset.Kernel(UnBeaching) + pset.Kernel(Ageing)

outfile = '/scratch/shared/delandmeter/northSea_plastic/test.nc'
pfile = ParticleFile(outfile, pset)
pfile.write(pset, pset[0].time)

tic = timelib.time()
ndays = 365*4+100
for d in range(ndays/2):
    day = 2 * d
    print('running %d / %d [time %g s]: %d particles ' % (day+1, ndays, timelib.time()-tic, len(pset)))
    pset.execute(kernel, runtime=delta(days=2), dt=900, verbose_progress=False, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
    pfile.write(pset, pset[0].time)

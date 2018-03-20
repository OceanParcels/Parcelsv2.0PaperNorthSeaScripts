from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle
from parcels import timer
import numpy as np

ptype = {'scipy': ScipyParticle, 'jit': JITParticle}


def fieldset(xdim=40, ydim=100):
    U = np.zeros((ydim, xdim), dtype=np.float32)
    V = np.zeros((ydim, xdim), dtype=np.float32)
    lon = np.linspace(0, 1, xdim, dtype=np.float32)
    lat = np.linspace(-60, 60, ydim, dtype=np.float32)
    depth = np.zeros(1, dtype=np.float32)
    time = np.zeros(1, dtype=np.float64)
    data = {'U': np.array(U, dtype=np.float32), 'V': np.array(V, dtype=np.float32)}
    dimensions = {'lat': lat, 'lon': lon, 'depth': depth, 'time': time}
    return FieldSet.from_data(data, dimensions)


def test_pset_create_lon_lat(fieldset, mode, npart=100):
    lon = np.linspace(0, 1, npart, dtype=np.float32)
    lat = np.linspace(1, 0, npart, dtype=np.float32)
    pset = ParticleSet(fieldset, lon=lon, lat=lat, pclass=ptype[mode])

npart = 10000

timer.root = timer.Timer('Main')
timer.fieldset = timer.Timer('FieldSet', parent=timer.root)
fset = fieldset()
timer.fieldset.stop()
timer.pset = timer.Timer('ParticleSet', parent=timer.root)
test_pset_create_lon_lat(fset, 'jit', npart=npart)
timer.pset.stop()
timer.root.stop()

print('Creating %d particles' % npart)
#timer.root.print_tree()
timer.pset.print_local()


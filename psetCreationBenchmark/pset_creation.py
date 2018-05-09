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


class Obj(object):
    def __init__(self):
        self.i = 4

#@profile
def test_pset_create_lon_lat(fieldset, mode, npart=100):
    lon = np.linspace(0, 1, npart, dtype=np.float32)
    lat = np.linspace(1, 0, npart, dtype=np.float32)
    #pset = ParticleSet(fieldset, lon=lon, lat=lat, pclass=ptype[mode])
    for i in range(npart):
        p = JITParticle(0,0,fieldset=fieldset, depth=0, cptr=None, time=0)
        #p = Obj()

npart = 80000

timer.root = timer.Timer('Main')
timer.fieldset = timer.Timer('FieldSet', parent=timer.root)
fset = fieldset()
timer.fieldset.stop()
timer.pset = timer.Timer('ParticleSet', parent=timer.root)
timer.JITP = timer.Timer('JITParticle', parent=timer.pset, start=False)
timer.cptr = timer.Timer('cptr', parent=timer.JITP, start=False)
timer.super = timer.Timer('super', parent=timer.JITP, start=False)
timer.gridIndex = timer.Timer('GridIndex', parent=timer.JITP, start=False)
timer.gobj = timer.Timer('GridIndex obj', parent=timer.gridIndex, start=False)
timer.gobjptr = timer.Timer('GridIndex obj ptr', parent=timer.gridIndex, start=False)
timer.gobjptr_1 = timer.Timer('GridIndex obj ptr 1', parent=timer.gobjptr, start=False)
timer.gobjptr_2 = timer.Timer('GridIndex obj ptr 2', parent=timer.gobjptr, start=False)
timer.gobjptr_3 = timer.Timer('GridIndex obj ptr 3', parent=timer.gobjptr, start=False)
timer.gobjptrval = timer.Timer('GridIndex obj ptr value', parent=timer.gridIndex, start=False)
test_pset_create_lon_lat(fset, 'jit', npart=npart)
timer.pset.stop()
timer.root.stop()

print('Creating %d particles' % npart)
#timer.root.print_tree()
timer.pset.print_tree()


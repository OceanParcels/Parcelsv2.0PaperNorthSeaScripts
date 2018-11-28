def AdvectionRK4(particle, fieldset, time, dt):
    if particle.beached == 0:
        (u1, v1) = fieldset.UV[time, particle.lon, particle.lat, particle.depth]
        lon1, lat1 = (particle.lon + u1*.5*dt, particle.lat + v1*.5*dt)

        (u2, v2) = fieldset.UV[time + .5 * dt, lon1, lat1, particle.depth]
        lon2, lat2 = (particle.lon + u2*.5*dt, particle.lat + v2*.5*dt)

        (u3, v3) = fieldset.UV[time + .5 * dt, lon2, lat2, particle.depth]
        lon3, lat3 = (particle.lon + u3*dt, particle.lat + v3*dt)

        (u4, v4) = fieldset.UV[time + dt, lon3, lat3, particle.depth]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * dt
        particle.beached = 2

def BeachTesting(particle, fieldset, time, dt):
    if particle.beached == 2 or particle.beached == 3:
        (u, v) = fieldset.UV[time, particle.lon, particle.lat, particle.depth]
        if u == 0 and v == 0:
            if particle.beached == 2:
                particle.beached = 4
            else:
                particle.beached = 1
        else:
            particle.beached = 0
           

def UnBeaching(particle, fieldset, time, dt):
    if particle.beached == 4:
        (ub, vb) = fieldset.UV_unbeach[time, particle.lon, particle.lat, particle.depth]
        particle.lon += ub * dt
        particle.lat += vb * dt

def Ageing(particle, fieldset, time, dt):
    particle.age += dt 

def DeleteParticle(particle, fieldset, time, dt):
    print("Particle lost !! (%g %g %g %g)" % (particle.lon, particle.lat, particle.depth, particle.time))
    particle.delete()

from body import *
from ssystem import *
import numpy as np
from astropy import constants as const
from astropy import units as u



# a=Planet()
# b=Planet(mass=0.2,density=1000)
# a.roche_limit(b)

# mass=10e6 * u.g
# radius=10e12 *u.cm
# density=1000
#
# B=Body(mass=mass,radius=radius,density=density)
#

S=Star(mass=1,radius=1)
#
# E=Planet(mass=1,density=5514)
# m=Planet(mass=0.0122,density=3340)
#
# roche=E.roche_limit(m)

SS=Ssystem(S)

print("planetary mass:"+str(SS.planetaryMass))
print("adding bodies")
SS.addBodies()
SS.orbital_dist_histogram()

#Earth=Planet()
#Moon=Planet(mass=0.0123,orbit_radius=1.0025,density=3300)

#SS.bodies.append(Earth)
#SS.bodies.append(Moon)

#SS.mass_histogram()
#SS.orbital_dist_histogram()

print(str(len(SS.bodies))+" in the system")
print ("consolidating")
#SS.consolidate()
SS.dump()
SS.orbital_dist_histogram()

pass

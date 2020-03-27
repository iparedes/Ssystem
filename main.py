from body import *
from ssystem import *
import numpy as np



# a=Planet()
# b=Planet(mass=0.2,density=1000)
# a.roche_limit(b)

E=Planet(mass=1,density=5514)

SS=Ssystem()

print("planetary mass:"+str(SS.planetaryMass))
print("adding bodies")
SS.addBodies()

#Earth=Planet()
#Moon=Planet(mass=0.0123,orbit_radius=1.0025,density=3300)

#SS.bodies.append(Earth)
#SS.bodies.append(Moon)

#SS.mass_histogram()
#SS.orbital_dist_histogram()
print(str(len(SS.bodies))+" in the system")
print ("consolidating")
SS.consolidateBodies()
SS.dump()

pass

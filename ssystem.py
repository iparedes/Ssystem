from astropy import constants as const
from scipy.stats import lognorm
import numpy as np
import matplotlib.pyplot as plt

import random
import bisect
from body import *

# todo: use quantities for the attributes
# todo: check the distribution of densities
# Let's use sun mass in Sun masses and planet mass in Earth masses

PLANETARY_MASS_PCTG=0.2  # Planetary mass as a percentage of the solar mass. Originally 0.2

S=0.985 # Sigma for the lognorm distribution
SCALE=1/10000 # Scale for the lognorm distribution. eq to exp(mu). Originally 1/10000

# Density in Kg/m3
MIN_DENSITY=500
MAX_DENSITY=10000

# parameters for the mass distribution function
SIGMA_MASS=1.0
MU_MASS=10.0

class Ssystem:

    # solarMass in Sun masses
    def __init__(self, solarMass=1, planetaryMassPctg=PLANETARY_MASS_PCTG):
        self.mass_dist_normal_std= np.sqrt(np.log(1 + (SIGMA_MASS/MU_MASS)**2))
        self.mass_dist_normal_mean= np.log(MU_MASS) - self.mass_dist_normal_std**2 / 2

        self.planetaryMassPctg=planetaryMassPctg # total mass of planets as a pctg of the sun mass
        self.Sun=Star(solarMass)

        # The sphere of influence is the maximum distance for bodies orbiting in the system
        self.SOI=self.Sun.SOI()
        #self.SOI=const.au*10

        # total planetary mass in Earth masses
        self.planetaryMass=self.Sun.mass_SM*(const.M_sun/const.M_earth)*self.planetaryMassPctg
        #self.planetaryMass=(solarMass*self.planetaryMassPctg)*(const.M_sun/const.M_earth)

        self.bodies=[]
        self.remainingPlanetaryMassPctg=1


    def addBodies(self):
        while (self.remainingPlanetaryMassPctg>0):
            mass=self.randMass()
            orbit=self.randOrbit()
            dens=self.randDensity()
            print("adding mass: "+str(mass))
            nb=Planet(mass=mass,density=dens,orbit_radius=orbit)
            bisect.insort_left(self.bodies, nb)



    # https://gamedev.stackexchange.com/questions/151721/procedural-generation-of-semi-correct-planetary-systems
    # Starting from most to least massive, calculate the hill sphere of each planet.
    # Any less massive planet in the hill sphere of a more massive planet becomes a moon of that planet.
    # Randomly-generate the orbital radius of the moon around the parent with a logarithmic distribution
    # between 0 and the sphere of influence of the parent.
    def consolidateBodies(self):
        for p in self.bodies:
            print("consolidating "+str(p))
            idx=self.bodies.index(p)
            for q in self.bodies[idx+1:]:
                if p!=q:
                    hills=p.hill_sphere(self.Sun)
                    dist=p.orbital_distance(q)
                    if dist<hills:
                        # make it a moon
                        print("mooning!")
                        self.bodies.remove(q)
                        q.orbital_radius=dist
                        p.satellites.append(q)




    # returns mass (in Earth masses for a new body from the remaining planetary mass
    # returns 0 if there is not remaining planetary mass
    # it uses a lognormal distribution with params set at the creation of the Ssystem
    def randMass(self):
        if self.remainingPlanetaryMassPctg>0:

            # SCREWED UP THIS. THE NEW FUNCTION TO GET A RANDOM SAMPLE IS TUNED TO RETURN A MASS IN EM FOR
            # THE NEW PLANET. IT SHOULD RETURN A PERCENTAGE!!! (OR NOT)


            #sample=lognorm.rvs(S, scale=SCALE) # sample is a pctg of the remaining planetary mass
            sample=np.random.lognormal(self.mass_dist_normal_mean, self.mass_dist_normal_std)

            # if sample is bigger than remaining, just take everything
            if sample>=self.planetaryMass:
                new_mass=self.planetaryMass
                self.planetaryMass=0
            else:
                new_mass=sample
                self.planetaryMass-=sample
            # if sample>=self.remainingPlanetaryMassPctg:
            #     new_mass=self.remainingPlanetaryMassPctg*self.planetaryMass
            #     self.remainingPlanetaryMassPctg=0
            # else:
            #     new_mass=sample*self.planetaryMass
            #     self.remainingPlanetaryMassPctg-=sample
        else:
            new_mass=0
        return new_mass

    # returns the radius of a orbit inside the SOI of the system
    # In Astronomic Units
    def randOrbit(self):
        orbit=random.randrange(0,self.SOI)
        orbit=orbit/const.au.value
        return orbit

    def randDensity(self):
        mean=MIN_DENSITY+((MAX_DENSITY-MIN_DENSITY)/2)
        stdev=mean/3
        d=random.gauss(mean,stdev)
        return d


    def dump(self):
        print (str(len(self.bodies))+" bodies")
        data=[(x.mass_EM,x.radius_ER,x.orbit_radius_AU) for x in self.bodies]
        print("smallest mass: "+str(min(data)[0]))
        print("largest mass: "+str(max(data)[0]))
        print("smallest size: "+str(min(data)[1]))
        print("largest size: "+str(max(data)[1]))
        print("closest: "+str(min(data)[2]))
        print("furthest: "+str(max(data)[2]))
        print
        for p in self.bodies:
            print(p.dump())

    def mass_histogram(self):
        masses=[x.mass_EM for x in self.bodies]
        plt.hist(masses)
        plt.show()

    def orbital_dist_histogram(self):
        dists=[x.orbit_radius/const.au for x in self.bodies]
        plt.hist(dists)
        plt.show()







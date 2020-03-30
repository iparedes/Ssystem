from astropy import constants as const
from astropy import units as u
from scipy.stats import *
import numpy as np
import matplotlib.pyplot as plt
import quantities as pq

import random
import bisect
from body import *

# todo: use quantities for the attributes. Right now is a freaking mess mixing values and quantities
# todo: check the distribution of densities
# todo: restrict original radius. The SOI edge is to far. Something between 10AU adn 50AU mabye...

PLANETARY_MASS_PCTG=0.0015  # Planetary mass as a percentage of the solar mass.

S=0.985 # Sigma for the lognorm distribution
SCALE=1/10000 # Scale for the lognorm distribution. eq to exp(mu). Originally 1/10000

# Density in Kg/m3
MIN_DENSITY=500
MAX_DENSITY=10000

# parameters for the mass distribution function
MU_MASS=1.0
SIGMA_MASS=10.0

class Ssystem:

    # solarMass in Sun masses
    def __init__(self, solarMass=1, planetaryMassPctg=PLANETARY_MASS_PCTG):
        self.mass_dist_normal_std= np.sqrt(np.log(1 + (SIGMA_MASS/MU_MASS)**2))
        self.mass_dist_normal_mean= np.log(MU_MASS) - self.mass_dist_normal_std**2 / 2

        self.planetaryMassPctg=planetaryMassPctg # total mass of planets as a pctg of the sun mass

        density=1410
        self.Sun=Star(mass=1,radius=1,density=density)

        # The sphere of influence is the maximum distance for bodies orbiting in the system
        self.SOI=self.Sun.SOI()
        # todo: SOI is maybe too much, putting a cap here
        self.SOI=random.randrange(10,30)*u.au
        #self.SOI=5*const.au.value
        #self.SOI=const.au*10

        # total planetary mass in Earth masses
        self.planetaryMass=self.Sun.mass*self.planetaryMassPctg
        self.planetaryMass=self.planetaryMass.to('M_earth')
        #self.planetaryMass=(solarMass*self.planetaryMassPctg)*(const.M_sun/const.M_earth)
        self.remainingPlanetaryMass=u.Quantity(self.planetaryMass)

        self.bodies=[]

        # create the function for the distribution of masses
        # todo: parametrize mean,desv,a,b
        self.mass_distr=self._get_truncated_normal(0,0.8,0,100)


    def _get_truncated_normal(self,mean=0, sd=1, low=0, upp=10):
        return truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)


    def addBodies(self):
        #while (self.remainingPlanetaryMassPctg>0):
        while (self.planetaryMass>0):
            mass=self.randMass()
            orbit=self.randOrbit()
            dens=self.randDensity()
            print("adding mass: "+str(mass))
            nb=Planet(mass=mass,density=dens,orbit_radius=orbit)
            #nb=Body(mass=mass,density=dens,orbit_radius=orbit)
            bisect.insort_left(self.bodies, nb)


    # Perform step 2 for all moon systems to resolve hill-sphere conflicts of moons.
    # Whether a moon can have a stable satellite is a matter of debate among the astronomy community
    # (no example is known in our solar system). When you dont want any moon-moons,
    # simply delete the smaller moon or put it on a different random orbit.
    def consolidate(self):
       bodies=0
       bodies_now=len(self.bodies)

       while (bodies!=bodies_now):
           self.stats()
           bodies=bodies_now
           # Creates moons
           self.consolidateBodies(self.bodies)

           # Look for impossible moons
           for b in self.bodies:
               if len(b.satellites)>1:
                   b.satellites.sort(key=lambda x: x.mass, reverse=True)
                   for p in b.satellites:
                       idx=b.satellites.index(p)
                       for q in b.satellites[idx+1:]:
                            if p!=q:
                                hills=p.hill_sphere(b)
                                dist=p.orbital_distance(q)
                                if dist<hills:
                                    print("breaking moon")
                                    # q is in the hill_sphere of p, so we remove it and add to the planet
                                    b.satellites.remove(q)
                                    b.mass+=q.mass
           self.roching()
           bodies_now=len(self.bodies)



    # https://gamedev.stackexchange.com/questions/151721/procedural-generation-of-semi-correct-planetary-systems
    # Starting from most to least massive, calculate the hill sphere of each planet.
    # Any less massive planet in the hill sphere of a more massive planet becomes a moon of that planet.
    # Randomly-generate the orbital radius of the moon around the parent with a logarithmic distribution
    # between 0 and the sphere of influence of the parent.
    def consolidateBodies(self,bodies):
        for p in self.bodies:
            print("consolidating "+str(p))
            idx=bodies.index(p)
            for q in bodies[idx+1:]:
                if p!=q:
                    hills=p.hill_sphere(self.Sun)
                    dist=p.orbital_distance(q)
                    if dist<hills:
                        # make it a moon
                        print("mooning!")
                        bodies.remove(q)
                        q.orbital_radius=dist
                        p.satellites.append(q)

    def roching(self):
        self.bodies.sort(key=lambda x: x.mass, reverse=True)
        for p in self.bodies:
            idx=self.bodies.index(p)
            for q in self.bodies[idx+1:]:
                if p!=q:
                    rochelimit=p.roche_limit(q)
                    dist=p.orbital_distance(q)
                    print("dist: "+str(dist)+" rocheL: "+str(rochelimit))
                    if dist<rochelimit:
                        print("Roching!!")
                        self.bodies.remove(q)
                        mass=q.mass+sum([x.mass for x in q.satellites])
                        p.mass+=mass


    # returns mass (in Earth masses for a new body from the remaining planetary mass
    # returns 0 if there is not remaining planetary mass
    # it uses a lognormal distribution with params set at the creation of the Ssystem
    def randMass(self):
        if self.remainingPlanetaryMass>0:


            #sample=lognorm.rvs(S, scale=SCALE) # sample is a pctg of the remaining planetary mass
            #sample=np.random.lognormal(self.mass_dist_normal_mean, self.mass_dist_normal_std)

            # sample is in earth masses
            sample=20*self.mass_distr.rvs() * u.M_earth

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
        a=self.SOI.to('au')
        orbit=random.uniform(0,a.value)
        #orbit=orbit/const.au.value
        return orbit

    def randDensity(self):
        mean=MIN_DENSITY+((MAX_DENSITY-MIN_DENSITY)/2)
        stdev=mean/3
        d=random.gauss(mean,stdev)
        return d



    def stats(self):
        print("SOI: "+str(self.SOI.to('lyr'))+" light years")
        print("# planets: "+str(len(self.bodies)))
        nmoons=sum([len(x.satellites) for x in self.bodies])
        print("# moons:"+str(nmoons))

    def dump(self):
        self.stats()

        self.bodies.sort(key=lambda x: x.orbit_radius, reverse=True)

        data=[(x.mass,x.radius,x.orbit_radius) for x in self.bodies]
        mins=list(map(min,zip(*data)))
        maxs=list(map(max,zip(*data)))

        print("smallest mass: "+str(mins[0].to('M_earth')))
        print("largest mass: "+str(maxs[0].to('M_earth')))
        print("smallest size: "+str(mins[1].to('R_earth')))
        print("largest size: "+str(maxs[1].to('R_earth')))
        print("closest: "+str(mins[2].to('au')))
        print("furthest: "+str(maxs[2].to('au')))
        print
        for p in self.bodies:
            print(p.dump())
        pass

    def mass_histogram(self):
        masses=[x.mass.to('M_earth').value for x in self.bodies]
        plt.hist(masses)
        plt.show()

    def orbital_dist_histogram(self):
        dists=[x.orbit_radius/const.au for x in self.bodies]
        plt.hist(dists)
        plt.show()







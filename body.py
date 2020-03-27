from astropy import constants as const
from math import pi

# Galaxy mass in sun masses
M_galaxy=2e12

# Sun distance to the galaxy center in meters
# Used to calculate the SOI of a star
# We are assuming that the generated solar systems are in remote corners of the galaxy, like ours
D_galaxy_center=2.5733e20

class Body:
    # mass and radius in Kgs and meters
    # orbit_radius [m]
    # density [Kg/m3]
    # If a density is given, then disregards the given radius (if any) and calculates it according to the mass
    def __init__(self,mass=0,radius=0,orbit_radius=0,density=0):
        self.mass=mass
        self.radius=radius

        # if density is given let's calculate the radius (even if a radius is given)
        if density>0:
            self.density=density
            self.radius=((3*self.mass.value)/(4*pi*self.density))**(1/3)
        else:
            # if radius is given let's calculate the density
            if self.radius>0:
                self.density=(3*self.mass)/(4*pi*(self.radius**3))

        # parent is the body at the centre of the orbit of the body
        self.parent=None
        self.orbit_radius=orbit_radius

        self.satellites=[]

    # bodies are compared acoording to their mass
    def __lt__(self, other):
        return self.mass < other.mass



    # returns the distance between the orbits of two bodies (m)
    def orbital_distance(self,other):
        return abs(self.orbit_radius-other.orbit_radius)

    # returns the radius of the Hill sphere of two bodies
    def hill_sphere(self,other):
        m,M=min(self.mass,other.mass),max(self.mass,other.mass)
        r=self.orbital_distance(other)*(m/(3*M))**(1/3)
        return r

    # returns the Roche limit of two bodies (primary and secondary) [m]
    def roche_limit(self,other):
        [pri,sec]=sorted([self,other],key=lambda x: x.mass, reverse=True)
        rl=pri.radius*(2*pri.density/sec.density)**(1/3)
        return rl

class Planet(Body):

    # mass and radius are entered in earth units, then converted to Kgs and meters
    # orbit_radius is entered in AU
    def __init__(self,mass=1,radius=1,orbit_radius=1,density=5514):

        masskg=mass*const.M_earth
        radiusm=radius*const.R_earth
        orbit_radiusm=orbit_radius*const.au
        Body.__init__(self,mass=masskg,radius=radiusm,orbit_radius=orbit_radiusm,density=density)

        # attributes for mass and radius in Earth units
        self.mass_EM=mass
        self.radius_ER=self.radius/const.R_earth
        self.orbit_radius_AU=orbit_radius

    def dump(self):
        cad="mass: "+str(self.mass_EM)+" radius: "+str(self.radius_ER)+" satellites: "+str(len(self.satellites))
        return cad


class Star(Body):
    # mass and radius are entered in sun units, then converted to Kgs and meters
    def __init__(self,mass=1,radius=1,orbit_radius=0,density=1410):

        masskg=mass*const.M_sun
        radiusm=radius*const.R_sun
        Body.__init__(self,mass=masskg,radius=radiusm,orbit_radius=orbit_radius,density=density)


        # attributes for mass and radius in Sun units
        self.mass_SM=mass
        self.radius_SR=self.radius/const.R_sun
        self.orbit_radius=0



    # returns the radius (m) of the sphere of influence of the Sun
    # Uses the mass of the galaxy as second body
    def SOI(self):
        soi=D_galaxy_center*(self.mass_SM/M_galaxy)**(2/5)
        return soi

from astropy import constants as const
from astropy import units as u
from math import pi


# Galaxy mass in sun masses
M_galaxy=2e12 * u.M_sun

# Sun distance to the galaxy center in meters
# Used to calculate the SOI of a star
# We are assuming that the generated solar systems are in remote corners of the galaxy, like ours
D_galaxy_center=2.5733e20 * u.m


class Body:
    # mass, radius, orbit_radius and density should be quantities (value and unit)
    # if no unit is given, converts as followsÑ
    # mass and radius in Kgs and meters
    # orbit_radius [m]
    # density [Kg/m3]
    # If a density is given, then disregards the given radius (if any) and calculates it according to the mass
    def __init__(self,mass,radius,orbit_radius=0,density=0,parent=None):

        try:
            self.mass=mass.to('kg')
        except AttributeError:
            self.mass=mass*u.kg

        try:
            self.radius=radius.to('m')
        except AttributeError:
            self.radius=radius*u.m

        try:
            self.orbit_radius=orbit_radius.to('m')
        except AttributeError:
            self.orbit_radius=orbit_radius*u.m

        # if density is given let's calculate the radius (even if a radius is given)
        if density>0:
            try:
                self.density=density.to(u.kg / u.m**3)
            except AttributeError:
                self.density=density*u.kg/u.m**3
            self.radius=((3*self.mass)/(4*pi*self.density))**(1/3)
        else:
            # if radius is given let's calculate the density
            if self.radius>0:
                self.density=(3*self.mass)/(4*pi*(self.radius**3))

        # parent is the body at the centre of the orbit of the body
        self.parent=parent
        self.orbit_radius=orbit_radius

        self.satellites=[]

    # bodies are compared acoording to their mass
    def __lt__(self, other):
        return self.mass < other.mass

    # returns the distance between the orbits of two bodies (m)
    def orbital_distance(self,other):
        return abs(self.orbit_radius-other.orbit_radius)


    # Hill Sphere: given a large mass (eg Sun) and a small mass (eg Earth),
    # can a tiny mass (eg Moon) find a stable orbit around the small mass?
    # (If the tiny mass goes outside the Hill Sphere of the small mass, no.)
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

    # SOI: given two large mass objects and a small object between them,
    # (eg sending a probe from Earth to Mars),
    # which massive object should we use as the origin of the frame of reference?
    # (The small object is within which massive object's SOI?)
    # returns the radius (m) of the sphere of influence of the body
    # Uses the mass of the parent as second body
    def SOI(self):
        soi=self.orbit_radius*(self.mass/self.parent.mass)**(2/5)
        return soi



class Planet(Body):

    # mass and radius are entered in earth units
    # orbit_radius is entered in AU
    # desity [kg/m^3]
    # expect floats and converts them to quantities
    def __init__(self,mass=1,radius=1,orbit_radius=1,density=5514,parent=None):

        try:
            mass=mass.to('M_earth')
        except AttributeError:
            mass=mass*u.M_earth

        try:
            radius=radius.to('R_earth')
        except AttributeError:
            radius=radius*u.R_earth

        try:
            orbit_radius=orbit_radius.to('au')
        except AttributeError:
            orbit_radius=orbit_radius*u.au

        try:
            density=density.to(u.kg / u.m**3)
        except AttributeError:
            density=density*u.kg/u.m**3

        Body.__init__(self,mass=mass,radius=radius,orbit_radius=orbit_radius,density=density,parent=parent)


    def dump(self):
        cad="mass: "+str(self.mass.to('M_earth'))+" radius: "+str(self.radius.to('R_earth'))+"distance: "+str(self.orbit_radius)+" satellites: "+str(len(self.satellites))
        return cad


class Star(Body):

    # mass and radius are entered in sun units
    # orbit_radius is entered in AU
    # desity [kg/m^3]
    # expect floats and converts them to quantities
    def __init__(self,mass=1,radius=1,orbit_radius=0,density=1410):


        try:
            mass=mass.to('M_sun')
        except AttributeError:
            mass=mass*u.M_sun

        try:
            radius=radius.to('R_sun')
        except AttributeError:
            radius=radius*u.R_sun

        try:
            orbit_radius=orbit_radius.to('au')
        except AttributeError:
            orbit_radius=orbit_radius*u.au

        try:
            density=density.to(u.kg / u.m**3)
        except AttributeError:
            density=density*u.kg/u.m**3

        Body.__init__(self,mass=mass,radius=radius,orbit_radius=orbit_radius,density=density)


    # returns the radius (m) of the sphere of influence of the Sun
    # Uses the mass of the galaxy as second body
    def SOI(self):
        soi=D_galaxy_center*(self.mass.to('M_sun')/M_galaxy)**(2/5)
        return soi.to('au')

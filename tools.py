
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

# Returns a normal distribution between two values
# It creates the normal with the values in the extremes
def normal_distro_AB(a,b):
    low,upp=min(a,b),max(a,b)
    mean=low+(upp-low)/2
    sd=(upp-low)/4
    #distro=scipy.stats.norm(mu,sd)
    distro=scipy.stats.truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

    return distro

# Creates a normal distribution to the right of the mean where a is the man and b/4 the sd
def truncated_normal(a,b):
    low,upp=min(a,b),max(a,b)
    mean=0
    sd=upp/4
    return scipy.stats.truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)


def log_normal(loc,s):
    #mu, sigma = 0.5, 4
    #normal_std = np.sqrt(np.log(1 + (sigma/mu)**2))
    #normal_mean = np.log(mu) - normal_std**2 / 2
    #distro = np.random.lognormal(normal_mean, normal_std)
    distro= scipy.stats.lognorm(loc=loc,s=s)
    #samples_log = np.random.lognormal(normal_mean, normal_std, 1000)
    return distro


#
# #distro=normal_distro_AB(50,100)
# distro=truncated_normal(8,12)
# samples=distro.rvs(10000)
#
# plt.hist(samples,bins=100)
# plt.show()
#
# pass


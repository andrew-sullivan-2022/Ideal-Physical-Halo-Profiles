import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as int
import h5py
import statistics as stat
import random

x=150

#Defining Variables
l = np.linspace(0, 3)
m = np.linspace(0, 4)
a = np.linspace(0, 1.5, num=x)
b = np.linspace(0, 0.5, num=x)
c = np.linspace(0, 100)
s = np.logspace(-4, 4, base=10, num=x)
t = np.logspace(-4, 4, base=10, num=x)
u = np.logspace(-4, 4, base=10, num=x)
S = np.logspace(-3, 1, base=10, num=x)
P = np.logspace(-3, 1, base=10, num=x)


mass = np.logspace(5, 20, base=10, num=1000)
sigma = np.linspace(0, 10**4, num=1000)


s50 = np.logspace(-4, 4, base=10, num=50)
s150 = np.logspace(-4, 4, base=10, num=150)
S25 = np.logspace(-3, 1, base=10, num=25)
S50 = np.logspace(-3, 1, base=10, num=50)


a25 = np.linspace(0, 1.5, num=25)
b25 = np.linspace(0, 0.5, num=25)

#Physical Constants
G = 4.3009*10**(-9)


def crit_density(h):
    return 2.7754*10**(11)*h**2


#Defining Functions

NFW_b0_c5_Radial_vel_disp = []
NFW_b0_c10_Radial_vel_disp = []
NFW_b0_c20_Radial_vel_disp = []
NFW_bquarter_c10_Radial_vel_disp = []
NFW_bhalf_c10_Radial_vel_disp = []
NFW_bquarter_c5_Radial_vel_disp = []
NFW_bhalf_c5_Radial_vel_disp = []
NFW_b0_c5_LOS_vel_disp = []
NFW_b0_c10_LOS_vel_disp = []
NFW_b0_c20_LOS_vel_disp = []
NFW_bquarter_c10_LOS_vel_disp = []
NFW_bhalf_c10_LOS_vel_disp = []
NFW_bquarter_c5_LOS_vel_disp = []
NFW_bhalf_c5_LOS_vel_disp = []
NFW_b0_c5_ap_vel_disp = []
NFW_b0_c10_ap_vel_disp = []
NFW_b0_c20_ap_vel_disp = []
NFW_bquarter_c10_ap_vel_disp = []
NFW_bhalf_c10_ap_vel_disp = []
NFW_bquarter_c5_ap_vel_disp = []
NFW_bhalf_c5_ap_vel_disp = []
a_0_b0_c10_Radial_vel_disp = []
a_05_b0_c10_Radial_vel_disp = []
a_15_b0_c10_Radial_vel_disp = []
a_0_bquarter_c10_Radial_vel_disp = []
a_05_bquarter_c10_Radial_vel_disp = []
a_15_bquarter_c10_Radial_vel_disp = []
a_0_bhalf_c10_Radial_vel_disp = []
a_05_bhalf_c10_Radial_vel_disp = []
a_15_bhalf_c10_Radial_vel_disp = []
a_0_b0_c5_Radial_vel_disp = []
a_05_b0_c5_Radial_vel_disp = []
a_15_b0_c5_Radial_vel_disp = []
a_0_bquarter_c5_Radial_vel_disp = []
a_05_bquarter_c5_Radial_vel_disp = []
a_15_bquarter_c5_Radial_vel_disp = []
a_0_bhalf_c5_Radial_vel_disp = []
a_05_bhalf_c5_Radial_vel_disp = []
a_15_bhalf_c5_Radial_vel_disp = []
a_0_b0_c10_LOS_vel_disp = []
a_05_b0_c10_LOS_vel_disp = []
a_15_b0_c10_LOS_vel_disp = []
a_0_bquarter_c10_LOS_vel_disp = []
a_05_bquarter_c10_LOS_vel_disp = []
a_15_bquarter_c10_LOS_vel_disp = []
a_0_bhalf_c10_LOS_vel_disp = []
a_05_bhalf_c10_LOS_vel_disp = []
a_15_bhalf_c10_LOS_vel_disp = []
a_0_b0_c5_LOS_vel_disp = []
a_05_b0_c5_LOS_vel_disp = []
a_15_b0_c5_LOS_vel_disp = []
a_0_bquarter_c5_LOS_vel_disp = []
a_05_bquarter_c5_LOS_vel_disp = []
a_15_bquarter_c5_LOS_vel_disp = []
a_0_bhalf_c5_LOS_vel_disp = []
a_05_bhalf_c5_LOS_vel_disp = []
a_15_bhalf_c5_LOS_vel_disp = []
a_0_b0_c10_ap_vel_disp = []
a_05_b0_c10_ap_vel_disp = []
a_15_b0_c10_ap_vel_disp = []
a_0_bquarter_c10_ap_vel_disp = []
a_05_bquarter_c10_ap_vel_disp = []
a_15_bquarter_c10_ap_vel_disp = []
a_0_bhalf_c10_ap_vel_disp = []
a_05_bhalf_c10_ap_vel_disp = []
a_15_bhalf_c10_ap_vel_disp = []
a_0_b0_c5_ap_vel_disp = []
a_05_b0_c5_ap_vel_disp = []
a_15_b0_c5_ap_vel_disp = []
a_0_bquarter_c5_ap_vel_disp = []
a_05_bquarter_c5_ap_vel_disp = []
a_15_bquarter_c5_ap_vel_disp = []
a_0_bhalf_c5_ap_vel_disp = []
a_05_bhalf_c5_ap_vel_disp = []
a_15_bhalf_c5_ap_vel_disp = []
xi_NFW_c10_rv_01 = []
xi_NFW_c5_rv_01 = []
xi_NFW_c10_rv_05 = []
xi_NFW_c5_rv_05 = []
xi_b0_c10_rv_01 = []
xi_b0_c5_rv_01 = []
xi_b0_c10_rv_05 = []
xi_b0_c5_rv_05 = []
xi_bquarter_c10_rv_01 = []
xi_bquarter_c5_rv_01 = []
xi_bquarter_c10_rv_05 = []
xi_bquarter_c5_rv_05 = []
xi_bhalf_c10_rv_01 = []
xi_bhalf_c5_rv_01 = []
xi_bhalf_c10_rv_05 = []
xi_bhalf_c5_rv_05 = []










#Defining Integrals
def general_conc_func(m, l, c):
    return 1/(int.quad(lambda s: (s**(2-m))/((1+c*s)**(l-m)), 0, 1)[0])


def conc_func(a, c):
    return 1/(int.quad(lambda s: (s**(2-a))/((1+c*s)**(3-a)), 0, 1)[0])


def mass_int(a, c, s):
    return int.quad(lambda t: t**(2-a)/((1+c*t)**(3-a)), 0, s)[0]


def potential_int(a, c, s):
    return int.quad(lambda s: s**(1-a)/((1+c*s)**(3-a)), s, np.infty)[0]


def surface_int(a, c, S):
    return int.quad(lambda s: (s**(1-a))/(((1+c*s)**(3-a))*((s**2 - S**2)**0.5)), S, np.infty)[0]


def radial_disp_int(a, b, c, s):
    return int.quad(lambda u: (u**(2*b-a-2)/((1+c*u)**(3-a))*mass_int(a, c, u)), s, np.infty)[0]


def los_disp_int(a, b, c, S):
    return int.quad(lambda s: (((1 - b*(S/s)**2)*s**(1-2*b))/((s**2 - S**2)**0.5))*radial_disp_int(a, b, c, s), S, np.infty)[0]


def iso_los_disp_int(a, c, S):
    return int.quad(lambda s: (((s**2 - S**2)**0.5)/((s**(a+2))*((1+c*s)**(3-a))))*(mass_int(a, c, s)), S, np.infty)[0]


def aperture_surface_int(a, c, S):
    return int.quad(lambda P: P*surface_int(a, c, P), 0, S)[0]


def aperture_disp_int(a, b, c, S):
    return int.quad(lambda P: P*los_disp_int(a, b, c, P), 0, S)[0]


def iso_aperture_disp_int(a, c, S):
    return int.quad(lambda P: P*iso_los_disp_int(a, c, P), 0, S)[0]










#Density Profiles
def General_Density(m, l, c, s):
    return general_conc_func(m, l, c)/(3*s**m*(1+c*s)**(l-m))


def Density(a, c, s):
    return conc_func(a, c)/(3*s**a*(1+c*s)**(3-a))


#Radial Velocity Dispersion
def Radial_vel_disp(a, b, c, s):
    return (conc_func(a, c)*s**(a - 2*b)*(1+c*s)**(3 - a)*radial_disp_int(a, b, c, s))**0.5


#Line of Sight Velocity Dispersion
def LOS_vel_disp(a, b, c, S):
    return (conc_func(a, c)*los_disp_int(a, b, c, S)/surface_int(a, c, S))**0.5


#Isotropic Line of Sight Velocity Dispersion
def Iso_LOS_vel_disp(a, c, S):
    return (conc_func(a, c)*iso_los_disp_int(a, c, S)/surface_int(a, c, S))**0.5


#Aperture Velocity Dispersion
def Aperture_vel_disp(a, b, c, S):
    return (conc_func(a, c)*aperture_disp_int(a, b, c, S)/aperture_surface_int(a, c, S))**0.5


#Isotropic Aperture Velocity Dispersion
def Iso_Aperture_vel_disp(a, c, S):
    return (conc_func(a, c)*iso_aperture_disp_int(a, c, S)/aperture_surface_int(a, c, S))**0.5










#Loading Curves

X_NFW_b0_c5_Radial_vel_disp = np.loadtxt('NFW_b0_c5_Radial_vel_disp.txt', delimiter = ',')
X_NFW_b0_c10_Radial_vel_disp = np.loadtxt('NFW_b0_c10_Radial_vel_disp.txt', delimiter = ',')
X_NFW_bquarter_c5_Radial_vel_disp = np.loadtxt('NFW_bquarter_c5_Radial_vel_disp.txt', delimiter = ',')
X_NFW_bquarter_c10_Radial_vel_disp = np.loadtxt('NFW_bquarter_c10_Radial_vel_disp.txt', delimiter = ',')
X_NFW_bhalf_c5_Radial_vel_disp = np.loadtxt('NFW_bhalf_c5_Radial_vel_disp.txt', delimiter = ',')
X_NFW_bhalf_c10_Radial_vel_disp = np.loadtxt('NFW_bhalf_c10_Radial_vel_disp.txt', delimiter = ',')
X_NFW_b0_c5_LOS_vel_disp = np.loadtxt('NFW_b0_c5_LOS_vel_disp.txt', delimiter = ',')
X_NFW_b0_c10_LOS_vel_disp = np.loadtxt('NFW_b0_c10_LOS_vel_disp.txt', delimiter = ',')
X_NFW_bquarter_c5_LOS_vel_disp = np.loadtxt('NFW_bquarter_c5_LOS_vel_disp.txt', delimiter = ',')
X_NFW_bquarter_c10_LOS_vel_disp = np.loadtxt('NFW_bquarter_c10_LOS_vel_disp.txt', delimiter = ',')
X_NFW_bhalf_c5_LOS_vel_disp = np.loadtxt('NFW_bhalf_c5_LOS_vel_disp.txt', delimiter = ',')
X_NFW_bhalf_c10_LOS_vel_disp = np.loadtxt('NFW_bhalf_c10_LOS_vel_disp.txt', delimiter = ',')
X_NFW_b0_c5_ap_vel_disp = np.loadtxt('NFW_b0_c5_ap_vel_disp.txt', delimiter = ',')
X_NFW_b0_c10_ap_vel_disp = np.loadtxt('NFW_b0_c10_ap_vel_disp.txt', delimiter = ',')
X_NFW_bquarter_c5_ap_vel_disp = np.loadtxt('NFW_bquarter_c5_ap_vel_disp.txt', delimiter = ',')
X_NFW_bquarter_c10_ap_vel_disp = np.loadtxt('NFW_bquarter_c10_ap_vel_disp.txt', delimiter = ',')
X_NFW_bhalf_c5_ap_vel_disp = np.loadtxt('NFW_bhalf_c5_ap_vel_disp.txt', delimiter = ',')
X_NFW_bhalf_c10_ap_vel_disp = np.loadtxt('NFW_bhalf_c10_ap_vel_disp.txt', delimiter = ',')
X_a_0_b0_c5_Radial_vel_disp = np.loadtxt('a_0_b0_c5_Radial_vel_disp.txt', delimiter = ',')
X_a_0_b0_c10_Radial_vel_disp = np.loadtxt('a_0_b0_c10_Radial_vel_disp.txt', delimiter = ',')
X_a_05_b0_c5_Radial_vel_disp = np.loadtxt('a_05_b0_c5_Radial_vel_disp.txt', delimiter = ',')
X_a_05_b0_c10_Radial_vel_disp = np.loadtxt('a_05_b0_c10_Radial_vel_disp.txt', delimiter = ',')
X_a_15_b0_c5_Radial_vel_disp = np.loadtxt('a_15_b0_c5_Radial_vel_disp.txt', delimiter = ',')
X_a_15_b0_c10_Radial_vel_disp = np.loadtxt('a_15_b0_c10_Radial_vel_disp.txt', delimiter = ',')
X_a_0_bquarter_c5_Radial_vel_disp = np.loadtxt('a_0_bquarter_c5_Radial_vel_disp.txt', delimiter = ',')
X_a_0_bquarter_c10_Radial_vel_disp = np.loadtxt('a_0_bquarter_c10_Radial_vel_disp.txt', delimiter = ',')
X_a_05_bquarter_c5_Radial_vel_disp = np.loadtxt('a_05_bquarter_c5_Radial_vel_disp.txt', delimiter = ',')
X_a_05_bquarter_c10_Radial_vel_disp = np.loadtxt('a_05_bquarter_c10_Radial_vel_disp.txt', delimiter = ',')
X_a_15_bquarter_c5_Radial_vel_disp = np.loadtxt('a_15_bquarter_c5_Radial_vel_disp.txt', delimiter = ',')
X_a_15_bquarter_c10_Radial_vel_disp = np.loadtxt('a_15_bquarter_c10_Radial_vel_disp.txt', delimiter = ',')
X_a_0_bhalf_c5_Radial_vel_disp = np.loadtxt('a_0_bhalf_c5_Radial_vel_disp.txt', delimiter = ',')
X_a_0_bhalf_c10_Radial_vel_disp = np.loadtxt('a_0_bhalf_c10_Radial_vel_disp.txt', delimiter = ',')
X_a_05_bhalf_c5_Radial_vel_disp = np.loadtxt('a_05_bhalf_c5_Radial_vel_disp.txt', delimiter = ',')
X_a_05_bhalf_c10_Radial_vel_disp = np.loadtxt('a_05_bhalf_c10_Radial_vel_disp.txt', delimiter = ',')
X_a_15_bhalf_c5_Radial_vel_disp = np.loadtxt('a_15_bhalf_c5_Radial_vel_disp.txt', delimiter = ',')
X_a_15_bhalf_c10_Radial_vel_disp = np.loadtxt('a_15_bhalf_c10_Radial_vel_disp.txt', delimiter = ',')
X_a_0_b0_c5_LOS_vel_disp = np.loadtxt('a_0_b0_c5_LOS_vel_disp.txt', delimiter = ',')
X_a_0_b0_c10_LOS_vel_disp = np.loadtxt('a_0_b0_c10_LOS_vel_disp.txt', delimiter = ',')
X_a_05_b0_c5_LOS_vel_disp = np.loadtxt('a_05_b0_c5_LOS_vel_disp.txt', delimiter = ',')
X_a_05_b0_c10_LOS_vel_disp = np.loadtxt('a_05_b0_c10_LOS_vel_disp.txt', delimiter = ',')
X_a_15_b0_c5_LOS_vel_disp = np.loadtxt('a_15_b0_c5_LOS_vel_disp.txt', delimiter = ',')
X_a_15_b0_c10_LOS_vel_disp = np.loadtxt('a_15_b0_c10_LOS_vel_disp.txt', delimiter = ',')
X_a_0_bquarter_c5_LOS_vel_disp = np.loadtxt('a_0_bquarter_c5_LOS_vel_disp.txt', delimiter = ',')
X_a_0_bquarter_c10_LOS_vel_disp = np.loadtxt('a_0_bquarter_c10_LOS_vel_disp.txt', delimiter = ',')
X_a_05_bquarter_c5_LOS_vel_disp = np.loadtxt('a_05_bquarter_c5_LOS_vel_disp.txt', delimiter = ',')
X_a_05_bquarter_c10_LOS_vel_disp = np.loadtxt('a_05_bquarter_c10_LOS_vel_disp.txt', delimiter = ',')
X_a_15_bquarter_c5_LOS_vel_disp = np.loadtxt('a_15_bquarter_c5_LOS_vel_disp.txt', delimiter = ',')
X_a_15_bquarter_c10_LOS_vel_disp = np.loadtxt('a_15_bquarter_c10_LOS_vel_disp.txt', delimiter = ',')
X_a_0_bhalf_c5_LOS_vel_disp = np.loadtxt('a_0_bhalf_c5_LOS_vel_disp.txt', delimiter = ',')
X_a_0_bhalf_c10_LOS_vel_disp = np.loadtxt('a_0_bhalf_c10_LOS_vel_disp.txt', delimiter = ',')
X_a_05_bhalf_c5_LOS_vel_disp = np.loadtxt('a_05_bhalf_c5_LOS_vel_disp.txt', delimiter = ',')
X_a_05_bhalf_c10_LOS_vel_disp = np.loadtxt('a_05_bhalf_c10_LOS_vel_disp.txt', delimiter = ',')
X_a_15_bhalf_c5_LOS_vel_disp = np.loadtxt('a_15_bhalf_c5_LOS_vel_disp.txt', delimiter = ',')
X_a_15_bhalf_c10_LOS_vel_disp = np.loadtxt('a_15_bhalf_c10_LOS_vel_disp.txt', delimiter = ',')
X_a_0_b0_c5_ap_vel_disp = np.loadtxt('a_0_b0_c5_ap_vel_disp.txt', delimiter = ',')
X_a_0_b0_c10_ap_vel_disp = np.loadtxt('a_0_b0_c10_ap_vel_disp.txt', delimiter = ',')
X_a_05_b0_c5_ap_vel_disp = np.loadtxt('a_05_b0_c5_ap_vel_disp.txt', delimiter = ',')
X_a_05_b0_c10_ap_vel_disp = np.loadtxt('a_05_b0_c10_ap_vel_disp.txt', delimiter = ',')
X_a_15_b0_c5_ap_vel_disp = np.loadtxt('a_15_b0_c5_ap_vel_disp.txt', delimiter = ',')
X_a_15_b0_c10_ap_vel_disp = np.loadtxt('a_15_b0_c10_ap_vel_disp.txt', delimiter = ',')
X_a_0_bquarter_c5_ap_vel_disp = np.loadtxt('a_0_bquarter_c5_ap_vel_disp.txt', delimiter = ',')
X_a_0_bquarter_c10_ap_vel_disp = np.loadtxt('a_0_bquarter_c10_ap_vel_disp.txt', delimiter = ',')
X_a_05_bquarter_c5_ap_vel_disp = np.loadtxt('a_05_bquarter_c5_ap_vel_disp.txt', delimiter = ',')
X_a_05_bquarter_c10_ap_vel_disp = np.loadtxt('a_05_bquarter_c10_ap_vel_disp.txt', delimiter = ',')
X_a_15_bquarter_c5_ap_vel_disp = np.loadtxt('a_15_bquarter_c5_ap_vel_disp.txt', delimiter = ',')
X_a_15_bquarter_c10_ap_vel_disp = np.loadtxt('a_15_bquarter_c10_ap_vel_disp.txt', delimiter = ',')
X_a_0_bhalf_c5_ap_vel_disp = np.loadtxt('a_0_bhalf_c5_ap_vel_disp.txt', delimiter = ',')
X_a_0_bhalf_c10_ap_vel_disp = np.loadtxt('a_0_bhalf_c10_ap_vel_disp.txt', delimiter = ',')
X_a_05_bhalf_c5_ap_vel_disp = np.loadtxt('a_05_bhalf_c5_ap_vel_disp.txt', delimiter = ',')
X_a_05_bhalf_c10_ap_vel_disp = np.loadtxt('a_05_bhalf_c10_ap_vel_disp.txt', delimiter = ',')
X_a_15_bhalf_c5_ap_vel_disp = np.loadtxt('a_15_bhalf_c5_ap_vel_disp.txt', delimiter = ',')
X_a_15_bhalf_c10_ap_vel_disp = np.loadtxt('a_15_bhalf_c10_ap_vel_disp.txt', delimiter = ',')
X_xi_b0_c10_rv_01 = np.loadtxt('xi_b0_c10_rv_01.txt', delimiter = ',')
X_xi_b0_c5_rv_01 = np.loadtxt('xi_b0_c5_rv_01.txt', delimiter = ',')
X_xi_b0_c10_rv_05 = np.loadtxt('xi_b0_c10_rv_05.txt', delimiter = ',')
X_xi_b0_c5_rv_05 = np.loadtxt('xi_b0_c5_rv_05.txt', delimiter = ',')
X_xi_bquarter_c10_rv_01 = np.loadtxt('xi_bquarter_c10_rv_01.txt', delimiter = ',')
X_xi_bquarter_c5_rv_01 = np.loadtxt('xi_bquarter_c5_rv_01.txt', delimiter = ',')
X_xi_bquarter_c10_rv_05 = np.loadtxt('xi_bquarter_c10_rv_05.txt', delimiter = ',')
X_xi_bquarter_c5_rv_05 = np.loadtxt('xi_bquarter_c5_rv_05.txt', delimiter = ',')
X_xi_bhalf_c10_rv_01 = np.loadtxt('xi_bhalf_c10_rv_01.txt', delimiter = ',')
X_xi_bhalf_c5_rv_01 = np.loadtxt('xi_bhalf_c5_rv_01.txt', delimiter = ',')
X_xi_bhalf_c10_rv_05 = np.loadtxt('xi_bhalf_c10_rv_05.txt', delimiter = ',')
X_xi_bhalf_c5_rv_05 = np.loadtxt('xi_bhalf_c5_rv_05.txt', delimiter = ',')
X_xi_NFW_c10_rv_01 = np.loadtxt('xi_NFW_c10_rv_01.txt', delimiter = ',')
X_xi_NFW_c5_rv_01 = np.loadtxt('xi_NFW_c5_rv_01.txt', delimiter = ',')
X_xi_NFW_c10_rv_05 = np.loadtxt('xi_NFW_c10_rv_05.txt', delimiter = ',')
X_xi_NFW_c5_rv_05 = np.loadtxt('xi_NFW_c5_rv_05.txt', delimiter = ',')










#NFW Density Profiles

# plt.plot(s, Density(1, 5, s), label='c=5', color='teal')
# plt.plot(s, Density(1, 10, s), label='c=10', color='blue')
# plt.plot(s, Density(1, 20, s), label='c=20', color='blueviolet')
# plt.plot(s, Density(1, 50, s), label='c=50', color='mediumvioletred')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10)
# plt.ylabel(r'$\rho/ v \rho_{c0}$')
# plt.yscale('log')
# plt.ylim(10**(-2), 10**5)
# plt.legend()
# plt.savefig('NFW_density.png', dpi=350)
# plt.show()










#Dehnen Density Profiles

# plt.figure(figsize=[10, 5])
#
#
# plt.subplot(1, 2, 1)
# plt.plot(s, General_Density(1, 4, 5, s), label='\u03b3 = 1', color='teal')
# plt.plot(s, General_Density(2, 4, 5, s), label='\u03b3 = 2', color='blue')
# plt.plot(s, Density(1, 5, s), label='NFW', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10)
# plt.ylabel(r'$\rho/ v \rho_{c0}$')
# plt.yscale('log')
# plt.ylim(10**(-2), 10**5)
# plt.legend(title='c=5')
#
#
# plt.subplot(1, 2, 2)
# plt.plot(s, General_Density(1, 4, 10, s), label='\u03b3 = 1', color='teal')
# plt.plot(s, General_Density(2, 4, 10, s), label='\u03b3 = 2', color='blue')
# plt.plot(s, Density(1, 10, s), label='NFW', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10)
# plt.ylabel(r'$\rho/ v \rho_{c0}$')
# plt.yscale('log')
# plt.ylim(10**(-2), 10**5)
# plt.legend(title='c=10')
#
# plt.savefig('Dehnen_density.png', dpi=350)
# plt.show()










#Alpha Density Profiles

# plt.figure(figsize=[10, 5])
#
#
# plt.subplot(1, 2, 1)
# plt.plot(s, Density(0, 5, s), label='\u03b1=0', color='teal')
# plt.plot(s, Density(0.5, 5, s), label='\u03b1=0.5', color='blue')
# plt.plot(s, Density(1, 5, s), label='\u03b1=1', color='blueviolet')
# plt.plot(s, Density(1.5, 5, s), label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10)
# plt.ylabel(r'$\rho/ v \rho_{c0}$')
# plt.yscale('log')
# plt.ylim(10**(-2), 10**5)
# plt.legend(title='c=5')
#
#
#
# plt.subplot(1, 2, 2)
# plt.plot(s, Density(0, 10, s), label='\u03b1=0', color='teal')
# plt.plot(s, Density(0.5, 10, s), label='\u03b1=0.5', color='blue')
# plt.plot(s, Density(1, 10, s), label='\u03b1=1', color='blueviolet')
# plt.plot(s, Density(1.5, 10, s), label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10)
# plt.ylabel(r'$\rho/ v \rho_{c0}$')
# plt.yscale('log')
# plt.ylim(10**(-2), 10**5)
# plt.legend(title='c=10')
#
#
# plt.savefig('Alpha_density.png', dpi=350)
# plt.show()










#NFW Radial Orbits

"""
# for i in range(len(s)):
#     NFW_b0_c5_Radial_vel_disp.append(Radial_vel_disp(1, 0, 5, s[i]))
#     NFW_b0_c10_Radial_vel_disp.append(Radial_vel_disp(1, 0, 10, s[i]))
#     NFW_b0_c20_Radial_vel_disp.append(Radial_vel_disp(1, 0, 20, s[i]))
#     print(i)
#
#
# for i in range(len(s)):
#     NFW_bquarter_c10_Radial_vel_disp.append(Radial_vel_disp(1, 0.25, 10, s[i]))
#     NFW_bhalf_c10_Radial_vel_disp.append(Radial_vel_disp(1, 0.5, 10, s[i]))
#     NFW_bquarter_c5_Radial_vel_disp.append(Radial_vel_disp(1, 0.25, 5, s[i]))
#     NFW_bhalf_c5_Radial_vel_disp.append(Radial_vel_disp(1, 0.5, 5, s[i]))
#     print(i)
#
#
# plt.figure(figsize=[10, 5])
#
#
# plt.subplot(1, 2, 1)
# plt.plot(s, NFW_b0_c5_Radial_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(s, NFW_bquarter_c5_Radial_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(s, NFW_bhalf_c5_Radial_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title = 'c=5')
#
#
# plt.subplot(1, 2, 2)
# plt.plot(s, NFW_b0_c10_Radial_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(s, NFW_bquarter_c10_Radial_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(s, NFW_bhalf_c10_Radial_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title = 'c=10')
#
#
# plt.savefig('NFW_radial_disp.png', dpi=350)
# plt.show()
"""


# plt.figure(figsize=[10, 5])
#
#
# plt.subplot(1, 2, 1)
# plt.plot(s150, X_NFW_b0_c5_Radial_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(s150, X_NFW_bquarter_c5_Radial_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(s150, X_NFW_bhalf_c5_Radial_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title = 'c=5')
#
#
# plt.subplot(1, 2, 2)
# plt.plot(s150, X_NFW_b0_c10_Radial_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(s150, X_NFW_bquarter_c10_Radial_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(s150, X_NFW_bhalf_c10_Radial_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title = 'c=10')
#
#
# plt.savefig('NFW_radial_disp.png', dpi=350)
# plt.show()










#NFW LOS Orbits

"""
# for i in range(len(S)):
#     NFW_b0_c5_LOS_vel_disp.append(Iso_LOS_vel_disp(1, 5, S[i]))
#     NFW_b0_c10_LOS_vel_disp.append(Iso_LOS_vel_disp(1, 10, S[i]))
#     NFW_b0_c20_LOS_vel_disp.append(Iso_LOS_vel_disp(1, 20, S[i]))
#     print(i)
#
#
# for i in range(len(S)):
#     NFW_bquarter_c10_LOS_vel_disp.append(LOS_vel_disp(1, 0.25, 10, S[i]))
#     NFW_bhalf_c10_LOS_vel_disp.append(LOS_vel_disp(1, 0.5, 10, S[i]))
#     NFW_bquarter_c5_LOS_vel_disp.append(LOS_vel_disp(1, 0.25, 5, S[i]))
#     NFW_bhalf_c5_LOS_vel_disp.append(LOS_vel_disp(1, 0.5, 5, S[i]))
#     print(i)
#
#
# plt.figure(figsize=[10, 5])
#
#
# plt.subplot(1, 2, 1)
# plt.plot(S, NFW_b0_c5_LOS_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(S, NFW_bquarter_c5_LOS_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(S, NFW_bhalf_c5_LOS_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title = 'c=5')
#
#
# plt.subplot(1, 2, 2)
# plt.plot(S, NFW_b0_c10_LOS_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(S, NFW_bquarter_c10_LOS_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(S, NFW_bhalf_c10_LOS_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title='c=10')
#
#
# plt.savefig('NFW_los_disp.png', dpi=350)
# plt.show()
"""


# plt.figure(figsize=[10, 5])
#
#
# plt.subplot(1, 2, 1)
# plt.plot(S50, X_NFW_b0_c5_LOS_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(S50, X_NFW_bquarter_c5_LOS_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(S50, X_NFW_bhalf_c5_LOS_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title = 'c=5')
#
#
# plt.subplot(1, 2, 2)
# plt.plot(S50, X_NFW_b0_c10_LOS_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(S50, X_NFW_bquarter_c10_LOS_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(S50, X_NFW_bhalf_c10_LOS_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title='c=10')
#
#
# plt.savefig('NFW_los_disp.png', dpi=350)
# plt.show()










#NFW Aperture Orbits

"""
# for i in range(len(S)):
#     NFW_b0_c5_ap_vel_disp.append(Iso_Aperture_vel_disp(1, 5, S[i]))
#     NFW_b0_c10_ap_vel_disp.append(Iso_Aperture_vel_disp(1, 10, S[i]))
#     NFW_b0_c20_ap_vel_disp.append(Iso_Aperture_vel_disp(1, 20, S[i]))
#     print(i)
#
#
# for i in range(len(S)):
#     NFW_bquarter_c10_ap_vel_disp.append(Aperture_vel_disp(1, 0.25, 10, S[i]))
#     NFW_bhalf_c10_ap_vel_disp.append(Aperture_vel_disp(1, 0.5, 10, S[i]))
#     NFW_bquarter_c5_ap_vel_disp.append(Aperture_vel_disp(1, 0.25, 5, S[i]))
#     NFW_bhalf_c5_ap_vel_disp.append(Aperture_vel_disp(1, 0.5, 5, S[i]))
#     print(i)
#
#
# plt.figure(figsize=[10, 5])
#
#
# plt.subplot(1, 2, 1)
# plt.plot(S, NFW_b0_c5_ap_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(S, NFW_bquarter_c5_ap_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(S, NFW_bhalf_c5_ap_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title ='c=5')
#
#
# plt.subplot(1, 2, 2)
# plt.plot(S, NFW_b0_c10_ap_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(S, NFW_bquarter_c10_ap_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(S, NFW_bhalf_c10_ap_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title ='c=10')
#
#
# plt.savefig('NFW_ap_disp.png', dpi=350)
# plt.show()
"""


# plt.figure(figsize=[10, 5])
#
#
# plt.subplot(1, 2, 1)
# plt.plot(S25, X_NFW_b0_c5_ap_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(S25, X_NFW_bquarter_c5_ap_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(S25, X_NFW_bhalf_c5_ap_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title ='c=5')
#
#
# plt.subplot(1, 2, 2)
# plt.plot(S25, X_NFW_b0_c10_ap_vel_disp, label='\u03b2=0', color='teal')
# plt.plot(S25, X_NFW_bquarter_c10_ap_vel_disp, label='\u03b2=0.25', color='blue')
# plt.plot(S25, X_NFW_bhalf_c10_ap_vel_disp, label='\u03b2=0.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.2)
# plt.legend(title ='c=10')
#
#
# plt.savefig('NFW_ap_disp.png', dpi=350)
# plt.show()










#Alpha Radial Orbits

"""
# for i in range(len(s)):
#     a_0_b0_c10_Radial_vel_disp.append(Radial_vel_disp(0, 0, 10, s[i]))
#     a_05_b0_c10_Radial_vel_disp.append(Radial_vel_disp(0.5, 0, 10, s[i]))
#     NFW_b0_c10_Radial_vel_disp.append(Radial_vel_disp(1, 0, 10, s[i]))
#     a_15_b0_c10_Radial_vel_disp.append(Radial_vel_disp(1.5, 0, 10, s[i]))
#     print(i)
#
#
# for i in range(len(s)):
#     a_0_bquarter_c10_Radial_vel_disp.append(Radial_vel_disp(0, 0.25, 10, s[i]))
#     a_05_bquarter_c10_Radial_vel_disp.append(Radial_vel_disp(0.5, 0.25, 10, s[i]))
#     NFW_bquarter_c10_Radial_vel_disp.append(Radial_vel_disp(1, 0.25, 10, s[i]))
#     a_15_bquarter_c10_Radial_vel_disp.append(Radial_vel_disp(1.5, 0.25, 10, s[i]))
#     print(i)
#
#
# for i in range(len(s)):
#     a_0_bhalf_c10_Radial_vel_disp.append(Radial_vel_disp(0, 0.5, 10, s[i]))
#     a_05_bhalf_c10_Radial_vel_disp.append(Radial_vel_disp(0.5, 0.5, 10, s[i]))
#     NFW_bhalf_c10_Radial_vel_disp.append(Radial_vel_disp(1, 0.5, 10, s[i]))
#     a_15_bhalf_c10_Radial_vel_disp.append(Radial_vel_disp(1.5, 0.5, 10, s[i]))
#     print(i)
#
#
# for i in range(len(s)):
#     a_0_b0_c5_Radial_vel_disp.append(Radial_vel_disp(0, 0, 5, s[i]))
#     a_05_b0_c5_Radial_vel_disp.append(Radial_vel_disp(0.5, 0, 5, s[i]))
#     NFW_b0_c5_Radial_vel_disp.append(Radial_vel_disp(1, 0, 5, s[i]))
#     a_15_b0_c5_Radial_vel_disp.append(Radial_vel_disp(1.5, 0, 5, s[i]))
#     print(i)
#
#
# for i in range(len(s)):
#     a_0_bquarter_c5_Radial_vel_disp.append(Radial_vel_disp(0, 0.25, 5, s[i]))
#     a_05_bquarter_c5_Radial_vel_disp.append(Radial_vel_disp(0.5, 0.25, 5, s[i]))
#     NFW_bquarter_c5_Radial_vel_disp.append(Radial_vel_disp(1, 0.25, 5, s[i]))
#     a_15_bquarter_c5_Radial_vel_disp.append(Radial_vel_disp(1.5, 0.25, 5, s[i]))
#     print(i)
#
#
# for i in range(len(s)):
#     a_0_bhalf_c5_Radial_vel_disp.append(Radial_vel_disp(0, 0.5, 5, s[i]))
#     a_05_bhalf_c5_Radial_vel_disp.append(Radial_vel_disp(0.5, 0.5, 5, s[i]))
#     NFW_bhalf_c5_Radial_vel_disp.append(Radial_vel_disp(1, 0.5, 5, s[i]))
#     a_15_bhalf_c5_Radial_vel_disp.append(Radial_vel_disp(1.5, 0.5, 5, s[i]))
#     print(i)
#
#
# plt.figure(figsize=[15, 10])
#
#
# plt.subplot(2, 3, 1)
# plt.plot(s, a_0_b0_c5_Radial_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(s, a_05_b0_c5_Radial_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(s, NFW_b0_c5_Radial_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(s, a_15_b0_c5_Radial_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 2)
# plt.legend(title='c=5, \u03b2=0')
#
#
# plt.subplot(2, 3, 2)
# plt.plot(s, a_0_bquarter_c5_Radial_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(s, a_05_bquarter_c5_Radial_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(s, NFW_bquarter_c5_Radial_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(s, a_15_bquarter_c5_Radial_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 2)
# plt.legend(title='c=5, \u03b2=0.25')
#
#
# plt.subplot(2, 3, 3)
# plt.plot(s, a_0_bhalf_c5_Radial_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(s, a_05_bhalf_c5_Radial_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(s, NFW_bhalf_c5_Radial_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(s, a_15_bhalf_c5_Radial_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 2)
# plt.legend(title='c=5, \u03b2=0.5')
#
#
# plt.subplot(2, 3, 4)
# plt.plot(s, a_0_b0_c10_Radial_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(s, a_05_b0_c10_Radial_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(s, NFW_b0_c10_Radial_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(s, a_15_b0_c10_Radial_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 2)
# plt.legend(title='c=10, \u03b2=0')
#
#
# plt.subplot(2, 3, 5)
# plt.plot(s, a_0_bquarter_c10_Radial_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(s, a_05_bquarter_c10_Radial_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(s, NFW_bquarter_c10_Radial_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(s, a_15_bquarter_c10_Radial_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 2)
# plt.legend(title='c=10, \u03b2=0.25')
#
#
# plt.subplot(2, 3, 6)
# plt.plot(s, a_0_bhalf_c10_Radial_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(s, a_05_bhalf_c10_Radial_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(s, NFW_bhalf_c10_Radial_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(s, a_15_bhalf_c10_Radial_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 2)
# plt.legend(title='c=10, \u03b2=0.5')
#
# plt.savefig('Alpha_radial_disp.png', dpi=350)
# plt.show()
"""

"""
# with open("a_0_b0_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_0_b0_c10_Radial_vel_disp))
# with open("a_05_b0_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_05_b0_c10_Radial_vel_disp))
# with open("NFW_b0_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(NFW_b0_c10_Radial_vel_disp))
# with open("a_15_b0_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_15_b0_c10_Radial_vel_disp))
# with open("a_0_b0_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_0_b0_c5_Radial_vel_disp))
# with open("a_05_b0_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_05_b0_c5_Radial_vel_disp))
# with open("NFW_b0_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(NFW_b0_c5_Radial_vel_disp))
# with open("a_15_b0_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_15_b0_c5_Radial_vel_disp))
#
#
# with open("a_0_bquarter_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_0_bquarter_c10_Radial_vel_disp))
# with open("a_05_bquarter_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_05_bquarter_c10_Radial_vel_disp))
# with open("NFW_bquarter_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(NFW_bquarter_c10_Radial_vel_disp))
# with open("a_15_bquarter_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_15_bquarter_c10_Radial_vel_disp))
# with open("a_0_bhalf_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_0_bhalf_c10_Radial_vel_disp))
# with open("a_05_bhalf_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_05_bhalf_c10_Radial_vel_disp))
# with open("NFW_bhalf_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(NFW_bhalf_c10_Radial_vel_disp))
# with open("a_15_bhalf_c10_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_15_bhalf_c10_Radial_vel_disp))
#
#
# with open("a_0_bquarter_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_0_bquarter_c5_Radial_vel_disp))
# with open("a_05_bquarter_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_05_bquarter_c5_Radial_vel_disp))
# with open("NFW_bquarter_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(NFW_bquarter_c5_Radial_vel_disp))
# with open("a_15_bquarter_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_15_bquarter_c5_Radial_vel_disp))
# with open("a_0_bhalf_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_0_bhalf_c5_Radial_vel_disp))
# with open("a_05_bhalf_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_05_bhalf_c5_Radial_vel_disp))
# with open("NFW_bhalf_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(NFW_bhalf_c5_Radial_vel_disp))
# with open("a_15_bhalf_c5_Radial_vel_disp.txt", "w") as output:
#     output.write(str(a_15_bhalf_c5_Radial_vel_disp))
"""


# plt.figure(figsize=[15, 10])
#
#
# plt.subplot(2, 3, 1)
# plt.plot(s150, X_a_0_b0_c5_Radial_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(s150, X_a_05_b0_c5_Radial_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(s150, X_NFW_b0_c5_Radial_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(s150, X_a_15_b0_c5_Radial_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0')
#
#
# plt.subplot(2, 3, 2)
# plt.plot(s150, X_a_0_bquarter_c5_Radial_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(s150, X_a_05_bquarter_c5_Radial_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(s150, X_NFW_bquarter_c5_Radial_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(s150, X_a_15_bquarter_c5_Radial_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0.25')
#
#
# plt.subplot(2, 3, 3)
# plt.plot(s150, X_a_0_bhalf_c5_Radial_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(s150, X_a_05_bhalf_c5_Radial_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(s150, X_NFW_bhalf_c5_Radial_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(s150, X_a_15_bhalf_c5_Radial_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0.5')
#
#
# plt.subplot(2, 3, 4)
# plt.plot(s150, X_a_0_b0_c10_Radial_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(s150, X_a_05_b0_c10_Radial_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(s150, X_NFW_b0_c10_Radial_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(s150, X_a_15_b0_c10_Radial_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0')
#
#
# plt.subplot(2, 3, 5)
# plt.plot(s150, X_a_0_bquarter_c10_Radial_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(s150, X_a_05_bquarter_c10_Radial_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(s150, X_NFW_bquarter_c10_Radial_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(s150, X_a_15_bquarter_c10_Radial_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0.25')
#
#
# plt.subplot(2, 3, 6)
# plt.plot(s150, X_a_0_bhalf_c10_Radial_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(s150, X_a_05_bhalf_c10_Radial_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(s150, X_NFW_bhalf_c10_Radial_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(s150, X_a_15_bhalf_c10_Radial_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$r/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_r/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0.5')
#
# plt.savefig('Alpha_radial_disp.png', dpi=350)
# plt.show()










#Alpha LOS Orbits

"""
# for i in range(len(S)):
#     a_0_b0_c10_LOS_vel_disp.append(Iso_LOS_vel_disp(0, 10, S[i]))
#     a_05_b0_c10_LOS_vel_disp.append(Iso_LOS_vel_disp(0.5, 10, S[i]))
#     NFW_b0_c10_LOS_vel_disp.append(Iso_LOS_vel_disp(1, 10, S[i]))
#     a_15_b0_c10_LOS_vel_disp.append(Iso_LOS_vel_disp(1.5, 10, S[i]))
#     print(i)
#
#
# for i in range(len(S)):
#     a_0_bquarter_c10_LOS_vel_disp.append(LOS_vel_disp(0, 0.25, 10, S[i]))
#     a_05_bquarter_c10_LOS_vel_disp.append(LOS_vel_disp(0.5, 0.25, 10, S[i]))
#     NFW_bquarter_c10_LOS_vel_disp.append(LOS_vel_disp(1, 0.25, 10, S[i]))
#     a_15_bquarter_c10_LOS_vel_disp.append(LOS_vel_disp(1.5, 0.25, 10, S[i]))
#     print(i)
#
#
# for i in range(len(S)):
#     a_0_bhalf_c10_LOS_vel_disp.append(LOS_vel_disp(0, 0.5, 10, S[i]))
#     a_05_bhalf_c10_LOS_vel_disp.append(LOS_vel_disp(0.5, 0.5, 10, S[i]))
#     NFW_bhalf_c10_LOS_vel_disp.append(LOS_vel_disp(1, 0.5, 10, S[i]))
#     a_15_bhalf_c10_LOS_vel_disp.append(LOS_vel_disp(1.5, 0.5, 10, S[i]))
#     print(i)
#
#
# for i in range(len(S)):
#     a_0_b0_c5_LOS_vel_disp.append(Iso_LOS_vel_disp(0, 5, S[i]))
#     a_05_b0_c5_LOS_vel_disp.append(Iso_LOS_vel_disp(0.5, 5, S[i]))
#     NFW_b0_c5_LOS_vel_disp.append(Iso_LOS_vel_disp(1, 5, S[i]))
#     a_15_b0_c5_LOS_vel_disp.append(Iso_LOS_vel_disp(1.5, 5, S[i]))
#     print(i)
#
#
# for i in range(len(S)):
#     a_0_bquarter_c5_LOS_vel_disp.append(LOS_vel_disp(0, 0.25, 5, S[i]))
#     a_05_bquarter_c5_LOS_vel_disp.append(LOS_vel_disp(0.5, 0.25, 5, S[i]))
#     NFW_bquarter_c5_LOS_vel_disp.append(LOS_vel_disp(1, 0.25, 5, S[i]))
#     a_15_bquarter_c5_LOS_vel_disp.append(LOS_vel_disp(1.5, 0.25, 5, S[i]))
#     print(i)
#
#
# for i in range(len(S)):
#     a_0_bhalf_c5_LOS_vel_disp.append(LOS_vel_disp(0, 0.5, 5, S[i]))
#     a_05_bhalf_c5_LOS_vel_disp.append(LOS_vel_disp(0.5, 0.5, 5, S[i]))
#     NFW_bhalf_c5_LOS_vel_disp.append(LOS_vel_disp(1, 0.5, 5, S[i]))
#     a_15_bhalf_c5_LOS_vel_disp.append(LOS_vel_disp(1.5, 0.5, 5, S[i]))
#     print(i)
#
#
# plt.figure(figsize=[15, 10])
#
#
# plt.subplot(2, 3, 1)
# plt.plot(S, a_0_b0_c5_LOS_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_b0_c5_LOS_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_b0_c5_LOS_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_b0_c5_LOS_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.6)
# plt.legend(title='c=5, \u03b2=0')
#
#
# plt.subplot(2, 3, 2)
# plt.plot(S, a_0_bquarter_c5_LOS_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_bquarter_c5_LOS_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_bquarter_c5_LOS_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_bquarter_c5_LOS_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.6)
# plt.legend(title='c=5, \u03b2=0.25')
#
#
# plt.subplot(2, 3, 3)
# plt.plot(S, a_0_bhalf_c5_LOS_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_bhalf_c5_LOS_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_bhalf_c5_LOS_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_bhalf_c5_LOS_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.6)
# plt.legend(title='c=5, \u03b2=0.5')
#
#
# plt.subplot(2, 3, 4)
# plt.plot(S, a_0_b0_c10_LOS_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_b0_c10_LOS_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_b0_c10_LOS_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_b0_c10_LOS_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.6)
# plt.legend(title='c=10, \u03b2=0')
#
#
# plt.subplot(2, 3, 5)
# plt.plot(S, a_0_bquarter_c10_LOS_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_bquarter_c10_LOS_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_bquarter_c10_LOS_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_bquarter_c10_LOS_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.6)
# plt.legend(title='c=10, \u03b2=0.25')
#
#
# plt.subplot(2, 3, 6)
# plt.plot(S, a_0_bhalf_c10_LOS_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_bhalf_c10_LOS_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_bhalf_c10_LOS_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_bhalf_c10_LOS_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.6)
# plt.legend(title='c=10, \u03b2=0.5')
#
#
# plt.savefig('Alpha_los_disp.png', dpi=350)
# plt.show()
"""

"""
# with open("a_0_b0_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_0_b0_c10_LOS_vel_disp))
# with open("a_05_b0_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_05_b0_c10_LOS_vel_disp))
# with open("NFW_b0_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(NFW_b0_c10_LOS_vel_disp))
# with open("a_15_b0_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_15_b0_c10_LOS_vel_disp))
# with open("a_0_b0_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_0_b0_c5_LOS_vel_disp))
# with open("a_05_b0_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_05_b0_c5_LOS_vel_disp))
# with open("NFW_b0_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(NFW_b0_c5_LOS_vel_disp))
# with open("a_15_b0_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_15_b0_c5_LOS_vel_disp))
#
#
# with open("a_0_bquarter_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_0_bquarter_c10_LOS_vel_disp))
# with open("a_05_bquarter_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_05_bquarter_c10_LOS_vel_disp))
# with open("NFW_bquarter_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(NFW_bquarter_c10_LOS_vel_disp))
# with open("a_15_bquarter_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_15_bquarter_c10_LOS_vel_disp))
# with open("a_0_bhalf_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_0_bhalf_c10_LOS_vel_disp))
# with open("a_05_bhalf_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_05_bhalf_c10_LOS_vel_disp))
# with open("NFW_bhalf_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(NFW_bhalf_c10_LOS_vel_disp))
# with open("a_15_bhalf_c10_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_15_bhalf_c10_LOS_vel_disp))
#
#
# with open("a_0_bquarter_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_0_bquarter_c5_LOS_vel_disp))
# with open("a_05_bquarter_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_05_bquarter_c5_LOS_vel_disp))
# with open("NFW_bquarter_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(NFW_bquarter_c5_LOS_vel_disp))
# with open("a_15_bquarter_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_15_bquarter_c5_LOS_vel_disp))
# with open("a_0_bhalf_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_0_bhalf_c5_LOS_vel_disp))
# with open("a_05_bhalf_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_05_bhalf_c5_LOS_vel_disp))
# with open("NFW_bhalf_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(NFW_bhalf_c5_LOS_vel_disp))
# with open("a_15_bhalf_c5_LOS_vel_disp.txt", "w") as output:
#     output.write(str(a_15_bhalf_c5_LOS_vel_disp))
"""


# plt.figure(figsize=[15, 10])
#
#
# plt.subplot(2, 3, 1)
# plt.plot(S50, X_a_0_b0_c5_LOS_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S50, X_a_05_b0_c5_LOS_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S50, X_NFW_b0_c5_LOS_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S50, X_a_15_b0_c5_LOS_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0')
#
#
# plt.subplot(2, 3, 2)
# plt.plot(S50, X_a_0_bquarter_c5_LOS_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S50, X_a_05_bquarter_c5_LOS_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S50, X_NFW_bquarter_c5_LOS_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S50, X_a_15_bquarter_c5_LOS_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0.25')
#
#
# plt.subplot(2, 3, 3)
# plt.plot(S50, X_a_0_bhalf_c5_LOS_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S50, X_a_05_bhalf_c5_LOS_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S50, X_NFW_bhalf_c5_LOS_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S50, X_a_15_bhalf_c5_LOS_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0.5')
#
#
# plt.subplot(2, 3, 4)
# plt.plot(S50, X_a_0_b0_c10_LOS_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S50, X_a_05_b0_c10_LOS_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S50, X_NFW_b0_c10_LOS_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S50, X_a_15_b0_c10_LOS_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0')
#
#
# plt.subplot(2, 3, 5)
# plt.plot(S50, X_a_0_bquarter_c10_LOS_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S50, X_a_05_bquarter_c10_LOS_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S50, X_NFW_bquarter_c10_LOS_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S50, X_a_15_bquarter_c10_LOS_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0.25')
#
#
# plt.subplot(2, 3, 6)
# plt.plot(S50, X_a_0_bhalf_c10_LOS_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S50, X_a_05_bhalf_c10_LOS_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S50, X_NFW_bhalf_c10_LOS_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S50, X_a_15_bhalf_c10_LOS_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_\ell/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0.5')
#
#
# plt.savefig('Alpha_los_disp.png', dpi=350)
# plt.show()










#Alpha Aperture Orbits

"""
# for i in range(len(S)):
    # a_0_b0_c10_ap_vel_disp.append(Iso_Aperture_vel_disp(0, 10, S[i]))
    # a_05_b0_c10_ap_vel_disp.append(Iso_Aperture_vel_disp(0.5, 10, S[i]))
    # NFW_b0_c10_ap_vel_disp.append(Iso_Aperture_vel_disp(1, 10, S[i]))
    # a_15_b0_c10_ap_vel_disp.append(Iso_Aperture_vel_disp(1.5, 10, S[i]))
    # print(i)


# for i in range(len(S)):
    # a_0_bquarter_c10_ap_vel_disp.append(Aperture_vel_disp(0, 0.25, 10, S[i]))
    # print(i)
    # a_05_bquarter_c10_ap_vel_disp.append(Aperture_vel_disp(0.5, 0.25, 10, S[i]))
    # print(i)
    # NFW_bquarter_c10_ap_vel_disp.append(Aperture_vel_disp(1, 0.25, 10, S[i]))
    # print(i)
    # a_15_bquarter_c10_ap_vel_disp.append(Aperture_vel_disp(1.5, 0.25, 10, S[i]))
    # print(i)


# for i in range(len(S)):
    # a_0_bhalf_c10_ap_vel_disp.append(Aperture_vel_disp(0, 0.5, 10, S[i]))
    # print(i)
    # a_05_bhalf_c10_ap_vel_disp.append(Aperture_vel_disp(0.5, 0.5, 10, S[i]))
    # print(i)
    # NFW_bhalf_c10_ap_vel_disp.append(Aperture_vel_disp(1, 0.5, 10, S[i]))
    # print(i)
    # a_15_bhalf_c10_ap_vel_disp.append(Aperture_vel_disp(1.5, 0.5, 10, S[i]))
    # print(i)


# for i in range(len(S)):
    # a_0_b0_c5_ap_vel_disp.append(Iso_Aperture_vel_disp(0, 5, S[i]))
    # a_05_b0_c5_ap_vel_disp.append(Iso_Aperture_vel_disp(0.5, 5, S[i]))
    # NFW_b0_c5_ap_vel_disp.append(Iso_Aperture_vel_disp(1, 5, S[i]))
    # a_15_b0_c5_ap_vel_disp.append(Iso_Aperture_vel_disp(1.5, 5, S[i]))


# for i in range(len(S)):
    # a_0_bquarter_c5_ap_vel_disp.append(Aperture_vel_disp(0, 0.25, 5, S[i]))
    # print(i)
    # a_05_bquarter_c5_ap_vel_disp.append(Aperture_vel_disp(0.5, 0.25, 5, S[i]))
    # print(i)
    # NFW_bquarter_c5_ap_vel_disp.append(Aperture_vel_disp(1, 0.25, 5, S[i]))
    # print(i)
    # a_15_bquarter_c5_ap_vel_disp.append(Aperture_vel_disp(1.5, 0.25, 5, S[i]))
    # print(i)


# for i in range(len(S)):
    # a_0_bhalf_c5_ap_vel_disp.append(Aperture_vel_disp(0, 0.5, 5, S[i]))
    # print(i)
    # a_05_bhalf_c5_ap_vel_disp.append(Aperture_vel_disp(0.5, 0.5, 5, S[i]))
    # print(i)
    # NFW_bhalf_c5_ap_vel_disp.append(Aperture_vel_disp(1, 0.5, 5, S[i]))
    # print(i)
    # a_15_bhalf_c5_ap_vel_disp.append(Aperture_vel_disp(1.5, 0.5, 5, S[i]))
    # print(i)
#    
#
# plt.figure(figsize=[15, 10])


# plt.subplot(2, 3, 1)
# plt.plot(S, a_0_b0_c5_ap_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_b0_c5_ap_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_b0_c5_ap_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_b0_c5_ap_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0')


# plt.subplot(2, 3, 2)
# plt.plot(S, a_0_bquarter_c5_ap_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_bquarter_c5_ap_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_bquarter_c5_ap_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_bquarter_c5_ap_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0.25')


# plt.subplot(2, 3, 3)
# plt.plot(S, a_0_bhalf_c5_ap_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_bhalf_c5_ap_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_bhalf_c5_ap_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_bhalf_c5_ap_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0.5')


# plt.subplot(2, 3, 4)
# plt.plot(S, a_0_b0_c10_ap_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_b0_c10_ap_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_b0_c10_ap_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_b0_c10_ap_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0')


# plt.subplot(2, 3, 5)
# plt.plot(S, a_0_bquarter_c10_ap_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_bquarter_c10_ap_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_bquarter_c10_ap_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_bquarter_c10_ap_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0.25')


# plt.subplot(2, 3, 6)
# plt.plot(S, a_0_bhalf_c10_ap_vel_disp, label='\u03b1=0', color='olive')
# plt.plot(S, a_05_bhalf_c10_ap_vel_disp, label='\u03b1=0.5', color='teal')
# plt.plot(S, NFW_bhalf_c10_ap_vel_disp, label='\u03b1=1', color='blue')
# plt.plot(S, a_15_bhalf_c10_ap_vel_disp, label='\u03b1=1.5', color='blueviolet')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0.5')
#
#
# plt.savefig('Alpha_2_ap_disp.png', dpi=350)
# plt.show()
"""

"""
with open("a_0_b0_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(a_0_b0_c10_ap_vel_disp))
with open("a_05_b0_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(a_05_b0_c10_ap_vel_disp))
with open("NFW_b0_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(NFW_b0_c10_ap_vel_disp))
with open("a_15_b0_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(a_15_b0_c10_ap_vel_disp))
with open("a_0_b0_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(a_0_b0_c5_ap_vel_disp))
with open("a_05_b0_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(a_05_b0_c5_ap_vel_disp))
with open("NFW_b0_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(NFW_b0_c5_ap_vel_disp))
with open("a_15_b0_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(a_15_b0_c5_ap_vel_disp))


with open("a_0_bquarter_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(a_0_bquarter_c10_ap_vel_disp))
with open("a_05_bquarter_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(a_05_bquarter_c10_ap_vel_disp))
with open("NFW_bquarter_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(NFW_bquarter_c10_ap_vel_disp))
with open("a_15_bquarter_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(a_15_bquarter_c10_ap_vel_disp))
with open("a_0_bhalf_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(a_0_bhalf_c10_ap_vel_disp))
with open("a_05_bhalf_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(a_05_bhalf_c10_ap_vel_disp))
with open("NFW_bhalf_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(NFW_bhalf_c10_ap_vel_disp))
with open("a_15_bhalf_c10_ap_vel_disp.txt", "w") as output:
    output.write(str(a_15_bhalf_c10_ap_vel_disp))


with open("a_0_bquarter_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(a_0_bquarter_c5_ap_vel_disp))
with open("a_05_bquarter_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(a_05_bquarter_c5_ap_vel_disp))
with open("NFW_bquarter_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(NFW_bquarter_c5_ap_vel_disp))
with open("a_15_bquarter_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(a_15_bquarter_c5_ap_vel_disp))
with open("a_0_bhalf_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(a_0_bhalf_c5_ap_vel_disp))
with open("a_05_bhalf_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(a_05_bhalf_c5_ap_vel_disp))
with open("NFW_bhalf_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(NFW_bhalf_c5_ap_vel_disp))
with open("a_15_bhalf_c5_ap_vel_disp.txt", "w") as output:
    output.write(str(a_15_bhalf_c5_ap_vel_disp))
"""


# plt.figure(figsize=[15, 10])
#
#
# plt.subplot(2, 3, 1)
# plt.plot(S25, X_a_0_b0_c5_ap_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S25, X_a_05_b0_c5_ap_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S25, X_NFW_b0_c5_ap_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S25, X_a_15_b0_c5_ap_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0')
#
#
# plt.subplot(2, 3, 2)
# plt.plot(S25, X_a_0_bquarter_c5_ap_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S25, X_a_05_bquarter_c5_ap_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S25, X_NFW_bquarter_c5_ap_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S25, X_a_15_bquarter_c5_ap_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0.25')
#
#
# plt.subplot(2, 3, 3)
# plt.plot(S25, X_a_0_bhalf_c5_ap_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S25, X_a_05_bhalf_c5_ap_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S25, X_NFW_bhalf_c5_ap_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S25, X_a_15_bhalf_c5_ap_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=5, \u03b2=0.5')
#
#
# plt.subplot(2, 3, 4)
# plt.plot(S25, X_a_0_b0_c10_ap_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S25, X_a_05_b0_c10_ap_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S25, X_NFW_b0_c10_ap_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S25, X_a_15_b0_c10_ap_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0')
#
#
# plt.subplot(2, 3, 5)
# plt.plot(S25, X_a_0_bquarter_c10_ap_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S25, X_a_05_bquarter_c10_ap_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S25, X_NFW_bquarter_c10_ap_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S25, X_a_15_bquarter_c10_ap_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0.25')
#
#
# plt.subplot(2, 3, 6)
# plt.plot(S25, X_a_0_bhalf_c10_ap_vel_disp, label='\u03b1=0', color='teal')
# plt.plot(S25, X_a_05_bhalf_c10_ap_vel_disp, label='\u03b1=0.5', color='blue')
# plt.plot(S25, X_NFW_bhalf_c10_ap_vel_disp, label='\u03b1=1', color='blueviolet')
# plt.plot(S25, X_a_15_bhalf_c10_ap_vel_disp, label='\u03b1=1.5', color='mediumvioletred')
# plt.xlabel(r'$R/r_v$')
# plt.xscale('log')
# plt.xlim(0.5*10**(-2), 0.5*10**1)
# plt.ylabel(r'$\sigma_a/ v_{cv}$')
# plt.ylim(0, 1.4)
# plt.legend(title='c=10, \u03b2=0.5')
#
#
# plt.savefig('Alpha_ap_disp.png', dpi=350)
# plt.show()










#Xi for NFW Orbits

"""
for i in range(len(b)):
    xi_NFW_c5_rv_01.append(Aperture_vel_disp(1, b[i], 5, 0.1))
    xi_NFW_c10_rv_01.append(Aperture_vel_disp(1, b[i], 10, 0.1))
    xi_NFW_c5_rv_05.append(Aperture_vel_disp(1, b[i], 5, 0.5))
    xi_NFW_c10_rv_05.append(Aperture_vel_disp(1, b[i], 10, 0.5))
    print(i)


plt.figure(figsize=[10, 5])


plt.subplot(1, 2, 1)
plt.plot(b, xi_NFW_c5_rv_01, label='c=5', color='teal')
plt.plot(b, xi_NFW_c10_rv_01, label='c=10', color='blue')
plt.xlim(0, 0.5)
plt.xlabel('\u03B2')
plt.ylim(0.5, 1)
plt.ylabel('$\sigma_a/ v_{cv}$')
plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
plt.legend(title='R=0.1$r_v$')


plt.subplot(1, 2, 2)
plt.plot(b, xi_NFW_c5_rv_05, label='c=5', color='teal')
plt.plot(b, xi_NFW_c10_rv_05, label='c=10', color='blue')
plt.xlim(0, 0.5)
plt.xlabel('\u03B2')
plt.ylim(0.5, 1)
plt.ylabel('$\sigma_a/ v_{cv}$')
plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
plt.legend(title='R=0.5$r_v$')


plt.savefig('NFW_xi.png', dpi=350)
plt.show()
"""

"""
with open("xi_NFW_c5_rv_01.txt", "w") as output:
    output.write(str(xi_NFW_c5_rv_01))
with open("xi_NFW_c10_rv_01.txt", "w") as output:
    output.write(str(xi_NFW_c10_rv_01))
with open("xi_NFW_c5_rv_05.txt", "w") as output:
    output.write(str(xi_NFW_c5_rv_05))
with open("xi_NFW_c10_rv_05.txt", "w") as output:
    output.write(str(xi_NFW_c10_rv_05))
"""


# plt.figure(figsize=[10, 3])
#
#
# plt.subplot(1, 2, 1)
# plt.plot(b25, X_xi_NFW_c5_rv_01, label='c=5', color='teal')
# plt.plot(b25, X_xi_NFW_c10_rv_01, label='c=10', color='blue')
# plt.xlim(0, 0.5)
# plt.xlabel('\u03B2')
# plt.ylim(0.6, 0.9)
# plt.ylabel('$\sigma_a/ v_{cv}$')
# plt.yticks([0.6, 0.7, 0.8, 0.9])
# plt.legend(title='R=0.1$r_v$')
#
#
# plt.subplot(1, 2, 2)
# plt.plot(b25, X_xi_NFW_c5_rv_05, label='c=5', color='teal')
# plt.plot(b25, X_xi_NFW_c10_rv_05, label='c=10', color='blue')
# plt.xlim(0, 0.5)
# plt.xlabel('\u03B2')
# plt.ylim(0.6, 0.9)
# plt.ylabel('$\sigma_a/ v_{cv}$')
# plt.yticks([0.6, 0.7, 0.8, 0.9])
# plt.legend(title='R=0.5$r_v$')
#
#
# plt.tight_layout()
# plt.savefig('NFW_xi.png', dpi=350)
# plt.show()










#Xi for Alpha Orbits

"""
for i in range(len(a)):
    xi_b0_c5_rv_01.append(Iso_Aperture_vel_disp(a[i], 5, 0.1))
    xi_b0_c10_rv_01.append(Iso_Aperture_vel_disp(a[i], 10, 0.1))
    xi_b0_c5_rv_05.append(Iso_Aperture_vel_disp(a[i], 5, 0.5))
    xi_b0_c10_rv_05.append(Iso_Aperture_vel_disp(a[i], 10, 0.5))
    print(i)


for i in range(len(a)):
    xi_bquarter_c5_rv_01.append(Aperture_vel_disp(a[i], 0.25, 5, 0.1))
    xi_bquarter_c10_rv_01.append(Aperture_vel_disp(a[i], 0.25, 10, 0.1))
    xi_bquarter_c5_rv_05.append(Aperture_vel_disp(a[i], 0.25, 5, 0.5))
    xi_bquarter_c10_rv_05.append(Aperture_vel_disp(a[i], 0.25, 10, 0.5))
    print(i)


for i in range(len(a)):
    xi_bhalf_c5_rv_01.append(Aperture_vel_disp(a[i], 0.5, 5, 0.1))
    xi_bhalf_c10_rv_01.append(Aperture_vel_disp(a[i], 0.5, 10, 0.1))
    xi_bhalf_c5_rv_05.append(Aperture_vel_disp(a[i], 0.5, 5, 0.5))
    xi_bhalf_c10_rv_05.append(Aperture_vel_disp(a[i], 0.5, 10, 0.5))
    print(i)


plt.figure(figsize=[15, 10])


plt.subplot(2, 3, 1)
plt.plot(a, xi_b0_c5_rv_01, label='c=5', color='teal')
plt.plot(a, xi_b0_c10_rv_01, label='c=10', color='blue')
plt.xlim(0, 1.5)
plt.xlabel('\u03B1')
plt.ylim(0.5, 1.2)
plt.ylabel('$\sigma_a/ v_{cv}$')
plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
plt.legend(title='R = 0.1$r_v$, \u03B2 = 0')


plt.subplot(2, 3, 2)
plt.plot(a, xi_bquarter_c5_rv_01, label='c=5', color='teal')
plt.plot(a, xi_bquarter_c10_rv_01, label='c=10', color='blue')
plt.xlim(0, 1.5)
plt.xlabel('\u03B1')
plt.ylim(0.5, 1.2)
plt.ylabel('$\sigma_a/ v_{cv}$')
plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
plt.legend(title='R = 0.1$r_v$, \u03B2 = 0.25')


plt.subplot(2, 3, 3)
plt.plot(a, xi_bhalf_c5_rv_01, label='c=5', color='teal')
plt.plot(a, xi_bhalf_c10_rv_01, label='c=10', color='blue')
plt.xlim(0, 1.5)
plt.xlabel('\u03B1')
plt.ylim(0.5, 1.2)
plt.ylabel('$\sigma_a/ v_{cv}$')
plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
plt.legend(title='R = 0.1$r_v$, \u03B2 = 0.5')


plt.subplot(2, 3, 4)
plt.plot(a, xi_b0_c5_rv_05, label='c=5', color='teal')
plt.plot(a, xi_b0_c10_rv_05, label='c=10', color='blue')
plt.xlim(0, 1.5)
plt.xlabel('\u03B1')
plt.ylim(0.5, 1.2)
plt.ylabel('$\sigma_a/ v_{cv}$')
plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
plt.legend(title='R = 0.5$r_v$, \u03B2 = 0')


plt.subplot(2, 3, 5)
plt.plot(a, xi_bquarter_c5_rv_05, label='c=5', color='teal')
plt.plot(a, xi_bquarter_c10_rv_05, label='c=10', color='blue')
plt.xlim(0, 1.5)
plt.xlabel('\u03B1')
plt.ylim(0.5, 1.2)
plt.ylabel('$\sigma_a/ v_{cv}$')
plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
plt.legend(title='R = 0.5$r_v$, \u03B2 = 0.25')


plt.subplot(2, 3, 6)
plt.plot(a, xi_bhalf_c5_rv_05, label='c=5', color='teal')
plt.plot(a, xi_bhalf_c10_rv_05, label='c=10', color='blue')
plt.xlim(0, 1.5)
plt.xlabel('\u03B1')
plt.ylim(0.5, 1.2)
plt.ylabel('$\sigma_a/ v_{cv}$')
plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
plt.legend(title='R = 0.5$r_v$, \u03B2 = 0.5')


plt.savefig('Alpha_xi.png', dpi=350)
plt.show()
"""

"""
with open("xi_b0_c5_rv_01.txt", "w") as output:
    output.write(str(xi_b0_c5_rv_01))
with open("xi_b0_c10_rv_01.txt", "w") as output:
    output.write(str(xi_b0_c10_rv_01))
with open("xi_b0_c5_rv_05.txt", "w") as output:
    output.write(str(xi_b0_c5_rv_05))
with open("xi_b0_c10_rv_05.txt", "w") as output:
    output.write(str(xi_b0_c10_rv_05))
with open("xi_bquarter_c5_rv_01.txt", "w") as output:
    output.write(str(xi_bquarter_c5_rv_01))
with open("xi_bquarter_c10_rv_01.txt", "w") as output:
    output.write(str(xi_bquarter_c10_rv_01))
with open("xi_bquarter_c5_rv_05.txt", "w") as output:
    output.write(str(xi_bquarter_c5_rv_05))
with open("xi_bquarter_c10_rv_05.txt", "w") as output:
    output.write(str(xi_bquarter_c10_rv_05))
with open("xi_bhalf_c5_rv_01.txt", "w") as output:
    output.write(str(xi_bhalf_c5_rv_01))
with open("xi_bhalf_c10_rv_01.txt", "w") as output:
    output.write(str(xi_bhalf_c10_rv_01))
with open("xi_bhalf_c5_rv_05.txt", "w") as output:
    output.write(str(xi_bhalf_c5_rv_05))
with open("xi_bhalf_c10_rv_05.txt", "w") as output:
    output.write(str(xi_bhalf_c10_rv_05))
"""


# plt.figure(figsize=[15, 6])
#
#
# plt.subplot(2, 3, 1)
# plt.plot(a25, X_xi_b0_c5_rv_01, label='c=5', color='teal')
# plt.plot(a25, X_xi_b0_c10_rv_01, label='c=10', color='blue')
# plt.xlim(0, 1.5)
# plt.xlabel('\u03B1')
# plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
# plt.ylim(0.6, 1)
# plt.ylabel('$\sigma_a/ v_{cv}$')
# plt.yticks([0.6, 0.7, 0.8, 0.9, 1.0])
# plt.legend(title='R = 0.1$r_v$, \u03B2 = 0')
#
#
# plt.subplot(2, 3, 2)
# plt.plot(a25, X_xi_bquarter_c5_rv_01, label='c=5', color='teal')
# plt.plot(a25, X_xi_bquarter_c10_rv_01, label='c=10', color='blue')
# plt.xlim(0, 1.5)
# plt.xlabel('\u03B1')
# plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
# plt.ylim(0.6, 1)
# plt.ylabel('$\sigma_a/ v_{cv}$')
# plt.yticks([0.6, 0.7, 0.8, 0.9, 1.0])
# plt.legend(title='R = 0.1$r_v$, \u03B2 = 0.25')
#
#
# plt.subplot(2, 3, 3)
# plt.plot(a25, X_xi_bhalf_c5_rv_01, label='c=5', color='teal')
# plt.plot(a25, X_xi_bhalf_c10_rv_01, label='c=10', color='blue')
# plt.xlim(0, 1.5)
# plt.xlabel('\u03B1')
# plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
# plt.ylim(0.6, 1)
# plt.ylabel('$\sigma_a/ v_{cv}$')
# plt.yticks([0.6, 0.7, 0.8, 0.9, 1.0])
# plt.legend(title='R = 0.1$r_v$, \u03B2 = 0.5')
#
#
# plt.subplot(2, 3, 4)
# plt.plot(a25, X_xi_b0_c5_rv_05, label='c=5', color='teal')
# plt.plot(a25, X_xi_b0_c10_rv_05, label='c=10', color='blue')
# plt.xlim(0, 1.5)
# plt.xlabel('\u03B1')
# plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
# plt.ylim(0.6, 1)
# plt.ylabel('$\sigma_a/ v_{cv}$')
# plt.yticks([0.6, 0.7, 0.8, 0.9, 1.0])
# plt.legend(title='R = 0.5$r_v$, \u03B2 = 0')
#
#
# plt.subplot(2, 3, 5)
# plt.plot(a25, X_xi_bquarter_c5_rv_05, label='c=5', color='teal')
# plt.plot(a25, X_xi_bquarter_c10_rv_05, label='c=10', color='blue')
# plt.xlim(0, 1.5)
# plt.xlabel('\u03B1')
# plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
# plt.ylim(0.6, 1)
# plt.ylabel('$\sigma_a/ v_{cv}$')
# plt.yticks([0.6, 0.7, 0.8, 0.9, 1.0])
# plt.legend(title='R = 0.5$r_v$, \u03B2 = 0.25')
#
#
# plt.subplot(2, 3, 6)
# plt.plot(a25, X_xi_bhalf_c5_rv_05, label='c=5', color='teal')
# plt.plot(a25, X_xi_bhalf_c10_rv_05, label='c=10', color='blue')
# plt.xlim(0, 1.5)
# plt.xlabel('\u03B1')
# plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
# plt.ylim(0.6, 1)
# plt.ylabel('$\sigma_a/ v_{cv}$')
# plt.yticks([0.6, 0.7, 0.8, 0.9, 1.0])
# plt.legend(title='R = 0.5$r_v$, \u03B2 = 0.5')
#
#
# plt.tight_layout()
# plt.savefig('Alpha_xi.png', dpi=350)
# plt.show()










#Virial Mass Equation

def virial_mass(vel, v, h, ximin, ximax):
    return [np.sqrt(3/(4*np.pi*G**3*v*crit_density(h)))*1/(ximax**3), np.sqrt(3/(4*np.pi*G**3*v*crit_density(h)))*1/(ximin**3)]*vel**3


def virial_mass_midrange(vel, v, h, ximin, ximax):
    return ((np.sqrt(3/(4*np.pi*G**3*v*crit_density(h)))*1/(ximax**3) + np.sqrt(3/(4*np.pi*G**3*v*crit_density(h)))*1/(ximin**3))/2)*vel**3










#Xi Minimum and Maximum Values for the Alpha-Profiles

xi_min_c5_rv_01 = 0.6052647998743508
xi_max_c5_rv_01 = 0.8770943442280277
xi_min_c10_rv_01 = 0.6626519156072777
xi_max_c10_rv_01 = 0.9658667401727521
xi_min_c5_rv_05 = 0.6261369273579641
xi_max_c5_rv_05 = 0.7398364856097094
xi_min_c10_rv_05 = 0.6499627305948097
xi_max_c10_rv_05 = 0.7966263636598186


#Xi Minimum and Maximum Values for the NFW Profile

xi_min_NFW_c5_rv_01 = 0.6672541640620692
xi_max_NFW_c5_rv_01 = 0.8251234042229985
xi_min_NFW_c10_rv_01 = 0.7659711268323034
xi_max_NFW_c10_rv_01 = 0.8828006133384602
xi_min_NFW_c5_rv_05 = 0.6550669850853098
xi_max_NFW_c5_rv_05 = 0.7075083191334789
xi_min_NFW_c10_rv_05 = 0.7062300320640201
xi_max_NFW_c10_rv_05 = 0.739939809226397










#Mass-200 & Mass-500 for the Alpha-Profiles

def M200_c5_rv_01(vel):
    return virial_mass(vel, 200, 0.6751, xi_min_c5_rv_01, xi_max_c5_rv_01)


def M200_c10_rv_01(vel):
    return virial_mass(vel, 200, 0.6751, xi_min_c10_rv_01, xi_max_c10_rv_01)


def M200_c5_rv_05(vel):
    return virial_mass(vel, 200, 0.6751, xi_min_c5_rv_05, xi_max_c5_rv_05)


def M200_c10_rv_05(vel):
    return virial_mass(vel, 200, 0.6751, xi_min_c10_rv_05, xi_max_c10_rv_05)


def M200_c5_rv_01_midrange(vel):
    return virial_mass_midrange(vel, 200, 0.6751, xi_min_c5_rv_01, xi_max_c5_rv_01)


def M200_c10_rv_01_midrange(vel):
    return virial_mass_midrange(vel, 200, 0.6751, xi_min_c10_rv_01, xi_max_c10_rv_01)


def M200_c5_rv_05_midrange(vel):
    return virial_mass_midrange(vel, 200, 0.6751, xi_min_c5_rv_05, xi_max_c5_rv_05)


def M200_c10_rv_05_midrange(vel):
    return virial_mass_midrange(vel, 200, 0.6751, xi_min_c10_rv_05, xi_max_c10_rv_05)


def M500_c5_rv_01(vel):
    return virial_mass(vel, 500, 0.6751, xi_min_c5_rv_01, xi_max_c5_rv_01)


def M500_c10_rv_01(vel):
    return virial_mass(vel, 500, 0.6751, xi_min_c10_rv_01, xi_max_c10_rv_01)


def M500_c5_rv_05(vel):
    return virial_mass(vel, 500, 0.6751, xi_min_c5_rv_05, xi_max_c5_rv_05)


def M500_c10_rv_05(vel):
    return virial_mass(vel, 500, 0.6751, xi_min_c10_rv_05, xi_max_c10_rv_05)


def M500_c5_rv_01_midrange(vel):
    return virial_mass_midrange(vel, 500, 0.6751, xi_min_c5_rv_01, xi_max_c5_rv_01)


def M500_c10_rv_01_midrange(vel):
    return virial_mass_midrange(vel, 500, 0.6751, xi_min_c10_rv_01, xi_max_c10_rv_01)


def M500_c5_rv_05_midrange(vel):
    return virial_mass_midrange(vel, 500, 0.6751, xi_min_c5_rv_05, xi_max_c5_rv_05)


def M500_c10_rv_05_midrange(vel):
    return virial_mass_midrange(vel, 500, 0.6751, xi_min_c10_rv_05, xi_max_c10_rv_05)


#Mass-200 & Mass-500 for the NFW Profile

def M200_NFW_c5_rv_01(vel):
    return virial_mass(vel, 200, 0.6751, xi_min_NFW_c5_rv_01, xi_max_NFW_c5_rv_01)


def M200_NFW_c10_rv_01(vel):
    return virial_mass(vel, 200, 0.6751, xi_min_NFW_c10_rv_01, xi_max_NFW_c10_rv_01)


def M200_NFW_c5_rv_05(vel):
    return virial_mass(vel, 200, 0.6751, xi_min_NFW_c5_rv_05, xi_max_NFW_c5_rv_05)


def M200_NFW_c10_rv_05(vel):
    return virial_mass(vel, 200, 0.6751, xi_min_NFW_c10_rv_05, xi_max_NFW_c10_rv_05)


def M200_NFW_c5_rv_01_midrange(vel):
    return virial_mass_midrange(vel, 200, 0.6751, xi_min_NFW_c5_rv_01, xi_max_NFW_c5_rv_01)


def M200_NFW_c10_rv_01_midrange(vel):
    return virial_mass_midrange(vel, 200, 0.6751, xi_min_NFW_c10_rv_01, xi_max_NFW_c10_rv_01)


def M200_NFW_c5_rv_05_midrange(vel):
    return virial_mass_midrange(vel, 200, 0.6751, xi_min_NFW_c5_rv_05, xi_max_NFW_c5_rv_05)


def M200_NFW_c10_rv_05_midrange(vel):
    return virial_mass_midrange(vel, 200, 0.6751, xi_min_NFW_c10_rv_05, xi_max_NFW_c10_rv_05)


def M500_NFW_c5_rv_01(vel):
    return virial_mass(vel, 500, 0.6751, xi_min_NFW_c5_rv_01, xi_max_NFW_c5_rv_01)


def M500_NFW_c10_rv_01(vel):
    return virial_mass(vel, 500, 0.6751, xi_min_NFW_c10_rv_01, xi_max_NFW_c10_rv_01)


def M500_NFW_c5_rv_05(vel):
    return virial_mass(vel, 500, 0.6751, xi_min_NFW_c5_rv_05, xi_max_NFW_c5_rv_05)


def M500_NFW_c10_rv_05(vel):
    return virial_mass(vel, 500, 0.6751, xi_min_NFW_c10_rv_05, xi_max_NFW_c10_rv_05)


def M500_NFW_c5_rv_01_midrange(vel):
    return virial_mass_midrange(vel, 500, 0.6751, xi_min_NFW_c5_rv_01, xi_max_NFW_c5_rv_01)


def M500_NFW_c10_rv_01_midrange(vel):
    return virial_mass_midrange(vel, 500, 0.6751, xi_min_NFW_c10_rv_01, xi_max_NFW_c10_rv_01)


def M500_NFW_c5_rv_05_midrange(vel):
    return virial_mass_midrange(vel, 500, 0.6751, xi_min_NFW_c5_rv_05, xi_max_NFW_c5_rv_05)


def M500_NFW_c10_rv_05_midrange(vel):
    return virial_mass_midrange(vel, 500, 0.6751, xi_min_NFW_c10_rv_05, xi_max_NFW_c10_rv_05)










#Alpha Profile Mass Estimate Curves

c5_lower_m200_r10 = 5.104*10**5*(sigma)**3
c5_upper_m200_r10 = 15.532*10**5*(sigma)**3
c10_lower_m200_r10 = 3.282*10**5*(sigma)**3
c10_upper_m200_r10 = 11.836*10**5*(sigma)**3

c5_lower_m200_r50 = 8.505*10**5*(sigma)**3
c5_upper_m200_r50 = 14.030*10**5*(sigma)**3
c10_lower_m200_r50 = 6.813*10**5*(sigma)**3
c10_upper_m200_r50 = 12.543*10**5*(sigma)**3

c5_lower_m500_r10 = 3.228*10**5*(sigma)**3
c5_upper_m500_r10 = 9.823*10**5*(sigma)**3
c10_lower_m500_r10 = 2.417*10**5*(sigma)**3
c10_upper_m500_r10 = 7.486*10**5*(sigma)**3

c5_lower_m500_r50 = 5.379*10**5*(sigma)**3
c5_upper_m500_r50 = 8.873*10**5*(sigma)**3
c10_lower_m500_r50 = 4.309*10**5*(sigma)**3
c10_upper_m500_r50 = 7.933*10**5*(sigma)**3


# plt.figure(figsize=[10, 10])
#
#
# plt.subplot(2, 2, 1)
# plt.plot(sigma, c5_lower_m200_r10, color = 'blue', label = 'c=5')
# plt.plot(sigma, c5_upper_m200_r10, color = 'blue')
# plt.plot(sigma, c10_lower_m200_r10, color = 'darkorange', label = 'c=10')
# plt.plot(sigma, c10_upper_m200_r10, color = 'darkorange')
# plt.fill_between(sigma, c5_lower_m200_r10, c5_upper_m200_r10, alpha=0.2)
# plt.fill_between(sigma, c10_lower_m200_r10, c10_upper_m200_r10, alpha=0.2)
# plt.xlim(100, 1000)
# plt.xlabel('$\sigma_a \quad [\mathrm{kms^{-1}}]$')
# plt.yscale('log')
# plt.ylim(0.5*10**12, 0.2*10**16)
# plt.ylabel('$M_{200} \quad [\mathrm{M}_\odot]$')
# plt.legend(title = '$R=0.1r_{200}$')
#
# plt.subplot(2, 2, 2)
# plt.plot(sigma, c5_lower_m200_r50, color = 'blue', label = 'c=5')
# plt.plot(sigma, c5_upper_m200_r50, color = 'blue')
# plt.plot(sigma, c10_lower_m200_r50, color = 'darkorange', label = 'c=10')
# plt.plot(sigma, c10_upper_m200_r50, color = 'darkorange')
# plt.fill_between(sigma, c5_lower_m200_r50, c5_upper_m200_r50, alpha=0.2)
# plt.fill_between(sigma, c10_lower_m200_r50, c10_upper_m200_r50, alpha=0.2)
# plt.xlim(100, 1000)
# plt.xlabel('$\sigma_a \quad [\mathrm{kms^{-1}}]$')
# plt.yscale('log')
# plt.ylim(0.5*10**12, 0.2*10**16)
# plt.ylabel('$M_{200} \quad [\mathrm{M}_\odot]$')
# plt.legend(title = '$R=0.5r_{200}$')
#
# plt.subplot(2, 2, 3)
# plt.plot(sigma, c5_lower_m500_r10, color = 'blue', label = 'c=5')
# plt.plot(sigma, c5_upper_m500_r10, color = 'blue')
# plt.plot(sigma, c10_lower_m500_r10, color = 'darkorange', label = 'c=10')
# plt.plot(sigma, c10_upper_m500_r10, color = 'darkorange')
# plt.fill_between(sigma, c5_lower_m500_r10, c5_upper_m500_r10, alpha=0.2)
# plt.fill_between(sigma, c10_lower_m500_r10, c10_upper_m500_r10, alpha=0.2)
# plt.xlim(100, 1000)
# plt.xlabel('$\sigma_a \quad [\mathrm{kms^{-1}}]$')
# plt.yscale('log')
# plt.ylim(0.5*10**12, 0.2*10**16)
# plt.ylabel('$M_{500} \quad [\mathrm{M}_\odot]$')
# plt.legend(title = '$R=0.1r_{500}$')
#
# plt.subplot(2, 2, 4)
# plt.plot(sigma, c5_lower_m500_r50, color = 'blue', label = 'c=5')
# plt.plot(sigma, c5_upper_m500_r50, color = 'blue')
# plt.plot(sigma, c10_lower_m500_r50, color = 'darkorange', label = 'c=10')
# plt.plot(sigma, c10_upper_m500_r50, color = 'darkorange')
# plt.fill_between(sigma, c5_lower_m500_r50, c5_upper_m500_r50, alpha=0.2)
# plt.fill_between(sigma, c10_lower_m500_r50, c10_upper_m500_r50, alpha=0.2)
# plt.xlim(100, 1000)
# plt.xlabel('$\sigma_a \quad [\mathrm{kms^{-1}}]$')
# plt.yscale('log')
# plt.ylim(0.5*10**12, 0.2*10**16)
# plt.ylabel('$M_{500} \quad [\mathrm{M}_\odot]$')
# plt.legend(title = '$R=0.5r_{500}$')
#
# plt.savefig('mass_estimates.png', dpi=350)
# plt.show()










#NFW Profile Mass Estimate Curves

NFW_c5_lower_m200_r10 = 6.131*10**5*(sigma)**3
NFW_c5_upper_m200_r10 = 11.593*10**5*(sigma)**3
NFW_c10_lower_m200_r10 = 5.006*10**5*(sigma)**3
NFW_c10_upper_m200_r10 = 7.664*10**5*(sigma)**3

NFW_c5_lower_m200_r50 = 9.725*10**5*(sigma)**3
NFW_c5_upper_m200_r50 = 12.252*10**5*(sigma)**3
NFW_c10_lower_m200_r50 = 8.501*10**5*(sigma)**3
NFW_c10_upper_m200_r50 = 9.778*10**5*(sigma)**3

NFW_c5_lower_m500_r10 = 3.877*10**5*(sigma)**3
NFW_c5_upper_m500_r10 = 7.332*10**5*(sigma)**3
NFW_c10_lower_m500_r10 = 3.166*10**5*(sigma)**3
NFW_c10_upper_m500_r10 = 4.847*10**5*(sigma)**3

NFW_c5_lower_m500_r50 = 6.150*10**5*(sigma)**3
NFW_c5_upper_m500_r50 = 7.749*10**5*(sigma)**3
NFW_c10_lower_m500_r50 = 5.376*10**5*(sigma)**3
NFW_c10_upper_m500_r50 = 6.184*10**5*(sigma)**3


# plt.figure(figsize=[10, 10])
#
#
# plt.subplot(2, 2, 1)
# plt.plot(sigma, NFW_c5_lower_m200_r10, color = 'blue', label = 'c=5')
# plt.plot(sigma, NFW_c5_upper_m200_r10, color = 'blue')
# plt.plot(sigma, NFW_c10_lower_m200_r10, color = 'darkorange', label = 'c=10')
# plt.plot(sigma, NFW_c10_upper_m200_r10, color = 'darkorange')
# plt.fill_between(sigma, NFW_c5_lower_m200_r10, NFW_c5_upper_m200_r10, alpha=0.2)
# plt.fill_between(sigma, NFW_c10_lower_m200_r10, NFW_c10_upper_m200_r10, alpha=0.2)
# plt.xlim(100, 1000)
# plt.xlabel('$\sigma_a \quad [\mathrm{kms^{-1}}]$')
# plt.yscale('log')
# plt.ylim(0.5*10**12, 0.2*10**16)
# plt.ylabel('$M_{200} \quad [\mathrm{M}_\odot]$')
# plt.legend(title = '$R=0.1r_{200}$')
#
# plt.subplot(2, 2, 2)
# plt.plot(sigma, NFW_c5_lower_m200_r50, color = 'blue', label = 'c=5')
# plt.plot(sigma, NFW_c5_upper_m200_r50, color = 'blue')
# plt.plot(sigma, NFW_c10_lower_m200_r50, color = 'darkorange', label = 'c=10')
# plt.plot(sigma, NFW_c10_upper_m200_r50, color = 'darkorange')
# plt.fill_between(sigma, NFW_c5_lower_m200_r50, NFW_c5_upper_m200_r50, alpha=0.2)
# plt.fill_between(sigma, NFW_c10_lower_m200_r50, NFW_c10_upper_m200_r50, alpha=0.2)
# plt.xlim(100, 1000)
# plt.xlabel('$\sigma_a \quad [\mathrm{kms^{-1}}]$')
# plt.yscale('log')
# plt.ylim(0.5*10**12, 0.2*10**16)
# plt.ylabel('$M_{200} \quad [\mathrm{M}_\odot]$')
# plt.legend(title = '$R=0.5r_{200}$')
#
# plt.subplot(2, 2, 3)
# plt.plot(sigma, NFW_c5_lower_m500_r10, color = 'blue', label = 'c=5')
# plt.plot(sigma, NFW_c5_upper_m500_r10, color = 'blue')
# plt.plot(sigma, NFW_c10_lower_m500_r10, color = 'darkorange', label = 'c=10')
# plt.plot(sigma, NFW_c10_upper_m500_r10, color = 'darkorange')
# plt.fill_between(sigma, NFW_c5_lower_m500_r10, NFW_c5_upper_m500_r10, alpha=0.2)
# plt.fill_between(sigma, NFW_c10_lower_m500_r10, NFW_c10_upper_m500_r10, alpha=0.2)
# plt.xlim(100, 1000)
# plt.xlabel('$\sigma_a \quad [\mathrm{kms^{-1}}]$')
# plt.yscale('log')
# plt.ylim(0.5*10**12, 0.2*10**16)
# plt.ylabel('$M_{500} \quad [\mathrm{M}_\odot]$')
# plt.legend(title = '$R=0.1r_{500}$')
#
# plt.subplot(2, 2, 4)
# plt.plot(sigma, NFW_c5_lower_m500_r50, color = 'blue', label = 'c=5')
# plt.plot(sigma, NFW_c5_upper_m500_r50, color = 'blue')
# plt.plot(sigma, NFW_c10_lower_m500_r50, color = 'darkorange', label = 'c=10')
# plt.plot(sigma, NFW_c10_upper_m500_r50, color = 'darkorange')
# plt.fill_between(sigma, NFW_c5_lower_m500_r50, NFW_c5_upper_m500_r50, alpha=0.2)
# plt.fill_between(sigma, NFW_c10_lower_m500_r50, NFW_c10_upper_m500_r50, alpha=0.2)
# plt.xlim(100, 1000)
# plt.xlabel('$\sigma_a \quad [\mathrm{kms^{-1}}]$')
# plt.yscale('log')
# plt.ylim(0.5*10**12, 0.2*10**16)
# plt.ylabel('$M_{500} \quad [\mathrm{M}_\odot]$')
# plt.legend(title = '$R=0.5r_{500}$')
#
# plt.savefig('NFW_mass_estimates.png', dpi=350)
# plt.show()










#Mass Error for the Alpha-Profiles

M200_c5_rv_01_percent_error = ((M200_c5_rv_01(1)[1])/M200_c5_rv_01_midrange(1) - 1)*100
M200_c5_rv_05_percent_error = ((M200_c5_rv_05(1)[1])/M200_c5_rv_05_midrange(1) - 1)*100
M200_c10_rv_01_percent_error = ((M200_c10_rv_01(1)[1])/M200_c10_rv_01_midrange(1) - 1)*100
M200_c10_rv_05_percent_error = ((M200_c10_rv_05(1)[1])/M200_c10_rv_05_midrange(1) - 1)*100
M500_c5_rv_01_percent_error = ((M500_c5_rv_01(1)[1])/M500_c5_rv_01_midrange(1) - 1)*100
M500_c5_rv_05_percent_error = ((M500_c5_rv_05(1)[1])/M500_c5_rv_05_midrange(1) - 1)*100
M500_c10_rv_01_percent_error = ((M500_c10_rv_01(1)[1])/M500_c10_rv_01_midrange(1) - 1)*100
M500_c10_rv_05_percent_error = ((M500_c10_rv_05(1)[1])/M500_c10_rv_05_midrange(1) - 1)*100


#Mass Error for the NFW Profile

M200_NFW_c5_rv_01_percent_error = ((M200_NFW_c5_rv_01(1)[1])/M200_NFW_c5_rv_01_midrange(1) - 1)*100
M200_NFW_c5_rv_05_percent_error = ((M200_NFW_c5_rv_05(1)[1])/M200_NFW_c5_rv_05_midrange(1) - 1)*100
M200_NFW_c10_rv_01_percent_error = ((M200_NFW_c10_rv_01(1)[1])/M200_NFW_c10_rv_01_midrange(1) - 1)*100
M200_NFW_c10_rv_05_percent_error = ((M200_NFW_c10_rv_05(1)[1])/M200_NFW_c10_rv_05_midrange(1) - 1)*100
M500_NFW_c5_rv_01_percent_error = ((M500_NFW_c5_rv_01(1)[1])/M500_NFW_c5_rv_01_midrange(1) - 1)*100
M500_NFW_c5_rv_05_percent_error = ((M500_NFW_c5_rv_05(1)[1])/M500_NFW_c5_rv_05_midrange(1) - 1)*100
M500_NFW_c10_rv_01_percent_error = ((M500_NFW_c10_rv_01(1)[1])/M500_NFW_c10_rv_01_midrange(1) - 1)*100
M500_NFW_c10_rv_05_percent_error = ((M500_NFW_c10_rv_05(1)[1])/M500_NFW_c10_rv_05_midrange(1) - 1)*100


#Mass Uncertainty for the Alpha-Profiles

M200_c5_rv_01_uncertainty = M200_c5_rv_01(1)[1] - M200_c5_rv_01_midrange(1)
M200_c5_rv_05_uncertainty = M200_c5_rv_05(1)[1] - M200_c5_rv_05_midrange(1)
M200_c10_rv_01_uncertainty = M200_c10_rv_01(1)[1] - M200_c10_rv_01_midrange(1)
M200_c10_rv_05_uncertainty = M200_c10_rv_05(1)[1] - M200_c10_rv_05_midrange(1)
M500_c5_rv_01_uncertainty = M500_c5_rv_01(1)[1] - M500_c5_rv_01_midrange(1)
M500_c5_rv_05_uncertainty = M500_c5_rv_05(1)[1] - M500_c5_rv_05_midrange(1)
M500_c10_rv_01_uncertainty = M500_c10_rv_01(1)[1] - M500_c10_rv_01_midrange(1)
M500_c10_rv_05_uncertainty = M500_c10_rv_05(1)[1] - M500_c10_rv_05_midrange(1)


#Mass Uncertainty for the NFW Profile

M200_NFW_c5_rv_01_uncertainty = M200_NFW_c5_rv_01(1)[1] - M200_NFW_c5_rv_01_midrange(1)
M200_NFW_c5_rv_05_uncertainty = M200_NFW_c5_rv_05(1)[1] - M200_NFW_c5_rv_05_midrange(1)
M200_NFW_c10_rv_01_uncertainty = M200_NFW_c10_rv_01(1)[1] - M200_NFW_c10_rv_01_midrange(1)
M200_NFW_c10_rv_05_uncertainty = M200_NFW_c10_rv_05(1)[1] - M200_NFW_c10_rv_05_midrange(1)
M500_NFW_c5_rv_01_uncertainty = M500_NFW_c5_rv_01(1)[1] - M500_NFW_c5_rv_01_midrange(1)
M500_NFW_c5_rv_05_uncertainty = M500_NFW_c5_rv_05(1)[1] - M500_NFW_c5_rv_05_midrange(1)
M500_NFW_c10_rv_01_uncertainty = M500_NFW_c10_rv_01(1)[1] - M500_NFW_c10_rv_01_midrange(1)
M500_NFW_c10_rv_05_uncertainty = M500_NFW_c10_rv_05(1)[1] - M500_NFW_c10_rv_05_midrange(1)










#Alpha Profile Xi & Mass Bounds

# print('R = 0.1rv, \u03B2 = 0, c=5 :   \u03BEmax = ' + str(np.max(X_xi_b0_c5_rv_01)) + ', \u03BEmin = ' + str(np.min(X_xi_b0_c5_rv_01)))
# print('R = 0.1rv, \u03B2 = 0.25, c=5 :   \u03BEmax = ' + str(np.max(X_xi_bquarter_c5_rv_01)) + ', \u03BEmin = ' + str(np.min(X_xi_bquarter_c5_rv_01)))
# print('R = 0.1rv, \u03B2 = 0.5, c=5 :   \u03BEmax = ' + str(np.max(X_xi_bhalf_c5_rv_01)) + ', \u03BEmin = ' + str(np.min(X_xi_bhalf_c5_rv_01)))
# print('R = 0.1rv, \u03B2 = 0, c=10 :   \u03BEmax = ' + str(np.max(X_xi_b0_c10_rv_01)) + ', \u03BEmin = ' + str(np.min(X_xi_b0_c10_rv_01)))
# print('R = 0.1rv, \u03B2 = 0.25, c=10 :   \u03BEmax = ' + str(np.max(X_xi_bquarter_c10_rv_01)) + ', \u03BEmin = ' + str(np.min(X_xi_bquarter_c10_rv_01)))
# print('R = 0.1rv, \u03B2 = 0.5, c=10 :   \u03BEmax = ' + str(np.max(X_xi_bhalf_c10_rv_01)) + ', \u03BEmin = ' + str(np.min(X_xi_bhalf_c10_rv_01)))
# print('R = 0.5rv, \u03B2 = 0, c=5 :   \u03BEmax = ' + str(np.max(X_xi_b0_c5_rv_05)) + ', \u03BEmin = ' + str(np.min(X_xi_b0_c5_rv_05)))
# print('R = 0.5rv, \u03B2 = 0.25, c=5 :   \u03BEmax = ' + str(np.max(X_xi_bquarter_c5_rv_05)) + ', \u03BEmin = ' + str(np.min(X_xi_bquarter_c5_rv_05)))
# print('R = 0.5rv, \u03B2 = 0.5, c=5 :   \u03BEmax = ' + str(np.max(X_xi_bhalf_c5_rv_05)) + ', \u03BEmin = ' + str(np.min(X_xi_bhalf_c5_rv_05)))
# print('R = 0.5rv, \u03B2 = 0, c=10 :   \u03BEmax = ' + str(np.max(X_xi_b0_c10_rv_05)) + ', \u03BEmin = ' + str(np.min(X_xi_b0_c10_rv_05)))
# print('R = 0.5rv, \u03B2 = 0.25, c=10 :   \u03BEmax = ' + str(np.max(X_xi_bquarter_c10_rv_05)) + ', \u03BEmin = ' + str(np.min(X_xi_bquarter_c10_rv_05)))
# print('R = 0.5rv, \u03B2 = 0.5, c=10 :   \u03BEmax = ' + str(np.max(X_xi_bhalf_c10_rv_05)) + ', \u03BEmin = ' + str(np.min(X_xi_bhalf_c10_rv_05)))
#
# print('R = 0.1rv, c=5 :   M200 = ' + str(M200_c5_rv_01(1)) + '*vel^3, ' + ' M500 = ' + str(M500_c5_rv_01(1)) + '*vel^3')
# print('R = 0.1rv, c=10 :   M200 = ' + str(M200_c10_rv_01(1)) + '*vel^3, ' + ' M500 = ' + str(M500_c10_rv_01(1)) + '*vel^3')
# print('R = 0.5rv, c=5 :   M200 = ' + str(M200_c5_rv_05(1)) + '*vel^3, ' + ' M500 = ' + str(M500_c5_rv_05(1)) + '*vel^3')
# print('R = 0.5rv, c=10 :   M200 = ' + str(M200_c10_rv_05(1)) + '*vel^3, ' + ' M500 = ' + str(M500_c10_rv_05(1)) + '*vel^3')
# print('R = 0.1rv, c=5 :   M200 midrange = ' + str(M200_c5_rv_01_midrange(1)) + '*vel^3, ' + ' M500 midrange = ' + str(M500_c5_rv_01_midrange(1)) + '*vel^3')
# print('R = 0.1rv, c=10 :   M200 midrange = ' + str(M200_c10_rv_01_midrange(1)) + '*vel^3, ' + ' M500 midrange = ' + str(M500_c10_rv_01_midrange(1)) + '*vel^3')
# print('R = 0.5rv, c=5 :   M200 midrange = ' + str(M200_c5_rv_05_midrange(1)) + '*vel^3, ' + ' M500 midrange = ' + str(M500_c5_rv_05_midrange(1)) + '*vel^3')
# print('R = 0.5rv, c=10 :   M200 midrange = ' + str(M200_c10_rv_05_midrange(1)) + '*vel^3, ' + ' M500 midrange = ' + str(M500_c10_rv_05_midrange(1)) + '*vel^3')
#
# print('R = 0.1rv, c=5 :   M200 uncertainty = ' + str(M200_c5_rv_01_uncertainty) + ',   M500 uncertainty = ' + str(M500_c5_rv_01_uncertainty))
# print('R = 0.1rv, c=10 :   M200 uncertainty = ' + str(M200_c10_rv_01_uncertainty) + ',   M500 uncertainty = ' + str(M500_c10_rv_01_uncertainty))
# print('R = 0.5rv, c=5 :   M200 uncertainty = ' + str(M200_c5_rv_05_uncertainty) + ',   M500 uncertainty = ' + str(M500_c5_rv_05_uncertainty))
# print('R = 0.5rv, c=10 :   M200 uncertainty = ' + str(M200_c10_rv_05_uncertainty) + ',   M500 uncertainty = ' + str(M500_c10_rv_05_uncertainty))
#
# print('R = 0.1rv, c=5 :   M200 error = ' + str(M200_c5_rv_01_percent_error) + ',   M500 error = ' + str(M500_c5_rv_01_percent_error))
# print('R = 0.1rv, c=10 :   M200 error = ' + str(M200_c10_rv_01_percent_error) + ',   M500 error = ' + str(M500_c10_rv_01_percent_error))
# print('R = 0.5rv, c=5 :   M200 error = ' + str(M200_c5_rv_05_percent_error) + ',   M500 error = ' + str(M500_c5_rv_05_percent_error))
# print('R = 0.5rv, c=10 :   M200 error = ' + str(M200_c10_rv_05_percent_error) + ',   M500 error = ' + str(M500_c10_rv_05_percent_error))


#NFW Profile Xi & Mass Bounds

# print('R = 0.1rv, \u03B1 = 1, c=5 :   \u03BEmax = ' + str(np.max(X_xi_NFW_c5_rv_01)) + ', \u03BEmin = ' + str(np.min(X_xi_NFW_c5_rv_01)))
# print('R = 0.1rv, \u03B1 = 1, c=10 :   \u03BEmax = ' + str(np.max(X_xi_NFW_c10_rv_01)) + ', \u03BEmin = ' + str(np.min(X_xi_NFW_c10_rv_01)))
# print('R = 0.5rv, \u03B1 = 1, c=5 :   \u03BEmax = ' + str(np.max(X_xi_NFW_c5_rv_05)) + ', \u03BEmin = ' + str(np.min(X_xi_NFW_c5_rv_05)))
# print('R = 0.5rv, \u03B1 = 1, c=10 :   \u03BEmax = ' + str(np.max(X_xi_NFW_c10_rv_05)) + ', \u03BEmin = ' + str(np.min(X_xi_NFW_c10_rv_05)))
#
# print('R = 0.1rv, c=5 :   M200 = ' + str(M200_NFW_c5_rv_01(1)) + '*vel^3, ' + ' M500 = ' + str(M500_NFW_c5_rv_01(1)) + '*vel^3')
# print('R = 0.1rv, c=10 :   M200 = ' + str(M200_NFW_c10_rv_01(1)) + '*vel^3, ' + ' M500 = ' + str(M500_NFW_c10_rv_01(1)) + '*vel^3')
# print('R = 0.5rv, c=5 :   M200 = ' + str(M200_NFW_c5_rv_05(1)) + '*vel^3, ' + ' M500 = ' + str(M500_NFW_c5_rv_05(1)) + '*vel^3')
# print('R = 0.5rv, c=10 :   M200 = ' + str(M200_NFW_c10_rv_05(1)) + '*vel^3, ' + ' M500 = ' + str(M500_NFW_c10_rv_05(1)) + '*vel^3')
# print('R = 0.1rv, c=5 :   M200 midrange = ' + str(M200_NFW_c5_rv_01_midrange(1)) + '*vel^3, ' + ' M500 midrange = ' + str(M500_NFW_c5_rv_01_midrange(1)) + '*vel^3')
# print('R = 0.1rv, c=10 :   M200 midrange = ' + str(M200_NFW_c10_rv_01_midrange(1)) + '*vel^3, ' + ' M500 midrange = ' + str(M500_NFW_c10_rv_01_midrange(1)) + '*vel^3')
# print('R = 0.5rv, c=5 :   M200 midrange = ' + str(M200_NFW_c5_rv_05_midrange(1)) + '*vel^3, ' + ' M500 midrange = ' + str(M500_NFW_c5_rv_05_midrange(1)) + '*vel^3')
# print('R = 0.5rv, c=10 :   M200 midrange = ' + str(M200_NFW_c10_rv_05_midrange(1)) + '*vel^3, ' + ' M500 midrange = ' + str(M500_NFW_c10_rv_05_midrange(1)) + '*vel^3')
#
# print('R = 0.1rv, c=5 :   M200 uncertainty = ' + str(M200_NFW_c5_rv_01_uncertainty) + ',   M500 uncertainty = ' + str(M500_NFW_c5_rv_01_uncertainty))
# print('R = 0.1rv, c=10 :   M200 uncertainty = ' + str(M200_NFW_c10_rv_01_uncertainty) + ',   M500 uncertainty = ' + str(M500_NFW_c10_rv_01_uncertainty))
# print('R = 0.5rv, c=5 :   M200 uncertainty = ' + str(M200_NFW_c5_rv_05_uncertainty) + ',   M500 uncertainty = ' + str(M500_NFW_c5_rv_05_uncertainty))
# print('R = 0.5rv, c=10 :   M200 uncertainty = ' + str(M200_NFW_c10_rv_05_uncertainty) + ',   M500 uncertainty = ' + str(M500_NFW_c10_rv_05_uncertainty))
#
# print('R = 0.1rv, c=5 :   M200 error = ' + str(M200_NFW_c5_rv_01_percent_error) + ',   M500 error = ' + str(M500_NFW_c5_rv_01_percent_error))
# print('R = 0.1rv, c=10 :   M200 error = ' + str(M200_NFW_c10_rv_01_percent_error) + ',   M500 error = ' + str(M500_NFW_c10_rv_01_percent_error))
# print('R = 0.5rv, c=5 :   M200 error = ' + str(M200_NFW_c5_rv_05_percent_error) + ',   M500 error = ' + str(M500_NFW_c5_rv_05_percent_error))
# print('R = 0.5rv, c=10 :   M200 error = ' + str(M200_NFW_c10_rv_05_percent_error) + ',   M500 error = ' + str(M500_NFW_c10_rv_05_percent_error))










#Idealised Simulations


r200 = 162.62235170749221
M200 = 9.9956*10**11
vc200 = 162.58940025258275
crit_density_Msol_kpc3 = 2.7754*10**2
m_factor = 0.000966*10**10
N = 126570


path= 'Users/andrewsullivan/PycharmProjects/IdealisedSims/'

cf = h5py.File('nfw_halo_m200_100_c_10.gdt.hdf5', "r")

#print(list(cf['PartType1']))
#print(list(cf['Header'].attrs.keys()))

velocity = cf['PartType1/Velocities'][()]
coordinates = cf['PartType1/Coordinates'][()]
ID = cf['PartType1/ParticleIDs'][()]

Number=cf['Header'].attrs['NumPart_Total'][()]
Mass=cf['Header'].attrs['MassTable'][()]

cf.close()










#Raw Data

x = coordinates[:, 0]
y = coordinates[:, 1]
z = coordinates[:, 2]

r = np.sqrt(x**2 + y**2 + z**2)
R_x = np.sqrt(y**2 + z**2)

v_x = velocity[:, 0]
v_y = velocity[:, 1]
v_z = velocity[:, 2]

v = np.sqrt(v_x**2 + v_y**2 + v_z**2)
v_r = (x*v_x + y*v_y + z*v_z)/(np.sqrt(x**2 + y**2 + z**2))

log_r = np.log10(np.sqrt(x**2 + y**2 + z**2))
log_R_x = np.log10(np.sqrt(y**2 + z**2))










#Radially-Binned Data

r_by_r_binned = []
v_r_by_r_binned = []
v_r_dispersion_by_r_binned = []
density_by_r_binned = []

zip_by_r = zip(r, v_r, ID)
sorted_by_r = sorted(zip_by_r, key = lambda x: x[0])
r_by_r = sorted(r)
log_r_by_r = sorted(log_r)
v_r_by_r = [i[1] for i in sorted_by_r]

vals = []
indices = []
vals.append(-0.6094955)
indices.append(0)

log_r_bins = 131
log_r_steps = 10**(np.arange(-0.20, 2.42, 0.02))


for i in np.arange(-0.20, 2.42, 0.02):  #131 steps
    vals.append(min((log_r_by_r), key=lambda x: abs(x-i)))


for i in range(1, log_r_bins+1):
    indices.append(log_r_by_r.index(vals[i]))


for i in range(0, log_r_bins):
    r_by_r_binned.append((log_r_by_r[indices[i] : indices[i+1]]))
    v_r_by_r_binned.append((v_r_by_r[indices[i]: indices[i+1]]))


for i in range(0, log_r_bins):
    v_r_dispersion_by_r_binned.append(stat.stdev(v_r_by_r_binned[i]))
    density_by_r_binned.append(m_factor * len(r_by_r_binned[i]) / ( (4 / 3) * np.pi * ((10 ** max(r_by_r_binned[i])) ** 3 - (10 ** min(r_by_r_binned[i])) ** 3)))


scaled_log_r_by_r = [i/r200 for i in log_r_steps]
scaled_density_by_r = [i/(200*crit_density_Msol_kpc3) for i in density_by_r_binned]
scaled_v_r_dispersion = [i/vc200 for i in v_r_dispersion_by_r_binned]










#Projected Radial Binned Data

R_x_by_R_binned = []
v_x_by_R_binned = []
v_x_dispersion_by_R_binned = []
v_x_ap_dispersion_by_R_binned = []

zip_proj_x = zip(R_x, v_x, ID)
sorted_proj_x = sorted(zip_proj_x, key=lambda x: x[0])
R_x_by_R = sorted(R_x)
log_R_x_by_R = sorted(log_R_x)
v_x_by_R = [i[1] for i in sorted_proj_x]

vals_proj_x = []
indices_proj_x = []
vals_proj_x.append(-1.9947113)
indices_proj_x.append(0)

log_R_x_bins = 147
log_R_x_steps = 10 ** (np.arange(-0.52, 2.42, 0.02))


for i in np.arange(-0.52, 2.42, 0.02):  # U=147 steps
    vals_proj_x.append(min((log_R_x_by_R), key=lambda x: abs(x - i)))


for i in range(1, log_R_x_bins + 1):
    indices_proj_x.append(log_R_x_by_R.index(vals_proj_x[i]))


for i in range(0, log_R_x_bins):
    R_x_by_R_binned.append((log_R_x_by_R[indices_proj_x[i]: indices_proj_x[i + 1]]))
    v_x_by_R_binned.append((v_x_by_R[indices_proj_x[i]: indices_proj_x[i + 1]]))


for i in range(0, log_R_x_bins):
    v_x_dispersion_by_R_binned.append(stat.stdev(v_x_by_R_binned[i]))


v_x_ap_dispersion_by_R_binned.append(v_x_dispersion_by_R_binned[0])


for i in range(1, log_R_x_bins):
    v_x_ap_dispersion_by_R_binned.append(stat.stdev(v_x_by_R[0:indices_proj_x[i]]))


scaled_log_R_x_by_R = [i / r200 for i in log_R_x_steps]
scaled_v_x_dispersion = [i / vc200 for i in v_x_dispersion_by_R_binned]
scaled_v_x_ap_dispersion = [i / vc200 for i in v_x_ap_dispersion_by_R_binned]










#Panel for Density, Radial Velocity Dispersion, LOS Velocity Dispersion & Aperture Velocity Dispersion for Idealised Simulation


# plt.figure(figsize=[10, 10])
#
#
# plt.subplot(2, 2, 1)
# plt.scatter(scaled_log_r_by_r, scaled_density_by_r, color='deepskyblue')
# plt.plot(s, Density(1, 10, s), label='NFW, c=10', color='darkblue')
# plt.xlabel('$r/r_{200}$')
# plt.xscale('log')
# plt.xlim(0.3*10**(-2), 0.2*10)
# plt.ylabel('$ \u03C1 / 200 \u03C1_{c0} $')
# plt.yscale('log')
# plt.ylim(10**(-2), 10**5)
#
#
# plt.subplot(2, 2, 2)
# plt.scatter(scaled_log_r_by_r, scaled_v_r_dispersion, color='deepskyblue')
# plt.plot(s150, X_NFW_b0_c10_Radial_vel_disp, label='NFW, c=10, \u03b2=0', color='darkblue')
# plt.xlabel('$r/r_{200}$')
# plt.xlim(0.3*10**(-2), 0.2*10**1)
# plt.ylabel('$\sigma_r /v_{c200}$')
# plt.xscale('log')
# plt.ylim(0, 1)
#
#
# plt.subplot(2, 2, 3)
# plt.scatter(scaled_log_R_x_by_R, scaled_v_x_dispersion, color='deepskyblue')
# plt.plot(S50, X_NFW_b0_c10_LOS_vel_disp, label='NFW, c=10, \u03b2=0', color='darkblue')
# plt.xlabel('$R/r_{200}$')
# plt.xlim(0.3*10**(-2), 0.2*10**1)
# plt.ylabel('$\sigma_{\ell} /v_{c200}$')
# plt.xscale('log')
# plt.ylim(0, 1)
#
#
# plt.subplot(2, 2, 4)
# plt.scatter(scaled_log_R_x_by_R, scaled_v_x_ap_dispersion, color='deepskyblue')
# plt.plot(S25, X_NFW_b0_c10_ap_vel_disp, label='NFW, c=10, \u03b2=0', color='darkblue')
# plt.xlabel('$R/r_{200}$')
# plt.xlim(0.3*10**(-2), 0.2*10**1)
# plt.ylabel('$\sigma_{a} /v_{c200}$')
# plt.xscale('log')
# plt.ylim(0, 1)
#
#
# plt.savefig('idealised_sim.png', dpi=350)
# plt.show()










#Sampling the Populations of Particles

#BLUE SAMPLE

rN10_blue = sorted(random.sample(range(788, N), 10))
rN100_blue = sorted(random.sample(range(788, N), 100))
rN1000_blue = sorted(random.sample(range(788, N), 1000))

N10_R_x_by_R_blue = [R_x_by_R[i] for i in rN10_blue]
N100_R_x_by_R_blue = [R_x_by_R[i] for i in rN100_blue]
N1000_R_x_by_R_blue = [R_x_by_R[i] for i in rN1000_blue]

N10_v_x_by_R_blue = [v_x_by_R[i] for i in rN10_blue]
N100_v_x_by_R_blue = [v_x_by_R[i] for i in rN100_blue]
N1000_v_x_by_R_blue = [v_x_by_R[i] for i in rN1000_blue]

N10_v_x_ap_dispersion_by_R_binned_blue = []
N100_v_x_ap_dispersion_by_R_binned_blue = []
N1000_v_x_ap_dispersion_by_R_binned_blue = []


for i in range(1, 10):
    N10_v_x_ap_dispersion_by_R_binned_blue.append(stat.stdev(N10_v_x_by_R_blue[:i+1]))


for i in range(1, 100):
    N100_v_x_ap_dispersion_by_R_binned_blue.append(stat.stdev(N100_v_x_by_R_blue[:i+1]))


for i in range(1, 1000):
    N1000_v_x_ap_dispersion_by_R_binned_blue.append(stat.stdev(N1000_v_x_by_R_blue[:i+1]))


N10_scaled_R_x_by_R_blue = [i/r200 for i in N10_R_x_by_R_blue]
N10_scaled_v_x_ap_dispersion_blue = [i/vc200 for i in N10_v_x_ap_dispersion_by_R_binned_blue]
N100_scaled_R_x_by_R_blue = [i/r200 for i in N100_R_x_by_R_blue]
N100_scaled_v_x_ap_dispersion_blue = [i/vc200 for i in N100_v_x_ap_dispersion_by_R_binned_blue]
N1000_scaled_R_x_by_R_blue = [i/r200 for i in N1000_R_x_by_R_blue]
N1000_scaled_v_x_ap_dispersion_blue = [i/vc200 for i in N1000_v_x_ap_dispersion_by_R_binned_blue]


#ORANGE SAMPLE

rN10_orange = sorted(random.sample(range(788, N), 10))
rN100_orange = sorted(random.sample(range(788, N), 100))
rN1000_orange = sorted(random.sample(range(788, N), 1000))

N10_R_x_by_R_orange = [R_x_by_R[i] for i in rN10_orange]
N100_R_x_by_R_orange = [R_x_by_R[i] for i in rN100_orange]
N1000_R_x_by_R_orange = [R_x_by_R[i] for i in rN1000_orange]

N10_v_x_by_R_orange = [v_x_by_R[i] for i in rN10_orange]
N100_v_x_by_R_orange = [v_x_by_R[i] for i in rN100_orange]
N1000_v_x_by_R_orange = [v_x_by_R[i] for i in rN1000_orange]

N10_v_x_ap_dispersion_by_R_binned_orange = []
N100_v_x_ap_dispersion_by_R_binned_orange = []
N1000_v_x_ap_dispersion_by_R_binned_orange = []


for i in range(1, 10):
    N10_v_x_ap_dispersion_by_R_binned_orange.append(stat.stdev(N10_v_x_by_R_orange[:i+1]))


for i in range(1, 100):
    N100_v_x_ap_dispersion_by_R_binned_orange.append(stat.stdev(N100_v_x_by_R_orange[:i+1]))


for i in range(1, 1000):
    N1000_v_x_ap_dispersion_by_R_binned_orange.append(stat.stdev(N1000_v_x_by_R_orange[:i+1]))


N10_scaled_R_x_by_R_orange = [i/r200 for i in N10_R_x_by_R_orange]
N10_scaled_v_x_ap_dispersion_orange = [i/vc200 for i in N10_v_x_ap_dispersion_by_R_binned_orange]
N100_scaled_R_x_by_R_orange = [i/r200 for i in N100_R_x_by_R_orange]
N100_scaled_v_x_ap_dispersion_orange = [i/vc200 for i in N100_v_x_ap_dispersion_by_R_binned_orange]
N1000_scaled_R_x_by_R_orange = [i/r200 for i in N1000_R_x_by_R_orange]
N1000_scaled_v_x_ap_dispersion_orange = [i/vc200 for i in N1000_v_x_ap_dispersion_by_R_binned_orange]


#MAGENTA SAMPLE

rN10_magenta = sorted(random.sample(range(788, N), 10))
rN100_magenta = sorted(random.sample(range(788, N), 100))
rN1000_magenta = sorted(random.sample(range(788, N), 1000))

N10_R_x_by_R_magenta = [R_x_by_R[i] for i in rN10_magenta]
N100_R_x_by_R_magenta = [R_x_by_R[i] for i in rN100_magenta]
N1000_R_x_by_R_magenta = [R_x_by_R[i] for i in rN1000_magenta]

N10_v_x_by_R_magenta = [v_x_by_R[i] for i in rN10_magenta]
N100_v_x_by_R_magenta = [v_x_by_R[i] for i in rN100_magenta]
N1000_v_x_by_R_magenta = [v_x_by_R[i] for i in rN1000_magenta]

N10_v_x_ap_dispersion_by_R_binned_magenta = []
N100_v_x_ap_dispersion_by_R_binned_magenta = []
N1000_v_x_ap_dispersion_by_R_binned_magenta = []


for i in range(1, 10):
    N10_v_x_ap_dispersion_by_R_binned_magenta.append(stat.stdev(N10_v_x_by_R_magenta[:i+1]))


for i in range(1, 100):
    N100_v_x_ap_dispersion_by_R_binned_magenta.append(stat.stdev(N100_v_x_by_R_magenta[:i+1]))


for i in range(1, 1000):
    N1000_v_x_ap_dispersion_by_R_binned_magenta.append(stat.stdev(N1000_v_x_by_R_magenta[:i+1]))


N10_scaled_R_x_by_R_magenta = [i/r200 for i in N10_R_x_by_R_magenta]
N10_scaled_v_x_ap_dispersion_magenta = [i/vc200 for i in N10_v_x_ap_dispersion_by_R_binned_magenta]
N100_scaled_R_x_by_R_magenta = [i/r200 for i in N100_R_x_by_R_magenta]
N100_scaled_v_x_ap_dispersion_magenta = [i/vc200 for i in N100_v_x_ap_dispersion_by_R_binned_magenta]
N1000_scaled_R_x_by_R_magenta = [i/r200 for i in N1000_R_x_by_R_magenta]
N1000_scaled_v_x_ap_dispersion_magenta = [i/vc200 for i in N1000_v_x_ap_dispersion_by_R_binned_magenta]


#PURPLE SAMPLE

rN10_purple = sorted(random.sample(range(788, N), 10))
rN100_purple = sorted(random.sample(range(788, N), 100))
rN1000_purple = sorted(random.sample(range(788, N), 1000))

N10_R_x_by_R_purple = [R_x_by_R[i] for i in rN10_purple]
N100_R_x_by_R_purple = [R_x_by_R[i] for i in rN100_purple]
N1000_R_x_by_R_purple = [R_x_by_R[i] for i in rN1000_purple]

N10_v_x_by_R_purple = [v_x_by_R[i] for i in rN10_purple]
N100_v_x_by_R_purple = [v_x_by_R[i] for i in rN100_purple]
N1000_v_x_by_R_purple = [v_x_by_R[i] for i in rN1000_purple]

N10_v_x_ap_dispersion_by_R_binned_purple = []
N100_v_x_ap_dispersion_by_R_binned_purple = []
N1000_v_x_ap_dispersion_by_R_binned_purple = []


for i in range(1, 10):
    N10_v_x_ap_dispersion_by_R_binned_purple.append(stat.stdev(N10_v_x_by_R_purple[:i+1]))


for i in range(1, 100):
    N100_v_x_ap_dispersion_by_R_binned_purple.append(stat.stdev(N100_v_x_by_R_purple[:i+1]))


for i in range(1, 1000):
    N1000_v_x_ap_dispersion_by_R_binned_purple.append(stat.stdev(N1000_v_x_by_R_purple[:i+1]))


N10_scaled_R_x_by_R_purple = [i/r200 for i in N10_R_x_by_R_purple]
N10_scaled_v_x_ap_dispersion_purple = [i/vc200 for i in N10_v_x_ap_dispersion_by_R_binned_purple]
N100_scaled_R_x_by_R_purple = [i/r200 for i in N100_R_x_by_R_purple]
N100_scaled_v_x_ap_dispersion_purple = [i/vc200 for i in N100_v_x_ap_dispersion_by_R_binned_purple]
N1000_scaled_R_x_by_R_purple = [i/r200 for i in N1000_R_x_by_R_purple]
N1000_scaled_v_x_ap_dispersion_purple = [i/vc200 for i in N1000_v_x_ap_dispersion_by_R_binned_purple]


#RED SAMPLE

rN10_red = sorted(random.sample(range(788, N), 10))
rN100_red = sorted(random.sample(range(788, N), 100))
rN1000_red = sorted(random.sample(range(788, N), 1000))

N10_R_x_by_R_red = [R_x_by_R[i] for i in rN10_red]
N100_R_x_by_R_red = [R_x_by_R[i] for i in rN100_red]
N1000_R_x_by_R_red = [R_x_by_R[i] for i in rN1000_red]

N10_v_x_by_R_red = [v_x_by_R[i] for i in rN10_red]
N100_v_x_by_R_red = [v_x_by_R[i] for i in rN100_red]
N1000_v_x_by_R_red = [v_x_by_R[i] for i in rN1000_red]

N10_v_x_ap_dispersion_by_R_binned_red = []
N100_v_x_ap_dispersion_by_R_binned_red = []
N1000_v_x_ap_dispersion_by_R_binned_red = []


for i in range(1, 10):
    N10_v_x_ap_dispersion_by_R_binned_red.append(stat.stdev(N10_v_x_by_R_red[:i+1]))


for i in range(1, 100):
    N100_v_x_ap_dispersion_by_R_binned_red.append(stat.stdev(N100_v_x_by_R_red[:i+1]))


for i in range(1, 1000):
    N1000_v_x_ap_dispersion_by_R_binned_red.append(stat.stdev(N1000_v_x_by_R_red[:i+1]))


N10_scaled_R_x_by_R_red = [i/r200 for i in N10_R_x_by_R_red]
N10_scaled_v_x_ap_dispersion_red = [i/vc200 for i in N10_v_x_ap_dispersion_by_R_binned_red]
N100_scaled_R_x_by_R_red = [i/r200 for i in N100_R_x_by_R_red]
N100_scaled_v_x_ap_dispersion_red = [i/vc200 for i in N100_v_x_ap_dispersion_by_R_binned_red]
N1000_scaled_R_x_by_R_red = [i/r200 for i in N1000_R_x_by_R_red]
N1000_scaled_v_x_ap_dispersion_red = [i/vc200 for i in N1000_v_x_ap_dispersion_by_R_binned_red]


# plt.figure(figsize=[10, 5])
#
#
# plt.subplot(1, 2, 1)
# plt.scatter(N100_scaled_R_x_by_R_blue[1:], N100_scaled_v_x_ap_dispersion_blue, color='deepskyblue', alpha=0.3)
# plt.scatter(N100_scaled_R_x_by_R_orange[1:], N100_scaled_v_x_ap_dispersion_orange, color='darkorange', alpha=0.3)
# plt.scatter(N100_scaled_R_x_by_R_magenta[1:], N100_scaled_v_x_ap_dispersion_magenta, color='magenta', alpha=0.3)
# plt.scatter(N100_scaled_R_x_by_R_purple[1:], N100_scaled_v_x_ap_dispersion_purple, color='indigo', alpha=0.3)
# plt.scatter(N100_scaled_R_x_by_R_red[1:], N100_scaled_v_x_ap_dispersion_red, color='red', alpha=0.3)
# plt.plot(S25, X_NFW_b0_c10_ap_vel_disp, color='darkblue')
# plt.xlabel('$R/r_{200}$')
# plt.xscale('log')
# plt.xlim(10**(-1.5), 10**0)
# plt.ylabel('$\sigma_{a} /v_{c200}$')
# plt.ylim(0.4, 1)
# plt.legend(title='N=100')
#
#
# plt.subplot(1, 2, 2)
# plt.scatter(N1000_scaled_R_x_by_R_blue[1:][::10], N1000_scaled_v_x_ap_dispersion_blue[::10], color='deepskyblue', alpha=0.3)
# plt.scatter(N1000_scaled_R_x_by_R_orange[1:][::10], N1000_scaled_v_x_ap_dispersion_orange[::10], color='darkorange', alpha=0.3)
# plt.scatter(N1000_scaled_R_x_by_R_magenta[1:][::10], N1000_scaled_v_x_ap_dispersion_magenta[::10], color='magenta', alpha=0.3)
# plt.scatter(N1000_scaled_R_x_by_R_purple[1:][::10], N1000_scaled_v_x_ap_dispersion_purple[::10], color='indigo', alpha=0.3)
# plt.scatter(N1000_scaled_R_x_by_R_red[1:][::10], N1000_scaled_v_x_ap_dispersion_red[::10], color='red', alpha=0.3)
# plt.plot(S25, X_NFW_b0_c10_ap_vel_disp, color='darkblue')
# plt.xlabel('$R/r_{200}$')
# plt.xscale('log')
# plt.xlim(10**(-1.5), 10**0)
# plt.ylabel('$\sigma_{a} /v_{c200}$')
# plt.ylim(0.4, 1)
# plt.legend(title='N=1000')
#
#
# plt.savefig('sampled_pops.png', dpi=350)
# plt.show()










#Random Scatter in Sampling Process

r01_v_x_ap_count = []
r05_v_x_ap_count = []


for i in range(1, 1000):
    n = sorted(random.sample(range(0, N), i))

    r01_velocity_count = []
    r05_velocity_count = []

    for j in n:
        if R_x[j] < 0.1*r200:
            r01_velocity_count.append(v_x[j])

        if R_x[j] < 0.5*r200:
            r05_velocity_count.append(v_x[j])


    if len(r01_velocity_count) > 1:
        r01_v_x_ap_count.append(stat.stdev(r01_velocity_count))

    if len(r01_velocity_count) <= 1:
        r01_v_x_ap_count.append(np.nan)

    if len(r05_velocity_count) > 1:
        r05_v_x_ap_count.append(stat.stdev(r05_velocity_count))

    if len(r05_velocity_count) <= 1:
        r05_v_x_ap_count.append(np.nan)


r01 = 0.7659711230449848
r05 = 0.7062300300214605

r01_scaled_v_x_ap_count = [i/vc200 for i in r01_v_x_ap_count]
r05_scaled_v_x_ap_count = [i/vc200 for i in r05_v_x_ap_count]

N_vals = np.arange(1, 1000, 1)


# plt.figure(figsize=[10, 5])
#
#
# plt.subplot(1, 2, 1)
# plt.scatter(N_vals, r01_scaled_v_x_ap_count, color='deepskyblue')
# plt.plot(N_vals, r01*N_vals**0, color='darkblue')
# plt.xlabel('$N$')
# plt.xlim(0, 1000)
# plt.ylabel('$\sigma_{a} /v_{c200}$')
# plt.ylim(0.2, 1.2)
# plt.legend(title='$R = 0.1r_{200}$')
#
#
# plt.subplot(1, 2, 2)
# plt.scatter(N_vals, r05_scaled_v_x_ap_count, color='deepskyblue')
# plt.plot(N_vals, r05*N_vals**0, color='darkblue')
# plt.xlabel('$N$')
# plt.xlim(0, 1000)
# plt.ylabel('$\sigma_{a} /v_{c200}$')
# plt.ylim(0.2, 1.2)
# plt.legend(title='$R = 0.5r_{200}$')
#
#
# plt.savefig('sample_stats.png', dpi=350)
# plt.show()
#!/usr/bin/env python3
import numpy as np
from scipy.interpolate import interp1d
import itertools

##SOME CONSTANTS##############################################
hbar = 6.582119514e-16       #eV s
c = 299792458                #m/s
pi = np.pi
###############################################################

##CALCULATES FLUORESCENCE LIFETIME IN S########################
def calc_lifetime(xd,yd):
    #Integrates the emission spectrum
    IntEmi = np.trapz(yd,xd)
    taxa   = (1/hbar)*IntEmi
    return 1/taxa 
###############################################################

##CALCULATES FORSTER RADIUS####################################
def radius(xa,ya,xd,yd,kappa):
    #Speed of light
    c = 299792458  #m/s
    
    #Finds the edges of interpolation
    minA = min(xa)
    minD = min(xd)
    maxA = max(xa)
    maxD = max(xd)
    MIN = max(minA,minD)
    MAX = min(maxA,maxD)

    X = np.linspace(MIN, MAX, 1000)
    f1 = interp1d(xa, ya, kind='cubic')
    f2 = interp1d(xd, yd, kind='cubic')
    YA = f1(X)
    YD = f2(X)

    #Calculates the overlap
    Overlap = YA*YD/(X**4)

    #Integrates overlap
    IntOver = np.trapz(Overlap, X)

    #Gets lifetime
    tau = calc_lifetime(xd,yd)	

    #Calculates radius
    c *= 1e10
    const   = (hbar**3)*(9*(c**4)*(kappa**2)*tau)/(8*pi)
    radius6 = const*IntOver
    radius  = radius6**(1/6)
    return radius
###############################################################

##CALCULATES THE MEAN AND ERROR################################
def avg_error(variable):
    mean   = (max(variable) + min(variable))/2
    error  = (max(variable) - min(variable))/2
    return mean, error
###############################################################

##AVG HOPPING DISTANCES########################################
def r_avg(alpha,moment,rmin,dim):
    neighbors = 100
    n1 = neighbors
    n2 = 0
    n3 = 0
    if dim > 1:
        n2 = neighbors
    if dim > 2:
        n3 = neighbors
    numerator = []
    denominator = []
    for i in range(0,n1+1):
        for j in range(0,n2+1):
            for k in range(0,n3+1):
                r = rmin*np.sqrt(i**2+j**2+k**2)
                if r != 0:
                    numerator.append(r/((alpha*moment + r)**6))
                    denominator.append(1.0/((alpha*moment + r)**6))
    average = np.sum(numerator)/np.sum(denominator)
    return average
###############################################################


##1D PROJECTIONS OF DIFFUSION LENGTH###########################
def difflen(radius, alpha, moment, r, dim, phi):
    ld = r*(radius**3)/((alpha*moment+r)**3)
    ld = np.sqrt(phi)*ld/np.sqrt(dim)
    #Conversion to nm
    ld = ld/10
    return ld
###############################################################

def run_ld(Abs, Emi, alpha, rmin, kappa, Phi):
    #Gets the Transition Dipole Moment in atomic units
    try:
        with open(Emi, 'r') as f:
            for line in f:
                moment = float(line.split('=')[1].split()[0])
                break
    except:
        moment = 0


    #Loads the data 
    data_abs = np.loadtxt(Abs)
    data_emi = np.loadtxt(Emi)

    #Gets energies, intensities and errors for the donor and acceptor
    xa, ya_max, ya_min = data_abs[:,0], data_abs[:,1] + data_abs[:,2], data_abs[:,1] - data_abs[:,2]
    xd, yd_max, yd_min = data_emi[:,0], data_emi[:,1] + data_emi[:,2], data_emi[:,1] - data_emi[:,2]


    #Lifetime calculations
    lifetimes = []
    for yd in [yd_min,yd_max]:
        lifetimes.append(calc_lifetime(xd,yd))
    mean_life, error_life = avg_error(lifetimes)

    #Radius calculations
    radii = []
    for ya, yd in itertools.product([ya_min,ya_max], [yd_min,yd_max]):
        radii.append(radius(xa,ya,xd,yd,kappa))

    mean_radius, error_radius = avg_error(radii)


    #Calculates average hopping distances for each case	
    dista = 0.25*alpha*moment + 1.25*rmin
    dist1 = r_avg(alpha,moment,rmin,1)
    dist2 = r_avg(alpha,moment,rmin,2)
    dist3 = r_avg(alpha,moment,rmin,3)

    LDA, LD1, LD2, LD3 = [], [], [], []
    for Rf in [mean_radius - error_radius,mean_radius + error_radius]:
        #Amorphous average hopping distance	and diffusion length
        LDA.append(difflen(Rf,alpha, moment, dista, 3, Phi))

        #1D crystal average hopping distance and diffusion length		
        LD1.append(difflen(Rf,alpha, moment, dist1, 1, Phi))

        #2D crystal average hopping distance and diffusion length	
        LD2.append(difflen(Rf,alpha, moment, dist2, 2, Phi))

        #3D crystal average hopping distance and diffusion length	 
        LD3.append(difflen(Rf,alpha, moment, dist3, 3, Phi))


    lda_mean, error_lda = avg_error(LDA)
    ld1_mean, error_ld1 = avg_error(LD1)
    ld2_mean, error_ld2 = avg_error(LD2)
    ld3_mean, error_ld3 = avg_error(LD3)


    with open("ld.lx", 'w') as f:
        f.write("Forster Radius:      {:.1f} +/- {:.1f} Å \n".format(mean_radius,error_radius))
        f.write("Radiative Lifetime:  {:.1f} +/- {:.1f} ps\n".format(mean_life*1e12,error_life*1e12))
        f.write("Avg. Dipole Moment:  {:.1f} a.u. \n".format(moment))
        f.write("Morphology   Avg_Hop_Distance(A)  Diffusion_Length(nm)\n")
        f.write("1D Crystal   {:<19.1f}  {:>.1f} +/- {:>.1f}\n".format(dist1, ld1_mean, error_ld1))
        f.write("2D Crystal   {:<19.1f}  {:>.1f} +/- {:>.1f}\n".format(dist2, ld2_mean, error_ld2))
        f.write("3D Crystal   {:<19.1f}  {:>.1f} +/- {:>.1f}\n".format(dist3, ld3_mean, error_ld3))
        f.write("Amorphous    {:<19.1f}  {:>.1f} +/- {:>.1f}\n".format(dista, lda_mean, error_lda))


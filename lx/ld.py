#!/usr/bin/env python3
import numpy as np
from scipy.interpolate import interp1d

##SOME CONSTANTS##############################################
hbar = 6.582119514e-16       #eV s
c    = 299792458             #m/s
pi = np.pi
###############################################################

##CALCULATES FLUORESCENCE LIFETIME IN S########################
def calc_lifetime(xd,yd,dyd):
    #Integrates the emission spectrum
    IntEmi = np.trapz(yd,xd)
    taxa   = (1/hbar)*IntEmi
    error  = (1/hbar)*np.sqrt(np.trapz((dyd**2),xd))
    dlife  = (1/taxa)*(error/taxa)
    return 1/taxa, dlife 
###############################################################

##CALCULATES FORSTER RADIUS####################################
def radius(xa,ya,dya,xd,yd,dyd,kappa):
    #Speed of light
    c = 299792458  #m/s
    
    #Finds the edges of interpolation
    minA = min(xa)
    minD = min(xd)
    maxA = max(xa)
    maxD = max(xd)
    MIN  = max(minA,minD)
    MAX  = min(maxA,maxD)

    X = np.linspace(MIN, MAX, 1000)
    f1 = interp1d(xa, ya, kind='cubic')
    f2 = interp1d(xd, yd, kind='cubic')
    f3 = interp1d(xa, dya, kind='cubic')
    f4 = interp1d(xd, dyd, kind='cubic')
    

    YA  = f1(X)
    YD  = f2(X)
    DYA = f3(X)
    DYD = f4(X)

    #Calculates the overlap
    Overlap = YA*YD/(X**4)

    #Overlap error
    OverError   = Overlap*np.sqrt((DYA/YA)**2 + (DYD/YD)**2)

    #Integrates overlap
    IntOver = np.trapz(Overlap, X)

    #Integrated Overlap Error
    DeltaOver = np.sqrt(np.trapz((OverError**2),X))       

    #Gets lifetime
    tau, delta_tau = calc_lifetime(xd,yd,dyd)	

    #Calculates radius sixth power
    c *= 1e10
    const   = (hbar**3)*(9*(c**4)*(kappa**2)*tau)/(8*pi)
    radius6 = const*IntOver

    #Relative error in radius6
    delta   = np.sqrt((DeltaOver/IntOver)**2 + (delta_tau/tau)**2)

    #Calculates radius
    radius  = radius6**(1/6)

    #Error in radius
    delta_radius = radius*delta/6

    return radius, delta_radius
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
    xa, ya, dya = data_abs[:,0], data_abs[:,1], data_abs[:,2]
    xd, yd, dyd = data_emi[:,0], data_emi[:,1], data_emi[:,2]


    #Lifetime calculations
    mean_life, error_life = calc_lifetime(xd,yd,dyd)
    

    #Radius calculations
    mean_radius, error_radius = radius(xa,ya,dya,xd,yd,dyd,kappa)


    #Calculates average hopping distances for each case	
    dista = 0.25*alpha*moment + 1.25*rmin
    dist1 = r_avg(alpha,moment,rmin,1)
    dist2 = r_avg(alpha,moment,rmin,2)
    dist3 = r_avg(alpha,moment,rmin,3)

    
    lda_mean  = difflen(mean_radius,alpha, moment, dista, 3, Phi)
    error_lda = 3*lda_mean*(error_radius/mean_radius)
    ld1_mean  = difflen(mean_radius,alpha, moment, dist1, 1, Phi)
    error_ld1 = 3*ld1_mean*(error_radius/mean_radius)
    ld2_mean  = difflen(mean_radius,alpha, moment, dist2, 2, Phi)
    error_ld2 = 3*ld2_mean*(error_radius/mean_radius)
    ld3_mean  = difflen(mean_radius,alpha, moment, dist3, 3, Phi)
    error_ld3 = 3*ld3_mean*(error_radius/mean_radius)


    with open("ld.lx", 'w') as f:
        f.write("Foerster Radius:      {:.1f} +/- {:.1f} AA \n".format(mean_radius,error_radius))
        f.write("Radiative Lifetime:  {:.3e} +/- {:.3e} s\n".format(mean_life,error_life))
        f.write("Avg. Dipole Moment:  {:.1f} a.u. \n".format(moment))
        f.write("Morphology   Avg_Hop_Distance(A)  Diffusion_Length(nm)\n")
        f.write("1D Crystal   {:<19.1f}  {:>.1f} +/- {:>.1f}\n".format(dist1, ld1_mean, error_ld1))
        f.write("2D Crystal   {:<19.1f}  {:>.1f} +/- {:>.1f}\n".format(dist2, ld2_mean, error_ld2))
        f.write("3D Crystal   {:<19.1f}  {:>.1f} +/- {:>.1f}\n".format(dist3, ld3_mean, error_ld3))
        f.write("Amorphous    {:<19.1f}  {:>.1f} +/- {:>.1f}\n".format(dista, lda_mean, error_lda))

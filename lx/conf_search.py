#!/usr/bin/env python3
import numpy as np
import os
import sys
import subprocess
import random
import shutil
import lx.tools
from scipy.stats import norm

c     = lx.tools.c
pi    = lx.tools.pi
hbar  = lx.tools.hbar
hbar2 = lx.tools.hbar2
kb    = lx.tools.kb


def distance_matrix(G):
    matrix = np.zeros((1,np.shape(G)[0]))
    for ind in range(np.shape(G)[0]):
            distances = G - G[ind,:]
            distances = np.sqrt(np.sum(np.square(distances),axis=1))    
            matrix = np.vstack((matrix,distances[np.newaxis,:]))
    matrix = matrix[1:,:]
    matrix[matrix==0] = np.inf
    return matrix

def bond(matrix,atoms):
    valence = {1:0.31, 2:0.28, 3:1.28, 4:0.96, 5:0.84, 6:0.76, 7:0.71, 8:0.66, 9:0.57, 10:0.58, 11:1.66, 12:1.41, 13:1.21, 14:1.11,
           15:1.07, 16:1.05, 17:1.02, 18:1.06, 19:2.03, 20:1.76, 21:1.7, 22:1.6, 23:1.53, 24:1.39, 25:1.61, 26:1.52, 27:1.50,
           28:1.24, 29:1.32, 30:1.22, 31:1.22, 32:1.2, 33:1.19, 34:1.20, 35:1.20, 36:1.16}
    CM = np.zeros(np.shape(matrix))
    while np.min(matrix) < 100:
        x,y = np.where(matrix == np.min(matrix))
        x, y = x[0], y[0]         
        if matrix[x,y] < 1.3*(valence[atoms[x]] + valence[atoms[y]]):
            CM[x,y] += 1
            CM[y,x] += 1      
        matrix[x,y] = np.inf
        matrix[y,x] = np.inf   
    return CM      

def fingerprint(file,folder):
    try:
        G, atoms = lx.tools.pega_geom(folder+'/'+file)
        atoms = np.array(atoms).astype(float)
        matrix = distance_matrix(G)
        cm = bond(matrix,atoms)
    except:
        cm = np.zeros((5,5))
    return cm              

##SAMPLES GEOMETRIES###########################################
def make_geoms(freqlog, num_geoms, T, header, bottom):
    lista = []
    counter = lx.tools.start_counter()
    _, atomos, A = lx.tools.sample_geometries(freqlog,num_geoms,T,3000)
    for n in range(0,np.shape(A)[1],3):
        Gfinal = A[:,n:n+3]
        lx.tools.write_input(atomos,Gfinal,header.replace("UUUUU",str((n+3)//3)),bottom.replace("UUUUU",str((n+3)//3)),"Geometry-"+str((n+3)//3+counter)+"-.com")
        lista.append("Geometry-"+str((n+3)//3+counter)+"-.com") 
    return lista      
############################################################### 

##GETS ENERGY FROM THE ORIGINAL FREQ LOG FILE##################
def get_energy_origin(freqlog):
    with open(freqlog, 'r') as f:
        for line in f:
            if 'SCF Done:' in line:
                line = line.split()
                scf  = float(line[4])*27.2114
            elif 'Rotational constants' in line:
                line = line.split()
                rot  = [float(line[3]),float(line[4]),float(line[5])]    
            elif 'Normal termination' in line:
                return scf, rot[0], rot[1], rot[2]
###############################################################

##GETS ENERGIES FROM OPT LOG FILES#############################
def get_energies(folder,original_molecule):
    nums,scfs,rotsx, rotsy, rotsz = [], [], [], [], []
    files = [i for i in os.listdir(folder) if '.log' in i and 'Geometry' in i]
    for file in files:
        if np.array_equal(original_molecule, fingerprint(file,folder)):
            with open(folder+'/'+file, 'r') as f:
                num = float(file.split('-')[1])
                for line in f:
                    if 'SCF Done:' in line:
                        line = line.split()
                        scf  = float(line[4])*27.2114
                    elif 'Rotational constants' in line:
                        line = line.split()
                        rotx = float(line[3])
                        roty = float(line[4])
                        rotz = float(line[5])
                    elif 'Normal termination' in line:
                        scfs.append(scf)
                        nums.append(num)
                        rotsx.append(rotx)
                        rotsy.append(roty)
                        rotsz.append(rotz)
    for file in files:
        try:
            shutil.move(file, 'Geometries/'+file)
            shutil.move(file[:-3]+'com', 'Geometries/'+file[:-3]+'com')
        except:
            pass
    nums = np.array(nums)
    scfs = np.array(scfs)
    rotsx = np.array(rotsx)
    rotsy = np.array(rotsy)
    rotsz = np.array(rotsz)
    return nums, scfs, rotsx, rotsy, rotsz
###############################################################

def measure(vec1,vec2,cr):
    vec1 = np.array(vec1)     
    vec2 = np.array(vec2)
    cr   = np.array(cr)
    dist = max(abs(vec1 -vec2) - cr)
    distance =  np.heaviside(dist,0)   
    return distance

##CLASSIFIES THE VARIOUS OPTIMIZED STRUCTURES##################
def classify(nums,scfs,rotsx, rotsy, rotsz,cr0,first):
    try:
        assert not first
        data = np.loadtxt('conformation.lx')
        if len(np.shape(data)) > 1:
            engs =  data[:,1].flatten()
            rotx =  data[:,4].flatten()
            roty =  data[:,5].flatten() 
            rotz =  data[:,6].flatten()
            last =  data[:,10].flatten()
            crix =  data[:,7].flatten()
            criy =  data[:,8].flatten()
            criz =  data[:,9].flatten()
            exam =  data[:,11].flatten()
        else:
            engs = np.array([data[1]])
            rotx = np.array([data[4]])
            roty = np.array([data[5]])
            rotz = np.array([data[6]]) 
            last = np.array([data[10]])
            crix = np.array([data[7]])
            criy = np.array([data[8]])
            criz = np.array([data[9]])
            exam = np.array([data[11]])
    except:
        rotx = np.array([])
        roty = np.array([])
        rotz = np.array([])
        engs = np.array([])
        last = np.array([])
        crix = np.array([])
        criy = np.array([])
        criz = np.array([])
        exam = np.array([])
    new = []
    for m in range(len(rotsx)):
        distances = []   
        ROTX = np.copy(rotx)
        ROTY = np.copy(roty)
        ROTZ = np.copy(rotz)
        for n in range(len(ROTX)):
            try:
                cr = [max(min(2*crix[n],cr0[0]),cr0[0]),max(min(2*crix[n],cr0[1]),cr0[1]),max(min(2*crix[n],cr0[2]),cr0[2])]
            except:
                cr = cr0    
            distance = measure([rotsx[m], rotsy[m], rotsz[m]], [ROTX[n],ROTY[n],ROTZ[n]],cr)
            distances.append(distance)
        try:
            a = distances.index(0)
            exam[a] = nums[m]
            engs[a] = (last[a]*engs[a] + scfs[m])/(last[a]+1)
            x2      = (last[a]*(crix[a]**2 +rotx[a]**2) + rotsx[m]**2)/(last[a]+1)
            y2      = (last[a]*(criy[a]**2 +roty[a]**2) + rotsy[m]**2)/(last[a]+1)
            z2      = (last[a]*(criz[a]**2 +rotz[a]**2) + rotsz[m]**2)/(last[a]+1)
            rotx[a] = (last[a]*rotx[a] + rotsx[m])/(last[a]+1)
            roty[a] = (last[a]*roty[a] + rotsy[m])/(last[a]+1)
            rotz[a] = (last[a]*rotz[a] + rotsz[m])/(last[a]+1)    
            crix[a] = np.sqrt(x2 -rotx[a]**2)
            criy[a] = np.sqrt(y2 -roty[a]**2)
            criz[a] = np.sqrt(z2 -rotz[a]**2)
            last[a] += 1
        except:
            new.append(nums[m])
            rotx  = np.append(rotx,rotsx[m])
            roty  = np.append(roty,rotsy[m])
            rotz  = np.append(rotz,rotsz[m])
            crix  = np.append(crix,0)
            criy  = np.append(criy,0)
            criz  = np.append(criz,0) 
            engs  = np.append(engs,scfs[m])
            last  = np.append(last,1)
            exam  = np.append(exam,nums[m])

    try:
        origin = random.choice(new)
    except:
        origin = 0

    probs  = np.exp(-1*(engs - min(engs))/0.026)
    probs  /= np.sum(probs)
    args   = np.argsort(engs)
    engs   = engs[args]
    probs  = probs[args]
    exam   = exam[args]
    rotx   = rotx[args]   
    roty   = roty[args]
    rotz   = rotz[args]
    crix   = crix[args]
    criy   = criy[args]
    criz   = criz[args]
    last   = last[args]
    
    with open('conformation.lx', 'w') as f: 
        f.write('{:6}  {:10}  {:10}  {:12}  {:10}  {:10}  {:10}  {:8}  {:8}  {:8}  {:6}  {:6}\n'.format('#Group','Energy(eV)','DeltaE(eV)','Prob@300K(%)','Rot1','Rot2','Rot3','Std1','Std2','Std3','Number','Last'))
        for i in range(len(probs)):
            f.write('{:<6}  {:<10.3f}  {:<10.3f}  {:<12.1f}  {:<10.7f}  {:<10.7f}  {:<10.7f}  {:<8.2e}  {:<8.2e}  {:<8.2e}  {:<6.0f}  {:<6.0f}\n'.format(i+1,engs[i],engs[i] -min(engs),100*probs[i],rotx[i],roty[i],rotz[i], crix[i], criy[i], criz[i], last[i],exam[i]))
    return int(origin), exam
###############################################################

##RUNS FREQ CALCULATION FOR NEW CONFORMATION###################
def rodar_freq(origin,nproc,mem,base,cm,batch_file,gaussian):
    geomlog = 'Geometries/Geometry-'+str(origin)+'-.log'
    G, atomos = lx.tools.pega_geom(geomlog) 
    header = "%nproc={}\n%mem={}\n# freq=(noraman) nosymm  {} \n\n{}\n\n{}\n".format(nproc,mem,base,'ABSSPCT',cm)
    file = "Freq-"+str(origin)+"-.com"
    lx.tools.write_input(atomos,G,header,'',file)
    lx.tools.rodar_lista([file],batch_file,gaussian,'conformation.lx')
    log = file[:-3]+'log'
    with open(log, 'r') as f:
        for line in f:
            if 'Normal termination' in line:
                return log
            elif 'Error termination' in line:
                return None
###############################################################


def main():
    freqlog   = sys.argv[1]
    base      = sys.argv[2]
    nproc     = sys.argv[3]
    mem       = sys.argv[4]
    T         = float(sys.argv[5])
    DT        = float(sys.argv[6])
    num_geoms = int(sys.argv[7])
    rounds    = int(sys.argv[8])
    script    = sys.argv[9]
    gaussian  = sys.argv[10]

    freq0 = freqlog
    try:
        os.mkdir('Geometries')
    except:
        pass    
    original_molecule = fingerprint(freqlog,'.')
    cm        = lx.tools.get_cm(freqlog) 
    header    = "%nproc={}\n%mem={}\n# opt nosymm  {} \n\n{}\n\n{}\n".format(nproc,mem,base,'ABSSPCT',cm)
    scf, rotx, roty, rotz  = get_energy_origin(freqlog)
    cr0 = [rotx/1000, roty/1000, rotz/1000]
    origin, conformation = classify(np.array([0]),np.array([scf]), np.array([rotx]),np.array([roty]),np.array([rotz]),cr0,True)
    files = [i for i in os.listdir('Geometries') if 'Geometry' in i and '.log' in i]
    if len(files) > 0:
        nums, scfs, rotsx, rotsy, rotsz = get_energies('Geometries',original_molecule)
        origin, conformation  = classify(nums,scfs,rotsx,rotsy,rotsz,cr0,False)
    else:
        pass    

    T0 = T

    for i in range(rounds):
        lista      = make_geoms(freqlog, num_geoms, T0, header, '')
        lx.tools.rodar_lista(lista,script, gaussian, 'conformation.lx')
        nums, scfs, rotsx,rotsy,rotsz = get_energies('.',original_molecule)
        origin, conformation  = classify(nums,scfs,rotsx,rotsy,rotsz,cr0,False)
        with open('conformation.lx', 'a') as f:
            f.write('\n#Round {}/{} Temperature: {} K'.format(i+1,rounds,T0))    
        if origin != 0:
            log  = rodar_freq(origin,nproc,mem,base,cm,script,gaussian)
            if log != None:
                freqlog = log
                T0 = T
        else:
            T0 += DT

    with open('conformation.lx', 'a') as f:
        f.write('\n#Search concluded!')

    try:
        os.mkdir('Conformers')
    except:
        pass    

    for i in range(len(conformation)):
        numero  = conformation[i]
        if numero == 0:
            freqlog = freq0    
        else:
            freqlog = 'Geometries/Geometry-{:.0f}-.log'.format(numero) 
        _, _, nproc, mem, scrf, _ = lx.tools.busca_input(freqlog)
        cm = lx.tools.get_cm(freqlog)
        header = '%nproc={}\n%mem={}\n%chk=Group_{}_.chk\n# {} {} opt\n\nTITLE\n\n{}\n'.format(nproc,mem,i+1,'pm6',scrf,cm)
        G, atomos = lx.tools.pega_geom(freqlog)
        lx.tools.write_input(atomos,G,header,'','Conformers/Group_{}_.com'.format(i+1))
    



if __name__ == "__main__":
    sys.exit(main())        


        
        
        


    
        


#!/usr/bin/env python3
import numpy as np
import os
import sys
import subprocess
import random
from lx.tools import *
import shutil
import time

##CHECKS WHETHER JOBS ARE DONE#################################
def hold_watch(files):
    rodando = files.copy()
    while len(rodando) > 0:
        rodando = watcher(rodando,1)
        if 'limit.lx' not in os.listdir('.'):
            with open('conformation.lx','a') as f:
                f.write('\n#Aborted!')
            sys.exit()
        time.sleep(30)    
###############################################################

##RUNS CALCULATIONS############################################
def rodar_opts(lista, batch_file): 
    for file in lista:
        subprocess.call(['bash', batch_file, file]) 
    hold_watch(lista)
###############################################################


##SAMPLES GEOMETRIES###########################################
def make_geoms(freqlog, num_geoms, T, header, bottom):
    lista = []
    F, M = pega_freq(freqlog)
    F[F < 0] *= -1    
    counter = start_counter()
    G, atomos = pega_geom(freqlog)
    NNC = pega_modos(G,freqlog)
    num_atom = np.shape(G)[0]   
    for n in range(1,num_geoms+1):
        A = np.zeros((3*num_atom,1))
        numbers = []
        for i in range(0,len(F)):
            scale = np.sqrt(hbar2/(2*M[i]*F[i]*np.tanh(hbar*F[i]/(2*kb*T))))
            normal = norm(scale=scale,loc=0)
            #Displacements in  Å
            q = normal.rvs()*1e10
            numbers.append(q)
            A += q*(np.expand_dims(NNC[:,i],axis=1))
        numbers = np.round(np.array(numbers)[np.newaxis,:],4)
        A = np.reshape(A,(num_atom,3))
        Gfinal = A + G  
        write_input(atomos,Gfinal,header.replace("UUUUU",str(n)),bottom.replace("UUUUU",str(n)),"Geometry-"+str(n+counter)+"-.com")
        lista.append("Geometry-"+str(n+counter)+"-.com") 
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
def get_energies():
    nums,scfs,rotsx, rotsy, rotsz = [], [], [], [], []
    files = [i for i in os.listdir('.') if '.log' in i and 'Geometry' in i]
    for file in files:
        with open(file, 'r') as f:
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
        shutil.move(file, 'Geometries/'+file)
        shutil.move(file[:-3]+'com', 'Geometries/'+file[:-3]+'com')
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
    dist = max(abs(vec1/vec2 -1) - cr)
    distance =  np.heaviside(dist,0)   #np.heaviside(1e3*np.sqrt(np.sum((vec1 - vec2)**2))-1,1) + np.heaviside(abs(e1-e2)- 0.01, 1)
    return distance

##CLASSIFIES THE VARIOUS OPTIMIZED STRUCTURES##################
def classify(nums,scfs,rotsx, rotsy, rotsz):
    try:
        data = np.loadtxt('conformation.lx')
        if len(np.shape(data)) > 1:
            engs =  data[:,1].flatten()
            rotx =  data[:,4].flatten()
            roty =  data[:,5].flatten() 
            rotz =  data[:,6].flatten()
            last =  data[:,10].flatten()
            crix =  data[:,11].flatten()
            criy =  data[:,12].flatten()
            criz =  data[:,13].flatten()
            exam =  data[:,14].flatten()
        else:
            engs = np.array([data[1]])
            rotx = np.array([data[4]])
            roty = np.array([data[5]])
            rotz = np.array([data[6]]) 
            last = np.array([data[10]])
            crix = np.array([data[11]])
            criy = np.array([data[12]])
            criz = np.array([data[13]])
            exam = np.array([data[14]])
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
                cr = [max(crix[n],0.01), max(criy[n],0.01), max(criz[n],0.01)]
            except:
                cr = [0.01, 0.01, 0.01]    
            distance = measure([rotsx[m], rotsy[m], rotsz[m]], [ROTX[n],ROTY[n],ROTZ[n]],cr)
            distances.append(distance)
        try:
            a = distances.index(0)
            exam[a] = nums[m]
            engs[a] = (last[a]*engs[a] + scfs[m])/(last[a]+1)
            rotx[a] = (last[a]*rotx[a] + rotsx[m])/(last[a]+1)
            roty[a] = (last[a]*roty[a] + rotsy[m])/(last[a]+1)
            rotz[a] = (last[a]*rotz[a] + rotsz[m])/(last[a]+1)    
            crix[a] = np.sqrt((last[a]*(crix[a]**2 +rotx[a]**2) + rotsx[m]**2)/(last[a]+1))
            criy[a] = np.sqrt((last[a]*(criy[a]**2 +roty[a]**2) + rotsy[m]**2)/(last[a]+1))
            criz[a] = np.sqrt((last[a]*(criz[a]**2 +rotz[a]**2) + rotsz[m]**2)/(last[a]+1))
            last[a] += 1 #nums[m]
        except:
            new.append(nums[m])
            rotx  = np.append(rotx,rotsx[m])
            roty  = np.append(roty,rotsy[m])
            rotz  = np.append(rotz,rotsz[m])
            crix  = np.append(crix,0.01)
            criy  = np.append(criy,0.01)
            criz  = np.append(criz,0.01) 
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

    with open('conformation.lx', 'w') as f: 
        f.write('{:6}\t{:10}\t{:10}\t{:12}\t{:10}\t{:10}\t{:10}\t{:10}\t{:10}\t{:10}\t{:10}\t{:10}'.format('#Group','Energy(eV)','DeltaE(eV)','Prob@300K(%)','Rot1','Rot2','Rot3','Std1','Std2','Std3','Number','Last\n'))
        for i in range(len(probs)):
            f.write('{:6}\t{:<10.3f}\t{:<10.3f}\t{:<12.1f}\t{:<10.7f}\t{:<10.7f}\t{:<10.7f}\t{:<10.3f}\t{:<10.3f}\t{:<10.3f}\t{:<10.0f}\t{:<10.0f}\n'.format(i+1,engs[i],engs[i] -min(engs),100*probs[i],rotx[i],roty[i],rotz[i], crix[i], criy[i], criz[i], last[i],exam[i]))
    return int(origin), last
###############################################################

##RUNS FREQ CALCULATION FOR NEW CONFORMATION###################
def rodar_freq(origin,nproc,mem,base,cm,batch_file):
    geomlog = 'Geometries/Geometry-'+str(origin)+'-.log'
    G, atomos = pega_geom(geomlog) 
    header = "%nproc={}\n%mem={}\n# freq=(noraman) nosymm  {} \n\n{}\n\n{}\n".format(nproc,mem,base,'ABSSPCT',cm)
    file = "Freq-"+str(origin)+"-.com"
    write_input(atomos,G,header,'',file)
    subprocess.call(['bash', batch_file, file]) 
    hold_watch([file])
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
    
    freq0 = freqlog
    try:
        os.mkdir('Geometries')
    except:
        pass    
    
    cm        = get_cm(freqlog) 
    header    = "%nproc={}\n%mem={}\n# opt nosymm  {} \n\n{}\n\n{}\n".format(nproc,mem,base,'ABSSPCT',cm)
    scf, rotx, roty, rotz  = get_energy_origin(freqlog)
    origin, conformation = classify(np.array([0]),np.array([scf]), np.array([rotx]),np.array([roty]),np.array([rotz]))
    files = [i for i in os.listdir('.') if 'Geometry' in i and '.log' in i]
    if len(files) > 0:
        nums, scfs, rotsx, rotsy, rotsz = get_energies()
        origin, conformation  = classify(nums,scfs,rotsx,rotsy,rotsz)
    else:
        pass    

    T0 = T

    for i in range(rounds):
        lista      = make_geoms(freqlog, num_geoms, T0, header, '')
        rodar_opts(lista,script)
        nums, scfs, rotsx,rotsy,rotsz = get_energies()
        origin, conformation  = classify(nums,scfs,rotsx,rotsy,rotsz)
        with open('conformation.lx', 'a') as f:
            f.write('\n#Round {}/{} Temperature: {} K'.format(i+1,rounds,T0))    
        if origin != 0:
            log  = rodar_freq(origin,nproc,mem,base,cm,script)
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
        _, _, nproc, mem, scrf, _ = busca_input(freqlog)
        cm = get_cm(freqlog)
        header = '%nproc={}\n%mem={}\n%chk=Group_{}_.chk\n# {} {} opt\n\nTITLE\n\n{}\n'.format(nproc,mem,i+1,'pm6',scrf,cm)
        G, atomos = pega_geom(freqlog)
        write_input(atomos,G,header,'','Conformers/Group_{}_.com'.format(i+1))
    



if __name__ == "__main__":
    sys.exit(main())        


        
        
        


    
        


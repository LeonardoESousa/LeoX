#!/usr/bin/env python3
import numpy as np
import os
import sys
import subprocess
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
                rot  = np.array([float(line[3]),float(line[4]),float(line[5])])    
            elif 'Normal termination' in line:
                return scf, rot
###############################################################

##GETS ENERGIES FROM OPT LOG FILES#############################
def get_energies(baseline):
    nums,scfs,rots = [], [], []
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
                    rot  = 1e3*np.round(np.sqrt(np.sum((np.array([float(line[3]),float(line[4]),float(line[5])]) - baseline)**2)),3)
                elif 'Normal termination' in line:
                    scfs.append(scf)
                    nums.append(num)
                    rots.append(rot)
    for file in files:
        shutil.move(file, 'Geometries/'+file)
        shutil.move(file[:-3]+'com', 'Geometries/'+file[:-3]+'com')
    nums = np.array(nums)
    scfs = np.array(scfs)
    rots = np.array(rots)
    return nums, scfs, rots
###############################################################

##CLASSIFIES THE VARIOUS OPTIMIZED STRUCTURES##################
def classify(nums,scfs,rots):
    try:
        data = np.loadtxt('conformation.lx')
        if len(np.shape(data)) > 1:
            old_rots = data[:,4] 
            nums     = np.append(data[:,5],nums)
            scfs     = np.append(data[:,1],scfs)
            rots     = np.append(data[:,4],rots)
        else:
            old_rots = np.array(data[4]) 
            nums     = np.append(data[5],nums)
            scfs     = np.append(data[1],scfs)
            rots     = np.append(data[4],rots)
    except:
        old_rots = [] 
    
    unrots = np.unique(rots)
    new = []
    for elem in unrots:
        if elem not in old_rots:
            new.append(elem)

    try:
        pos = np.where(rots == max(new))[0][0]
        origin = nums[pos]
    except:
        origin = 0

    groups = np.unique(rots)
    conformation = [[] for _ in groups]
    engs  = np.zeros(len(groups))
    for i in range(len(nums)):
        indice = np.where(groups == rots[i])[0][0]
        conformation[indice].append(str(int(nums[i])))
        engs[indice]  = scfs[i]

    probs  = np.exp(-1*(engs - min(engs))/0.026)
    probs  /= np.sum(probs)
    args   = np.argsort(engs)
    engs   = engs[args]
    probs  = probs[args]
    groups = groups[args]
    conformation = [conformation[i] for i in args]

    with open('conformation.lx', 'w') as f: 
        f.write('#Group    Energy(eV)    DeltaE(eV)    Prob@300K(%)    ObjFunction    First\n')
        for i in range(len(probs)):
            f.write('{:5}     {:<10.3f}    {:<10.3f}    {:<5.1f}           {:<5.1f}           {:5}\n'.format(i+1,engs[i],engs[i] -min(engs),100*probs[i],groups[i],conformation[i][-1]))
    return int(origin), conformation
###############################################################

##RUNS FREQ CALCULATION FOR NEW CONFORMATION###################
def rodar_freq(origin,nproc,mem,base,cm,batch_file):
    geomlog = 'Geometries/Geometry-'+str(origin)+'-.log'
    G, atomos = pega_geom(geomlog) 
    header = "%nproc={}\n%mem={}\n# freq=(noraman)  {} \n\n{}\n\n{}\n".format(nproc,mem,base,'ABSSPCT',cm)
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
    num_geoms = int(sys.argv[6])
    rounds    = int(sys.argv[7])
    script    = sys.argv[8]
    
    freq0 = freqlog
    try:
        os.mkdir('Geometries')
    except:
        pass    
    
    cm        = get_cm(freqlog) 
    header    = "%nproc={}\n%mem={}\n# opt  {} \n\n{}\n\n{}\n".format(nproc,mem,base,'ABSSPCT',cm)
    scf, baseline  = get_energy_origin(freqlog)
    origin, conformation = classify(np.array([0]),np.array([scf]), np.array([0]))
    try:
        nums, scfs, rots = get_energies(baseline)
        origin, conformation  = classify(nums,scfs,rots)
    except:
        pass    

    T0 = T

    for i in range(rounds):
        lista      = make_geoms(freqlog, num_geoms, T0, header, '')
        rodar_opts(lista,script)
        nums, scfs, rots = get_energies(baseline)
        origin, conformation  = classify(nums,scfs,rots)
        with open('conformation.lx', 'a') as f:
            f.write('\n#Round {}/{} Temperature: {} K'.format(i+1,rounds,T0))    
        if origin != 0:
            log  = rodar_freq(origin,nproc,mem,base,cm,script)
            if log != None:
                freqlog = log
                T0 = T
        else:
            T0 += 100

    with open('conformation.lx', 'a') as f:
        f.write('\n#Search concluded!')

    try:
        os.mkdir('Conformers')
    except:
        pass    

    for i in range(len(conformation)):
        numero  = conformation[i][0]
        if numero == '0':
            freqlog = freq0    
        else:
            freqlog = 'Geometries/Geometry-{}-.log'.format(numero) 
        _, _, nproc, mem, scrf, _ = busca_input(freqlog)
        cm = get_cm(freqlog)
        header = '%nproc={}\n%mem={}\n%chk=Group_{}_.chk\n# {} {} opt\n\nTITLE\n\n{}\n'.format(nproc,mem,i+1,'pm6',scrf,cm)
        G, atomos = pega_geom(freqlog)
        write_input(atomos,G,header,'','Conformers/Group_{}_.com'.format(i+1))
    



if __name__ == "__main__":
    sys.exit(main())        


        
        
        


    
        


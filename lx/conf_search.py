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
            with open('omega.lx','a') as f:
                f.write('#Aborted!')
            sys.exit()
        time.sleep(60)    
###############################################################

##RUNS CALCULATIONS############################################
def rodar_opts(lista, batch_file): 
    for file in lista:
        subprocess.call(['bash', batch_file, file]) 
    hold_watch(lista)
###############################################################

##CHECKS FOR EXISTING GEOMETRIES###############################
def my_counter():
    files = [file for file in os.listdir('.') if ".com" in file and "Geometry" in file]
    return len(files)
###############################################################


##SAMPLES GEOMETRIES###########################################
def make_geoms(freqlog, num_geoms, T, header, bottom):
    lista = []
    F, M = pega_freq(freqlog)
    F[F < 0] *= -1    
    counter = my_counter()
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

def get_energy_origin(freqlog):
    with open(freqlog, 'r') as f:
        for line in f:
            if 'SCF Done:' in line:
                line = line.split()
                scf  = float(line[4])*27.2114
            elif 'Normal termination' in line:
                return scf
    

def get_energies():
    nums,scfs = [], []
    files = [i for i in os.listdir('.') if '.log' in i and 'Geometry' in i]
    for file in files:
        with open(file, 'r') as f:
            num = float(file.split('-')[1])
            for line in f:
                if 'SCF Done:' in line:
                    line = line.split()
                    scf  = float(line[4])*27.2114
                elif 'Normal termination' in line:
                    scfs.append(scf)
                    nums.append(num)
    for file in files:
        shutil.move(file, 'Geometries/'+file)
        shutil.move(file[:-3]+'com', 'Geometries/'+file[:-3]+'com')
    nums = np.array(nums)
    scfs = np.round(np.array(scfs),1)
    return nums, scfs


def classify(nums,scfs):
    try:
        data = np.loadtxt('conformation.lx')
        old_engs = data[:,1] 
        nums     = np.append(data[:,4],nums)
        scfs     = np.append(data[:,1],scfs)
    except:
        old_engs = [] 
        

    engs = np.unique(scfs)
    new = []
    for elem in engs:
        if elem not in old_engs:
            new.append(elem)

    try:
        pos = np.where(scfs == max(new))[0][0]
        origin = nums[pos]
    except:
        origin = 0

    scfs -= min(scfs)
    scfs = np.round(scfs,1)
    groups = np.unique(scfs)
    boltz = np.exp(-1*groups/0.026)
    total  = np.sum(boltz)
    conformation = [[] for _ in groups]
    probs = 100*boltz/total
    for i in range(len(nums)):
        indice = np.where(groups == scfs[i])[0][0]
        conformation[indice].append(str(int(nums[i])))

    with open('conformation.lx', 'w') as f: 
        f.write('#Group    Energy(eV)    DeltaE(eV)    Prob@300K(%)    First\n')
        for i in range(len(probs)):
            f.write('{:5}     {:<10.1f}    {:<10.1f}    {:<5.1f}       {:5}\n'.format(i+1,engs[i],groups[i],probs[i],conformation[i][0]))
            #f.write('#Geometries: {}\n\n'.format(', '.join(conformation[i])))
    return int(origin)

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



def main():
    freqlog   = sys.argv[1]
    base      = sys.argv[2]
    nproc     = sys.argv[3]
    mem       = sys.argv[4]
    T         = float(sys.argv[5])
    num_geoms = int(sys.argv[6])
    rounds    = int(sys.argv[7])
    script    = sys.argv[8]
    
    try:
        os.mkdir('Geometries')
    except:
        pass    
    
    cm = get_cm(freqlog) 
    header = "%nproc={}\n%mem={}\n# opt  {} \n\n{}\n\n{}\n".format(nproc,mem,base,'ABSSPCT',cm)
    scf = get_energy_origin(freqlog)
    origin = classify(np.array([0]),np.array([scf]))
    try:
        nums, scfs = get_energies()
        origin = classify(nums,scfs)
    except:
        pass    

    T0 = T

    for _ in range(rounds):
        lista      = make_geoms(freqlog, num_geoms, T0, header, '')
        rodar_opts(lista,script)
        nums, scfs = get_energies()
        origin = classify(nums,scfs)
        if origin != 0:
            log = rodar_freq(origin,nproc,mem,base,cm,script)
            if log != None:
                freqlog = log
                T0 = T
        else:
            T0 += 100




if __name__ == "__main__":
    sys.exit(main())        


        
        
        


    
        


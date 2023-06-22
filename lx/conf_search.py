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
    _, atomos, A = lx.tools.sample_geometries(freqlog,num_geoms,T,3000,warning=False)
    for n in range(0,np.shape(A)[1],3):
        Gfinal = A[:,n:n+3]
        lx.tools.write_input(atomos,Gfinal,header.replace("UUUUU",str((n+3)//3)),bottom.replace("UUUUU",str((n+3)//3)),"Geometry-"+str((n+3)//3+counter)+"-.com")
        lista.append("Geometry-"+str((n+3)//3+counter)+"-.com") 
    return lista      
############################################################### 

##GETS ENERGY FROM THE ORIGINAL FREQ LOG FILE##################
def get_energy_origin(freqlog):
    exc = 0
    with open(freqlog, 'r') as f:
        for line in f:
            if 'SCF Done:' in line:
                line = line.split()
                scf  = float(line[4])*27.2114
            elif 'Total Energy,' in line:
                line = line.split()
                exc  = float(line[4])*27.2114    
            elif 'Rotational constants' in line:
                line = line.split()
                rot  = [float(line[3]),float(line[4]),float(line[5])]    
            elif 'Normal termination' in line:
                if exc != 0:
                    scf = exc
                return scf, np.array([rot[0], rot[1], rot[2]])
###############################################################

##GETS ENERGIES FROM OPT LOG FILES#############################
def get_energies(folder,original_molecule):
    nums,scfs,rots = [], [], []
    files = [i for i in os.listdir(folder) if '.log' in i and 'Geometry' in i]
    for file in files:
        exc = 0
        if np.array_equal(original_molecule, fingerprint(file,folder)):
            with open(folder+'/'+file, 'r') as f:
                num = float(file.split('-')[1])
                for line in f:
                    if 'SCF Done:' in line:
                        line = line.split()
                        scf  = float(line[4])*27.2114
                    elif 'Total Energy,' in line:
                        line = line.split()
                        exc  = float(line[4])*27.2114    
                    elif 'Rotational constants' in line:
                        line = line.split()
                        rot  = [float(line[3]),float(line[4]),float(line[5])]
                    elif 'Normal termination' in line:
                        if exc != 0:
                            scf = exc
                        scfs.append(scf)
                        nums.append(num)
                        rots.append(rot)
                        
    for file in files:
        try:
            shutil.move(file, 'Geometries/'+file)
            shutil.move(file[:-3]+'com', 'Geometries/'+file[:-3]+'com')
        except:
            pass
    return np.array(nums), np.array(scfs), np.array(rots)
###############################################################

def measure(vec1,vec2,cr):
    vec1 = np.array(vec1)     
    vec2 = np.array(vec2)
    cr   = np.array(cr)
    dist = max(abs(vec1 -vec2) - cr)
    distance =  np.heaviside(dist,0)   
    return distance

class Conformation:
    def __init__(self,rot,energy,identity,num) -> None:
        self.rot = rot[np.newaxis,:]
        self.energy = [energy]
        self.identity = identity
        self.std = rot/1000
        self.num = [num]

    def add_rot(self,rot,num,energy):
        self.rot = np.vstack((self.rot,rot[np.newaxis,:]))
        newstd = np.std(self.rot,axis=0)    
        #get higher std
        self.std = np.where(newstd > self.std, newstd, self.std)
        self.num.append(num)
        self.energy.append(energy)

    def merge_rot(self,rot):
        self.rot = np.vstack((self.rot,rot))
        newstd = np.std(self.rot,axis=0)
        #get higher std
        self.std = np.where(newstd > self.std, newstd, self.std)


    def get_avg(self):
        return np.mean(self.rot,axis=0)  



def internal_comparison(conformations):
    remove = []
    for i in range(len(conformations)-1):
        for j in range(i+1,len(conformations)):
            distance = measure(conformations[i].get_avg(),conformations[j].get_avg(),conformations[j].std+conformations[i].std)
            if distance == 0:
                conformations[i].num += conformations[j].num
                conformations[i].energy += conformations[j].energy
                conformations[i].merge_rot(conformations[j].rot)
                remove.append(j)
                break
    # remove elements from list whose index is in remove
    conformations = [i for j, i in enumerate(conformations) if j not in remove]
    return conformations        
     

def classify(conformations,folder):
    nums, scfs, rots = get_energies(folder,conformations[0].identity)
    for i in range(rots.shape[0]):
        for conformation in conformations:
            distance = measure(rots[i,:],conformation.get_avg(),conformation.std)
            if distance == 0:
                conformation.add_rot(rots[i,:],nums[i],scfs[i])
                break
            conformations.append(Conformation(rots[i,:],scfs[i],conformations[0].identity,nums[i])) 
    if len(conformations) > 1:
        conformations  = internal_comparison(conformations)           
    return conformations

def write_report(conformations,round,total_rounds,temp):
    engs = []
    for conformation in conformations:
        engs.append(np.mean(conformation.energy))
    engs = np.array(engs)    
    argsort = np.argsort(engs)
    engs = engs[argsort]
    #sort conformations as list
    conformations = [conformations[i] for i in argsort]
    deltae = engs - min(engs)
    probs = np.exp(-(deltae)/(0.026))
    probs = probs/sum(probs)

    with open('conformation.lx','w') as f:
        f.write('{:6}  {:10}  {:10}  {:12}  {:10}  {:10}  {:10}  {:10}  {:10}  {:10}  {:6}  {:6}\n'.format('#Group','Energy(eV)','DeltaE(eV)','Prob@300K(%)','Rot1','Rot2','Rot3','Std1','Std2','Std3','Number','Last'))
        for i in range(len(conformations)):
            rot = conformations[i].get_avg()
            std = conformations[i].std
            last = conformations[i].num[-1]
            total= len(conformations[i].num)
            f.write(f'{i+1:<6}  {engs[i]:<10.3f}  {deltae[i]:<10.3f}  {100*probs[i]:<12.1f}  {rot[0]:<10.7f}  {rot[1]:<10.7f}  {rot[2]:<10.7f}  {std[0]:<10.7f}  {std[1]:<10.7f}  {std[2]:<10.7f}  {total:<6.0f}  {last:<6.0f}\n')
        f.write(f'\n#Round {round}/{total_rounds} Temperature: {temp} K')    


##RUNS FREQ CALCULATION FOR NEW CONFORMATION###################
def rodar_freq(origin,nproc,mem,base,cm,batch_file,gaussian):
    geomlog = f'Geometries/Geometry-{origin:.0f}-.log'
    G, atomos = lx.tools.pega_geom(geomlog) 
    header = "%nproc={}\n%mem={}\n# freq=(noraman) nosymm  {} \n\n{}\n\n{}\n".format(nproc,mem,base,'ABSSPCT',cm)
    file = f"Freq-{origin:.0f}-.com"
    lx.tools.write_input(atomos,G,header,'',file)
    lx.tools.rodar_lista([file],batch_file,gaussian,'conformation.lx',1)
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
    numjobs   = int(sys.argv[9])
    script    = sys.argv[10]
    gaussian  = sys.argv[11]
    T0 = T
    freq0 = freqlog
    if 'td' in base.lower():
        opt = '=loose'
    else:
        opt = ''    
    
    try:
        os.mkdir('Geometries')
    except:
        pass    
    cm        = lx.tools.get_cm(freqlog) 
    header    = f"%nproc={nproc}\n%mem={mem}\n# opt{opt} nosymm  {base} \n\nTitle\n\n{cm}\n"
    scf, rot = get_energy_origin(freqlog)
    conformations = [Conformation(rot,scf,fingerprint(freqlog,'.'),0)]
    files = [i for i in os.listdir('Geometries') if 'Geometry' in i and '.log' in i]
    if len(files) > 0:
        conformations = classify(conformations,'Geometries')
        write_report(conformations,0,rounds,T0)
    else:
        pass    

    
    groups = len(conformations)
    for i in range(rounds):
        lista      = make_geoms(freqlog, num_geoms, T0, header, '')
        lx.tools.rodar_lista(lista,script, gaussian, 'conformation.lx',numjobs)
        conformations  = classify(conformations,'.')
        write_report(conformations,i+1,rounds,T0)
        
        if len(conformations) != groups:
            log  = rodar_freq(conformations[-1].num[-1],nproc,mem,base,cm,script,gaussian)
            if log != None:
                freqlog = log
                T0 = T
            groups = len(conformations)    
        else:
            T0 += DT

    with open('conformation.lx', 'a') as f:
        f.write('\n#Search concluded!')

    try:
        os.mkdir('Conformers')
    except:
        pass    

    for i in range(len(conformations)):
        numero  = conformations[i].num[-1]
        if numero == 0:
            freqlog = freq0    
        else:
            freqlog = 'Geometries/Geometry-{:.0f}-.log'.format(numero) 
        _, _, nproc, mem, scrf, _ = lx.tools.busca_input(freqlog)
        cm = lx.tools.get_cm(freqlog)
        header = '%nproc={}\n%mem={}\n%chk=Group_{}_.chk\n# {} {} opt\n\nTITLE\n\n{}\n'.format(nproc,mem,i+1,'pm6',scrf,cm)
        G, atomos = lx.tools.pega_geom(freqlog)
        lx.tools.write_input(atomos,G,header,'','Conformers/Geometry-{}-.com'.format(i+1))
    



if __name__ == "__main__":
    sys.exit(main())        


        
        
        


    
        


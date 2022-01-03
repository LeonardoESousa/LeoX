#!/usr/bin/env python3
import numpy as np
import os
import sys
import subprocess
from lx.tools import *
import shutil
import time

##GENERATES NEUTRAL INPUT######################################
def gera_optcom(atomos,G,base,nproc,mem,omega,op):
    header = "%nproc=JJJ\n%mem=MEM\n# BASE iop(3/108=MMMMM00000) iop(3/107=MMMMM00000) {}\n\nTITLE\n\n0 1\n".format(op)
    header = header.replace("JJJ",nproc).replace("MEM", mem).replace("BASE", base).replace("MMMMM",omega)
    file = "OPT_"+omega+"_.com"
    write_input(atomos,G,header,'',file)
    return file
###############################################################

##GENERATES ION INPUTS#########################################
def gera_ioncom(atomos,G,base,nproc,mem,omega):
    header = "%nproc=JJJ\n%mem=MEM\n# BASE iop(3/108=MMMMM00000) iop(3/107=MMMMM00000)\n\nTITLE\n\n1 2\n"
    header = header.replace("JJJ",nproc).replace("MEM", mem).replace("BASE", base).replace("MMMMM",omega)
    file1 = "pos_"+omega+"_.com"
    file2 = "neg_"+omega+"_.com"
    write_input(atomos,G,header,'',file1)
    header = "%nproc=JJJ\n%mem=MEM\n# BASE iop(3/108=MMMMM00000) iop(3/107=MMMMM00000) \n\nTITLE\n\n-1 2\n"
    header = header.replace("JJJ",nproc).replace("MEM", mem).replace("BASE", base).replace("MMMMM",omega)
    write_input(atomos,G,header,'',file2)
    return [file1,file2]
###############################################################

##GETS ENERGY FROM LOG#########################################
def pega_energia(file):
    energias = []
    with open(file, 'r') as f:
        for line in f:
            if "SCF Done:" in line:
                line = line.split()
                energias.append(float(line[4]))
    return min(energias)
###############################################################

##MAP TO FLOAT WITH EXCEPTION##################################
def to_float(num):
    try:
        num = float(num)
    except:
        num = np.nan
    return num        
###############################################################

##GETS HOMO ENERGY FROM LOG####################################
def pega_homo(file):
    HOMOS = []
    with open(file, 'r') as f:
        for line in f:
            if ('OPT' in file and "Optimized Parameters" in line) :
                HOMOS = [] 
            if "occ. eigenvalues" in line:
                line = line.split()
                homos = line[4:]
                HOMOS.extend(homos)
    if len(HOMOS) == 0:
        with open(file, 'r') as f:
            for line in f:
                if "occ. eigenvalues" in line:
                    line = line.split()
                    homos = line[4:]
                    HOMOS.extend(homos)
    #HOMOS = list(map(to_float,HOMOS))
    return float(HOMOS[-1])  #max(HOMOS)
###############################################################

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
def rodar_omega(atomos,G,base,nproc,mem,omega,op,batch_file): 
    omega = "{:05.0f}".format(omega)
    file  = gera_optcom(atomos,G,base,nproc,mem,omega,op)
    remover = [file]
    rodando = watcher([file],1)
    for file in rodando:
        subprocess.call(['bash', batch_file, file]) 
    hold_watch(rodando)
    if op == 'opt':
        G, atomos = pega_geom(file[:-3]+"log")
    neutro      = pega_energia(file[:-3]+"log")
    homo_neutro = pega_homo(file[:-3]+"log")  
    files = gera_ioncom(atomos,G,base,nproc,mem,omega)
    rodando = watcher(files,1)
    for file in rodando:
        subprocess.call(['bash', batch_file, file]) 
    hold_watch(rodando)
    for file in files:
        remover.append(file)
        if "pos" in file:
            cation     = pega_energia(file[:-3]+"log")
        elif "neg" in file:
            anion      = pega_energia(file[:-3]+"log")
            homo_anion = pega_homo(file[:-3]+"log")        
    try:
        os.mkdir("Logs")
    except:
        pass
    for file in remover:
        shutil.move(file, 'Logs/'+file)
        shutil.move(file[:-3]+'log', 'Logs/'+file[:-3]+'log')
    J = np.sqrt(((homo_neutro + cation - neutro)**2 + (homo_anion + neutro - anion)**2))*(27.2114)
    return J, G, atomos
###############################################################

##WRITES LOG WITH RESULTS######################################
def write_tolog(omegas,Js,frase):    
    with open("omega.lx", 'w') as f:
        f.write('#{}    {}\n'.format('w(10^4 bohr^-1)','J(eV)'))
        list1, list2 = zip(*sorted(zip(omegas, Js)))
        for i in range(len(list1)):
            f.write("{:05.0f}               {:.4f}\n".format(list1[i],list2[i])) 
        f.write("\n{} {:05.0f}\n".format(frase,list1[list2.index(min(list2))])) 
###############################################################

def main():
    geomlog = sys.argv[1]
    base    = sys.argv[2]
    nproc   = sys.argv[3]
    mem     = sys.argv[4]
    omega1  = sys.argv[5]
    passo   = sys.argv[6]
    relax   = sys.argv[7]
    script  = sys.argv[8]

    try:
        int(nproc)
        passo  = float(passo)*10000
        omega1 = float(omega1)*10000
    except:
        fatal_error('nproc, omega and step must be numbers. Goodbye!')
    if relax.lower() == 'y':
        op = 'opt'
    elif relax.lower() == 'n':
        op = ''
    else:
        fatal_error('It must be either y or n. Goodbye!')    

    omegas, Js = [], []
    oms, jotas = [], []
    try:
        with open("omega.lx", 'r') as f:
            for line in f:
                line = line.split()
                if len(line) == 2 and '#' not in line:
                    om = float(line[0])
                    omegas.append(om)
                    Js.append(float(line[1]))
        menor = omegas[Js.index(min(Js))]
        G, atomos = pega_geom('Logs/OPT_{:05.0f}_.log'.format(menor))            
    except:
        G, atomos  = pega_geom(geomlog)

    while passo > 25:
        if omega1 in omegas:
            ind = omegas.index(omega1)
            J = Js[ind]
        else:
            J, G, atomos = rodar_omega(atomos,G,base,nproc,mem,omega1,op,script)              
            omegas.append(omega1)
            Js.append(J)  
        oms.append(omega1)
        jotas.append(J)    
        try:
            if jotas[-1] - jotas[-2] > 0:     
                passo = int(passo/2)
                sign = -1*np.sign(int(oms[-1]) - int(oms[-2])) 
            else:           
                sign = +1*np.sign(int(oms[-1]) - int(oms[-2]))
            omega1 += sign*passo 

        except:
            omega1 += passo
        write_tolog(omegas,Js,'#Best value so far:')


    write_tolog(omegas,Js,'#Done! Optimized value:')
    menor = omegas[Js.index(min(Js))]
    log = 'Logs/OPT_{:05.0f}_.log'.format(menor)
    G, atomos = pega_geom(log)    
    base, _, nproc, mem, scrf, _ = busca_input(log)
    cm = get_cm(log)
    header = '%nproc={}\n%mem={}\n# {} {}\n\nTITLE\n\n{}\n'.format(nproc,mem,base,scrf,cm)
    write_input(atomos,G, header,'', 'tuned_w.com')



if __name__ == "__main__":
    sys.exit(main())        


        
        
        


    
        


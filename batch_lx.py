import os
import time
import subprocess
import numpy as np
import sys
import shutil

def watcher(rodando,counter):
    done = []
    for input in rodando: 
        term = 0
        try:
            with open('Geometry-'+str(input)+'-.log', 'r') as f:
                for line in f:
                    if 'Normal termination' in line:
                        term += 1
            if term == counter:
                done.append(input)
        except:
            pass        
    for elem in done:
        del rodando[rodando.index(elem)]                                
    return rodando

def set_factor(file):
    factor = 1
    with open(file, 'r') as f:
        for line in f:
            if 'Link1' in line:
                factor = 2
    return factor

def limite():
    numero = np.loadtxt('../limit.lx')
    return numero

try:
    batch_file = sys.argv[1]
    shutil.copy(batch_file,'Geometries')

    os.chdir('Geometries')
    inputs = [i for i in os.listdir('.') if 'Geometr' in i and '.com' in i]
    inputs = sorted(inputs, key=lambda pair: float(pair.split('-')[1]))

    factor = set_factor(inputs[0])

    inputs = watcher(inputs,factor)

    rodando = []
    for i in range(len(inputs)):
        rodando = watcher(rodando,factor)
        nlim = limite()
        a = subprocess.call(['bash',batch_file, inputs[i]])
        rodando.append(i)
        while len(rodando) >= nlim:
            time.sleep(20)
            rodando = watcher(rodando,factor)
            nlim = limite()
except:
    print('Something went wrong! Abort.')           
    sys.exit()
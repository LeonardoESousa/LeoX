#!/usr/bin/env python3
import os
import time
import subprocess
import numpy as np
import sys
import shutil
import lx.tools         

##CHECKS WHETHER THE JOB IS TWO STEP###########################
def set_factor(file):
    factor = 1
    with open(file, 'r') as f:
        for line in f:
            if 'Link1' in line:
                factor = 2
                break
    return factor
###############################################################

#SETS NUMBER OF SIMULTANEOUS JOBS##############################
def limite():
    numero = np.loadtxt('../limit.lx')
    return numero
###############################################################

##MAIN LOOP####################################################
try:
    batch_file = sys.argv[1]
    num        = int(sys.argv[2])
    command    = sys.argv[3]
    scripts    = [i for i in os.listdir('.') if '.sh' in i]
    for file in scripts:
        shutil.copy(file,'Geometries')
    os.chdir('Geometries')
    inputs = [i for i in os.listdir('.') if 'Geometr' in i and '.com' in i]
    inputs = sorted(inputs, key=lambda pair: float(pair.split('-')[1]))
    factor = set_factor(inputs[0])
    inputs = lx.tools.watcher(inputs,factor)
    if len(inputs) == 0:
        print('No jobs left to run! Goodbye!')
        sys.exit()
    rodando = []
    queue, batch_num = 0, 0
    newcommand = ''
    leftover = len(inputs)%num
    for input in inputs:
        rodando = lx.tools.watcher(rodando,factor)
        nlim    = limite()
        newcommand += '{} {} \n'.format(command, input)
        queue += 1
        if queue == num or (queue == leftover and batch_num >= len(inputs) - leftover):
            newcommand += 'wait'
            with open('cmd_{}_.sh'.format(batch_num), 'w') as q:
                q.write(newcommand)
            a = subprocess.call(['bash',batch_file, 'cmd_{}_.sh'.format(batch_num)])
            queue = 0
            newcommand = ''
        rodando.append(input)
        batch_num += 1
        while len(rodando)/num >= nlim:
            time.sleep(20)
            rodando = lx.tools.watcher(rodando,factor)
            nlim = limite()
except:
    print('Something went wrong! Abort.')           
###############################################################
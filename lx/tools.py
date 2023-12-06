#!/usr/bin/env python3
import os
import sys
import time
import subprocess
from scipy.stats import norm
import numpy as np
import pandas as pd
from lx.ld import run_ld
import lx.parser

##SOME CONSTANTS##############################################
EPSILON_0 =   lx.parser.EPSILON_0    #F/m
HBAR_EV =     lx.parser.HBAR_EV      #eV s
HBAR_J =      lx.parser.HBAR_J       #J s
MASS_E =      lx.parser.MASS_E       #kg
LIGHT_SPEED = lx.parser.LIGHT_SPEED  #m/s
E_CHARGE =    lx.parser.E_CHARGE     #C
BOLTZ_EV =    lx.parser.BOLTZ_EV     #eV/K
AMU =         lx.parser.AMU          #kg
###############################################################

def distance_matrix(geom):
    matrix = np.zeros((1,np.shape(geom)[0]))
    for ind in range(np.shape(geom)[0]):
        distances = geom - geom[ind,:]
        distances = np.sqrt(np.sum(np.square(distances),axis=1))
        matrix = np.vstack((matrix,distances[np.newaxis,:]))
    matrix = matrix[1:,:]
    return matrix

def adjacency(geom,atoms):
    covalent_radii = {1:0.31, 2:0.28, 3:1.28, 4:0.96, 5:0.84, 6:0.76, 7:0.71, 8:0.66, 9:0.57, 10:0.58, 11:1.66, 12:1.41, 13:1.21, 14:1.11,
           15:1.07, 16:1.05, 17:1.02, 18:1.06, 19:2.03, 20:1.76, 21:1.7, 22:1.6, 23:1.53, 24:1.39, 25:1.61, 26:1.52, 27:1.50,
           28:1.24, 29:1.32, 30:1.22, 31:1.22, 32:1.2, 33:1.19, 34:1.20, 35:1.20, 36:1.16}
    atoms = np.array(atoms).astype(float)
    dist_matrix = distance_matrix(geom)
    adj_matrix = np.zeros(np.shape(dist_matrix))
    #connectivity matrix
    for i in range(np.shape(dist_matrix)[0]):
        for j in range(i, np.shape(dist_matrix)[1]):
            r_e = (covalent_radii[atoms[i]] + covalent_radii[atoms[j]]) + 0.4
            if 0.8 < dist_matrix[i,j] < r_e:
                adj_matrix[i,j] = adj_matrix[j,i] = 1
    return adj_matrix


def fingerprint(file,folder):
    geom, atoms = lx.parser.pega_geom(folder+'/'+file)
    cm = adjacency(geom,atoms)
    return cm

##SAVES OPT GEOMETRY###########################################
def salva_geom(G,atomos):
    atomos = np.array([atomos]).astype(float)
    atomos = atomos.T
    G = np.hstack((atomos,G))
    np.savetxt('opt_geom.lx', G, delimiter='\t', fmt=['%1.1u','%+1.5f','%+1.5f','%+1.5f'])
    print("The optimized geometry that is used is saved in the opt_geom.lx file!")
###############################################################

##WRITES ATOMS AND XYZ COORDS TO FILE##########################
def write_input(atomos,geom,header,bottom,file):
    with open(file, 'w', encoding='utf-8') as f:
        f.write(header)
        for i, atomo in enumerate(atomos):
            texto = f"{atomo:2s}  {geom[i,0]:.6f}  {geom[i,1]:.6f}  {geom[i,2]:.6f}\n"
            f.write(texto)
        f.write("\n"+bottom+'\n')
###############################################################

##CHECKS FOR EXISTING GEOMETRIES###############################
def start_counter():
    files = [file for file in os.listdir('Geometries') if ".com" in file and "Geometr" in file]
    return len(files)
###############################################################

##DISTORTS GEOMETRIES IN DIRECTION OF IMAG FREQ################
def distort(freqlog):
    num_geoms = 1
    T = 300
    G, atomos = lx.parser.pega_geom(freqlog)
    F, M      = lx.parser.pega_freq(freqlog)
    NNC       = lx.parser.pega_modos(G,freqlog)
    num_atom  = np.shape(G)[0]
    A = np.zeros((3*num_atom,num_geoms))
    # check for imaginary frequencies
    if np.all(F > 0):
        lx.parser.fatal_error("No imaginary frequencies found. Exiting...")
    for i in range(len(F)):
        if F[i] < 0:
            f = -1*F[i]
            q = 3*np.sqrt(HBAR_J/(2*M[i]*f*np.tanh(HBAR_EV*f/(2*BOLTZ_EV*T))))
            q = np.array(q)
            A += np.outer(NNC[:,i],q)
    for n in range(np.shape(A)[1]):
        A1 = np.reshape(A[:,n],(num_atom,3))
        try:
            Gfinal = np.hstack((Gfinal,A1 + G))
        except:
            Gfinal = A1 + G
    base, _, nproc, mem, _, _ = lx.parser.busca_input(freqlog)
    cm = lx.parser.get_cm(freqlog)
    header = f"%nproc={nproc}\n%mem={mem}\n# {base} opt freq=noraman\n\nDISTORTED GEOM\n\n{cm}\n"
    bottom = '\n'
    for n in range(0,np.shape(A)[1],3):
        Gfinal = Gfinal[:,n:n+3]
        write_input(atomos,Gfinal,header,bottom,"distorted.com")
    print("New input file written to distorted.com")
###############################################################


def sample_geometries(freqlog,num_geoms,temp,limit=np.inf,warning=True):
    geom, atomos = lx.parser.pega_geom(freqlog)
    old = adjacency(geom,atomos)
    freqs, masses = lx.parser.pega_freq(freqlog)
    normal_coord = lx.parser.pega_modosLP(geom,freqlog)
    # check for negative frequencies
    if warning:
        lx.parser.double_check(freqlog)
    else:
        freqs[freqs < 0] *= -1
        mask = freqs < limit*(LIGHT_SPEED*100*2*np.pi)
        freqs = freqs[mask]
        normal_coord = normal_coord[:,mask]
    structures = np.zeros((geom.shape[0],geom.shape[1],num_geoms))
    scales = 1e10*np.sqrt(HBAR_J/(2*masses*freqs*np.tanh(HBAR_EV*freqs/(2*BOLTZ_EV*temp))))
    for j in range(num_geoms):
        ok = False
        while not ok:
            start_geom = geom.copy()
            nums = []
            for i, scale in enumerate(scales):
                normal = norm(scale=scale,loc=0)
                q = normal.rvs(size=1)
                start_geom += q*normal_coord[:,:,i]
                nums.append(q)
            new = adjacency(start_geom,atomos)
            if 0.5*np.sum(np.abs(old-new)) < 2:
                ok = True
                structures[:,:,j] = start_geom
            else:
                print('rejected',0.5*np.sum(np.abs(old-new)))
        nums = np.array(nums).T
        try:
            numbers = np.vstack((numbers,nums))
        except UnboundLocalError:
            numbers = nums
    numbers = np.round(numbers,4)
    return numbers, atomos, structures
###############################################################

##MAKES ENSEMBLE###############################################
def make_ensemble(freqlog, num_geoms, temp, header, bottom):
    try:
        os.mkdir('Geometries')
    except FileExistsError:
        pass
    counter = start_counter()
    print("\nGenerating geometries...\n")
    numbers, atomos, A = sample_geometries(freqlog,num_geoms,temp)
    F, M      = lx.parser.pega_freq(freqlog)
    #convert numbers to dataframe
    numbers = pd.DataFrame(numbers,columns=[f"mode_{i+1}" for i in range(np.shape(numbers)[1])])
    #check if file exists
    if os.path.isfile(f'Magnitudes_{temp:.0f}K_.lx'):
        data = pd.read_csv(f'Magnitudes_{temp:.0f}K_.lx')
        # get only columns with mode_ in the name
        data = data.filter(regex='mode_')
        #remove nan values
        data = data.dropna()
        # join data and numbers on axis 0
        numbers = pd.concat([data,numbers],axis=0,ignore_index=True)
    # concatenate frequencies and masses to numbers
    numbers = pd.concat([pd.DataFrame(F,columns=['freq']),pd.DataFrame(M,columns=['mass']),numbers],axis=1)
    numbers.to_csv(f'Magnitudes_{temp:.0f}K_.lx',index=False)
    for n in range(np.shape(A)[2]):
        Gfinal = A[:,:,n]
        write_input(atomos,Gfinal,header.replace("UUUUU",str(n+1)),bottom.replace("UUUUU",str(n+1)),f"Geometries/Geometry-{n+1+counter}-.com")
        progress = 100*(n+1)/num_geoms
        text = f"{progress:2.1f}%"
        print(' ', text, "of the geometries done.",end="\r", flush=True)
    print("\n\nDone! Ready to run.")
################################################################

##COLLECTS RESULTS##############################################
def gather_data(opc, tipo):
    files = [file for file in os.listdir('Geometries') if ".log" in file and "Geometr" in file ]
    files = [i for i in files if 'Normal termination' in open('Geometries/'+i, 'r').read()]
    files = sorted(files, key=lambda file: float(file.split("-")[1]))
    with open("Samples.lx", 'w') as f:
        for file in files:
            num = file.split("-")[1]
            broadening = opc
            numeros, energies, fs, scfs = [], [], [], []
            corrected, total_corrected = -1, -1
            with open('Geometries/'+file, 'r') as g:
                for line in g: 
                    if "Excited State" in line:
                        line = line.split()
                        numeros.append(line[2])
                        energies.append(line[4])
                        fs.append(line[8][2:])
                    elif "Corrected transition energy" in line:
                        line = line.split()
                        corrected = line[4]
                    elif "Total energy after correction" in line:
                        line = line.split()
                        total_corrected = 27.2114*float(line[5])
                    elif "SCF Done:" in line:
                        line = line.split()
                        scfs.append(27.2114*float(line[4]))
                if len(numeros) > 0:
                    f.write("Geometry "+num+":  Vertical transition (eV) Oscillator strength Broadening Factor (eV) \n")
                    if corrected != -1 and tipo == 'abs': #abspcm
                        f.write(f"Excited State {numeros[0]}\t{corrected}\t{fs[0]}\t{broadening}\n")
                    elif corrected != -1 and tipo == 'emi': #emipcm
                        energy = total_corrected - scfs[-1]
                        f.write(f"Excited State {numeros[0]}\t{energy:.3f}\t{fs[0]}\t{broadening}\n")
                    elif corrected == -1 and tipo == 'emi':
                        f.write(f"Excited State {numeros[0]}\t{energies[0]}\t{fs[0]}\t{broadening}\n")
                    else:
                        for i in range(len(energies)):
                            f.write(f"Excited State {numeros[i]}\t{energies[i]}\t{fs[i]}\t{broadening}\n")
                    f.write("\n")
###############################################################


##NORMALIZED GAUSSIAN##########################################
def gauss(x,v,s):
    y =  (1/(np.sqrt(2*np.pi)*s))*np.exp(-0.5*((x-v)/s)**2)
    return y
###############################################################


##COMPUTES AVG TRANSITION DIPOLE MOMENT########################
def calc_tdm(O,V):
    #Energy terms converted to J
    term = E_CHARGE*(HBAR_J**2)/V
    dipoles = np.sqrt(3*term*O/(2*MASS_E))
    #Conversion in au
    dipoles *= 1.179474389E29
    return np.mean(dipoles)
###############################################################

##PREVENTS OVERWRITING#########################################
def naming(arquivo):
    new_arquivo = arquivo
    if arquivo in os.listdir('.'):
        duplo = True
        vers = 2
        while duplo:
            new_arquivo = str(vers)+arquivo
            if new_arquivo in os.listdir('.'):
                vers += 1
            else:
                duplo = False
    return new_arquivo
###############################################################

##CALCULATES FLUORESCENCE LIFETIME IN S########################
def calc_emi_rate(xd,yd,dyd):
    #Integrates the emission spectrum
    IntEmi = np.trapz(yd,xd)
    taxa   = (1/HBAR_EV)*IntEmi
    error  = (1/HBAR_EV)*np.sqrt(np.trapz((dyd**2),xd))
    return taxa, error
###############################################################

##COMPUTES SPECTRA#############################################
def spectra(tipo, num_ex, nr):
    if tipo == "abs":
        constante = (np.pi*(E_CHARGE**2)*HBAR_EV)/(2*nr*MASS_E*LIGHT_SPEED*EPSILON_0)*10**(20)
    elif tipo == 'emi':
        constante = ((nr**2)*(E_CHARGE**2)/(2*np.pi*HBAR_EV*MASS_E*(LIGHT_SPEED**3)*EPSILON_0))
    V, O, S = [], [], []
    N = 0
    with open("Samples.lx", 'r') as f:
        for line in f:
            if "Geometry" in line:
                N += 1
            elif "Excited State" in line and int(line.split()[2][:-1]) in num_ex:
                line = line.split()
                if float(line[3]) <= 0:
                    print('Ignoring geom with negative vertical transition!')
                    N -= 1
                else: 
                    V.append(float(line[3]))
                    O.append(float(line[4]))
                    S.append(float(line[5]))
    coms = start_counter()
    if len(V) == 0 or len(O) == 0:
        lx.parser.fatal_error("You need to run steps 1 and 2 first! Goodbye!")
    elif len(V) != coms*max(num_ex):
        print("Number of log files is less than the number of inputs. Something is not right! Computing the spectrum anyway...")
    V = np.array(V)
    O = np.array(O)
    S = np.array(S)
    if tipo == 'abs':
        espectro = (constante*O)
    else:
        espectro = (constante*(V**2)*O)
        tdm = calc_tdm(O,V)
    left = max(min(V)-3*max(S),0.001)
    right = max(V)+ 3*max(S)
    x  = np.linspace(left,right, int((right-left)/0.01))
    if tipo == 'abs':
        arquivo = 'cross_section.lx'
        primeira = "{:8s} {:8s} {:8s}\n".format("#Energy(ev)", "cross_section(A^2)", "error")
    else:
        arquivo = 'differential_rate.lx'
        primeira = "{:4s} {:4s} {:4s} TDM={:.3f} au\n".format("#Energy(ev)", "diff_rate", "error",tdm)
    arquivo = naming(arquivo)
    y = espectro[:,np.newaxis]*gauss(x,V[:,np.newaxis],S[:,np.newaxis])
    mean_y =   np.sum(y,axis=0)/N 
    #Error estimate
    sigma  =   np.sqrt(np.sum((y-mean_y)**2,axis=0)/(N*(N-1))) 

    if tipo == 'emi':
        #Emission rate calculations
        mean_rate, error_rate = calc_emi_rate(x, mean_y,sigma) 
        segunda = '# Total Rate S1 -> S0: {:5.2e} +/- {:5.2e} s^-1\n'.format(mean_rate,error_rate)
    else:
        segunda = '# Absorption from State: S0\n'

    print(N, "geometries considered.")
    with open(arquivo, 'w') as f:
        f.write(primeira)
        f.write(segunda)
        for i in range(0,len(x)):
            text = f"{x[i]:.6f} {mean_y[i]:.6e} {sigma[i]:.6e}\n"
            f.write(text)
    print(f'Spectrum printed in the {arquivo} file')
###############################################################

##CHECKS PROGRESS##############################################
def andamento():
    try:
        coms = [file for file in os.listdir("Geometries") if 'Geometr' in file and '.com' in file]
        logs = [file for file in os.listdir("Geometries") if 'Geometr' in file and '.log' in file]
        factor = 1
        with open('Geometries/'+coms[0], 'r') as f:
            for line in f:
                if 'Link1' in line:
                    factor = 2
        count = 0
        error = 0
        for file in logs:
            with open('Geometries/'+file, 'r') as f:
                for line in f:
                    if "Normal termination" in line:
                        count += 1
                    elif "Error termination" in line:
                        error += 1
        print("\n\nThere are", int(count/factor), "successfully completed calculations out of", len(coms), "inputs")
        if error > 0:
            print(f"There are {error} failed jobs. If you used option 2, check the nohup.out file for details.")                
        print(np.round(100*(count+error)/(factor*len(coms)),1), "% of the calculations have been run.")
    except:
        print('No files found! Check the folder!')
###############################################################


##FETCHES  FILES###############################################
def fetch_file(frase,ends):
    files = []
    for file in [i for i in os.listdir('.')]:
        for end in ends:
            if end in file:
                files.append(file)
    if len(files) == 0:
        lx.parser.fatal_error(f"No {frase} file found. Goodbye!")
    freqlog = 'nada0022'
    for file in sorted(files):
        print("\n"+file)
        resp = input(f'Is this the {frase} file? y ou n?\n')
        if resp.lower() == 'y':
            freqlog = file
            break
    if freqlog == 'nada0022':
        lx.parser.fatal_error(f"No {frase} file found. Goodbye!")
    return freqlog
###############################################################

##RUNS TASK MANAGER############################################
def batch():
    script   = fetch_file('batch script?',['.sh'])    
    num      = input("Number of jobs in each batch?\n")
    limite   = input("Maximum number of jobs to be submitted simultaneously?\n")
    gaussian = input('g16 or g09?\n')
    try:
        int(limite)
        int(num)
    except:
        lx.parser.fatal_error("These must be integers. Goodbye!")
    folder = os.path.dirname(os.path.realpath(__file__))
    with open('limit.lx','w') as f:
        f.write(limite)
    subprocess.Popen(['nohup', 'python3', folder+'/batch_lx.py', script, num, gaussian, '&'])
###############################################################


##RUNS W TUNING################################################
def omega_tuning():
    geomlog = fetch_file('input or log',['.com','.log'])
    base, _, nproc, mem, _, _ = lx.parser.busca_input(geomlog)
    if 'IOP' in base.upper() and ('108' in base or '107' in base):
        base2 = base.split()
        for elem in base2:
            if '/' in elem:
                base = elem
                break
    omega1 = '0.1'
    passo  = '0.025'
    relax  = 'y'
    print('This is the configuration taken from the file:\n')
    print(f'Functional/basis: {base}')
    print('%nproc='+nproc)
    print('%mem='+mem)
    print(f'Initial Omega: {omega1} bohr^-1')
    print(f'Step: {passo} bohr^-1')
    print('Optimize at each step: yes')
    change = input('Are you satisfied with these parameters? y or n?\n')
    if change == 'n':
        base   = default(base,f"Functional/basis is {base}. If ok, Enter. Otherwise, type functional/basis.\n")
        nproc  = default(nproc,f'nproc={nproc}. If ok, Enter. Otherwise, type it.\n')
        mem    = default(mem,f"mem={mem}. If ok, Enter. Otherwise, type it.\n")
        omega1 = default(omega1,f"Initial omega is {omega1} bohr^-1. If ok, Enter. Otherwise, type it.\n")
        passo  = default(passo,f"Initial step is {passo} bohr^-1. If ok, Enter. Otherwise, type it.\n")
        relax  = default(relax,"Optimize at each step: yes. If ok, Enter. Otherwise, type n\n")
    
    script = fetch_file('batch script',['.sh'])
    gaussian = input('g16 or g09?\n')
    folder = os.path.dirname(os.path.realpath(__file__))
    with open('limit.lx','w') as f:
        f.write('Running')
    subprocess.Popen(['nohup', 'python3', folder+'/omega.py', geomlog, base, nproc, mem, omega1, passo, relax, script, gaussian, '&'])
###############################################################

##RUNS CONFORMATIONAL SEARCH###################################
def conformational():
    freqlog = fetch_file('frequency',['.log'])
    F, _ = lx.parser.pega_freq(freqlog)
    F_active = F[:40]
    T  = int(HBAR_EV*F_active[-1]/BOLTZ_EV)
    DT = int(T/10)
    T, DT = str(T), str(DT)
    base, _, nproc, mem, _, _ = lx.parser.busca_input(freqlog)
    print('This is the configuration taken from the file:\n')
    print(f'Functional/basis: {base}')
    print(f'%nproc={nproc}')
    print(f'%mem={mem}')
    print(f'Initial Temperature: {T} K')
    print(f'Temperature step: {DT} K')
    change = input('Are you satisfied with these parameters? y or n?\n')
    if change == 'n':
        base   = default(base,f"Functional/basis is {base}. If ok, Enter. Otherwise, type functional/basis.\n")
        nproc  = default(nproc,f'nproc={nproc}. If ok, Enter. Otherwise, type it.\n')
        mem    = default(mem,f"mem={mem}. If ok, Enter. Otherwise, type it.\n")
        T      = default(T,f"Initial temperature is {T} K. If ok, Enter. Otherwise, type it.\n")
        DT     = default(DT,f"Temperature step is {DT} K. If ok, Enter. Otherwise, type it.\n")
    script    = fetch_file('batch script',['.sh'])
    num_geoms = input("Number of geometries sampled at each round?\n")
    rounds    = input("Number of rounds?\n")
    numjobs   = input("Number of jobs in each batch?\n")
    gaussian  = input('g16 or g09?\n')
    try:
        int(num_geoms)
        int(rounds)
        int(numjobs)
    except:
        lx.parser.fatal_error("These must be integers. Goodbye!")
    with open('limit.lx','w') as f:
        f.write('Running')
    folder = os.path.dirname(os.path.realpath(__file__))
    subprocess.Popen(['nohup', 'python3', folder+'/conf_search.py', freqlog, base, nproc, mem, T, DT, num_geoms, rounds,numjobs, script, gaussian, '&'])
###############################################################


##FINDS SUITABLE VALUE FOR STD#################################
def detect_sigma():
    try:
        files = [i for i in os.listdir('.') if 'Magnitudes' in i and '.lx' in i]
        file  = files[0]
        temp = float(file.split('_')[1].strip('K'))
        sigma =  np.round(BOLTZ_EV*temp,3)
    except:
        sigma = 0.000
    return sigma
###############################################################

##CHECKS SPECTRUM TYPE#########################################
def get_spec():
    coms = [file for file in os.listdir("Geometries") if 'Geometr' in file and '.com' in file]
    with open('Geometries/'+coms[0],'r') as f:
        for line in f:
            if 'ABSSPCT' in line:
                tipo = 'absorption'
                break
            elif 'EMISPCT' in line:
                tipo = 'emission'
                break
    return tipo
###############################################################

##QUERY FUNCTION###############################################
def default(a,frase):
    b = input(frase)
    if b == '':
        return a
    return b
###############################################################

##SETS DIELECTRIC CONSTANTS####################################
def set_eps(scrf):
    if 'READ' in scrf.upper():
        eps1 = input("Type the static dielectric constant.\n")
        eps2 = input("Type the dynamic dielectric constant (n^2).\n")
        try:
            float(eps1)
            float(eps2)
        except:
            lx.parser.fatal_error("The constants must be numbers. Goodbye!")
        epss = "Eps="+eps1+"\nEpsInf="+eps2+"\n\n"
    else:
        epss = '\n'
    return epss
###############################################################

##STOP SUBMISSION OF JOBS######################################
def abort_batch():
    choice = input('Are you sure you want to prevent new jobs from being submitted? y or n?\n')
    if choice == 'y':
        try:
            os.remove('limit.lx')
            print('Done!')
        except:
            print('Could not find the files. Maybe you are in the wrong folder.')
    else:
        print('OK, nevermind')
###############################################################

##DELETES CHK FILES############################################
def delchk(input_file,term):
    num = input_file.split('-')[1]
    if term == 1:
        a = ''
    elif term == 2:
        a = '2'
    try:
        os.remove(f'step{a}_{num}.chk')
    except:
        pass
###############################################################

##CHECKS WHETHER JOBS ARE DONE#################################
def watcher(files,counter):
    rodando = files.copy()
    done = []
    for input_file in rodando:
        term = 0
        error = False
        try:
            with open(input_file[:-3]+'log', 'r') as f:
                for line in f:
                    if 'Normal termination' in line:
                        term += 1
                        if counter == 2:
                            delchk(input,term)
                    elif 'Error termination' in line:
                        error = True
                        print(f'The following job returned an error: {input_file}')
                        print('Please check the file for any syntax errors.')
            if term == counter or error:
                done.append(input)
        except:
            pass
    for elem in done:
        del rodando[rodando.index(elem)]
    return rodando
###############################################################


##GETS SPECTRA#################################################
def search_spectra():
    Abs, Emi = 'None', 'None'
    candidates = [i for i in os.listdir('.') if '.lx' in i]
    for candidate in candidates:
        with open(candidate, 'r') as f:
            for line in f:
                if 'cross_section' in line:
                    Abs = candidate
                elif 'diff_rate' in line:
                    Emi = candidate
                break
    return Abs, Emi
###############################################################

##RUNS EXCITON ANALYSIS########################################
def ld():
    Abs, Emi = search_spectra()
    print(f'Absorption file: {Abs}')
    print(f'Emission file: {Emi}')
    check = input('Are these correct? y or n?\n')
    if check == 'n':
        Abs = input('Type name of the absorption spectrum file\n')
        Emi = input('Type name of the emission spectrum file\n')

    kappa = input('Orientation Factor (k^2):\n')
    rmin  = input("Average intermolecular distance in Å:\n")
    Phi   = input("Fluorescence quantum yield (from 0 to 1):\n")
    try:
        rmin  = float(rmin)
        kappa = np.sqrt(float(kappa))
        Phi   = float(Phi)
    except:
        lx.parser.fatal_error('These features must be numbers. Goodbye!')    
    if Phi > 1 or Phi < 0:
        lx.parser.fatal_error('Quantum yield must be between 0 and 1. Goodbye!')

    correct = input('Include correction for short distances? y or n?\n')
    if correct == 'y':
        alpha = 1.15*0.53
        print('Employing correction!')
    else:
        alpha = 0
        print('Not employing correction!')

    print('Computing...')
    try:
        run_ld(Abs, Emi, alpha, rmin, kappa, Phi)
        print('Results can be found in the ld.lx file')
    except:
        print('Something went wrong. Check if the name of the files are correct.')
###############################################################

##CHECKS WHETHER JOBS ARE DONE#################################
def hold_watch(files, log):
    rodando = files.copy()
    while len(rodando) > 0:
        rodando = watcher(rodando,1)
        if 'limit.lx' not in os.listdir('.'):
            with open(log,'a') as f:
                f.write('\n#Aborted!')
            sys.exit()
        time.sleep(30)
###############################################################

##RUNS CALCULATIONS############################################
def rodar_lista(lista, batch_file, gaussian, log, num=1): 
    #number of scripts is integer division of number of files by num
    n = len(lista)//num
    for i in range(n):
        with open(f'cmd_{i}.sh', 'w') as f:
            for file in lista[i*num:(i+1)*num]:
                f.write(f'{gaussian} {file}\n')
        subprocess.call(['bash', batch_file, f'cmd_{i}.sh'])
    if len(lista)%num != 0:
        with open(f'cmd_{n+1}.sh', 'w') as f:
            for file in lista[n*num:]:
                f.write(f'{gaussian} {file}\n')
        subprocess.call(['bash', batch_file, f'cmd_{n+1}.sh'])
    hold_watch(lista, log)
###############################################################

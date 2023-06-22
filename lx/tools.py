#!/usr/bin/env python3
import numpy as np
import os
import sys
from scipy.stats import norm
import time
import subprocess
import pandas as pd

##SOME CONSTANTS##############################################
epsilon0 = 8.854187817e-12   #F/m
hbar = 6.582119514e-16       #eV s
hbar2 = 1.054571800e-34      #J s
mass = 9.10938356e-31        #kg
c = 299792458                #m/s
e = 1.60217662e-19           #C
kb = 8.6173303e-5            #eV/K
amu = 1.660539040e-27        #kg
pi = np.pi
###############################################################

##ERROR FUNCTION###############################################
def fatal_error(msg):
    print(msg)
    sys.exit()
###############################################################

##GETS FREQUENCIES AND REDUCED MASSES##########################
def pega_freq(freqlog):
    F, M = [], []
    with open(freqlog, 'r') as f:
        for line in f:
            if "Frequencies --" in line:
                line = line.split()
                for j in range(2,len(line)):
                    if float(line[j]) in F:
                        pass
                    F.append(float(line[j]))
            elif "Red. masses --" in line:
                line = line.split()
                for j in range(3,len(line)):
                    M.append(float(line[j]))
            elif 'Thermochemistry' in line:
                break        
    #conversion in angular frequency
    F = np.array(F)*(c*100*2*pi) 
    try:
        f = F[0]
    except:
        fatal_error("No frequencies in the log file! Goodbye!")
    #conversion from amu to kg
    M = np.asarray(M)*amu
    return F, M
###############################################################

##GETS ATOMS AND LAST GEOMETRY IN FILE#########################
def pega_geom(freqlog):
    if ".log" in freqlog:
        status = 0
        busca = "orientation:"
        n = -1
        with open(freqlog, 'r') as f:
            for line in f:
                if busca in line and 'Dipole' not in line:
                    n = 0
                    G = np.zeros((1,3))
                    atomos = []
                elif n >= 0 and n < 4:
                    n += 1
                elif n >= 4 and "---------------------------------------------------------------------" not in line:    
                    line = line.split()
                    NG = []
                    for j in range(3,len(line)):
                        NG.append(float(line[j]))
                    atomos.append(line[1])
                    G = np.vstack((G,NG))       
                    n += 1  
                elif "---------------------------------------------------------------------" in line and n>1:
                    n = -1       
    else:
        G = np.zeros((1,3))
        atomos = []
        with open(freqlog, 'r') as f:
            for line in f:
                line = line.split()
                try:
                    vetor = np.array([float(line[1]),float(line[2]), float(line[3])])
                    atomos.append(line[0])
                    G = np.vstack((G,vetor))
                except:
                    pass
    G = G[1:,:]                 
    return G, atomos
###############################################################

##SAVES OPT GEOMETRY###########################################
def salva_geom(G,atomos):
    atomos = np.array([atomos]).astype(float)
    atomos = atomos.T
    G = np.hstack((atomos,G))
    np.savetxt('opt_geom.lx', G, delimiter='\t', fmt=['%1.1u','%+1.5f','%+1.5f','%+1.5f'])
    print("The optimized geometry that is used is saved in the opt_geom.lx file!")
###############################################################

##GETS NORMAL COORDINATES IN HIGH PRECISION####################
def pega_modosHP(G, freqlog):
    F, M = pega_freq(freqlog)
    n = -1
    num_atom = np.shape(G)[0]
    NC = np.zeros((3*num_atom,1))
    with open(freqlog, 'r') as f:
        for line in f:
            if n == 0:
                line = line.split()[3:]
                C = np.asarray([float(i) for i in line])
                n += 1
            elif n < 0 or n > 3*num_atom:
                if "Coord Atom Element:" in line:
                    n = 0
            elif n > 0 and n < 3*num_atom:
                line = line.split()[3:]
                line = np.asarray([float(i) for i in line])
                C = np.vstack((C,line))
                n += 1  
            elif n == 3*num_atom:
                NC = np.hstack((NC,C))
                n += 1
    NC = NC[:,1:]
    MM = np.zeros((1,len(F)))
    M = np.expand_dims(M,axis=0)
    for _ in range(0,3*num_atom):
        MM = np.vstack((MM,M))
    M = MM[1:,:]
    return NC
###############################################################

##GETS NORMAL COORDINATES IN REGULAR PRECISION#################
def pega_modosLP(G,freqlog):
    F, M = pega_freq(freqlog)
    C = []
    n = -1
    num_atom = np.shape(G)[0]
    with open(freqlog, 'r') as f:
        for line in f:
            if n < 0 or n >= num_atom:
                if "Atom  AN      X      Y      Z" in line:
                    n = 0
                else:
                    pass
            elif n >= 0 and n < num_atom:
                line = line.split()
                for j in range(2,len(line)):
                    C.append(float(line[j]))
                n += 1  
                
    num_modos = len(F)
    
    l = 0
    p = 0
    NNC = np.zeros((num_atom,1))
    while l < num_modos:
        NC = np.zeros((1,3))
        k =0
        while k < num_atom:     
            B = np.asarray(C[3*(l+3*k)+p:3*(l+3*k)+3+p])
            NC = np.vstack((NC,B))
            k += 1      
        NNC = np.hstack((NNC,NC[1:,:]))
        l += 1
        if l%3 == 0 and l != 0:
            p = p + (num_atom-1)*9  
    NNC = NNC[:,1:] #matriz com as coordenadas normais de cada modo
    D = np.zeros((3*num_atom,1))
    for i in range(0,len(F)):
        normal = NNC[:,3*i:3*i+3].flatten()
        normal = np.expand_dims(normal,axis=1)
        D = np.hstack((D,normal))
    D = D[:,1:]
    MM = np.zeros((1,len(F)))
    M = np.expand_dims(M,axis=0)
    for i in range(0,3*num_atom):
        MM = np.vstack((MM,M))
    M = MM[1:,:]
    return D
###############################################################

##DETECTS WHETHER HIGH PRECISION IS USED#######################
def pega_modos(G,freqlog):
    x = 'LP'
    with open(freqlog, 'r') as f:
        for line in f:
            if "Coord Atom Element:" in line:
                x = 'HP'
                break
    if x == 'LP':
        return pega_modosLP(G,freqlog)
    else:
        return pega_modosHP(G,freqlog)
###############################################################

##WRITES ATOMS AND XYZ COORDS TO FILE##########################
def write_input(atomos,G,header,bottom,file):
    with open(file, 'w') as f:
        f.write(header)
        for i in range(0,len(atomos)):
            texto = "{:2s}  {:.14f}  {:.14f}  {:.14f}\n".format(atomos[i],G[i,0],G[i,1],G[i,2])
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
    G, atomos = pega_geom(freqlog)
    F, M      = pega_freq(freqlog)
    NNC       = pega_modos(G,freqlog)
    num_atom  = np.shape(G)[0]
    A = np.zeros((3*num_atom,num_geoms))
    # check for imaginary frequencies
    if np.all(F > 0):
        fatal_error("No imaginary frequencies found. Exiting...")
    for i in range(0,len(F)):
        if F[i] < 0:
            f = -1*F[i]
            q = 3*np.sqrt(hbar2/(2*M[i]*f*np.tanh(hbar*f/(2*kb*T))))
            q = np.array(q)
            A += np.outer(NNC[:,i],q)
    for n in range(np.shape(A)[1]):
        A1 = np.reshape(A[:,n],(num_atom,3))
        try:
            Gfinal = np.hstack((Gfinal,A1 + G))
        except:
            Gfinal = A1 + G
    base, _, nproc, mem, scrf, spec = busca_input(freqlog)
    cm = get_cm(freqlog)
    header = "%nproc={}\n%mem={}\n# {} opt freq=noraman\n\nDISTORTED GEOM\n\n{}\n".format(nproc,mem,base,cm)
    bottom = '\n'
    for n in range(0,np.shape(A)[1],3):
        Gfinal = Gfinal[:,n:n+3]
        write_input(atomos,Gfinal,header,bottom,"distorted.com")
    print("New input file written to distorted.com")    
###############################################################

##CHECKS FREQ FILES############################################
def double_check(freqlog):
    nproc, mem, header = get_input_params(freqlog)
    header = header.lower()
    if "opt" in header and 'freq' in header and 'iop(' in header:
        optfreqissue = True
    with open(freqlog, 'r') as f:
        for line in f:
            if "Non-Optimized Parameters" in line:
                print('*'*50)
                print("WARNING: Non-optimized parameters detected in your frequency file.")
                print('Even though the frequencies may be all real, your structure is still not fully optimized.')
                print('This may lead to inaccurate results.')
                if optfreqissue:
                    print('In your case, this may be due to running an opt freq calculation using a single input file and IOP options.')
                    print('Gaussian does not carry the IOP options to the frequency calculation when using a single input file.')
                    print('To avoid this issue, run the optimization and frequency calculations separately.')
                else:
                    print('To learn more about this issue, check https://gaussian.com/faq3/ .')        
                print('Proceed at your own risk.')
                print('*'*50)
                print('\n')
###############################################################              

##SAMPLES GEOMETRIES###########################################
def sample_geometries(freqlog,num_geoms,T, limit=np.inf, warning=True):
    G, atomos = pega_geom(freqlog)
    F, M      = pega_freq(freqlog)
    NNC       = pega_modos(G,freqlog)
    # check for negative frequencies
    if warning:
        if np.any(F < 0):
            fatal_error("Imaginary frequencies detected. Check your frequency file. Goodbye.")
        double_check(freqlog)
    else:
        F[F < 0] *= -1
        mask = F < limit*(c*100*2*pi)
        F = F[mask]
        NNC = NNC[:,mask]
    num_atom  = np.shape(G)[0]
    A = np.zeros((3*num_atom,num_geoms))
    for i in range(0,len(F)):
        scale = np.sqrt(hbar2/(2*M[i]*F[i]*np.tanh(hbar*F[i]/(2*kb*T))))
        normal = norm(scale=scale,loc=0)
        #Displacements in  Å
        q = normal.rvs(size=num_geoms)*1e10
        try:
            numbers = np.hstack((numbers,q[:,np.newaxis]))
        except:
            numbers = q[:,np.newaxis]
        A += np.outer(NNC[:,i],q)
    for n in range(np.shape(A)[1]):
        A1 = np.reshape(A[:,n],(num_atom,3))
        try:
            Gfinal = np.hstack((Gfinal,A1 + G))
        except:
            Gfinal = A1 + G     
    numbers = np.round(numbers,4)
    return numbers, atomos, Gfinal
###############################################################

##MAKES ENSEMBLE###############################################
def make_ensemble(freqlog, num_geoms, T, header, bottom):
    try:
        os.mkdir('Geometries')
    except:
        pass        
    counter = start_counter()   
    print("\nGenerating geometries...\n")
    numbers, atomos, A = sample_geometries(freqlog,num_geoms,T)
    F, M      = pega_freq(freqlog)
    #convert numbers to dataframe
    numbers = pd.DataFrame(numbers,columns=[f"mode_{i+1}" for i in range(np.shape(numbers)[1])])
    #check if file exists
    if os.path.isfile(f'Magnitudes_{T:.0f}K_.lx'):
        data = pd.read_csv(f'Magnitudes_{T:.0f}K_.lx')
        # get only columns with mode_ in the name
        data = data.filter(regex='mode_')
        #remove nan values
        data = data.dropna()
        # join data and numbers on axis 0
        numbers = pd.concat([data,numbers],axis=0,ignore_index=True)
    # concatenate frequencies and masses to numbers
    numbers = pd.concat([pd.DataFrame(F,columns=['freq']),pd.DataFrame(M,columns=['mass']),numbers],axis=1)
    numbers.to_csv(f'Magnitudes_{T:.0f}K_.lx',index=False)
    for n in range(0,np.shape(A)[1],3):
        Gfinal = A[:,n:n+3]  
        write_input(atomos,Gfinal,header.replace("UUUUU",str((n+3)//3)),bottom.replace("UUUUU",str((n+3)//3)),"Geometries/Geometry-"+str((n+3)//3+counter)+"-.com")
        progress = 100*((n+3)//3)/num_geoms
        text = "{:2.1f}%".format(progress)
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
                        f.write("Excited State {}\t{}\t{}\t{}\n".format(numeros[0],corrected, fs[0], broadening))
                    elif corrected != -1 and tipo == 'emi': #emipcm     
                        energy = total_corrected - scfs[-1]
                        f.write("Excited State {}\t{:.3f}\t{}\t{}\n".format(numeros[0],energy,fs[0],broadening))
                    elif corrected == -1 and tipo == 'emi':
                        f.write("Excited State {}\t{}\t{}\t{}\n".format(numeros[0],energies[0],fs[0],broadening))
                    else:
                        for i in range(len(energies)):
                            f.write("Excited State {}\t{}\t{}\t{}\n".format(numeros[i],energies[i],fs[i],broadening))
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
    term = e*(hbar2**2)/V
    dipoles = np.sqrt(3*term*O/(2*mass))
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
    taxa   = (1/hbar)*IntEmi
    error  = (1/hbar)*np.sqrt(np.trapz((dyd**2),xd))
    return taxa, error 
###############################################################

##COMPUTES SPECTRA############################################# 
def spectra(tipo, num_ex, nr):
    if tipo == "abs":
        constante = (np.pi*(e**2)*hbar)/(2*nr*mass*c*epsilon0)*10**(20)
    elif tipo == 'emi':
        constante = ((nr**2)*(e**2)/(2*np.pi*hbar*mass*(c**3)*epsilon0))
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
        fatal_error("You need to run steps 1 and 2 first! Goodbye!")
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
            text = "{:.6f} {:.6e} {:.6e}\n".format(x[i],mean_y[i], sigma[i])
            f.write(text)
    print('Spectrum printed in the {} file'.format(arquivo))                
############################################################### 

##GETS INPUT PARAMS FROM LOG FILES#############################
def get_input_params(freqlog):
    nproc, mem, header = '', '', ''
    with open(freqlog, 'r') as f:
        search = False
        for line in f:
            if '%nproc' in line.lower():
                line = line.split('=')
                nproc = line[-1].replace('\n','')
            elif '%mem' in line.lower():
                line = line.split('=')
                mem = line[-1].replace('\n','')
            elif "#" in line and not search and header == '':
                search = True
                header += line.lstrip().replace('\n','')
            elif search and '----------' not in line:
                header += line.lstrip().replace('\n','')
            elif search and '----------' in line:
                search = False
                break
    return nproc, mem, header        
###############################################################

##CHECKS THE FREQUENCY LOG'S LEVEL OF THEORY###################
def busca_input(freqlog):
    base = 'lalala'
    exc = ''
    header = ''
    nproc = '4'
    mem   = '1GB'
    scrf  = ''
    nproc, mem, header = get_input_params(freqlog)
    
    if 'TDA' in header.upper():
        exc = 'tda'
        spec = 'EMISPCT'
    elif 'TD' in header.upper():
        exc = 'td'
        spec = 'EMISPCT'
    else:
        spec = 'ABSSPCT'

    if 'SCRF' in header.upper():
        new = header.split()
        for elem in new:
            if 'SCRF' in elem:
                scrf = elem
                break     
        
    header = header.split()
    base = ''
    for elem in header:
        if "/" in elem and 'IOP' not in elem.upper():
            base += elem.replace('#','')
        elif 'IOP' in elem.upper() and ('108' in elem or '107' in elem):
            base += ' '+elem
    return base, exc, nproc, mem, scrf, spec                
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
            print("There are {} failed jobs. If you used option 2, check the nohup.out file for details.".format(error))                
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
        fatal_error("No {} file found. Goodbye!".format(frase))
    freqlog = 'nada0022'    
    for file in sorted(files):
        print("\n"+file)
        resp = input('Is this the {} file? y ou n?\n'.format(frase))
        if resp.lower() == 'y':
            freqlog = file
            break
    if freqlog == 'nada0022':
        fatal_error("No {} file found. Goodbye!".format(frase))
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
        fatal_error("These must be integers. Goodbye!")
    
    import subprocess
    folder = os.path.dirname(os.path.realpath(__file__)) 
    with open('limit.lx','w') as f:
        f.write(limite)
    subprocess.Popen(['nohup', 'python3', folder+'/batch_lx.py', script, num, gaussian, '&'])
###############################################################


##RUNS W TUNING################################################
def omega_tuning():
    geomlog = fetch_file('input or log',['.com','.log'])
    base, _, nproc, mem, _, _ = busca_input(geomlog)
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
    print('Functional/basis: {}'.format(base))
    print('%nproc='+nproc)    
    print('%mem='+mem)
    print('Initial Omega: {} bohr^-1'.format(omega1))
    print('Step: {} bohr^-1'.format(passo))
    print('Optimize at each step: yes')
    change = input('Are you satisfied with these parameters? y or n?\n')
    if change == 'n':
        base   = default(base,"Functional/basis is {}. If ok, Enter. Otherwise, type functional/basis.\n".format(base))
        nproc  = default(nproc,'nproc={}. If ok, Enter. Otherwise, type it.\n'.format(nproc))
        mem    = default(mem,"mem={}. If ok, Enter. Otherwise, type it.\n".format(mem))
        omega1 = default(omega1,"Initial omega is {} bohr^-1. If ok, Enter. Otherwise, type it.\n".format(omega1))       
        passo  = default(passo,"Initial step is {} bohr^-1. If ok, Enter. Otherwise, type it.\n".format(passo))       
        relax  = default(relax,"Optimize at each step: yes. If ok, Enter. Otherwise, type n\n")
    
    script = fetch_file('batch script',['.sh'])    
    gaussian = input('g16 or g09?\n')
    import subprocess
    folder = os.path.dirname(os.path.realpath(__file__)) 
    with open('limit.lx','w') as f:
        f.write('Running')
    subprocess.Popen(['nohup', 'python3', folder+'/omega.py', geomlog, base, nproc, mem, omega1, passo, relax, script, gaussian, '&'])
###############################################################

##RUNS CONFORMATIONAL SEARCH###################################
def conformational():
    freqlog = fetch_file('frequency',['.log'])
    F, _ = pega_freq(freqlog)
    F_active = F[:40]
    T  = int(hbar*F_active[-1]/kb)
    DT = int(T/10)
    T, DT = str(T), str(DT)
    base, _, nproc, mem, _, _ = busca_input(freqlog)
    print('This is the configuration taken from the file:\n')
    print('Functional/basis: {}'.format(base))
    print('%nproc='+nproc)    
    print('%mem='+mem)
    print('Initial Temperature: {} K'.format(T))
    print('Temperature step: {} K'.format(DT))
    change = input('Are you satisfied with these parameters? y or n?\n')
    if change == 'n':
        base   = default(base,"Functional/basis is {}. If ok, Enter. Otherwise, type functional/basis.\n".format(base))
        nproc  = default(nproc,'nproc={}. If ok, Enter. Otherwise, type it.\n'.format(nproc))
        mem    = default(mem,"mem={}. If ok, Enter. Otherwise, type it.\n".format(mem))
        T      = default(T,"Initial temperature is {} K. If ok, Enter. Otherwise, type it.\n".format(T))
        DT     = default(DT,"Temperature step is {} K. If ok, Enter. Otherwise, type it.\n".format(DT))
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
        fatal_error("These must be integers. Goodbye!")
    with open('limit.lx','w') as f:
        f.write('Running')
    import subprocess
    folder = os.path.dirname(os.path.realpath(__file__)) 
    subprocess.Popen(['nohup', 'python3', folder+'/conf_search.py', freqlog, base, nproc, mem, T, DT, num_geoms, rounds,numjobs, script, gaussian, '&'])
###############################################################


##FINDS SUITABLE VALUE FOR STD#################################    
def detect_sigma():
    try:
        files = [i for i in os.listdir('.') if 'Magnitudes' in i and '.lx' in i]
        file  = files[0]
        temp = float(file.split('_')[1].strip('K'))
        sigma =  np.round(kb*temp,3)
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

##FETCHES REFRACTIVE INDEX##################################### 
def get_nr():
    buscar = False
    coms = [file for file in os.listdir("Geometries") if 'Geometr' in file and '.com' in file]
    with open('Geometries/'+coms[0],'r') as f:
        for line in f:
            if 'SCRF' in line.upper():
                buscar = True
                break
    if buscar:
        logs = [file for file in os.listdir("Geometries") if 'Geometr' in file and '.log' in file]
        for log in logs:
            with open('Geometries/'+log,'r') as f:
                for line in f:
                    if 'Solvent' in line and 'Eps' in line:
                        line = line.split()
                        nr = np.sqrt(float(line[6]))
                        return nr
    else:
        return 1                
###############################################################

##FETCHES CHARGE AND MULTIPLICITY##############################
def get_cm(freqlog):
    with open(freqlog,'r') as f:
        for line in f:
            if 'Charge' in line and 'Multiplicity' in line:
                line = line.split()
                charge = line[2]
                mult   = line[5]
                break
    return charge+' '+mult
###############################################################

##QUERY FUNCTION###############################################
def default(a,frase):
    b = input(frase)
    if b == '':
        return a
    else:
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
            fatal_error("The constants must be numbers. Goodbye!")
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
def delchk(input,term):
    num = input.split('-')[1]
    if term == 1:
        a = ''
    elif term == 2:
        a = '2'
    try:        
        os.remove('step{}_{}.chk'.format(a,num))
    except:
        pass      
###############################################################

##CHECKS WHETHER JOBS ARE DONE#################################
def watcher(files,counter):
    rodando = files.copy()
    done = []
    for input in rodando: 
        term = 0
        error = False
        try:
            with open(input[:-3]+'log', 'r') as f:
                for line in f:
                    if 'Normal termination' in line:
                        term += 1
                        if counter == 2:
                            delchk(input,term)
                    elif 'Error termination' in line:
                        error = True
                        print('The following job returned an error: {}'.format(input))
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
    print('Absorption file: {}'.format(Abs))
    print('Emission file: {}'.format(Emi))
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
        fatal_error('These features must be numbers. Goodbye!')    
    if Phi > 1 or Phi < 0:
        fatal_error('Quantum yield must be between 0 and 1. Goodbye!')

    correct = input('Include correction for short distances? y or n?\n')
    if correct == 'y':
        alpha = 1.15*0.53 
        print('Employing correction!')
    else:
        alpha = 0
        print('Not employing correction!')
    
    print('Computing...')
    from lx.ld import run_ld 
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
                f.write('{} {}\n'.format(gaussian,file))
        subprocess.call(['bash', batch_file, f'cmd_{i}.sh']) 
    if len(lista)%num != 0:
        with open(f'cmd_{n+1}.sh', 'w') as f:
            for file in lista[n*num:]:
                f.write('{} {}\n'.format(gaussian,file))
        subprocess.call(['bash', batch_file, f'cmd_{n+1}.sh'])
    hold_watch(lista, log)
###############################################################

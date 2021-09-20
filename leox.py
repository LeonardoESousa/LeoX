#!/usr/bin/env python3
import numpy as np
import os
import random
import sys


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
        G = np.zeros((1,4))
        with open(freqlog, 'r') as f:
            for line in f:
                line = line.split()
                try:
                    vetor = np.array([float(line[0]),float(line[1]),float(line[2]), float(line[3])])
                    G = np.vstack((G,vetor))
                except:
                    pass
    try:
        G = G[1:,:]                 
    except:
        fatal_error("No geometry in the log file! Goodbye!")
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
        return pega_modosLP(G,freqlog)
###############################################################

##DISPLACE GEOMETRY IN DIRECTIONS WITH IMAGINARY FREQ##########
def shake(freqlog, T):
    F, M = pega_freq(freqlog)
    G, atomos = pega_geom(freqlog)
    NNC = pega_modos(G,freqlog)
    num_atom = np.shape(G)[0]
    A1 = np.zeros((3*num_atom,1))
    A2 = np.zeros((3*num_atom,1))
    F = F[F < 0]
    if len(F) == 0:
        fatal_error("No imaginary frquencies in the log file. Goodbye!")
    F = -1*F
    for i in range(len(F)): # LL:
        q = [-1*T,T]
        A1 += q[0]*(np.expand_dims(NNC[:,i],axis=1))
        A2 += q[1]*(np.expand_dims(NNC[:,i],axis=1))
    A1 = np.reshape(A1,(num_atom,3))
    A2 = np.reshape(A2,(num_atom,3))
    Gfinal  = A1 + G
    Gfinal2 = A2 + G
    with open("shaken.lx", 'w') as f:
        f.write('#Geometry with displacement of '+str(T)+" A:\n" )
        for k in range(0, np.shape(Gfinal)[0]):
            text = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[k],Gfinal[k,0],Gfinal[k,1],Gfinal[k,2])
            f.write(text+"\n")
        f.write('\n#Geometry with displacement of '+str(-T)+" A:\n" )
        for k in range(0, np.shape(Gfinal2)[0]):
            text = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[k],Gfinal2[k,0],Gfinal2[k,1],Gfinal2[k,2])
            f.write(text+"\n")
    print("There are 2 geometries saved on shaken.lx!")
###############################################################

##CHECKS FOR EXISTING GEOMETRIES###############################
def start_counter():
    files = [file for file in os.listdir('Geometries') if ".com" in file and "Geometry-" in file]
    return len(files)
###############################################################

##SAMPLES GEOMETRIES###########################################
def sample_geom(freqlog, num_geoms, T, header, bottom):
    F, M = pega_freq(freqlog)
    if F[0] < 0:
        fatal_error("Imaginary frequency! Goodbye!")
    try:
        os.mkdir('Geometries')
    except:
        pass        
    counter = start_counter()
    G, atomos = pega_geom(freqlog)
    salva_geom(G,atomos)
    NNC = pega_modos(G,freqlog)
    num_atom = np.shape(G)[0]   
    print("\nGenerating geometries...\n")
    with open('Magnitudes.lx', 'w') as file:
        for n in range(1,num_geoms+1):
            A = np.zeros((3*num_atom,1))
            numbers = []
            for i in range(0,len(F)):
                x = np.linspace(-5, 5, 10000) #ja em angstrom
                boltz = np.tanh(hbar*F[i]/(2*kb*T))
                prob = np.sqrt((M[i]*F[i]*(boltz))/(np.pi*hbar2))*np.exp(-M[i]*F[i]*((x*(10**(-10)))**2)*(boltz)/hbar2)*(abs(x[1]-x[0])*10**(-10)) #com temperatura
                q = random.choices(x, prob)
                numbers.append(q[0])
                A += q[0]*(np.expand_dims(NNC[:,i],axis=1))
            numbers = np.round(np.array(numbers)[np.newaxis,:],4)
            np.savetxt(file, numbers, delimiter='\t', fmt='%s')
            A = np.reshape(A,(num_atom,3))
            Gfinal = A + G  
            with open("Geometries/Geometry-"+str(n+counter)+"-.com", 'w') as f:
                    f.write(header.replace("UUUUU",str(n)))
                    for k in range(0, np.shape(Gfinal)[0]):
                        text = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[k],Gfinal[k,0],Gfinal[k,1],Gfinal[k,2])
                        f.write(text+"\n")
                    f.write("\n"+bottom.replace("UUUUU",str(n)))
            progress = 100*n/num_geoms
            text = "{:2.1f}%".format(progress)
            print(' ', text, "of the geometries done.",end="\r", flush=True)
    
    print("\n\nDone! Ready to run.")   
###############################################################    

            
##COLLECTS RESULTS############################################## 
def gather_data(opc, tipo):
    files = [file for file in os.listdir('Geometries') if ".log" in file and "Geometry-" in file ]    
    files = sorted(files, key=lambda file: float(file.split("-")[1])) 
    with open("Samples.lx", 'w') as f:
        for file in files:
            num = file.split("-")[1]
            broadening = opc
            f.write("Geometry "+num+":  Vertical transition (eV) Oscillator strength Vibronic Shift (eV) Broadening Factor (eV) \n")
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
                if corrected != -1 and tipo == 'abs': #abspcm
                    vibronic = 0
                    f.write("Excited State 1:\t"+corrected+"\t"+fs[0]+"\t"+str(vibronic)+"\t"+str(broadening)+"\n")
                elif corrected != -1 and tipo == 'emi': #emipcm     
                    energy = str(np.round(total_corrected - scfs[-1],3))
                    vibronic = 0
                    f.write("Excited State 1:\t"+energy+"\t"+fs[0]+"\t"+str(vibronic)+"\t"+str(broadening)+"\n")
                elif corrected == -1 and tipo == 'emi':
                    vibronic = 0
                    f.write("Excited State "+numeros[0]+"\t"+energies[0]+"\t"+fs[0]+"\t"+str(vibronic)+"\t"+str(broadening)+"\n")
                else:
                    for i in range(len(energies)):
                        vibronic = 0
                        f.write("Excited State "+numeros[i]+"\t"+energies[i]+"\t"+fs[i]+"\t"+str(vibronic)+"\t"+str(broadening)+"\n")
                f.write("\n")   
############################################################### 


##NORMALIZED GAUSSIAN##########################################
def gauss(x,v,s):
    y =  (1/(np.sqrt(2*np.pi)*s))*np.exp(-0.5*((x-v)/s)**2)
    return y
###############################################################

##COMPUTES SPECTRA############################################# 
def spectra(tipo, num_ex, nr):
    if tipo == "abs":
        constante = (np.pi*(e**2)*hbar)/(2*nr*mass*c*epsilon0)*10**(20)
    elif tipo == 'emi':
        constante = ((nr**2)*(e**2)/(2*np.pi*hbar*mass*(c**3)*epsilon0))
    V, O, D, S = [], [], [], []
    N = 0
    with open("Samples.lx", 'r') as f:
        for line in f:
            if "Geometry" in line:
                N += 1
            elif "Excited State" in line and int(line.split()[2][:-1]) in num_ex:
                line = line.split()
                V.append(float(line[3]))
                O.append(float(line[4]))
                D.append(float(line[5]))
                S.append(float(line[6]))
    coms = [file for file in os.listdir(".") if 'Geometry-' in file and '.com' in file]
    if len(V) == 0 or len(O) == 0:
        fatal_error("You need to run steps 1 and 2 first! Goodbye!")
    elif len(V) != len(coms)*max(num_ex):
        print("Number of log files is less than the number of inputs. Something is not right! Computing the spectrum anyway...")
    V = np.asarray(V)
    O = np.asarray(O)
    D = np.asarray(D)
    S = np.asarray(S)
    if tipo == 'abs':
        espectro = (constante*O)
    else:
        espectro = (constante*(V**2)*O)
    x  = np.linspace(min(V)-3*max(np.sqrt(S)), max(V)+ 3*max(np.sqrt(S)), 200)
    y  = np.zeros((1,len(x)))
    if tipo == 'abs':
        arquivo = 'cross_section.lx'
        primeira = "{:4s} {:4s} {:4s}\n".format("#Energy(ev)", "cross_section(A^2)", "error")
    else:
        arquivo = 'differential_rate.lx'
        primeira = "{:4s} {:4s} {:4s}\n".format("#Energy(ev)", "diff_rate", "error")
    for i in range(0,len(espectro)):
        contribution = espectro[i]*gauss(x,V[i]+D[i],S[i])
        y  = np.vstack((y,contribution[np.newaxis,:]))

    y = y[1:,:]
    mean_y = np.mean(y,axis=0)
    #Error estimate
    sigma = np.std(y - mean_y,axis=0,ddof=1)/np.sqrt(len(mean_y))
    
    print(N, "geometries considered.")     
    with open(arquivo, 'w') as f:
        f.write(primeira)
        for i in range(0,len(x)):
            text = "{:.6f} {:.6e} {:.6e}\n".format(x[i],mean_y[i], sigma[i])
            f.write(text)
############################################################### 


def busca_input(freqlog):
    base = 'lalala'
    exc = False
    header = ''
    nproc = '4'
    mem   = '1GB'
    with open(freqlog, 'r') as f:
        search = False
        for line in f:
            if '%nproc' in line:
                line = line.split('=')
                nproc = line[-1].replace('\n','')
            elif '%mem' in line:
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
    if 'TD' in header.upper():
        exc = True
    header = header.split()
    base = ''
    for elem in header:
        if "/" in elem and 'IOP' not in elem.upper():
            base += elem.replace('#','')
        elif 'IOP' in elem.upper() and ('108' in elem or '107' in elem):
            base += ' '+elem
    return base, exc, nproc, mem                

def batch(gauss):
    files = [file for file in os.listdir(".") if 'Geometry-' in file and '.com' in file]
    logs  = [file for file in os.listdir(".") if 'Geometry-' in file and '.log' in file]
    prontos = []
    for file in logs:
        with open(file, 'r') as f:
            for line in f:
                if "Normal termination" in line:
                    prontos.append(file.split("-")[1])
    for file in sorted(files, key=lambda file: float(file.split("-")[1])):
        if file.split("-")[1] not in prontos:
            try:
                call(['ts', gauss, file])
            except:
                call(['tsp', gauss, file])

##CHECKS PROGRESS##############################################
def andamento():
    coms = [file for file in os.listdir("Geometries") if 'Geometry-' in file and '.com' in file]
    logs = [file for file in os.listdir("Geometries") if 'Geometry-' in file and '.log' in file]
    factor = 1
    with open('Geomtries/'+coms[0], 'r') as f:
        for line in f:
            if 'Link1' in line:
                factor = 2
    count = 0
    for file in logs:
        with open('Geometries/'+file, 'r') as f:
            for line in f:
                if "Normal termination" in line:
                    count += 1
    print("\n\nThere are", int(count/factor), "completed calculations out of", len(coms), "inputs")                
    print("It is", np.round(100*count/(factor*len(coms)),1), "% done.")                
###############################################################

##FETCHES LOG FILE#############################################
def busca_log(frase):                   
    files = [file for file in os.listdir('.') if ".log" in file and "Geometry-" not in file]
    if len(files) == 0:
        fatal_error("No frequency log found. Goodbye!")
    freqlog = 'nada0022'    
    for file in files:
        print("\n"+file)
        resp = input(frase+" y ou n?\n")
        if resp.lower() == 'y':
            freqlog = file
            break
    if freqlog == 'nada0022':
        fatal_error("No frequency log found. Goodbye!")
    return freqlog  
###############################################################    



print("#                       #     #")
print("#        ######   ####   #   # ")
print("#        #       #    #   # #  ")
print("#        #####   #    #    #   ")
print("#        #       #    #   # #  ")
print("#        #       #    #  #   # ")
print("#######  ######   ####  #     #")
print("----SPECTRA FOR THE PEOPLE!----\n")
print("Your options:\n")
print("I have frequency calculations done. I want to generate the inputs for the spectrum calculation - type 1")
print("My inputs are set, I want to run the spectrum calculations - type 2")
print("Calculations are done, I want to generate the spectrum - type 3")
print("I want to check the progess of the calculations - type 4")
print("I want to shake a molecule to help me get rid of imaginary frequencies - type 5")
op = input()
if op == '1':
    freqlog = busca_log("Is this the log file for the frequency calculation?")
    base, temtd, nproc, mem = busca_input(freqlog)
    print("\n"+base)
    resp = input("Are basis and functional correct? If so, pres Enter. Otherwise, type functional/basis.\n")
    if resp != "":
        base = resp 
    adicional = input("If there are extra keywords, type them. Otherwise, press Enter.\n")
    base += " "+adicional
    num_ex = input("How many excited states?\n")
    try:
        num_ex = int(num_ex)
    except:
        fatal_error("This must be a number! Goodbye!")
    print('%nproc='+nproc)    
    print('%mem='+mem)
    procmem = input('Are Nproc and Mem correct? y or n?\n')
    if procmem.lower() != 'y':
        nproc = input('nproc?\n')
        mem = input("mem?\n")
    num_geoms = int(input("How many geometries to be sampled?\n"))
    tda = 'TD'
    tamm = input('Use TDA (Tamm-Dancoff Approximation)? y or n?\n')
    if tamm.lower() == 's':
        tda = 'TDA'
    pcm = input("Include state specific solvent approach? y or n?\n")
    if pcm.lower() == 's':
        solv = input("What is the  solvent? If you want to specify the dielectric constants yourself, type read.\n")
        if solv.lower() == "read":
            eps1 = input("Type the static dielectric constant.\n")
            eps2 = input("Type the dynamic dielectric constant (n^2).\n")
            try:
                float(eps1)
                float(eps2)
            except:
                fatal_error("The constants must be numbers. Goodbye!")
            epss = "Eps="+eps1+"\nEpsInf="+eps2+"\n\n"
        else:
            solv = "SOLVENT="+solv
            epss = "\n"
        if temtd:
            print("Inputs suitable for emission spectra!\n")
            header = "%chk=step_UUUUU.chk\n%nproc="+nproc+"\n%mem="+mem+"\n# "+base+" "+tda+"=(NSTATES="+str(num_ex)+") SCRF=(CorrectedLR,NonEquilibrium=Save,"+solv+")\n\nTITLE\n\n0 1\n"
            bottom = epss+"\n--Link1--\n%nproc="+nproc+"\n%mem="+mem+"\n%oldchk=step_UUUUU.chk\n%chk=step2_UUUUU.chk\n# "+base+" GUESS=READ GEOM=CHECKPOINT SCRF(NonEquilibrium=Read,"+solv+")\n\nTITLE\n\n0 1\n\n"+epss
        else:
            print("Inputs suitable for absortion spectra!!\n")
            header = "%nproc="+nproc+"\n%mem="+mem+"\n# "+base+" SCRF=(CorrectedLR,"+solv+") "+tda+"=(NSTATES=3)\n\nTITLE\n\n0 1\n"
            bottom = epss
    elif pcm == 'n':
        header = "%nproc="+nproc+"\n%Mem="+mem+"\n# "+tda+"=(NStates="+str(num_ex)+") "+base+" \n\nTITLE\n\n0 1\n"
        bottom ="\n\n"
    else:
        fatal_error("It should be y or n. Goodbye!")
    T = float(input("Temperature in Kelvin?\n")) #K
    if T <= 0:
        fatal_error("Have you heard about absolute zero? Goodbye!")
    sample_geom(freqlog, num_geoms, T, header, bottom)    
elif op == '3':
    opc = input("What is the standard deviation of the gaussians? Should be typically around kT.\n")
    try:
        opc = float(opc)
    except: 
        fatal_error("It must be a number. Goodbye!")  
    print("What kind of spectrum?")
    tipo = input("Type abs (absorption) or emi (emission).\n")
    if tipo != 'abs' and tipo != 'emi':
        fatal_error("It must be either one. Goodbye!")
    estados = input("How many excited states?\n")
    try:
        estados = int(estados)
    except:
        fatal_error("It must be an integer! Goodbye!")
    nr = input("What is the refractive index?\n")
    try:
        nr = float(nr)
    except:
        fatal_error("It must be a number. Goodbye!")
    num_ex = range(0,estados+1)
    num_ex = list(map(int,num_ex))
    gather_data(opc, tipo)
    spectra(tipo, num_ex, nr)
elif op == '2':
    op = input("O ts está pronto já? s ou n?\n")
    if op != 's':
        fatal_error("Então vá aprontar o ts, animal!")
    gaussian = input("Which Gaussian? g09 ou g16?\n")
    batch(gaussian) 
elif op == '4':
    andamento()
elif op == '5':
    freqlog = busca_log("Is this the frequency calculation log file?")
    T = float(input("Magnitude of the displacement in Å? \n")) #K
    shake(freqlog,T)
else:
    fatal_error("It must be one of the options... Goodbye!")


    
        


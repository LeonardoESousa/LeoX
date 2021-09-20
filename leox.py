#!/usr/bin/env python3
import numpy as np
import os
import random
import sys
from decimal import Decimal


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
        print("No frequencies in the log file! Goodbye!") 
        sys.exit()
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
<<<<<<< HEAD:leox.py
        print("No geometry in the log file! Goodbye!")
=======
        print("Sem geometria no log de frequência! Adeus!")
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
        sys.exit()
    return G, atomos
###############################################################

##SAVES OPT GEOMETRY###########################################
def salva_geom(G,atomos):
    atomos = np.array([atomos]).astype(float)
    atomos = atomos.T
    G = np.hstack((atomos,G))
    np.savetxt('opt_geom.txt', G, delimiter='\t', fmt=['%1.1u','%+1.5f','%+1.5f','%+1.5f'])
    print("The optimized geometry that is used is saved in the opt_geom.txt file!")
###############################################################

<<<<<<< HEAD:leox.py
=======
def salva_geom(G,atomos):
    atomos = np.array([atomos]).astype(float)
    atomos = atomos.T
    G = np.hstack((atomos,G))
    np.savetxt('opt_geom.txt', G, delimiter='\t', fmt=['%1.1u','%+1.5f','%+1.5f','%+1.5f'])
    print("A geometria otimizada que vai ser usada está salva no arquivo opt_geom.txt!")


def pega_massas(freqlog,G):
    _ , M = pega_freq(freqlog)
    num_atom = np.shape(G)[0]
    massas = []
    with open(freqlog, 'r') as f:
        for line in f:
            if "has atomic number" in line:
                line = line.split()
                massas.append(float(line[-1]))
                massas.append(float(line[-1]))
                massas.append(float(line[-1]))
    massas = np.expand_dims(np.asarray(massas),axis=1)
    atomos = np.zeros((3*num_atom,1))
    for _ in range(0,len(M)):
        atomos = np.hstack((atomos,massas))
    massas = atomos[:,1:]
    return massas*amu
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py

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
    #NC = NC*np.sqrt(massas/M)
    #print(np.round(np.matmul(NC.T,NC),2))
    return NC

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
    #D = D*np.sqrt(massas/M)
    return D

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

def shake(freqlog, T):
    F, M = pega_freq(freqlog)
    G, atomos = pega_geom(freqlog)
    NNC = pega_modos(G,freqlog)
    num_atom = np.shape(G)[0]
    A1 = np.zeros((3*num_atom,1))
    A2 = np.zeros((3*num_atom,1))
    F = F[F < 0]
    if len(F) == 0:
<<<<<<< HEAD:leox.py
        print("No imaginary frquencies in the log file. Goodbye!")
=======
        print("Não há frequências imaginárias! Adeus!")
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
        sys.exit()
    F = -1*F
    for i in range(len(F)): # LL:
        q = [-1*T,T]
        A1 += q[0]*(np.expand_dims(NNC[:,i],axis=1))
        A2 += q[1]*(np.expand_dims(NNC[:,i],axis=1))
    A1 = np.reshape(A1,(num_atom,3))
    A2 = np.reshape(A2,(num_atom,3))
    Gfinal  = A1 + G
    Gfinal2 = A2 + G
    with open("shaken.xyz", 'w') as f:
<<<<<<< HEAD:leox.py
        f.write('#Geometry with displacement of '+str(T)+" A:\n" )
        for k in range(0, np.shape(Gfinal)[0]):
            text = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[k],Gfinal[k,0],Gfinal[k,1],Gfinal[k,2])
            f.write(text+"\n")
        f.write('\n#Geometry with displacement of '+str(-T)+" A:\n" )
        for k in range(0, np.shape(Gfinal2)[0]):
            text = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[k],Gfinal2[k,0],Gfinal2[k,1],Gfinal2[k,2])
            f.write(text+"\n")
    print("There are 2 geometries saved on shaken.xyz!")
=======
        f.write('#Geometria com deformação de '+str(T)+" A:\n" )
        for k in range(0, np.shape(Gfinal)[0]):
            text = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[k],Gfinal[k,0],Gfinal[k,1],Gfinal[k,2])
            f.write(text+"\n")
        f.write('\n#Geometria com deformação de '+str(-T)+" A:\n" )
        for k in range(0, np.shape(Gfinal2)[0]):
            text = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[k],Gfinal2[k,0],Gfinal2[k,1],Gfinal2[k,2])
            f.write(text+"\n")
    print("Há 2 geometrias salvas no arquivo shaken.xyz!")
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py


def sample_geom(freqlog, num_geoms, T, header, bottom):
    F, M = pega_freq(freqlog)
    if F[0] < 0:
<<<<<<< HEAD:leox.py
        print("Imaginary frequency! Goodbye!")
        sys.exit()
    try:
        os.mkdir('Geometries')
    except:
        pass        
=======
        print("Frequência negativa! Volte para casa 1")
        sys.exit()
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
    G, atomos = pega_geom(freqlog)
    salva_geom(G,atomos)
    NNC = pega_modos(G,freqlog)
    num_atom = np.shape(G)[0]   
<<<<<<< HEAD:leox.py
    print("\nGenerating geometries...\n")
    with open('Magnitudes.lx', 'w') as file:
=======
    print("\nGerando geometrias...\n")
    with open('Magnitudes.txt', 'w') as file:
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
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
<<<<<<< HEAD:leox.py
            with open("Geometries/Geometry-"+str(n)+"-.com", 'w') as f:
=======
            with open("Geometria-"+str(n)+"-.com", 'w') as f:
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
                    f.write(header.replace("UUUUU",str(n)))
                    for k in range(0, np.shape(Gfinal)[0]):
                        text = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[k],Gfinal[k,0],Gfinal[k,1],Gfinal[k,2])
                        f.write(text+"\n")
                    f.write("\n"+bottom.replace("UUUUU",str(n)))
            progress = 100*n/num_geoms
            text = "%2.1f" % progress
<<<<<<< HEAD:leox.py
            print(' ', text, "% of the geometries done.",end="\r", flush=True)
=======
            print(text, "% das geometrias feitas.",end="\r", flush=True)
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
    
    print("\n\nC'est fini! Ready to run.")   
    
    #Plots the histograms of normal coordinates
    #random.seed(10)
    #ind =0
    #for T in [0, 300]:
    #   for i in range(0,len(F)):
    #       x = np.linspace(-5, 5, 10000) #ja em angstrom
    #       boltz = np.tanh(hbar*F[i]/(2*kb*T))
    #       print(hbar*F[i])
    #       #print(hbar, F[i]/(2*np.pi*(10**12)), boltz, np.sinh(2*boltz), np.tanh(boltz))
    #       #prob1 = np.sqrt((M[i]*F[i])/(pi*hbar2))*np.exp(-M[i]*F[i]*((x*(10**(-10)))**2)/hbar2)*(abs(x[1]-x[0])*10**(-10)) #sem temperatura
    #       #prob2 = (((2*np.pi*kb*T*e)/(M[i]*(F[i]**2)))**(-1/2))*np.exp(-M[i]*(F[i]**2)*((x*(10**(-10)))**2)/(2*kb*T*e))*(abs(x[1]-x[0])*10**(-10))#    np.sqrt((M[i]*F[i])/(2*np.pi*hbar2*np.sinh(2*boltz)))*np.exp(-M[i]*F[i]*((x*(10**(-10)))**2)*np.tanh(boltz)/hbar2)*(abs(x[1]-x[0])*10**(-10))
    #       prob2 = np.sqrt((M[i]*F[i]*(boltz))/(np.pi*hbar2))*np.exp(-M[i]*F[i]*((x*(10**(-10)))**2)*(boltz)/hbar2)*(abs(x[1]-x[0])*10**(-10))
    #       #q = random.choices(x, prob1, k=1000)
    #       q2 = random.choices(x, prob2, k=1000)
    #       #plt.hist(q, bins=40)
    #       plt.subplot(2,3,ind+1)
    #       ind += 1
    #       plt.xlim([-0.3,0.3])
    #       plt.ylim([0,80])
    #       #plt.plot(x,prob2)
    #       plt.hist(q2, bins=40)
    #plt.show()
            

<<<<<<< HEAD:leox.py
def gather_data(opc, tipo):
    files = [file for file in os.listdir('Geometries') if ".log" in file and "Geometry-" in file ]    
    files = sorted(files, key=lambda file: float(file.split("-")[1])) 
=======
def gather_data(G, freqlog, opc, tipo):
    F, M = pega_freq(freqlog)
    NNC = pega_modos(G,freqlog)
    massas = pega_massas(freqlog,G)
    files = [file for file in os.listdir('.') if ".log" in file and "Geometria-" in file ]    
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
    with open("Samples.lx", 'w') as f:
        for file in files:
            num = file.split("-")[1]
            broadening = opc
<<<<<<< HEAD:leox.py
            f.write("Geometry "+num+":  Vertical transition (eV) Oscillator strength Vibronic Shift (eV) Broadening Factor (eV) \n")
=======
            f.write("Geometria "+num+":  Vertical transition (eV) Oscillator strength Vibronic Shift (eV) Broadening Factor (eV) \n")
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
            numeros, energies, fs, scfs = [], [], [], []
            corrected, total_corrected = -1, -1
            with open(file, 'r') as g:
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
                else:
                    for i in range(len(energies)):
                        vibronic = 0
                        f.write("Excited State "+numeros[i]+"\t"+energies[i]+"\t"+fs[i]+"\t"+str(vibronic)+"\t"+str(broadening)+"\n")
                f.write("\n")   


def gauss(x,v,s):
    y =  (1/(np.sqrt(2*np.pi)*s))*np.exp(-0.5*((x-v)/s)**2)
    return y


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
        print("You need to run steps 1 and 2 first! Goodbye!")
        sys.exit()
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
    x = np.linspace(min(V)-3*max(np.sqrt(S)), max(V)+ 3*max(np.sqrt(S)), 200)
    y = np.zeros((len(x),))
    if tipo == 'abs':
        arquivo = 'cross_section.lx'
        primeira = "%4s %4s" % ("#Energy(ev)", "cross_section(A^2)\n")
        for i in range(0,len(espectro)):
            y = y + espectro[i]*gauss(x,V[i]+D[i],S[i])
    else:
        arquivo = 'differential_rate.lx'
        primeira = "%4s %4s" % ("#Energy(ev)", "diff_rate\n")
        for i in range(0,len(espectro)):
            y = y + espectro[i]*gauss(x,V[i]+D[i],S[i])
    print(N, "geometries considered.")     
    y = y/float(N)
    with open(arquivo, 'w') as f:
        f.write(primeira)
        for i in range(0,len(x)):
            text = "%4.6f %8.6E\n" % (x[i], Decimal(y[i]))
            f.write(text)

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

def andamento():
    coms = [file for file in os.listdir("Geometries") if 'Geometry-' in file and '.com' in file]
    logs = [file for file in os.listdir("Geometries") if 'Geometry-' in file and '.log' in file]
    count = 0
    for file in logs:
        with open(file, 'r') as f:
            for line in f:
                if "Normal termination" in line:
                    count += 1
    print("\n\nThere are", count, "completed calculations out of", len(coms), "inputs")                
    print("It is", np.round(100*count/len(coms),1), "% done.")                

def busca_log(frase):                   
    files = [file for file in os.listdir('.') if ".log" in file and "Geometry-" not in file]
    if len(files) == 0:
        print("No frequency log found. Goodbye!")
        sys.exit()
    freqlog = 'nada0022'    
    for file in files:
        print("\n"+file)
        resp = input(frase+" y ou n?\n")
        if resp.lower() == 'y':
            freqlog = file
            break
    if freqlog == 'nada0022':
        print("No frequency log found. Goodbye!")
        sys.exit()
    return freqlog  
    
    
print("#                       #     #")
print("#        ######   ####   #   # ")
print("#        #       #    #   # #  ")
print("#        #####   #    #    #   ")
print("#        #       #    #   # #  ")
print("#        #       #    #  #   # ")
print("#######  ######   ####  #     #")
print("----SPECTRA FOR THE PEOPLE!----\n")
<<<<<<< HEAD:leox.py
print("Your options:\n")
print("I have frequency calculations done. I want to generate the inputs for the spectrum calculation - type 1")
print("My inputs are set, I want to run the spectrum calculations - type 2")
print("Calculations are done, I want to generate the spectrum - type 3")
print("I want to check the progess of the calculations - type 4")
print("I want to shake a molecule to help me get rid of imaginary frequencies - type 5")
op = input()
if op == '1':
    freqlog = busca_log("Is this the log file for the frequency calculation?")
=======
print("Quer fazer o que??\n")
print("Tenho só o log de frequência. Quero gerar geometrias - digite 1")
print("Geometrias prontas, quero botar para rodar com o ts - digite 2")
print("Tudo pronto, quero gerar o espectro - digite 3")
print("Quero saber a quantas anda essa joça - digite 4")
print("Quero sacudir uma molécula para me livrar de frequências imaginárias - digite 5")
op = input()
if op == '1':
    freqlog = busca_log("É esse o log de frequência?")
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
    base, temtd, nproc, mem = busca_input(freqlog)
    print("\n"+base)
    resp = input("Are basis and functional correct? If so, pres Enter. Otherwise, type functional/basis.\n")
    if resp != "":
        base = resp 
<<<<<<< HEAD:leox.py
    adicional = input("If there are extra keywords, type them. Otherwise, press Enter.\n")
=======
    adicional = input("Se houver comandos adicionais, digite-os. Caso contrário, aperte Enter.\n")
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
    base += " "+adicional
    num_ex = input("How many excited states?\n")
    try:
        num_ex = int(num_ex)
    except:
        print("This must be a number! Goodbye!")
        sys.exit()
    print('%nproc='+nproc)    
    print('%mem='+mem)
<<<<<<< HEAD:leox.py
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
=======
    procmem = input('Nproc e Mem estão corretos? s ou n?\n  ')
    if procmem != 's':
        nproc = input('nproc?\n')
        mem = input("mem?\n")
    num_geoms = int(input("Quantas geometrias?\n")) #numero de gometrias para gerar
    tda = 'TD'
    tamm = input('Usar TDA (Tamm-Dancoff Approximation)? s ou n?\n')
    if tamm == 's':
        tda = 'TDA'
    pcm = input("Incluir solvente (SS-PCM)? s ou n?\n")
    if pcm == 's':
        solv = input("Qual o solvente? Se quiser especificar as constantes dielétricas, digite read.\n")
        if solv == "read":
            eps1 = input("Digite a constante dielétrica estática.\n")
            eps2 = input("Digite a constante dielétrica dinâmica (n^2).\n")
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
            try:
                float(eps1)
                float(eps2)
            except:
                print("The constants must be numbers. Goodbye!")
                sys.exit()
            epss = "Eps="+eps1+"\nEpsInf="+eps2+"\n\n"
        else:
            solv = "SOLVENT="+solv
            epss = "\n"
        if temtd:
<<<<<<< HEAD:leox.py
            print("Inputs suitable for emission spectra!\n")
            header = "%chk=step_UUUUU.chk\n%nproc="+nproc+"\n%mem="+mem+"\n# "+base+" "+tda+"=(NSTATES=3) SCRF=(CorrectedLR,NonEquilibrium=Save,"+solv+")\n\nTITLE\n\n0 1\n"
            bottom = epss+"\n--Link1--\n%nproc="+nproc+"\n%mem="+mem+"\n%oldchk=step_UUUUU.chk\n%chk=step2_UUUUU.chk\n# "+base+" GUESS=READ GEOM=CHECKPOINT SCRF(NonEquilibrium=Read,"+solv+")\n\nTITLE\n\n0 1\n\n"+epss
        else:
            print("Inputs suitable for absortion spectra!!\n")
            header = "%nproc="+nproc+"\n%mem="+mem+"\n# "+base+" SCRF=(CorrectedLR,"+solv+") "+tda+"=(NSTATES=3)\n\nTITLE\n\n0 1\n"
=======
            print("Preparando inputs para espectro de emissão!\n")
            header = "%chk=stepUUUUU.chk\n%nproc="+nproc+"\n%mem="+mem+"\n# "+base+" "+tda+"=(NSTATES=1,EQSOLV) SCRF=(CorrectedLR,"+solv+") IOp(10/74=20)\n\nTITLE\n\n0 1\n"
            bottom = "NonEq=Write\n\n--Link1--\n%oldchk=stepUUUUU.chk\n%chk=step2UUUUU.chk\n# "+base+" GUESS=READ GEOM=CHECKPOINT SCRF("+solv+")\n\nTITLE\n\n0 1\n\nNONEQ=Read\n"+epss
        else:
            print("Preparando inputs para espectro de absorção!\n")
            header = "%nproc="+nproc+"\n%mem="+mem+"\n# "+base+" SCRF=(CorrectedLR,"+solv+") "+tda+"=(NSTATES=1,NonEqSolv) IOP(10/74=10)\n\nTITLE\n\n0 1\n"
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
            bottom = epss
    elif pcm == 'n':
        header = "%nproc="+nproc+"\n%Mem="+mem+"\n# "+tda+"=(NStates="+str(num_ex)+") "+base+" \n\nTITLE\n\n0 1\n"
        bottom ="\n\n"
    else:
        print("It should be y or n. Goodbye!.")
        sys.exit()
<<<<<<< HEAD:leox.py
    T = float(input("Temperature in Kelvin?\n")) #K
    if T <= 0:
        print("Have you heard about absolute zero? Goodbye!")
=======
    T = float(input("Temperatura em Kelvin?\n")) #K
    if T <= 0:
        print("Você está em violação da Terceira Lei da Termodinâmica... Aí não dá.")
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
        sys.exit()
    sample_geom(freqlog, num_geoms, T, header, bottom)    
elif op == '3':
    opc = input("What is the standard deviation of the gaussians? Should be typically around kT.\n")
    try:
        opc = float(opc)
    except: 
        print("It must be a number. Goodbye!!")
        sys.exit()  
    print("What kind of spectrum?")
    tipo = input("Type abs (absorption) or emi (emission).\n")
    if tipo != 'abs' and tipo != 'emi':
        print("It must be either one. Goodbye!")
        sys.exit()
    estados = input("How many excited states?\n")
    try:
        estados = int(estados)
    except:
        print("It must be an integer! Goodbye!")
        sys.exit()  
    nr = input("What is the refractive index?\n")
    try:
        nr = float(nr)
    except:
        print("It must be a number. Goodbye!")
        sys.exit()
    num_ex = range(0,estados+1)
    num_ex = list(map(int,num_ex))
<<<<<<< HEAD:leox.py
    gather_data(opc, tipo)
=======
    gather_data(G,freqlog, opc, tipo)
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py
    spectra(tipo, num_ex, nr)
elif op == '2':
    op = input("O ts está pronto já? s ou n?\n")
    if op != 's':
        print("Então vá aprontar o ts, animal!")
        sys.exit()
    gaussian = input("Which Gaussian? g09 ou g16?\n")
    batch(gaussian) 
elif op == '4':
    andamento()
elif op == '5':
<<<<<<< HEAD:leox.py
    freqlog = busca_log("Is this the frequency calculation log file?")
    T = float(input("Magnitude of the displacement in Å? \n")) #K
    shake(freqlog,T)
else:
    print("It must be one of the options... Goodbye!")
    sys.exit()


    
        

=======
    freqlog = busca_log("É esse o log de frequência?")
    T = float(input("Magnitude da deformação? (Algo entre 0.1 e 0.5 Angstrom)\n")) #K
    shake(freqlog,T)
else:
    print("Tem que ser um dos cinco, animal!")
    sys.exit()
>>>>>>> f6e4fc07a97286085a4954e159b92cb477ea906c:leox2.py

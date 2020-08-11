#!/usr/bin/env python3
#import matplotlib.pyplot as plt
import numpy as np
import os
import random
import sys
from subprocess import call
from decimal import Decimal
from scipy.linalg import orth

epsilon0 = 8.854187817*10**(-12) #F/m
hbar = 6.582119514*10**(-16) #eV s
hbar2 = 1.054571800*10**(-34) #J s
mass = 9.10938356*10**(-31) # kg
c = 299792458 #m/s
e = 1.60217662*10**(-19) #C
pi = np.pi
kb = 8.6173303*(10**(-5)) #eV/K
amu = 1.660539040*10**(-27) #kg


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
    
    F = np.asarray(F)*(c*100*2*pi) #converte em frequencia angular
    try:
        if F[0] < 0:
            print("Frequência negativa! Volte para casa 1")
            sys.exit()
    except:
        print("Log sem frequência! Volte para casa 1") 
        sys.exit()
    M = np.asarray(M)*amu # converte amu em kg
    return F, M

def pega_geom(freqlog):
    atomos = []
    if ".log" in freqlog:
        status = 0
        status2 = 0
        busca = "Input orientation:"
        with open(freqlog, 'r') as f:
            for line in f:
                if "Standard orientation" in line:
                    busca = "Standard orientation:"
                if 'geom=allcheckpoint' in line:
                    status2 = 1 
                if "Optimized Parameters" in line and status2 == 0:
                    status = 1
        print("\nEstou usando o", busca,"\n")           
        G = np.zeros((1,3))
        n = -1
        with open(freqlog, 'r') as f:
            for line in f:
                if status == 1 and "Optimized Parameters" in line:
                    status = 0
                if n < 0 and status == 0:
                    if busca in line:
                        n = 0
                    else:
                        pass
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
                    break       
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
    G = G[1:,:] #input geometry                 
    return G, atomos

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

def pega_modosHP(G, freqlog):
    #massas = pega_massas(freqlog,G)
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
    #massas = pega_massas(freqlog,G)
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

def sample_geom(freqlog, num_geoms, T, header, bottom):
    F, M = pega_freq(freqlog)
    G, atomos = pega_geom(freqlog)
    NNC = pega_modos(G,freqlog)
    num_atom = np.shape(G)[0]   
    print("\nGerando geometrias...\n")
    for n in range(1,num_geoms+1):
        A = np.zeros((3*num_atom,1))
        for i in range(0,len(F)):
            x = np.linspace(-5, 5, 10000) #ja em angstrom
            boltz = np.tanh(hbar*F[i]/(2*kb*T))
            prob = np.sqrt((M[i]*F[i]*(boltz))/(np.pi*hbar2))*np.exp(-M[i]*F[i]*((x*(10**(-10)))**2)*(boltz)/hbar2)*(abs(x[1]-x[0])*10**(-10)) #com temperatura
            q = random.choices(x, prob)
            A += q[0]*(np.expand_dims(NNC[:,i],axis=1))
        A = np.reshape(A,(num_atom,3))
        Gfinal = A + G  
        with open("Geometria-"+str(n)+"-.com", 'w') as f:
                f.write(header.replace("UUUUU",str(n)))
                for k in range(0, np.shape(Gfinal)[0]):
                    text = "%2s % 2.14f % 2.14f % 2.14f" % (atomos[k],Gfinal[k,0],Gfinal[k,1],Gfinal[k,2])
                    f.write(text+"\n")
                f.write("\n"+bottom.replace("UUUUU",str(n)))
        progress = 100*n/num_geoms
        text = "%2.1f" % progress
        print(text, "% das geometrias feitas.",end="\r", flush=True)
    
    print("\n\nC'est fini! Já pode botar pra rodar.")   
    
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
    

def rotaciona(G1,G,massas):
    atoms = [massas[:,1][3*i] for i in range(int(len(massas[:,1])/3))]
    G1 = G1[:,1:]
    cm_G1 = np.array([np.average(G1[:,0],weights=atoms), np.average(G1[:,1],weights=atoms), np.average(G1[:,2],weights=atoms)])
    cm_G = np.array([np.average(G[:,0],weights=atoms), np.average(G[:,1],weights=atoms), np.average(G[:,2],weights=atoms)])
    trans11 = cm_G1 
    trans22 = cm_G
    trans1, trans2 = cm_G1, cm_G 
    for i in range(np.shape(G)[0]):
        trans1 = np.vstack((trans1,trans11))
        trans2 = np.vstack((trans2,trans22))
    trans2 = trans2[1:,:]
    trans1 = trans1[1:,:]
    G = G - trans2
    G1 = G1 - trans1
    shifts = [1000000]
    for i in range(np.shape(G1)[0]):
        a = G[i,:] #vetor da geom final
        b = G1[i,:] #vetor da geom sampleada
        v = np.cross(a,b)/(np.sqrt(np.inner(a,a))*np.sqrt(np.inner(b,b)))
        s = np.sqrt(np.inner(v,v))
        c = np.inner(a,b)/(np.sqrt(np.inner(a,a))*np.sqrt(np.inner(b,b)))
        V = np.array([[0, -v[2], v[1]],[v[2], 0, -v[0]],[-v[1], v[0], 0]])
        if s != 0:
            R = np.identity(3) + V + np.matmul(V,V)*((1-c)/s**2)
        else:
            R = np.identity(3)
        NG = np.zeros([1,3])
        for j in range(np.shape(G1)[0]):
            vetor = np.matmul(R,G[j,:])
            NG = np.vstack((NG,vetor))
        NG = NG[1:,:]
        NG2 = (NG - G1)**2
        shift = (np.average(NG2.flatten('F')))#,weights=atoms1))
        shifts.append(shift)
        if shift <= min(shifts):
            GeomFinal = NG - G1 
    #print(min(shifts))
    #print(G1)
    #print(GeomFinal,"\n\n")
    #input("pausa")
    return GeomFinal
    
            

def gather_data(G, freqlog, opc):
    F, M = pega_freq(freqlog)
    NNC = pega_modos(G,freqlog)
    massas = pega_massas(freqlog,G)
    files = [file for file in os.listdir('.') if ".log" in file and "Geometria-" in file ]    
    with open("Samples.lx", 'w') as f:
        for file in files:
            num = file.split("-")[1]
            vibronic, broadening = 0, opc
            f.write("Geometria "+num+":  Vertical transition (eV) Oscillator strength Vibronic Shift (eV) Broadening Factor (eV) \n")
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
                if corrected != -1 and len(scfs) == 1: #abspcm
                    f.write("Excited State 1\t"+corrected+"\t"+fs[0]+"\t"+str(vibronic)+"\t"+str(broadening)+"\n")
                elif corrected != -1 and len(scfs) == 2: #emipcm     
                    energy = str(np.round(total_corrected - scfs[-1],3))
                    f.write("Excited State 1\t"+energy+"\t"+fs[0]+"\t"+str(vibronic)+"\t"+str(broadening)+"\n")
                else:
                    for i in range(len(energies)):
                        f.write("Excited State "+numeros[i]+"\t"+energies[i]+"\t"+fs[i]+"\t"+str(vibronic)+"\t"+str(broadening)+"\n")
                f.write("\n")   


def gauss(x,v,s):
    y =  (1/(np.sqrt(2*np.pi)*s))*np.exp(-0.5*((x-v)/s)**2)
    return y



def spectra(tipo, num_ex):
    if tipo == "abs":
        constante = (np.pi*(e**2)*hbar)/(2*mass*c*epsilon0)*10**(20)
    elif tipo == 'emi':
        constante = ((e**2)/(2*np.pi*hbar*mass*(c**3)*epsilon0))
    V, O, D, S = [], [], [], []
    N = 0
    with open("Samples.lx", 'r') as f:
        for line in f:
            if "Geometria" in line:
                N += 1
            elif "Excited State" in line and int(line.split()[2][:-1]) in num_ex:
                line = line.split()
                V.append(float(line[3]))
                O.append(float(line[4]))
                D.append(float(line[5]))
                S.append(float(line[6]))
    coms = [file for file in os.listdir(".") if 'Geometria-' in file and '.com' in file]
    if len(V) == 0 or len(O) == 0:
        print("Primeiro tem que fazer as opções 1 e 2 do programa, sua besta!")
        sys.exit()
    elif len(V) != len(coms)*max(num_ex):
        print("Número de logs é inferior ao número de inputs. Algo errado não está certo! Vou fazer o espectro mesmo assim.")
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
        primeira = "%4s %4s" % ("#Energy (ev)", "cross_section (A^2)\n")
        for i in range(0,len(espectro)):
            y = y + espectro[i]*gauss(x,V[i]+D[i],S[i])
    else:
        arquivo = 'differential_rate.lx'
        primeira = "%4s %4s" % ("#Energy (ev)", "diff_rate\n")
        for i in range(0,len(espectro)):
            y = y + espectro[i]*gauss(x,V[i]+D[i],S[i])
    print(N, "geometrias consideradas")     
    y = y/float(N)
    with open(arquivo, 'w') as f:
        f.write(primeira)
        for i in range(0,len(x)):
            text = "%4.6f %8.6E\n" % (x[i], Decimal(y[i]))
            f.write(text)

def busca_input(freqlog):
    base = 'lalala'
    exc = False
    with open(freqlog, 'r') as f:
        for line in f:
            if "#" in line:
                if 'TD' in line.upper():
                    exc = True
                line = line.split()
                for elem in line:
                    if "/" in elem:
                        base = elem
                        break
    return base, exc                 

def batch(gauss):
    files = [file for file in os.listdir(".") if 'Geometria-' in file and '.com' in file]
    logs = [file for file in os.listdir(".") if 'Geometria-' in file and '.log' in file]
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
    coms = [file for file in os.listdir(".") if 'Geometria-' in file and '.com' in file]
    logs = [file for file in os.listdir(".") if 'Geometria-' in file and '.log' in file]
    count = 0
    for file in logs:
        with open(file, 'r') as f:
            for line in f:
                if "Normal termination" in line:
                    count += 1
    print("\n\nHá", count, "TDs já prontos de um total de", len(coms), "inputs")                
    print("Está", 100*count/len(coms), "% completo")                

def busca_log(frase):                   
    files = [file for file in os.listdir('.') if ".log" in file and "Geometria-" not in file]
    if len(files) == 0:
        print("Sem log de Frequência.")
        sys.exit()
    freqlog = 'nada0022'    
    for file in files:
        print("\n"+file)
        resp = input(frase+" s ou n?\n")
        if resp == 's':
            freqlog = file
            break
    if freqlog == 'nada0022':
        print("Sem log de Frequência")
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
print("Quer fazer o que??\n")
print("Tenho só o log de frequência. Quero gerar geometrias - digite 1")
print("Geometrias prontas, quero botar para rodar com o ts - digite 2")
print("Tudo pronto, quero gerar o espectro - digite 3")
print("Quero saber a quantas anda essa joça - digite 4")
op = input()
if op == '1':
    freqlog = busca_log("É esse o log de frequência?")
    base, temtd = busca_input(freqlog)
    print("\n"+base)
    resp = input("Funcional e base estão corretos? Se sim, Enter. Se não, escreva (funcional/base).\n")
    if resp != "":
        base = resp 
    adicional = input("Se houver comandos adicionais (iops, por exemplo), digite-os. Caso contrário, aperte Enter.\n")
    base += " "+adicional
    num_ex = input("Quantos estados excitados?\n")
    try:
        num_ex = int(num_ex)
    except:
        print("Deu ruim!")
        sys.exit()
    nproc = input('nproc?\n')
    mem = input("mem?\n")
    num_geoms = int(input("Quantas geometrias?\n")) #numero de gometrias para gerar
    pcm = input("Incluir solvente (SS-PCM)? s ou n?\n")
    if pcm == 's':
        solv = input("Qual o solvente? Se quiser especificar as constantes dielétricas, digite read.\n")
        if solv == "read":
            eps1 = input("Digite a constante dielétrica estática.\n")
            eps2 = input("Digite a constante dielétrica dinâmica (n^2).\n")
            try:
                float(eps1)
                float(eps2)
            except:
                print("As constantes devem ser números.")
                sys.exit()
            epss = "Eps="+eps1+"\nEpsInf="+eps2+"\n\n"
        else:
            epss = "\n"
        if temtd:
            print("Preparando inputs para espectro de emissão!\n")
            header = "%chk=stepUUUUU.chk\n%nproc="+nproc+"\n%mem="+mem+"\n# "+base+" TD(NSTATES=1,EQSOLV) SCRF=(CorrectedLR,"+solv+") IOp(10/74=20)\n\nTITLE\n\n0 1\n"
            bottom = "NonEq=Write\n\n--Link1--\n%oldchk=stepUUUUU.chk\n%chk=step2UUUUU.chk\n# "+base+" GUESS=READ GEOM=CHECKPOINT SCRF("+solv+")\n\nTITLE\n\n0 1\n\nNONEQ=Read\n"+epss
        else:
            print("Preparando inputs para espectro de absorção!\n")
            header = "%nproc="+nproc+"\n%mem="+mem+"GB\n# "+base+" SCRF=(CorrectedLR,"+solv+") TD=(NSTATES=1,NonEqSolv) IOP(10/74=10)\n\nTITLE\n\n0 1\n"
            bottom = epss
    elif pcm == 'n':
        header = "%nproc="+nproc+"\n%Mem="+mem+"\n# td=(NStates="+str(num_ex)+") "+base+" \n\nTITLE\n\n0 1\n"
        bottom ="\n\n"
    else:
        print("É s ou n. Adeus.")
        sys.exit()
    T = float(input("Temperatura em Kelvin?\n")) #K
    if T == 0:
        print("Você está em violação da Terceira Lei da Termodinâmica... Aí não dá.")
        sys.exit()
    sample_geom(freqlog, num_geoms, T, header, bottom)    
elif op == '3':
    opc = input("Qual o desvio padrão das gaussianas? O padrão é usar kT.\n")
    try:
        opc = float(opc)
    except: 
        print("Tem que ser um número, diabo!")
        sys.exit()  
    freqlog = busca_log("É esse o log de frequência?")
    G, atm = pega_geom(freqlog)
    print("Espectro de que?")
    tipo = input("Digite abs (absorção) ou emi (emissão)\n")
    if tipo != 'abs' and tipo != 'emi':
        print("Tem que ser um dos dois, animal!")
        sys.exit()
    estados = input("Quantos estados excitados?\n")
    try:
        estados = int(estados)
    except:
        print("Fez alguma cagada!")
        sys.exit()  
    num_ex = range(0,estados+1)
    num_ex = list(map(int,num_ex))
    gather_data(G,freqlog, opc)
    spectra(tipo, num_ex)
elif op == '2':
    op = input("O ts está pronto já? s ou n?\n")
    if op != 's':
        print("Então vá aprontar o ts, animal!")
        sys.exit()
    gaussian = input("Qual o gaussian? g09 ou g16?\n")
    batch(gaussian) 
elif op == '4':
    andamento()
else:
    print("Tem que ser um dos quatro, animal!")
    sys.exit()


    
        

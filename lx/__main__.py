#!/usr/bin/env python3
import sys
import lx.tools
from lx.conf_search import classify_only

def interface():
    print("#                       #     #")
    print("#        ######   ####   #   # ")
    print("#        #       #    #   # #  ")
    print("#        #####   #    #    #   ")
    print("#        #       #    #   # #  ")
    print("#        #       #    #  #   # ")
    print("#######  ######   ####  #     #")
    print("--SPECTRA - EXCITON - TUNING--\n")
    print("Choose your option:\n")
    print('SPECTRUM SIMULATIONS:')
    print("\t1 - Generate the inputs for the spectrum calculation")
    print("\t2 - Run the spectrum calculations")
    print("\t3 - Generate the spectrum")
    print("\t4 - Check the progress of the calculations")
    print('EXCITON ANALYSIS:')
    print("\t5 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths")
    print('CONFORMATIONAL ANALYSIS:')
    print("\t6 - Perform conformational search")
    print('OTHER FEATURES:')
    print("\t7 - Perform long-range parameter tuning")
    print("\t8 - Get rid of imaginary frequencies")  
    print("\t9 - Retrieve last geometry from log file") 
    print("\t10 - Abort my calculations")
    op = input()
    if op == '1':
        freqlog = lx.tools.fetch_file("frequency",['.log'])
        F, _  = lx.tools.pega_freq(freqlog)
        if F[0] < 0:
            lx.tools.fatal_error("Imaginary frequency! Goodbye!")
        cm = lx.tools.get_cm(freqlog)
        base, temtd, nproc, mem, scrf, spec = lx.tools.busca_input(freqlog)
        if temtd == '':
            tda = 'TD'
        else:
            tda = temtd.upper()
        print('\nThe suggested configurations for you are:\n')
        print('Functional/basis: {}'.format(base))
        print('Solvent Method: {}'.format(scrf))
        print('Charge and Multiplicity: {}'.format(cm))
        print('TD-DFT or TDA-DFT: {}-DFT'.format(tda))
        print('%nproc='+nproc)    
        print('%mem='+mem)
        change = input('Are you satisfied with these parameters? y or n?\n')
        if change.lower() == 'n':     
            base  = lx.tools.default(base,"Functional/basis is {}. If ok, Enter. Otherwise, type functional/basis.\n".format(base))
            scrf  = lx.tools.default(scrf,"SCRF keyword is {}. If ok, Enter. Otherwise, type the desired one.\n".format(scrf))
            cm    = lx.tools.default(cm,'Charge and multiplicity is {}. If ok, Enter. Otherwise, type charge and multiplicity Ex.: 0 1\n'.format(cm))
            nproc = lx.tools.default(nproc,'nproc is {}. If ok, Enter. Otherwise, type it.\n'.format(nproc))
            mem   = lx.tools.default(mem,"mem is {}. If ok, Enter. Otherwise, type it.\n".format(mem))
            tamm  = input('Use TDA (Tamm-Dancoff Approximation)? y or n?\n')
            if tamm.lower() == 'y':
                tda = 'TDA'
            else:
                tda = 'TD'

        num_ex = input("How many excited states?\n")
        try:
            num_ex = int(num_ex)
        except:
            lx.tools.fatal_error("This must be a number! Goodbye!")
        num_geoms = int(input("How many geometries to be sampled?\n"))
        pcm = input("Include state specific solvent approach? y or n?\n")
        if pcm.lower() == 'y':
            solv = input("What is the solvent? If you want to specify the dielectric constants yourself, type read.\n")
            epss = lx.tools.set_eps(solv)
            if epss == '\n':
                solv = "SOLVENT="+solv
            if temtd:
                print("Inputs suitable for emission spectra!\n")    
                header = f"%chk=step_UUUUU.chk\n%nproc={nproc}\n%mem={mem}\n# {base} {tda}=(NSTATES={num_ex}) SCRF=(CorrectedLR,NonEquilibrium=Save,{solv})\n\n{spec}\n\n{cm}\n" 
                bottom = f"{epss}\n--Link1--\n%nproc={nproc}\n%mem={mem}\n%oldchk=step_UUUUU.chk\n%chk=step2_UUUUU.chk\n# {base} GUESS=READ GEOM=CHECKPOINT SCRF(NonEquilibrium=Read,{solv})\n\nTITLE\n\n{cm}\n\n{epss}" 
            else:
                print("Inputs suitable for absortion spectra!!\n")
                header = f"%nproc={nproc}\n%mem={mem}\n# {base} SCRF=(CorrectedLR,{solv}) {tda}=(NSTATES={num_ex})\n\n{spec}\n\n{cm}\n"
                bottom = epss
        elif pcm == 'n':
            epss = lx.tools.set_eps(scrf)
            header = f"%nproc={nproc}\n%Mem={mem}\n# {tda}=(NStates={num_ex}) {base} {scrf} \n\n{spec}\n\n{cm}\n"
            bottom = epss+'\n\n'
        else:
            lx.tools.fatal_error("It should be y or n. Goodbye!")
        T = float(input("Temperature in Kelvin?\n"))
        if T <= 0:
            lx.tools.fatal_error("Have you heard about absolute zero? Goodbye!")
        lx.tools.make_ensemble(freqlog, num_geoms, T, header, bottom)    
    elif op == '2':
        lx.tools.batch()
    elif op == '3':
        opc  = lx.tools.detect_sigma()
        tipo = lx.tools.get_spec()
        nr   = lx.tools.get_nr()
        print('The spectrum will be run with the following parameters:\n')
        print(f'Spectrum type: {tipo.title()}')
        print(f'Standard deviation of: {opc:.3f} eV')
        print(f'Refractive index: {nr:.3f}\n')
        change = input('Are you satisfied with these parameters? y or n?\n')
        if change.lower() == 'n':
            opc = input("What is the standard deviation of the gaussians?\n")
            try:
                opc = float(opc)
            except: 
                lx.tools.fatal_error("It must be a number. Goodbye!")
            nr = input("What is the refractive index?\n")
            try:
                nr = float(nr)
            except: 
                lx.tools.fatal_error("It must be a number. Goodbye!")
            tipo = input("What kind of spectrum? Type abs (absorption) or emi (emission)\n")
            if tipo != 'abs' and tipo != 'emi':
                lx.tools.fatal_error('It must be either one. Goodbye!')
        tipo = tipo[:3]
        if tipo == 'abs':
            estados = input("How many excited states?\n")
            try:
                estados = int(estados)
            except:
                lx.tools.fatal_error("It must be an integer! Goodbye!")
        else:
            estados = 1
        num_ex = range(0,estados+1)
        num_ex = list(map(int,num_ex))
        lx.tools.gather_data(opc, tipo)
        lx.tools.spectra(tipo, num_ex, nr)
    elif op == '4':
        lx.tools.andamento()
    elif op == '5':
        lx.tools.ld()
    elif op == '6':
        question = input('Classify only? y or n?\n')
        if question.lower() == 'y':
            try:
                classify_only()
            except:
                lx.tools.fatal_error('Something went wrong. Your folder may not contain Geomtetry log files. Goodbye!')
        else:
            lx.tools.conformational()
    elif op == '7':
        lx.tools.omega_tuning()
    elif op == '8':
        freqlog = lx.tools.fetch_file("Frequency log with imaginary frequencies",['.log'])
        lx.tools.distort(freqlog)
    elif op == '9':
        freqlog = lx.tools.fetch_file("log",['.log'])
        base, _, nproc, mem, scrf, _ = lx.tools.busca_input(freqlog)
        cm = lx.tools.get_cm(freqlog)
        header = f'%nproc={nproc}\n%mem={mem}\n# {base} {scrf}\n\nTITLE\n\n{cm}\n'
        G, atomos = lx.tools.pega_geom(freqlog)
        lx.tools.write_input(atomos,G,header,'','geom.lx')
        print('Geometry saved in the geom.lx file.')
    elif op == '10':
        lx.tools.abort_batch()
    else:
        lx.tools.fatal_error("It must be one of the options... Goodbye!")

def main():
    try:
        freqlog = sys.argv[1]
        G, atomos = lx.tools.pega_geom(freqlog)
        print(len(atomos))
        print('\n')
        for i in range(len(atomos)):
            print(f"{atomos[i]:2s}  {G[i,0]:.7f}  {G[i,1]:.7f}  {G[i,2]:.7f}")
    except:
        interface()

if __name__ == "__main__":
    sys.exit(main())
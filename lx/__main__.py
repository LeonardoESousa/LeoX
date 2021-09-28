#!/usr/bin/env python3
import sys
from lx.tools import *


def main():
    print("#                       #     #")
    print("#        ######   ####   #   # ")
    print("#        #       #    #   # #  ")
    print("#        #####   #    #    #   ")
    print("#        #       #    #   # #  ")
    print("#        #       #    #  #   # ")
    print("#######  ######   ####  #     #")
    print("----SPECTRA FOR THE PEOPLE!----\n")
    print("Choose your option:\n")
    print('SPECTRUM SIMULATIONS:')
    print("\t1 - Generate the inputs for the spectrum calculation")
    print("\t2 - Run the spectrum calculations")
    print("\t3 - Generate the spectrum")
    print("\t4 - Check the progress of the calculations")
    print('EXCITON ANALYSIS:')
    print("\t5 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths")
    print('OTHER FEATURES:')
    print("\t6 - Perform long-range parameter tuning") 
    print("\t7 - Retrieve last geometry from log file") 
    print("\t8 - Distort a molecule in the direction of imaginary normal modes")
    print("\t9 - Abort my calculations")
    op = input()
    if op == '1':
        freqlog = fetch_file("frequency",['.log'])
        cm = get_cm(freqlog)
        base, temtd, nproc, mem, scrf, spec = busca_input(freqlog)
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
            base  = default(base,"Functional/basis is {}. If ok, Enter. Otherwise, type functional/basis.\n".format(base))
            scrf  = default(scrf,"SCRF keyword is {}. If ok, Enter. Otherwise, type the desired one.\n".format(scrf))
            cm    = default(cm,'Charge and multiplicity is {}. If ok, Enter. Otherwise, type charge and multiplicity Ex.: 0 1\n'.format(cm))
            nproc = default(nproc,'nproc is {}. If ok, Enter. Otherwise, type it.\n'.format(nproc))
            mem   = default(mem,"mem is {}. If ok, Enter. Otherwise, type it.\n".format(mem))
            tamm  = input('Use TDA (Tamm-Dancoff Approximation)? y or n?\n')
            if tamm.lower() == 'y':
                tda = 'TDA'
            else:
                tda = 'TD'

        num_ex = input("How many excited states?\n")
        try:
            num_ex = int(num_ex)
        except:
            fatal_error("This must be a number! Goodbye!")
        num_geoms = int(input("How many geometries to be sampled?\n"))
        pcm = input("Include state specific solvent approach? y or n?\n")
        if pcm.lower() == 'y':
            solv = input("What is the solvent? If you want to specify the dielectric constants yourself, type read.\n")
            epss = set_eps(solv)
            if epss == '\n':
                solv = "SOLVENT="+solv
            if temtd:
                print("Inputs suitable for emission spectra!\n")    
                header = "%chk=step_UUUUU.chk\n%nproc={}\n%mem={}\n# {} {}=(NSTATES={}) SCRF=(CorrectedLR,NonEquilibrium=Save,{})\n\n{}\n\n{}\n".format(nproc,mem,base,tda,num_ex,solv,spec,cm) 
                bottom = "{}\n--Link1--\n%nproc={}\n%mem={}\n%oldchk=step_UUUUU.chk\n%chk=step2_UUUUU.chk\n# {} GUESS=READ GEOM=CHECKPOINT SCRF(NonEquilibrium=Read,{})\n\nTITLE\n\n{}\n\n{}".format(epss,nproc,mem,base,solv,cm,epss) 
            else:
                print("Inputs suitable for absortion spectra!!\n")
                header = "%nproc={}\n%mem={}\n# {} SCRF=(CorrectedLR,{}) {}=(NSTATES={})\n\n{}\n\n{}\n".format(nproc,mem,base,solv,tda,num_ex,spec,cm)
                bottom = epss
        elif pcm == 'n':
            epss = set_eps(scrf)
            header = "%nproc={}\n%Mem={}\n# {}=(NStates={}) {} {} \n\n{}\n\n{}\n".format(nproc,mem,tda,num_ex,base,scrf,spec,cm)
            bottom = epss+'\n\n'
        else:
            fatal_error("It should be y or n. Goodbye!")
        T = float(input("Temperature in Kelvin?\n"))
        if T <= 0:
            fatal_error("Have you heard about absolute zero? Goodbye!")
        sample_geom(freqlog, num_geoms, T, header, bottom)    
    elif op == '2':
        batch() 
    elif op == '3':
        opc = detect_sigma()
        tipo = get_spec()
        nr = get_nr() 
        print('The spectrum will be run with the following parameters:\n')
        print('Spectrum type: {}'.format(tipo.title()))
        print('Standard deviation of: {:.3f} eV'.format(opc))
        print('Refractive index: {:.3f}\n'.format(nr))
        change = input('Are you satisfied with these parameters? y or n?\n')
        if change.lower() == 'n':
            opc = input("What is the standard deviation of the gaussians?\n")
            try:
                opc = float(opc)
            except: 
                fatal_error("It must be a number. Goodbye!")  
            tipo = input("What kind of spectrum? Type abs (absorption) or emi (emission)\n")
            if tipo != 'abs' and tipo != 'emi':
                fatal_error('It must be either one. Goodbye!')
        else:
            tipo = tipo[:3]
            if tipo == 'abs':
                estados = input("How many excited states?\n")
                try:
                    estados = int(estados)
                except:
                    fatal_error("It must be an integer! Goodbye!")
            else:
                estados = 1
        num_ex = range(0,estados+1)
        num_ex = list(map(int,num_ex))
        gather_data(opc, tipo)
        spectra(tipo, num_ex, nr)
    elif op == '4':
        andamento()
    elif op == '5':
        ld()
    elif op == '6':
        omega_tuning()
    elif op == '7':
        freqlog = fetch_file("log",['.log'])
        base, _, nproc, mem, scrf, _ = busca_input(freqlog)
        cm = get_cm(freqlog)
        header = '%nproc={}\n%mem={}\n# {} {}\n\nTITLE\n\n{}\n'.format(nproc,mem,base,scrf,cm)
        G, atomos = pega_geom(freqlog)
        write_input(atomos,G,header,'geom.lx')
        print('Geometry saved in the geom.lx file.')    
    elif op == '8':
        freqlog = fetch_file("frequency",['.log'])
        base, temtd, nproc, mem, scrf, _ = busca_input(freqlog)
        cm = get_cm(freqlog)
        header = '%nproc={}\n%mem={}\n# {} FREQ=noraman {} {}\n\nTITLE\n\n{}\n'.format(nproc,mem,base,temtd,scrf,cm)
        T = float(input("Magnitude of the displacement in Å? \n")) #K
        shake(freqlog,T,header)
    elif op == '9':
        abort_batch()
    else:
        fatal_error("It must be one of the options... Goodbye!")


    
if __name__ == "__main__":
    sys.exit(main())        


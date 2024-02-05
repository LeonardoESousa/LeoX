# LeoX - Light emission and exciton diffusion in organic molecules 

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![license](https://img.shields.io/github/license/LeonardoESousa/LeoX?style=plastic)]()
[![down](https://img.shields.io/pypi/dm/LeoX)]()
[![maint](https://img.shields.io/maintenance/yes/2024)]()
[![commit](https://img.shields.io/github/last-commit/LeonardoESousa/LeoX?style=plastic)]()



A package for absorption and fluorescence spectrum simulations using the nuclear ensemble method along with TD(A)-DFT. Estimation of singlet exciton properties (Förster radius, lifetime, diffusion length). Long-range separation parameter tuning. Stochastic conformational search.  Interfaces with the Gaussian (09 or 16) package.


Table of Contents
=================
<!--ts-->
* [Cite as:](#cite-as)
* [What does this program do?](#what-does-this-program-do)
* [What is necessary to use it?](#what-is-necessary-to-use-it)
* [How to install it?](#how-to-install-it)
* [How to use it?](#how-to-use-it)
   
<!--te-->

## Cite as:

> de Sousa, L. E., Bueno, F. T., e Silva, G. M., da Silva Filho, D. A., & de Oliveira Neto, P. H. (2019). Fast predictions of exciton diffusion length in organic materials. Journal of Materials Chemistry C, 7(14), 4066-4071.


## What does this program do?

1.  Spectrum simulation:
    - Calculates Absorption and Fluorescence spectrum simulations using TD(A)-DFT.
    - Calculations include vibrational contributions to the spectra. 
    - Optionally, they may also include solvent effects either by PCM or by a state specific solvation model.
2.  Exciton properties:   
    - Calculates Förster radius for transfers between two molecules of equal or different type.
    - Calculates fluorescence lifetimes.
    - Calculates singlet exciton diffusion lengths.
3.  Extra features:
    - Tunes the w parameter of long-range separated functionals.
    - Extracts last geometry from Gaussian log file.
    - Runs a stochastic coformational search algorithm.


## What is necessary to use it?

 -  The program requires that the Gaussian quantum chemistry software (G09 or G16) be installed, as it interfaces with it.

-   The first step for running spectrum calculations or conformational searches is providing a Gaussian log file for a frequency calculation in the S0 state, if the goal is computing an absorption spectrum, or S1 state, if the objective is calculating a fluoresence spectrum. All frequencies must be real.  

-   To obtain estimates for Förster radius, fluorescence lifetimes and singlet exciton diffusion lengths, it is necessary to first perform both absorption and emission spectra calculations for the molecule of interest.

## How to install it?

Run:

`pip install LeoX`

Alternatively, clone the repository to your computer. Inside the LeoX folder, run:

`pip install .`

Depending on the operating system, the commands above may need to be replace by

`pip3 install LeoX` or `pip3 install .`

Once installed, you should be able to run the program from any folder by just using the `lx` command.

## How to use it?

1. For spectrum simulations:

    - Create a folder for your project. Add the log file for the frequency calculation to your folder. Run the lx command. Choose option 1 and follow the instructions to select the parameters of the calculation.
    - Add a bash script file to the folder to control execution. This file depends on which batch system you use. Examples for users of slurm or task spooler (ts) are presented [here](https://github.com/LeonardoESousa/LeoX/tree/master/batch_examples).
    - Run the `lx` command again, choose option 2 and follow the instructions. Once the calculations are running, you may use option 4 to check the progress or option 9 to abort.
    <img width="575" alt="LeoX1" src="https://user-images.githubusercontent.com/94139072/144780278-ef3b8ced-af82-4a9d-a8a5-667ebdcbdd71.png">

    - Once all calculations are done, run the `lx` command and choose option 3. Follow the instructions to set the parameters and the spectrum will be generated.
    <img width="436" alt="LeoX2" src="https://user-images.githubusercontent.com/94139072/144780487-7d0a5800-c925-4dbc-8b88-041b6d7f35b3.png">
    
    - After the spectrum is generated, a file with the extension .lx will be created. For absorption spectra, the file will be named "cross_section.lx", whereas for fluorescence spectra, the file will be named "differential_rate.lx"
    - It is possible to generate a spectrum with partially concluded calculations. Option 3 can be used multiple times without overwriting the previously generated spectrum. Each time a spectrum is generated, it will have a number at the beginning of the file name, such as "2cross_section.lx".
    - Each spectrum file is expected to have three columns as indicated below, where the third column can be used to estimate if the amount of sampled geometries are enough for accurate spectrum simulation.
    <img width="380" alt="LeoX3" src="https://user-images.githubusercontent.com/94139072/144781496-a7c1e1cc-56ea-4de4-85bd-6a0e7ed37d7e.png">


2. For exciton properties:

    - For exciton properties, you must first calculate the fluorescence and absorption spectra of the donor and acceptor molecules of interest to you. 
    - Once this is done, copy the spectra to a folder and inside this folder run the `lx` command. Choose option 5. Follow the instructions to set the calculation parameters, exemplified below.  
    <img width="450" alt="LeoX4" src="https://user-images.githubusercontent.com/94139072/144785795-70dec39b-c63d-41a7-9613-fcf4cd1c244f.png">
   
    - The correction for short distances takes into account the transition dipole moment of the donor molecule, extracted directly from the Gaussian log. Details can be obtained in the indicated published paper.
    - A file called `ld.lx` will be generated with all the information.
    <img width="450" alt="Captura de Tela 2021-12-06 às 01 13 13" src="https://user-images.githubusercontent.com/94139072/144786165-547a9a2e-cca3-434a-90a7-82c610d97d7e.png">


    - Importantly, diffusion length estimates are only sensible if donor and acceptor molecules are of the same kind. To calculate the Förster radius for transfers between different molecules, you must provide the fluorescence spectrum of the donor molecule and the absorption spectrum of the acceptor molecule. The fluorescence lifetime shown will correspond to that of the donor molecule.  

4. For conformational searches:
    
    - Create a folder. Add a frequency Gaussian .log file. Include also a bash script file according to your batch system (follow examples for slurm and task spooler [here](https://github.com/LeonardoESousa/LeoX/tree/master/batch_examples)). 
    - Run the `lx` command in the folder and choose option 6. Follow the instructions to set the calculation parameters. 
    - Answer the queries to determine the parameters of the conformational search (level of theory, initial temperature, temperature step, number of rounds, number of sampled geometries at each round). 
    - Once the search starts running, the file conformation.lx will show the conformations found thus far, along with the average energies and Boltzmann populations at 300 K of each. 
    - After the search is over, a folder named Conformers is created with the geometries of each conformer written to a gaussian input file that can be used for further optimization. 

3. For range separation parameter tuning:

    - Create a folder. Add either a Gaussian .log file or .com file (for any kind of calculation). Include also a bash script file according to your batch system (follow examples for slurm and task spooler [here](https://github.com/LeonardoESousa/LeoX/tree/master/batch_examples)). 
    - Run the `lx` command in the folder and choose option 7. Follow the instructions to set the calculation parameters. 
    - You may choose between a relaxed on unrelaxed tuning procedure. In the case of the former, geometry optimizations are run for each range separation parameter value. In the case of unrelaxed tuning, the geometry provided in the .log or .com file will be used for all calculations. 
 

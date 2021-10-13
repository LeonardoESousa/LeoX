# LeoX - Light-excited organic molecules: eXciton analysis package 

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![license](https://img.shields.io/github/license/LeonardoESousa/LeoX?style=plastic)]()
[![down](https://img.shields.io/github/downloads/LeonardoESousa/LeoX/total?style=plastic)]()
[![maint](https://img.shields.io/maintenance/yes/2021)]()
[![commit](https://img.shields.io/github/last-commit/LeonardoESousa/LeoX?style=plastic)]()


A package for absorption and fluorescence spectrum simulations using the nuclear ensemble method along with TD(A)-DFT. Estimation of singlet exciton properties (Förster radius, lifetime, diffusion length) and long-range separation parameter tuning. Interfaces with the Gaussian (09 or 16) package.


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
    - Extract slast geometry from Gaussian log file.
    - Distorts a molecule's geometry in the direction of imaginary normal modes.


## What is necessary to use it?

 -  The program requires that the Gaussian quantum chemistry software (G09 or G16) be installed, as it interfaces with it.

-   The first step for running spectrum calculations is providing a Gaussian log file for a frequency calculation in the S0 state, if the goal is computing an absorption spectrum, or S1 state, if the objective is calculating a fluoresence spectrum. All frequencies must be real.  

-   To obtain the estimates of Förster radius, fluorescence lifetimes and singlet exciton diffusion lengths, it is necessary to first perform both absorption and emission spectra calculations for the molecule of interest.

## How to install it?

Run:

`pip install LeoX`

Alternatively, clone the repository to your computer. Inside the LeoX folder, run:

`pip install .`

Once installed, you should be able to run the program from any folder by just using the `lx` command.

## How to use it?

1. For spectrum simulations:

    - Create a folder for your project. Add the log file for the frequency calculation to your folder. Run the lx command. Choose option 1 and follow the instructions to select the parameters of the calculation.
    - Add a bash script file to the folder to control execution. This file depends on which batch system you use. Examples for users of slurm or task spooler (ts) are presented [here](https://github.com/LeonardoESousa/LeoX/tree/master/batch_examples).
    - Run the `lx` command again, choose option 2 and follow the instructions. Once the calculations are running, you may use option 4 to check the progress or option 9 to abort.
    - Once all calculations are done, run the `lx` command and choose option 3. Follow the instructions to set the parameters and the spectrum will be generated.

2. For exciton properties:

    - For exciton properties, you must first calculate the fluorescence and absorption spectra of the donor and acceptor molecules of interest to you. 
    - Once this is done, copy the spectra to a folder and inside this folder run the `lx` command. Choose option 5. Follow the instructions to set the calculation parameters. A file called `ld.lx` will be generated with all the information. 
    - Importantly, diffusion length estimates are only sensible if donor and acceptor molecules are of the same kind. To calculate the Förster radius for transfers between different molecules, you must provide the fluorescence spectrum of the donor molecule and the absorption spectrum of the acceptor molecule. The fluorescence lifetime shown will correspond to that of the donor molecule.  

3. For range separation tuning:

    - Create a folder. Add either a Gaussian .log file or .com file (for any kind of calculation). Include also a bash script file according to your batch system (follow examples for slurm and task spooler [here](https://github.com/LeonardoESousa/LeoX/tree/master/batch_examples)). 
    - Run the `lx` command in the folder and choose option 6. Follow the instructions to set the calculation parameters. 
    - You may choose between a relaxed on unrelaxed tuning procedure. In the case of the former, geometry optimizations are run for each range separation parameter value. In the case of unrelaxed tuning, the geometry provided in the .log or .com file will be used for all calculations. 
 
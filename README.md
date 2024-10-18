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
   
<!--te-->

## Cite as:

> de Sousa, L. E., Bueno, F. T., e Silva, G. M., da Silva Filho, D. A., & de Oliveira Neto, P. H. (2019). Fast predictions of exciton diffusion length in organic materials. Journal of Materials Chemistry C, 7(14), 4066-4071.


## What does this program do?

1.  Spectrum simulation (DEPRECATED - USE [NEMO](https://github.com/LeonardoESousa/NEMO) INSTEAD):
    - Calculates Absorption and Fluorescence spectrum simulations using TD(A)-DFT.
    - Calculations include vibrational contributions to the spectra. 
    - Optionally, they may also include solvent effects either by PCM or by a state specific solvation model.
2.  Exciton properties (DEPRECATED - USE [NEMO](https://github.com/LeonardoESousa/NEMO) INSTEAD):   
    - Calculates Förster radius for transfers between two molecules of equal or different type.
    - Calculates fluorescence lifetimes.
    - Calculates singlet exciton diffusion lengths.
3. Conformational Search:
    - Runs a stochastic coformational search algorithm. See [Tutorial](https://github.com/LeonardoESousa/LeoX/blob/master/Tutorial/Tutorial.md)
4.  Extra features:
    - Tunes the w parameter of long-range separated functionals.
    - Extracts last geometry from Gaussian log file.
    


## What is necessary to use it?

 -  The program requires that the Gaussian quantum chemistry software (G09 or G16) be installed, as it interfaces with it.

## How to install it?

The easiest way to install is to use pip:

`pip install LeoX`

This will install the latest released version.

To install the version with the latest commit, run:

`pip install git+https://github.com/LeonardoESousa/LeoX`

Alternatively, clone the repository to your computer. 

`git clone https://github.com/LeonardoESousa/LeoX`

A new folder named **LeoX** will appear. Move to this folder and install using pip:

```
cd LeoX
pip install .
```

Once installed, you should be able to run the program from any folder in your computer by just using the `lx` command.

## How to use it?

See [Tutorial](https://github.com/LeonardoESousa/LeoX/blob/master/Tutorial/Tutorial.md) for conformational search.
 

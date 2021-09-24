# LeoX

## What does this program do?

 - Absorption and Fluorescence spectrum simulations using TD(A)-DFT including vibrational contributions.
 - Estimates Förster radius, fluorescence lifetimes and singlet exciton diffusion lengths.

  coupled to the Gaussian quantum chemistry package. 

## What is necessary to use it?

The program that the Gaussian quantum chemistry software (G09 or G16) be installed, as it interfaces with it.
The first step for spectrum calculations is providing a Gaussian log file for a frequency calculation in the S0 state, if the goal is computing an absorption spectrum, or S1 state, if the objective is calculating a fluoresence spectrum. All frequencies must be real.  

To obtain the estimates of Förster radius, fluorescence lifetimes and singlet exciton diffusion lengths, it is necessary to first perform both absorption and emission spectra calculations for the molecule of interest.

## How to install?

Clone the repository to any folder in your computer. In your terminal, navigate to the LeoX folder and add permission to execute by doing:

`chmod +x leox`

Now, open your .bashrc file and include the following:

`export PATH="/path/to/folder/LeoX:$PATH"`

In your terminal type:

`source ~/.bashrc`

Now you can run the program from any folder by just using the leox.py command.

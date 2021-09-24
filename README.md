# LeoX

## What does this program do?

 - Absorption and Fluorescence spectrum simulations using TD(A)-DFT.
 - Calculations include vibrational contributions to the spectra. Optionally, it may also include solvent effects either by PCM or by a state specific solvation model.
 - Estimates Förster radius, fluorescence lifetimes and singlet exciton diffusion lengths.

## What is necessary to use it?

The program requires that the Gaussian quantum chemistry software (G09 or G16) be installed, as it interfaces with it.
The first step for running spectrum calculations is providing a Gaussian log file for a frequency calculation in the S0 state, if the goal is computing an absorption spectrum, or S1 state, if the objective is calculating a fluoresence spectrum. All frequencies must be real.  

To obtain the estimates of Förster radius, fluorescence lifetimes and singlet exciton diffusion lengths, it is necessary to first perform both absorption and emission spectra calculations for the molecule of interest.

## How to install?

Run:

`pip install LeoX`

Now you can run the program from any folder by just using the `lx` command.

## How to use it?

Create a folder for your project. Add the log file for the frequency calculation in your folder. Run the lx command. Choose option 1 and follow the instructions.



# LeoX

## What does this program do?
1.  Spectrum simulation:
    - Absorption and Fluorescence spectrum simulations using TD(A)-DFT.
    - Calculations include vibrational contributions to the spectra. Optionally, it may also include solvent effects either by PCM or by a state specific solvation model.
2.  Exciton properties:   
    - Estimate Förster radius for transfers between two molecules of equal or different type.
    - Estimate fluorescence lifetimes.
    - Estimate singlet exciton diffusion lengths.
3.  Extra features:
    - Tune the parameter of long-range separated functionals.
    - Extract last geometry from Gaussian log file.
    - Distort a molecule's geometry in the direction of imaginary normal modes.

## What is necessary to use it?

 -  The program requires that the Gaussian quantum chemistry software (G09 or G16) be installed, as it interfaces with it.

-   The first step for running spectrum calculations is providing a Gaussian log file for a frequency calculation in the S0 state, if the goal is computing an absorption spectrum, or S1 state, if the objective is calculating a fluoresence spectrum. All frequencies must be real.  

-   To obtain the estimates of Förster radius, fluorescence lifetimes and singlet exciton diffusion lengths, it is necessary to first perform both absorption and emission spectra calculations for the molecule of interest.

## How to install?

Run:

`pip install LeoX`

Alternatively, clone the repository to your computer. Inside the LeoX folder, run:

`pip install .`

Onde installed, you should be able to run the program from any folder by just using the `lx` command.

## How to use it?

Create a folder for your project. Add the log file for the frequency calculation in your folder. Run the lx command. Choose option 1 and follow the instructions.



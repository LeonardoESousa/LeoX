# Tutorial for Conformational Search

## Installing the package

The easiest way to install is to use pip:

`pip install LeoX`

This will install the lates released version.

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

# Preliminary Steps

Start by creating a folder. Let's call this folder **Search**.

The first step before running a conformation search is running an opt freq calculation using an initial geometry of your choice. Since the search requires running many calculation, it is advisable to use a cheap computational method. In the exampl below we will use the pm6 method for the initial search. Later, we will use more accurated methods for refining the results.

The opt freq input file, which we will call **freq.com**, should then look something like the input below:

```
%nproc=40
%mem=100GB
# pm6 opt freq=noraman 

TITLE

0 1
#ADD XYZ COORDINATES HERE
```

Run this opt freq calculation and make sure that the optimized geometry is converged.

Next, we need to create a submission script, which we name **g16.sh**. The following is an example that first sets up some options and loads the appropriate modules. Adjust it according to the cluster you use. The essential line is the `bash $1` line.

```
#!/bin/bash
#SBATCH --mem=100GB
#SBATCH --time=1-0
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --partition=xeon40el8

module load Gaussian/16


bash $1
```

Now we create another script in a fille named **batch.sh**. It should contain the following in the case of a slurm queue system:

```
sbatch g16.sh $1
```

Your folder structure should look like the following:


```
    Search
    ├── freq.com
    ├── freq.log
    ├── g16.sh
    └── batch.sh
```

```
 ▄█          ▄████████  ▄██████▄  ▀████    ▐████▀
███         ███    ███ ███    ███   ███▌   ████▀
███         ███    █▀  ███    ███    ███  ▐███
███        ▄███▄▄▄     ███    ███    ▀███▄███▀
███       ▀▀███▀▀▀     ███    ███    ████▀██▄
███         ███    █▄  ███    ███   ▐███  ▀███
███▌    ▄   ███    ███ ███    ███  ▄███     ███▄
█████▄▄██   ██████████  ▀██████▀  ████       ███▄
▀
Version: 0.8.1

Choose your option:

SPECTRUM SIMULATIONS:
        1 - Generate the inputs for the spectrum calculation
        2 - Run the spectrum calculations
        3 - Generate the spectrum
        4 - Check the progress of the calculations
EXCITON ANALYSIS:
        5 - Estimate Förster radius, fluorescence lifetime and exciton diffusion lengths
CONFORMATIONAL ANALYSIS:
        6 - Perform conformational search
OTHER FEATURES:
        7 - Perform long-range parameter tuning
        8 - Get rid of imaginary frequencies
        9 - Retrieve last geometry from log file
        10 - Abort my calculations
```

Select option **6** and enter.

```
6
Classify only? y or n?
```

Classify only will be used after the search is concluded. So now, answer n.

```
n
This is the configuration taken from the file:

Functional/basis: pm6
%nproc=40
%mem=100GB
Initial Temperature: 352 K
Temperature step: 35 K
Are you satisfied with these parameters? y or n?
```

The program fetches the level of theory, nproc and mem values from the **freq.log** file present in the folder. If you agree with the parameters, type y and enter. Otherwise, you type n and make the changes you want.


```
y
Number of geometries sampled at each round?
```

We will go here with 10 geometries, but feel free to choose as many as you want. The more geometries, the better. However, the whole procedure will also take longer.

```
10
Number of rounds?
```

Same issue when it comes to the number of rounds. Again we go with 10.

```
10
Number of jobs in each batch?
```

Here we select how many jobs in each batch. For example, if you select 1, then 10 jobs will be submitted at each round. If you select 10, a single job with 10 calculations will be submitted with calculations run in sequence.  

```
10
g16 or g09?
```

Select your gaussian version

```
g16
```

Now the conformational search is on its way.

## Checking Results


As the conformational search proceeds, a file named **conformation.lx** will be generated. It will look something as follows:

```
#Group  Energy(eV)  DeltaE(eV)  Prob@300K(%)  Rot1        Rot2        Rot3        Std1        Std2        Std3        Number  Last
1       1.864       0.000       20.9          0.0246947   0.0097615   0.0093145   0.0000247   0.0000098   0.0000093   1       38
2       1.876       0.012       13.4          0.0225811   0.0096333   0.0078288   0.0000226   0.0000096   0.0000078   1       37
3       1.876       0.012       13.1          0.0228896   0.0091604   0.0075435   0.0000229   0.0000092   0.0000076   5       25

#Round 1/5 Temperature: 350.0 K
```


The first column **Group** lists the conformations found so far. The column **Energy(eV)** gives the average energy of the geometries classified in that conformation. The columns **DeltaE(eV)** shows the energy difference with respect to the lowest energy conformation found so far. 

The next column **Prob@300K(%)** gives the Boltzmann population at 300K for each conformation. The next 3 columns, **Rot1**,**Rot2** and **Rot3** give the average values of the rotational constants. Similarly, **Std1**, **Std2** and **Std3** give the corresponding standard deviation.

Finally, the column **Number** shows how many geometries were classified in that category and the colum **Last** shows the number of the last geometry classified as belonging to that conformation.

At the end of the file, information about the round progression and temperature used can be found.


## Refining Results

Once the conformational search is completed. A folder named **Conformers** will be generated. Inside it, a number of Gaussian input files will be present, one for each conformation identified. These are input files for optimization. At this point, you should select a more accurate level of theory, add it to each of the input files and run these optimizations. Once they are completed, make sure the geometries have converged and run the **lx** command again. Once more, select option **6**. Now, when asked "Classify only?", answer y.

A new **conformation.lx** file will be generated, with parameters taken from the more accurate level of theory. In this procedure, conformations that were originally different may collapse to the same category as a result. Now, the procedure is completed.


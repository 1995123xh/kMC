# kMC
A kinetic Monte-Carlo model for predicting isotopic structures of small organic molecules in catagensis

### Prerequisite
Install Matlab. Version needs to be newer than R2015b. Add the 'parellel computing toolbox' package when prompt during installation. You can type the following Matlab command to check if your installation has it:
```bash
parpool()
```

### Input
Every input (e.g. 'kIA.mat') is a binary data file that contains two types of datas:
1. 'atoms' is a character array of all the atoms in the molecule.
2. 'BoM' that is the bond order matrix showing the bondings between atoms. If there is a double bond between atom#6 and atom#7, BoM(6,7)==BoM(7,6)==2. 1.5 is used for aromatic bonds. Note that this matrix has to be complete instead of an upper triangle. 

### Run
To run, change the input file name in 'MC_isotopologue_seperategraphs_multibin_newcollector.m' and execute it.

### Data processing
Process the saved results file with:
```bash
parseoutput(FILENAME)
```

### Contact
hxie@caltech.edu

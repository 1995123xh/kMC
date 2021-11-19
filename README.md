# kMC
A kinetic Monte-Carlo model for predicting isotopic structures of small organic molecules in catagensis, presented in two papers:           <em>
Xie, H., Ponton, C., Formolo, M. J., Lawson, M., Ellis, G. S., Lewan, M. D., ... & Eiler, J. M. (2020). Position-specific distribution of hydrogen isotopes in natural propane: Effects of thermal cracking, equilibration and biodegradation. Geochimica et Cosmochimica Acta, 290, 235-256.</em>  
<em>Xie, H., Formolo, M. J., & Eiler, J. M. (2022). Position-specific distribution of hydrogen isotopes in natural propane: Effects of thermal cracking, equilibration and biodegradation. Submitted to Geochimica et Cosmochimica Acta, 290, 235-256.</em>  
Two versions of the model are stored in different folders, for each reference respectively. 

### Prerequisite
Install Matlab. Version needs to be newer than R2015b. Add the 'parellel computing toolbox' package when prompt during installation. You can type the following Matlab command to check if your installation has it:
```bash
parpool()
```

### Input
Every input (e.g. 'kIA.mat') is a binary data file that contains two types of datas:
1. 'atoms' is a character array of all the atoms in the molecule.
2. 'BoM' that is the bond order matrix showing the bondings between atoms. If there is a double bond between atom#6 and atom#7, BoM(6,7)==BoM(7,6)==2. 1.5 is used for aromatic bonds. Note that this matrix has to be complete instead of an upper triangle. 
3. 'ringinfo' is an array that describes the smallest memeber of rings each atom is in.

### Run
To run, change the input file name in 'MC_isotopologue_seperategraphs_multibin_newcollector.m' (or 'MC_2020ThermalProgram.m' for scheme A of the 2022 paper or 'MC_2020ThermalProgram_bondM.m' for scheme B of the 2022 paper) and execute it. Results will be automatically saved.

### Data processing
Process the saved results file with:
```bash
parseoutput(FILENAME)
```

### Contact
hxie@caltech.edu

**"Supervised Machine Learning to predict the Energy and Forces of Small Organic Molecular Systems"**


**OBJECTIVE**
We aim to predict the total energy and the forces for molecular systems utilizing as training data output files 
from dynamics run with Qchem Software. _(ref QChem)_


**Data Processing Step by Step**

Dynamics data from QChem:
We run _Ab Initio_ Molecular dynamics using QChem, we study the evolution of a molecular system in time. 
In our case, we are giving initial thermal energy which value is chosen accordingly to the problem we are interested on.
e.g. in the Interstellar Medium (ISM) the temperatures are low, say 50K. 

**1) Raw data (QChem output) to extended xyz**

The conversion of a set of outputs to a single extendedxyz file that merges all the dynamics for a particular system is done executing:

**python xyz_maker.py list_of_outputs.out name_of_merged_file_extended.xyz**

the **list_of_outputs.out** are QChem outputs from Molecular Dynamics that contain for each time step:
1)  "Total energy in the final basis set": Total energy of the system in Hartree. 
2)  "Standard Nuclear Orientation": Cartesian coordinates of the system in angstroms.
3)  "Gradient of SCF Energy": The gradient in hartree per Bohr.

The units are convertedto get the extended xyz file **name_of_merged_file_extended.xyz** in units relevant to work on sGDML _(refsGDML)_

1) Energy: Hartree to eV (factor:  27.211386245988)
2) Cartesian coordinates: no modification done
3) Gradient: Hartree per Bohr (factor: 51.422086190832). We need the Forces, thus we additionally multiply by -1.

The resulting extended xyz file format is as follows:  <br/>

number of atoms <br/>
energy <br/>
atom<sub>1</sub>id    x<sub>1</sub> y<sub>1</sub> z<sub>1</sub> force_x<sub>1</sub> force_y<sub>1</sub>   force_z<sub>1</sub>  <br/>
... <br/>
...  <br/>
atom<sub>n</sub>id    x<sub>n</sub> y<sub>n</sub> z<sub>n</sub> force_x<sub>n</sub> force_y<sub>n</sub>   force_z<sub>n</sub>  <br/>

**2) From extended xyz to npz binary file that is used to train sGDML model** 
This step is done using a py script from sGMDL: 

sgdml_dataset_from_extxyz.py **name_of_merged_file_extended.xyz** 

This step creates a file  **name_of_merged_file_extended.npz**


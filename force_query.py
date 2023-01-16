#########################################################################################
#                                                                                       #
# Author: Estefania Rossich Molina                                                      #
# Jan 2023                                                                              #
#                                                                                       #
# This script uses sGDML to predict the forces reading the coordinates from MD in Qchem #
# It is executed by force_query.sh                                                      #
#                                                                                       #
#########################################################################################

import numpy as np
from sgdml.predict import GDMLPredict
from sgdml.utils import io
import os #needed to import variables from bash!  

trained_model = 'h2o_10k_20trajs-unknown-train2000-sym2.npz' #this name needs to be introduced manually
model = np.load(trained_model) #this needs to be done only once

print ("Model name", trained_model) #estefi
gdml = GDMLPredict(model)

#print("***")
#print(str(os.environ["TESTF"]))
#print("***")

r,_ = io.read_xyz(str(os.environ["TESTF"])) #3 atoms #this file is generated from the QChem MD output at each time step 
e,f = gdml.predict(r)

f = f.reshape (3,3)  
print("Predicted sGDML Forces:\n") 
print (f) # (1,9) if n_atoms = 3.

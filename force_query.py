#to execute this script, the environment file needs to be first sourced, so, 
#we don't execute this directly, we do it by the force_query.sh file 
#doing ./force_query.sh 

import numpy as np
from sgdml.predict import GDMLPredict
from sgdml.utils import io

trained_model = 'h2o_10k_20trajs-unknown-train2000-sym2.npz' #this name needs to be introduced manually
model = np.load(trained_model) #this needs to be done only once

print ("Model name", model) #estefi
gdml = GDMLPredict(model)

r,_ = io.read_xyz('h2o.xyz') #3 atoms #this file is generated from the QChem MD output at each time step 
e,f = gdml.predict(r)

f = f.reshape (3,3)  
print("ML predicted forces") 
print (f) # (1,9) if n_atoms = 3.

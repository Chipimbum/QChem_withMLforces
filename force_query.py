import numpy as np
from sgdml.predict import GDMLPredict
from sgdml.utils import io

trained_model = 'h2o_10k_20trajs-unknown-train2000-sym2.npz'
model = np.load(trained_model) #this needs to be done only once

print ("Model name", model) #estefi
gdml = GDMLPredict(model)

r,_ = io.read_xyz('h2o.xyz') #3 atoms
e,f = gdml.predict(r)

f = f.reshape (3,3)  
print("ML predicted forces") 
print (f) # (1,9) if n_atoms = 3.

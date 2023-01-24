#!/bin/bash

###################################################################################
#                                                                                 #
# Author: Estefania Rossich Molina                                                #
# Jan 2023                                                                        #
#                                                                                 #
# This bash interpreter allows the execution of sGDML (in Python) along with      #
# the execution of QChem (in C and fortran)                                       #
#                                                                                 #
#                                                                                 #
###################################################################################

source /usr/people/tamars/fanirm/anaconda_spike/setup_anaconda_spike.sh #Python environment needed for sGDML
source /home/SCR/setup/qcsetup_6.0.sh                                   #Qchem path in Spike

HOME_DIR=`pwd`
export QCSCRATCH=/home/scr/$USER

QCIN=$1 #name of the input file needs to be given by the keyboard 
TESTF=$1.xyz #system variable that will be exported. It contains the geom of the molecule to predict forces using a trained model on sGDML

if [ -z "$QCIN" ]; then
    echo "please enter QChem input name (without its extension .in)"
#      read -p "Type qchem input name : " QCIN
        exit
        fi

#Running QChem dynamics 
qchem -save $QCIN.in $QCIN.out scr_$QCIN

cp -r $QCSCRATCH/scr_$QCIN $HOME_DIR
rm -rf $QCSCRATCH/scr_$QCIN

#Read coords and Energy from QChem output to a file.xyz (not extended) to predict its forces.
python  $HOME_DIR/xyz_maker.py $1.out $1.xyz  

export QCIN 
export TESTF  #"test file name" for sGDML testing
python  $HOME_DIR/force_query.py #executing the force calculation using a trained model 

#Read the forces generated in Python (sGDML) using c (routine to be introduced in Qchem) 
./force_reader.exe

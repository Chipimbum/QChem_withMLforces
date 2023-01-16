#!/bin/bash

#this script calls sgdml from a bash interpreter to be able to do the 
#source environment file setup_anaconda_spike.sh  

source /usr/people/tamars/fanirm/anaconda_spike/setup_anaconda_spike.sh

python  /home/scr/fanirm/h20_dyns/force_query.py #executing the force calculation using a trained model 

echo -e "\nForce has been calculated using sGDML"

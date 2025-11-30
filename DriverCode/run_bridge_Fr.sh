#!/bin/bash

# Author: Radu Cimpeanu
# Date: 23/08/2025

# Additional resolution levels and prescribed contact angles can be added below

for LEVEL in 8; do
	for ANGLES in 30; do
	
		# create new output folder
		mkdir Results-Angle$ANGLES-Level$LEVEL-CastorCopperFr-ModelComparison-Fig12  
		
		# everything from -L$BASILISK/gl onward is just for visualisation (but needs a correct installation)
		qcc -O2 -w -fopenmp -Wall liquid_bridge_Fr_Fig12_Set2.c -lm -o Results-Angle$ANGLES-Level$LEVEL-CastorCopperFr-ModelComparison-Fig12/2D_Bridge_Fig12.x -L$BASILISK/gl -lglutils -lfb_tiny
		
		# enter folder
		cd Results-Angle$ANGLES-Level$LEVEL-CastorCopperFr-ModelComparison-Fig12 
		mkdir Slices
		mkdir Interfaces
		mkdir Animations

		# execute code given resources
		export OMP_NUM_THREADS=4
		# currently the parameters are:
		# 1: contact angle value
		# 2: maximum resolution level
		./2D_Bridge_Fig12.x $ANGLES $LEVEL
		
		# exit folder      
		cd ..  
	
	done
done         

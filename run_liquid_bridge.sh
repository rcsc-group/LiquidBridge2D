#!/bin/bash

# Author: Radu Cimpeanu
# Date: 3/11/2024

# Additional resolution levels and contact angles can be added below

for LEVEL in 8; do
	for ANGLES in 45; do
	
		# create new output folder
		mkdir Results-Angle$ANGLES-Level$LEVEL  
		
		# everything from -L$BASILISK/gl onward is just for visualisation (but needs a correct installation)
		qcc -O2 -w -fopenmp -Wall liquid_bridge.c -lm -o Results-Angle$ANGLES-Level$LEVEL/liquid_bridge.x -L$BASILISK/gl -lglutils -lfb_tiny
		
		# enter folder
		cd Results-Angle$ANGLES-Level$LEVEL
		mkdir Slices
		mkdir Interfaces
		mkdir Animations

		# execute code given resources
		export OMP_NUM_THREADS=1
		./liquid_bridge.x $ANGLES $LEVEL
		
		# exit folder      
		cd ..  
	
	done
done         

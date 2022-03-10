#!/bin/bash

########################################
############# CSCI 2951-O ##############
########################################
E_BADARGS=65
if [ $# -ne 1 ]
then
	echo "Usage: `basename $0` <input>"
	exit $E_BADARGS
fi
	
input=$1

# change this for your own installations!
export CP_SOLVER_EXEC=/Applications/CPLEX_Studio129/cpoptimizer/bin/x86-64_osx/cpoptimizer

# add the solver jar to the classpath and run
python3.9 src/main.py $input

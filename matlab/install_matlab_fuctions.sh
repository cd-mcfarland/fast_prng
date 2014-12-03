#!/bin/sh

# This script compiles and copies all PRNGs into your MATLAB path. It compiles 
# using mex, which must be installed and paired with a compatible C compiler. 
# You must also ensure that $USERPATH is a searched directory of your MATLAB
# installation. 

USERPATH="$HOME/Documents/MATLAB/"   

for file in "cdm_rand" "cdm_exprnd" "cdm_randn"
do
    echo "Compiling $file ...";
    mex $file.c;
    cp $file.mex* $USERPATH;
done


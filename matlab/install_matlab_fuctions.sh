#!/bin/sh

:<<'COMMENT'

This script compiles and copies all PRNGs into your MATLAB path. It compiles 
via mex, which must be installed and paired with a compatible C compiler. Use: 

 $ mex -setup

to see you current mex setup. Finally, before running this script ensure that 
$USERPATH is a searched directory of your MATLAB installation and $MEXCOMMAND 
is the correct path.

COMMENT

USERPATH="$HOME/Documents/MATLAB/"   
MEXCOMMAND="mex" 

set -e
for file in "fast_rand" "fast_exprnd" "fast_randn"
do
    echo "Compiling $file ...";
    $MEXCOMMAND -outdir $USERPATH $file.c;
done


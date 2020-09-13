#!/bin/bash
# -------------------------------------------------
# --- script for running test examples of ONSAS ---
# -------------------------------------------------

# enter examples folder
cd examples

# run test examples 
octave runTestProblems.m

# read boolean and delete file
num=( $(<auxVerifBoolean.dat) )
rm auxVerifBoolean.dat

# go back to root folder
cd ..

# sets output status and exits
exit $num


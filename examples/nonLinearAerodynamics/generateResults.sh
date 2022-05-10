# Code to generete results of ONSAS example: nonLinearAerodynamic 
#
# Generate Julia codes
# before adding BoundaryValueDiffEq, Plots, FileIO, DataFrames, CSV libraries
julia DiffEq.jl
# Add julia folder into the .basrc file using e.g: 'export PATH=$PATH:'/home/user/tools/julia/bin/'
# Generate ONSAS example results
octave --eval nonLinearAerodyamics

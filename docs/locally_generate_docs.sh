#Bash to generate the ONSAS.m documents authomatically

#Run the script to trsansform .m into .md
cd src
octave --eval bringONSASmFilesToONSASdocs
cd ..

# Add julia folder into the .basrc file using e.g: 'export PATH=$PATH:'/home/user/tools/julia/bin/'

# Run julia to create html with documenter
julia --project=. make.jl       

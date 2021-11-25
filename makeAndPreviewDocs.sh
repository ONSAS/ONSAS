#!/bin/bash

cd src
octave --eval "bringONSASmFilesToONSASdocs"
cd ..

cp "examples/staticVonMisesTruss/output/vonMisesTrussCheck.png" docs/src/
cp "examples/uniaxialExtension/output/verifUniaxial.png" docs/src/
cp "examples/uniformCurvatureCantilever/output/verifCantileverBeam.png" docs/src/

# make documention
julia docs/make.jl $1

# preview build
if [ $1 = pdf ]
then
  evince docs/build/ONSAS.m.pdf
elif [ $1 = html ]
then
  firefox docs/build/index.html
else
  echo "ERROR 'pdf' or 'html' argument must be provided to this script"
  exit 1
fi

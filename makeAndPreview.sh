#!/bin/bash

if [ -z "$ONSAS_PATH" ]; then
    echo "Need to set ONSAS_PATH"
    exit 1
fi

cd src
octave --eval "bringONSASmFilesToONSASdocs('$ONSAS_PATH')"
cd ..

cp "$ONSAS_PATH/examples/staticVonMisesTruss/output/vonMisesTrussCheck.png" docs/src/
cp "$ONSAS_PATH/examples/uniaxialExtension/output/verifUniaxial.png" docs/src/
cp "$ONSAS_PATH/examples/uniformCurvatureCantilever/output/verifCantileverBeam.png" docs/src/

# make documention
julia docs/make.jl $1

# preview build
if [ $1 = pdf ]
then
  evince docs/build/ONSAS.m.pdf
elif [ $1 = html ]
then
  firefox docs/build/index.html
fi

#!/bin/bash

GSLINC="-I/opt/local/include/gsl"
GSLLIB="-L/opt/local/lib/ -lgsl"

g++-mp-6 main.cpp $GSLINC $GSLLIB -o JBD

# Set parameters
wBD="1000.0"
OmegaM0="0.269"
OmeagaR0="8.4e-5"
zini="1e6"
npts="10000"
outname="output.txt"

# Run code
./JBD $outname $wBD $OmegaM0 $OmeagaR0 $zini $npts


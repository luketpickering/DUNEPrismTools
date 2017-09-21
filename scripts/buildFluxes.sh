#!/bin/sh
#./BuildFluxes.exe -i ../forMike/ -d 0.1,5 -vb 0_8:0.05
#./BuildFluxes.exe -i ../forMike/ -d 0.1,5 -b 200,0,10 -o D00NFluxes.root
./BuildFluxes.exe -i ../forMike/ -d 0.1,3.5 -BO 3000 -r 200,0,10 -I -C -o D00NFluxes_3.5.root
./BuildFluxes.exe -i ../forMike/ -d 0.1,3.5 -b 200,0,10 -o D00NFluxes_3.5_nobo.root -I -C

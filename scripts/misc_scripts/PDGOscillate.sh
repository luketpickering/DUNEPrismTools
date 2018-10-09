#!/bin/sh

#PDG2014
#./OscillateFlux.exe -i D00NFluxes.root,LBNF_numu_mrad_0 -o D00NFluxes_osc.root,LBNF_numu_mrad_0_osc -n 14,14 -p 0.297,0.0214,0.437,7.37E-5,2.5E-3,0 -d 5.8

#Osc CDR
#./OscillateFlux.exe -i D00NFluxes_3.5.root,LBNF_numu_mrad_0 -o D00NFluxes_CDROsc_3.5.root,LBNF_numu_mrad_0_osc -n 14,14 -p 0.304,0.0217,0.452,7.5E-5,2.457E-3,0 -d 5.8
#./OscillateFlux.exe -i D00NFluxes_3.5_BO4000.root,LBNF_numu_mrad_0 -o D00NFluxes_CDROsc_3.5_BO4000.root,LBNF_numu_mrad_0_osc -n 14,14 -p 0.304,0.0217,0.452,7.5E-5,2.457E-3,0 -d 5.8
#./OscillateFlux.exe -i D00NFluxes_3.5_nobo.root,LBNF_numu_mrad_0 -o D00NFluxes_CDROsc_3.5_nobo.root,LBNF_numu_mrad_0_osc -n 14,14 -p 0.304,0.0217,0.452,7.5E-5,2.457E-3,0 -d 5.8


./OscillateFlux.exe -i D00NFluxes_3.5_BO4000.root,LBNF_numu_mrad_0 -o talk_plots/D00NFluxes_NOvA_Osc_3.5_BO4000.root,LBNF_numu_mrad_0_osc -n 14,14 -p 0.297,0.0214,0.404,7.37E-5,2.67E-3,0 -d 5.8
./OscillateFlux.exe -i D00NFluxes_3.5_BO4000.root,LBNF_numu_mrad_0 -o talk_plots/D00NFluxes_T2K_Osc_3.5_BO4000.root,LBNF_numu_mrad_0_osc -n 14,14 -p 0.297,0.0214,0.534,7.37E-5,2.54E-3,0 -d 5.8

./OscillateFlux.exe -i D00NFluxes_3.5_nobo.root,LBNF_numu_mrad_0 -o talk_plots/D00NFluxes_NOvA_Osc_3.5_nobo.root,LBNF_numu_mrad_0_osc -n 14,14 -p 0.297,0.0214,0.404,7.37E-5,2.67E-3,0 -d 5.8
./OscillateFlux.exe -i D00NFluxes_3.5_nobo.root,LBNF_numu_mrad_0 -o talk_plots/D00NFluxes_T2K_Osc_3.5_nobo.root,LBNF_numu_mrad_0_osc -n 14,14 -p 0.297,0.0214,0.534,7.37E-5,2.54E-3,0 -d 5.8

#!/bin/sh

# ./FitFluxes.exe -i D00NFluxes_CDROsc.root,LBNF_numu_mrad_0_osc -o D00NPrism_CDRParams_FirstOsc.root -f D00NFluxes.root -c 20 -n 100000 -m 0 -E -l 1.5,3.8
# ./FitFluxes.exe -i D00NFluxes_CDROsc.root,LBNF_numu_mrad_0_osc -o D00NPrism_CDRParams_FirstSecondOsc.root -f D00NFluxes.root -c 20 -n 100000 -m 0 -E -l 0.5,3.8
# ./FitFluxes.exe -i D00NFluxes_CDROsc.root,LBNF_numu_mrad_0_osc -o D00NPrism_CDRParams_SecondOsc.root -f D00NFluxes.root -c 20 -n 100000 -m 0 -E -l 0.5,1.5
# ./FitFluxes.exe -i D00NFluxes_CDROsc.root,LBNF_numu_mrad_0_osc -o D00NPrism_CDRParams_All.root -f D00NFluxes.root -c 20 -n 100000 -m 0 -E

#./FitFluxes.exe -i D00NFluxes_CDROsc_3.5_BO4000.root,LBNF_numu_mrad_0_osc -o D00NPrism_CDRParams_3.5_BO4000_all.root -f D00NFluxes_3.5_BO4000.root -c 20 -n 100000
#./FitFluxes.exe -i D00NFluxes_CDROsc_3.5_BO4000.root,LBNF_numu_mrad_0_osc -o D00NPrism_CDRParams_3.5_BO4000.root -f D00NFluxes_3.5_BO4000.root -c 20 -n 100000 -m 0 -l 0.35,3.8

#./FitFluxes.exe -i D00NFluxes_CDROsc_3.5.root,LBNF_numu_mrad_0_osc -o D00NPrism_CDRParams_3.5.root -f D00NFluxes_3.5.root -c 20 -n 100000 -m 0 -l 0.5,3.8
#./FitFluxes.exe -i D00NFluxes_CDROsc_3.5.root,LBNF_numu_mrad_0_osc -o D00NPrism_CDRParams_3.5_lowfit.root -f D00NFluxes_3.5.root -c 20 -n 100000 -m 0 -l 0,3.8
#./FitFluxes.exe -i D00NFluxes_CDROsc_3.5_nobo.root,LBNF_numu_mrad_0_osc -o D00NPrism_CDRParams_3.5_nobo.root -f D00NFluxes_3.5_nobo.root -c 20 -n 100000 -m 0 -l 0.5,3.8
#./FitFluxes.exe -i D00NFluxes_CDROsc_3.5_nobo.root,LBNF_numu_mrad_0_osc -o D00NPrism_CDRParams_3.5_nobo_lowfit.root -f D00NFluxes_3.5_nobo.root -c 20 -n 100000 -m 0 -l 0.35,3.8
#./FitFluxes.exe -i D00NFluxes_CDROsc_3.5_nobo.root,LBNF_numu_mrad_0_osc -o D00NPrism_CDRParams_3.5_nobo_highfit.root -f D00NFluxes_3.5_nobo.root -c 20 -n 100000 -m 0 -l 0.5,4.2

./FitFluxes.exe -i talk_plots/D00NFluxes_NOvA_Osc_3.5_BO4000.root,LBNF_numu_mrad_0_osc -o talk_plots/D00NPrism_NOvA_Osc_3.5_BO4000.root -f D00NFluxes_3.5_BO4000.root -c 30 -n 100000 -m 0 -l 0.35,3.8
./FitFluxes.exe -i talk_plots/D00NFluxes_T2K_Osc_3.5_BO4000.root,LBNF_numu_mrad_0_osc -o talk_plots/D00NPrism_T2K_Osc_3.5_BO4000.root -f D00NFluxes_3.5_BO4000.root -c 30 -n 100000 -m 0 -l 0.35,3.8

./FitFluxes.exe -i talk_plots/D00NFluxes_NOvA_Osc_3.5_nobo.root,LBNF_numu_mrad_0_osc -o talk_plots/D00NPrism_NOvA_Osc_3.5_nobo.root -f D00NFluxes_3.5_nobo.root -c 30 -n 100000 -m 0 -l 0.35,3.8
./FitFluxes.exe -i talk_plots/D00NFluxes_T2K_Osc_3.5_nobo.root,LBNF_numu_mrad_0_osc -o talk_plots/D00NPrism_T2K_Osc_3.5_nobo.root -f D00NFluxes_3.5_nobo.root -c 30 -n 100000 -m 0 -l 0.35,3.8

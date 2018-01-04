# ${DUNEPRISMTOOLSROOT}/scripts/FarmThrowGENIEEventsJobs.sh -p nominal_2.5E8/DUNEPrismFluxes/FHC -i /pnfs/dune/persistent/users/picker24/nominal_2.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/

# ${DUNEPRISMTOOLSROOT}/scripts/FarmThrowGENIEEventsJobs.sh -p nominal_2.5E8/DUNEPrismFluxes/RHC -i /pnfs/dune/persistent/users/picker24/nominal_2.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/flux/

${DUNEPRISMTOOLSROOT}/scripts/FarmThrowGENIEEventsJobs.sh -p nominal_7.5E8/DUNEPrismFluxes/FHC_FD -f LArDUNEFV_25.8mx23.8mx56.6m_1287km -g ${DUNEPRISMTOOLSROOT}/configs/OnAxisLArBox_FD_25.8mx23.8mx56.6m.xml -e 1.28E19 -N 100 -i /pnfs/dune/persistent/users/picker24/nominal_7.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/

${DUNEPRISMTOOLSROOT}/scripts/FarmThrowGENIEEventsJobs.sh -p nominal_7.5E8/DUNEPrismFluxes/RHC_FD -f LArDUNEFV_25.8mx23.8mx56.6m_1287km -g ${DUNEPRISMTOOLSROOT}/configs/OnAxisLArBox_FD_25.8mx23.8mx56.6m.xml -e 1.28E19 -N 100 -i /pnfs/dune/persistent/users/picker24/nominal_7.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/flux/

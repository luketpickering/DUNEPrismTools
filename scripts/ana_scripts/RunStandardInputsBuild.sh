#!/bin/bash

####### ND 4m wide, No FV, No acceptance
dp_RunSelection \
 -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30/*.root" \
 -v 20 \
 -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30.Selected.NoAcceptance.root \

if [[ -L Selected.4mwide.0cmFV.NoAcceptance.root ]]; then
  rm Selected.4mwide.0cmFV.NoAcceptance.root
fi
ln -s /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30.Selected.NoAcceptance.root \
  Selected.4mwide.0cmFV.NoAcceptance.root

dp_ProduceEfficiencyCorrector \
  -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30/*.root" \
  -v 20 \
  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30.EffCorrector.NoAcceptance.root \
  -b 0_10:0.1

if [[ -L EffCorrector.4mwide.0cmFV.NoAcceptance.root ]]; then
  rm EffCorrector.4mwide.0cmFV.NoAcceptance.root
fi
ln -s /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30.EffCorrector.NoAcceptance.root \
  EffCorrector.4mwide.0cmFV.NoAcceptance.root

dp_EfficiencyWeightFriendTreeBuilder \
  -i Selected.4mwide.0cmFV.NoAcceptance.root \
  -E EffCorrector.4mwide.0cmFV.NoAcceptance.root \
  -o EffFriend.4mwide.0cmFV.NoAcceptance.detpos_corr.root

dp_EfficiencyWeightFriendTreeBuilder \
  -i Selected.4mwide.0cmFV.NoAcceptance.root \
  -E EffCorrector.4mwide.0cmFV.NoAcceptance.root \
  -M 5 \
  -o EffFriend.4mwide.0cmFV.NoAcceptance.abspos_corr.root


dp_XSecVarFriendTreeBuilder \
   -i Selected.4mwide.0cmFV.NoAcceptance.root \
   -N 500 \
   -o XSecThrowFriend.4mwide.0cmFV.NoAcceptance.root \
   -S 2207


####### ND 4m wide, No FV, Continuous
# dp_RunSelection \
#  -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30/*.root" \
#  -v 20 \
#  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30.Selected.root \
#  -A 2
#
# if [[ -L Selected.4mwide.0cmFV.root ]]; then
#   rm Selected.4mwide.0cmFV.root
# fi
# ln -s /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30.Selected.root \
#   Selected.4mwide.0cmFV.root
#
# dp_ProduceEfficiencyCorrector \
#   -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30/*.root" \
#   -v 20 \
#   -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30.EffCorrector.root \
#   -b 0_10:0.1 -A 2
#
# if [[ -L EffCorrector.4mwide.0cmFV.root ]]; then
#   rm EffCorrector.4mwide.0cmFV.root
# fi
# ln -s /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30.EffCorrector.root \
#   EffCorrector.4mwide.0cmFV.root
#
# dp_EfficiencyWeightFriendTreeBuilder \
#   -i Selected.4mwide.0cmFV.root \
#   -E EffCorrector.4mwide.0cmFV.root \
#   -o EffFriend.4mwide.0cmFV.root

# dp_XSecVarFriendTreeBuilder \
#    -i Selected.4mwide.0cmFV.root \
#    -N 500 \
#    -o XSecThrowFriend.4mwide.0cmFV.root \
#    -S 2207

####### ND 7m wide, 100 cm FV, Continuous
# dp_RunSelection \
#  -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07/*.root" \
#  -v 20 \
#  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07.Selected.root \
#  -A 2 -FV 100,0,0
#
# if [[ -L Selected.7mwide.100cmFV.root ]]; then
#   rm Selected.7mwide.100cmFV.root
# fi
# ln -s /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07.Selected.root \
#   Selected.7mwide.100cmFV.root

# dp_ProduceEfficiencyCorrector \
#   -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07/*.root" \
#   -v 20 \
#   -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07.EffCorrector.root \
#   -b 0_10:0.1 -A 2
#
# if [[ -L EffCorrector.7mwide.100cmFV.root ]]; then
#   rm EffCorrector.7mwide.100cmFV.root
# fi
# ln -s /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07.EffCorrector.root \
#   EffCorrector.7mwide.100cmFV.root
#
# dp_EfficiencyWeightFriendTreeBuilder \
#   -i Selected.7mwide.100cmFV.root \
#   -E EffCorrector.7mwide.100cmFV.root \
#   -o EffFriend.7mwide.100cmFV.root

# dp_XSecVarFriendTreeBuilder \
#    -i Selected.7mwide.100cmFV.root \
#    -N 500 \
#    -o XSecThrowFriend.7mwide.100cmFV.root \
#    -S 2207

###### ND 7m wide, 100 cm FV, NoAcceptance
dp_RunSelection \
 -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07/*.root" \
 -v 20 \
 -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07.Selected.NoAcceptance.root \
 -FV 100,0,0

if [[ -L Selected.7mwide.100cmFV.NoAcceptance.root ]]; then
  rm Selected.7mwide.100cmFV.NoAcceptance.root
fi
ln -s /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07.Selected.NoAcceptance.root \
  Selected.7mwide.100cmFV.NoAcceptance.root

dp_ProduceEfficiencyCorrector \
  -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07/*.root" \
  -v 20 \
  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07.EffCorrector.NoAcceptance.root \
  -b 0_10:0.1

if [[ -L EffCorrector.7mwide.100cmFV.NoAcceptance.root ]]; then
  rm EffCorrector.7mwide.100cmFV.NoAcceptance.root
fi
ln -s /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07.EffCorrector.NoAcceptance.root \
  EffCorrector.7mwide.100cmFV.NoAcceptance.root

dp_EfficiencyWeightFriendTreeBuilder \
  -i Selected.7mwide.100cmFV.NoAcceptance.root \
  -E EffCorrector.7mwide.100cmFV.NoAcceptance.root \
  -o EffFriend.7mwide.100cmFV.NoAcceptance.detpos_corr.root

dp_EfficiencyWeightFriendTreeBuilder \
  -i Selected.7mwide.100cmFV.NoAcceptance.root \
  -E EffCorrector.7mwide.100cmFV.NoAcceptance.root \
  -M 5 \
  -o EffFriend.7mwide.100cmFV.NoAcceptance.abspos_corr.root

dp_XSecVarFriendTreeBuilder \
   -i Selected.7mwide.100cmFV.NoAcceptance.root \
   -N 500 \
   -o XSecThrowFriend.7mwide.100cmFV.NoAcceptance.root \
   -S 2207

#### FD

dp_RunSelection \
 -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_FD/Processed.2018-04-03/*.root" \
 -v 20 \
 -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_FD/Processed.2018-04-03.Selected.root \
 -FV 100,100,100

 if [[ -L Selected.FD.root ]]; then
   rm Selected.FD.root
 fi
 ln -s /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_FD/Processed.2018-04-03.Selected.root \
   Selected.FD.root

dp_ProduceEfficiencyCorrector \
  -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_FD/Processed.2018-04-03/*.root" \
  -v 20 \
  -b 0_10:0.1 \
  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_FD/Processed.2018-04-03.EffCorrector.root

if [[ -L EffCorrector.FD.root ]]; then
  rm EffCorrector.FD.root
fi
ln -s /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_FD/Processed.2018-04-03.EffCorrector.root \
  EffCorrector.FD.root

dp_EfficiencyWeightFriendTreeBuilder \
 -i Selected.FD.root \
 -E EffCorrector.FD.root \
 -o EffFriend.FD.root

dp_XSecVarFriendTreeBuilder \
  -i Selected.FD.root \
  -N 500 \
  -o XSecThrowFriend.FD.root \
  -S 2207

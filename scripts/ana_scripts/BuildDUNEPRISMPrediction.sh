# #!/bin/bash
# dp_PRISMAnalysis  \
#   -NI Selected.4mwide.0cmFV.NoAcceptance.root -NE EffFriend.4mwide.0cmFV.NoAcceptance.root \
#   -FI Selected.FD.root -FE EffFriend.FD.root \
#   -F FluxFit.root \
#   -o PRISMPrediction.4mwide.0cmFV.NoAcceptance.root \
#   -b 0_10:0.25 \
#   # -FX XSecThrowFriend.FD.root \
#   # -NX XSecThrowFriend.4mwide.0cmFV.NoAcceptance.root \
#

for i in dm230.0022.sin2theta230.5 dm230.0025.sin2theta230.5 \
  dm230.0028.sin2theta230.5 dm230.0022.sin2theta230.65 \
  dm230.0025.sin2theta230.65 dm230.0028.sin2theta230.65; do
  dp_PRISMAnalysis  \
    -NI Selected.4mwide.0cmFV.NoAcceptance.root -NE EffFriend.4mwide.0cmFV.NoAcceptance.abspos_corr.root \
    -FI Selected.FD.root -FE EffFriend.FD.root \
    -F nominal.${i}.root    \
    -o PRISMPrediction.4mwide.0cmFV.NoAcceptance.abspos_corr.ERec.${i}.root \
    -b 0_10:0.25 -PV 3 \
    # -FX XSecThrowFriend.FD.root \
    # -NX XSecThrowFriend.4mwide.0cmFV.root \
  dp_PRISMAnalysis  \
    -NI Selected.4mwide.0cmFV.NoAcceptance.root -NE EffFriend.4mwide.0cmFV.NoAcceptance.detpos_corr.root \
    -FI Selected.FD.root -FE EffFriend.FD.root \
    -F nominal.${i}.root    \
    -o PRISMPrediction.4mwide.0cmFV.NoAcceptance.detpos_corr.ERec.${i}.root \
    -b 0_10:0.25 -PV 3 \
    # -FX XSecThrowFriend.FD.root \
    # -NX XSecThrowFriend.4mwide.0cmFV.root \
  dp_PRISMAnalysis  \
    -NI Selected.4mwide.0cmFV.NoAcceptance.root -NE EffFriend.4mwide.0cmFV.NoAcceptance.abspos_corr.root \
    -FI Selected.FD.root -FE EffFriend.FD.root \
    -F nominal.${i}.root    \
    -o PRISMPrediction.4mwide.0cmFV.NoAcceptance.abspos_corr.EAvail.${i}.root \
    -b 0_10:0.25 -PV 2 \
    # -FX XSecThrowFriend.FD.root \
    # -NX XSecThrowFriend.4mwide.0cmFV.root \
  dp_PRISMAnalysis  \
    -NI Selected.4mwide.0cmFV.NoAcceptance.root -NE EffFriend.4mwide.0cmFV.NoAcceptance.detpos_corr.root \
    -FI Selected.FD.root -FE EffFriend.FD.root \
    -F nominal.${i}.root    \
    -o PRISMPrediction.4mwide.0cmFV.NoAcceptance.detpos_corr.EAvail.${i}.root \
    -b 0_10:0.25 -PV 2 \
    # -FX XSecThrowFriend.FD.root \
    # -NX XSecThrowFriend.7mwide.100cmFV.root \
done

for i in dm230.0022.sin2theta230.5 dm230.0025.sin2theta230.5 \
  dm230.0028.sin2theta230.5 dm230.0022.sin2theta230.65 \
  dm230.0025.sin2theta230.65 dm230.0028.sin2theta230.65; do
  dp_PRISMAnalysis  \
    -NI Selected.7mwide.100cmFV.NoAcceptance.root -NE EffFriend.7mwide.100cmFV.NoAcceptance.detpos_corr.root \
    -FI Selected.FD.root -FE EffFriend.FD.root \
    -F nominal.${i}.root    \
    -o PRISMPrediction.7mwide.100cmFV.NoAcceptance.detpos_corr.ERec.${i}.root \
    -b 0_10:0.25 -PV 3 \
    # -FX XSecThrowFriend.FD.root \
    # -NX XSecThrowFriend.7mwide.100cmFV.root \
  dp_PRISMAnalysis  \
    -NI Selected.7mwide.100cmFV.NoAcceptance.root -NE EffFriend.7mwide.100cmFV.NoAcceptance.abspos_corr.root \
    -FI Selected.FD.root -FE EffFriend.FD.root \
    -F nominal.${i}.root    \
    -o PRISMPrediction.7mwide.100cmFV.NoAcceptance.abspos_corr.ERec.${i}.root \
    -b 0_10:0.25 -PV 3 \
    # -FX XSecThrowFriend.FD.root \
    # -NX XSecThrowFriend.7mwide.100cmFV.root \
  dp_PRISMAnalysis  \
    -NI Selected.7mwide.100cmFV.NoAcceptance.root -NE EffFriend.7mwide.100cmFV.NoAcceptance.detpos_corr.root \
    -FI Selected.FD.root -FE EffFriend.FD.root \
    -F nominal.${i}.root    \
    -o PRISMPrediction.7mwide.100cmFV.NoAcceptance.detpos_corr.EAvail.${i}.root \
    -b 0_10:0.25 -PV 2 \
    # -FX XSecThrowFriend.FD.root \
    # -NX XSecThrowFriend.7mwide.100cmFV.root \
  dp_PRISMAnalysis  \
    -NI Selected.7mwide.100cmFV.NoAcceptance.root -NE EffFriend.7mwide.100cmFV.NoAcceptance.abspos_corr.root \
    -FI Selected.FD.root -FE EffFriend.FD.root \
    -F nominal.${i}.root    \
    -o PRISMPrediction.7mwide.100cmFV.NoAcceptance.abspos_corr.EAvail.${i}.root \
    -b 0_10:0.25 -PV 2  \
    # -FX XSecThrowFriend.FD.root \
    # -NX XSecThrowFriend.7mwide.100cmFV.root \
done

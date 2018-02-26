#include <math.h>
class DepoParticle {
 public:
  int PDG;
  int trackID;
  double energy;

  double xBound[2];
  double yBound;
  double zBound;

 private:
};

class DepoLepton : public DepoParticle {
 public:
  int flagExitBack = 0;
  int flagExitFront = 0;
  int flagExitY = 0;
  int flagExitXHigh = 0;
  int flagExitXLow = 0;
  int flagLepContained = 1;

  // Final position and momentum of the muon track
  // Used for determining where it exits
  double xf = 0.;
  double yf = 0.;
  double zf = 0.;

  double pxf = 0.;
  double pyf = 0.;
  double pzf = 0.;

  double eDepPrimary;    // From lepton itself
  double eDepSecondary;  // From any secondary/tertiary...

  DepoLepton(int inPDG, int inTrackID, double inE, double *inXBound,
             double inYBound, double inZBound) {
    PDG = inPDG;
    trackID = inTrackID;
    energy = inE;
    xBound[0] = inXBound[0];
    xBound[1] = inXBound[1];
    yBound = inYBound;
    zBound = inZBound;
  };

  void CheckContained() {
    if (xf <= xBound[0])
      flagExitXLow = 1;
    else if (xf >= xBound[1])
      flagExitXHigh = 1;
    else if (fabs(yf) >= yBound)
      flagExitY = 1;
    else if (zf >= zBound)
      flagExitBack = 1;
    else if (zf <= -1. * zBound)
      flagExitFront = 1;

    flagLepContained = !(flagExitXLow || flagExitXHigh || flagExitY ||
                         flagExitBack || flagExitFront);
  };

 private:
};

class DepoHadron : public DepoParticle {
 public:
  // In Bounds
  double eDepPrimaryIn;    // From hadron itself
  double eDepSecondaryIn;  // From any secondary/tertiary...

  // Out of Bounds
  double eDepPrimaryOut;    // From hadron itself
  double eDepSecondaryOut;  // From any secondary/tertiary...

  //    int flagNoEHadOut;

  DepoHadron(int inPDG, int inTrackID, double inE, double *inXBound,
             double inYBound, double inZBound) {
    PDG = inPDG;
    trackID = inTrackID;
    energy = inE;
    xBound[0] = inXBound[0];
    xBound[1] = inXBound[1];
    yBound = inYBound;
    zBound = inZBound;
  };
};

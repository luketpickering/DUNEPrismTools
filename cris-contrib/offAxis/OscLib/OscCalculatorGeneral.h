#ifndef OSCCALCULATORGENERAL_H
#define OSCCALCULATORGENERAL_H

#include "OscLib/func/IOscCalculator.h"

namespace osc
{

  class OscCalculatorGeneral: public IOscCalculatorAdjustable
  {
  public:
    OscCalculatorGeneral();
    virtual ~OscCalculatorGeneral();

    // Baseline in km
    virtual void SetL(double L) {fL = L;}
    // Density in g/cm^3
    virtual void SetRho(double rho) {fRho = rho;}
    // in eV^2
    virtual void SetDmsq21(double dmsq21) {fDmsq21 = dmsq21;}
    // This is a signed quantity, use a negative value for inverted hierarchy
    virtual void SetDmsq32(double dmsq32) {fDmsq32 = dmsq32;}
    // In radians
    virtual void SetTh12(double th12);
    virtual void SetTh13(double th13);
    virtual void SetTh23(double th23);
    virtual void SetdCP(double dCP);

    void SetNSIEpsilonMuTau(double emutau) {fEMuTau = emutau;}
    double GetNSIEpsilonMuTau() const {return fEMuTau;}

    virtual double P(int from, int to, double E);

    virtual TMD5* GetParamsHash() const override
    {
      // Default isn't good enough if we need to consider NSI
      if(fEMuTau) return 0;
      return IOscCalculatorAdjustable::GetParamsHashDefault("General");
    }

    struct Priv;
  protected:
    Priv* const d;

    double fEMuTau;

  private:
    OscCalculatorGeneral(const OscCalculatorGeneral&);
    OscCalculatorGeneral& operator=(const OscCalculatorGeneral&);
  };

} // end namespace

#endif

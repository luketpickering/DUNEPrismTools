#ifndef OSC_OSCCALCULATORPMNSOPT_H
#define OSC_OSCCALCULATORPMNSOPT_H

//                                                                      //
// \file OscCalculatorPMNSOpt.h                                         //
//                                                                      //
// Adapt the PMNSOpt calculator to standard interface                   //
// <bckhouse@caltech.edu>                                               //
//                                                                      //

#include "IOscCalculator.h"
#include "PMNSOpt.h"

namespace osc
{
  class OscCalculatorPMNSOpt: public IOscCalculatorAdjustable
  {
  public:
    OscCalculatorPMNSOpt();
    virtual ~OscCalculatorPMNSOpt();

    virtual double P(int flavBefore, int flavAfter, double E);

    virtual void SetL     (double L     ){fPropDirty = true; fL      = L;}
    virtual void SetRho   (double rho   ){fPropDirty = true; fRho    = rho;}
    virtual void SetDmsq21(double dmsq21){fDmDirty   = true; fDmsq21 = dmsq21;}
    virtual void SetDmsq32(double dmsq32){fDmDirty   = true; fDmsq32 = dmsq32;}
    virtual void SetTh12  (double th12  ){fMixDirty  = true; fTh12   = th12;}
    virtual void SetTh13  (double th13  ){fMixDirty  = true; fTh13   = th13;}
    virtual void SetTh23  (double th23  ){fMixDirty  = true; fTh23   = th23;}
    virtual void SetdCP   (double dCP   ){fMixDirty  = true; fdCP    = dCP;}

    virtual TMD5* GetParamsHash() const override
    {
      return IOscCalculatorAdjustable::GetParamsHashDefault("PMNSOpt");
    }
  protected:
    PMNSOpt fPMNSOpt;

    bool fMixDirty;
    bool fDmDirty;
    bool fPropDirty;
    double fPrevE;
    int fPrevAnti;
  };

} // namespace

#endif

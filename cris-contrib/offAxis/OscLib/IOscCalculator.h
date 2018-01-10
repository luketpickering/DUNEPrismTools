#ifndef IOSCCALCULATOR_H
#define IOSCCALCULATOR_H

#include "TMD5.h"

//                                                                      //
// \file    IOscCalculator.h                                            //
//                                                                      //
// \brief   General interface to oscillation calculators                //
// \author  Christopher Backhouse - bckhouse@caltech.edu                //
//                                                                      //

namespace osc
{
  class IOscCalculator
  {
  public:
    virtual ~IOscCalculator() {}

    virtual double P(int flavBefore, int flavAfter, double E) = 0;

    virtual TMD5* GetParamsHash() const {return 0;}
  };

  class NoOscillations: public IOscCalculator
  {
  public:
    virtual double P(int from, int to, double /*E*/)
    {
      if(from == to) return 1;
      return 0;
    }

    //virtual TMD5* GetParamsHash() const override
    //{
    //  TMD5* ret = new TMD5;
    //  const char* txt = "NoOscillations";
    //  ret->Update((unsigned char*)txt, strlen(txt));
    //  ret->Final();
    //  return ret;
    //}
  };

  class IOscCalculatorAdjustable : public IOscCalculator
  {
  public:
    // These setters are left unimplemented here, since calculators may want
    // to compute additional values when these are set.
    virtual void SetL     (double L     ) = 0;
    virtual void SetRho   (double rho   ) = 0;
    virtual void SetDmsq21(double dmsq21) = 0;
    virtual void SetDmsq32(double dmsq32) = 0;
    virtual void SetTh12  (double th12  ) = 0;
    virtual void SetTh13  (double th13  ) = 0;
    virtual void SetTh23  (double th23  ) = 0;
    virtual void SetdCP   (double dCP   ) = 0;

    virtual double GetL     () const { return fL      ; }
    virtual double GetRho   () const { return fRho    ; }
    virtual double GetDmsq21() const { return fDmsq21 ; }
    virtual double GetDmsq32() const { return fDmsq32 ; }
    virtual double GetTh12  () const { return fTh12   ; }
    virtual double GetTh13  () const { return fTh13   ; }
    virtual double GetTh23  () const { return fTh23   ; }
    virtual double GetdCP   () const { return fdCP    ; }

  protected:
    TMD5* GetParamsHashDefault(const std::string& txt) const;

    // Set by the user. Generally useful to derived classes
    double fRho; // density (g/cm^3)
    double fL; // baseline (km)
    double fDmsq21;
    double fDmsq32;
    double fTh12;
    double fTh13;
    double fTh23;
    double fdCP;
  };

} // namespace

#endif

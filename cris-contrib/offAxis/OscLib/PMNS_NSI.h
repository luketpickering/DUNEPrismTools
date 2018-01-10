











#ifndef PMNS_NSI_H
#define PMNS_NSI_H
#include <list>
#include <complex>

#include "OscLib/func/PMNSOpt.h"

namespace osc {
  class PMNS_NSI : public PMNSOpt {
  public:
    PMNS_NSI();
    virtual ~PMNS_NSI();
    
    void SetNSI(double eps_ee,      double eps_emu,      double eps_etau,
                double eps_mumu,    double eps_mutau,    double eps_tautau,
                double delta_emu=0, double delta_etau=0, double delta_mutau=0);

  protected:
    virtual void SolveHam(double E, double Ne, int anti);
    
    double  fEps_ee;        
    double  fEps_mumu;      
    double  fEps_tautau;    
    complex fEps_emu;       
    complex fEps_etau;      
    complex fEps_mutau;     
    bool    fResetNSI;      
  };
}
#endif


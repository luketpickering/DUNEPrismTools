






























#ifndef PMNSOPT_H
#define PMNSOPT_H
#include <list>
#include <complex>

// Some useful complex numbers
static std::complex<double> zero(0.0,0.0);
static std::complex<double> one (1.0,0.0);

// Unit conversion constants
static const double kKm2eV  = 5.06773103202e+09; 
static const double kK2     = 4.62711492217e-09; 
static const double kGeV2eV = 1.0e+09;           

//G_F in units of GeV^-2
static const double kGf     = 1.166371e-5;

namespace osc {
  class PMNSOpt {
  public:
    PMNSOpt();
    virtual ~PMNSOpt();

    virtual void SetMix(double th12, double th23, double th13, double deltacp);

    virtual void SetDeltaMsqrs(double dm21, double dm32);

    virtual void PropMatter(double L, double E, double Ne, int anti=1);
    virtual void PropMatter(const std::list<double>& L,
                    double                   E,
                    const std::list<double>& Ne,
                    int anti);

    virtual void PropVacuum(double L, double E, int anti=1);

    virtual double P(int flv) const;

    virtual void ResetToFlavour(int flv=1);

  protected:
    // A shorthand...
    typedef std::complex<double> complex;

    virtual void BuildHlv();

    virtual void SolveHam(double E, double Ne, int anti);

    virtual void SetVacuumEigensystem(double E, int anti);

    double  fDm21;          
    double  fDm31;          
    double  fTheta12;       
    double  fTheta23;       
    double  fTheta13;       
    double  fDeltaCP;       
    complex fHlv[3][3];     
    complex fEvec[3][3];    
    double  fEval[3];       
    complex fNuState[3];    
    double  fCachedNe;      
    double  fCachedE;       
    int     fCachedAnti;    
    bool    fBuiltHlv;      
  };
}
#endif

































#ifndef PMNS_H
#define PMNS_H
#include <list>
#include <complex>

namespace osc {
  class PMNS {
  public:
    PMNS();
    
    PMNS(double th12,
         double th23,
         double th13,
         double deltacp,
         double dms12,
         double dms23);
    
    void PrintMix() const;
    
    void PrintDeltaMsqrs() const;
    
    double P(int i, int j) const;
    
    void SetMix(double th12, double th23, double th13, double deltacp);
    
    void SetMixBWCP(double th1, double th2, double th3, double deltacp);
    
    void SetDeltaMsqrs(double dm21, double dm32);
    
    void PropVacuum(double L, double E, int anti);
    
    void PropMatter(double L, double E, double Ne, int anti);
    void PropMatter(const std::list<double>& L,
                    double                   E,
                    const std::list<double>& Ne,
                    int anti);
    
    void Reset();
    
  private:
    // A shorthand...
    typedef std::complex<double> complex;
    
  private:
    void Multi(complex A[][3], const complex B[][3], const complex C[][3]);
    
    void EvalEqn2(complex A[][3],
                  const complex U[][3],
                  const complex Udagg[][3],
                  const double  dmsqr[][3],
                  double L,
                  double E); 
    
    void EvalEqn5(complex       twoEH[][3],
                  const complex U[][3],
                  const complex Udagg[][3],
                  const double  dmsqr[][3],
                  double        E,
                  double        Gf,
                  double        Ne);
    
    void EvalEqn10(complex       A[][3],
                   const complex U[][3],
                   const complex X[][3],
                   const complex Udagg[][3]);
    
    void EvalEqn11(complex X[][3],
                   double L, double E, 
                   const complex twoEH[][3],
                   const double  Msqr[],
                   const double  dMsqr[][3]);
    
    void EvalEqn21(double Msqr[],
                   double alpha,
                   double beta,
                   double gamma);
    
    void EvalEqn22(double& alpha,
                   double& beta,
                   double& gamma,
                   double  E,
                   double  Gf,
                   double  Ne,
                   const double dmsqr[][3],
                   const complex U[][3]);
    
    void SortEigenvalues(double       dMsqr[][3],
                         const double dmsqr[][3],
                         const double MsqrVac[],
                         double       Msqr[]);
    
    void DumpMatrix(const complex M[][3]) const;
    
  private:
    double  fdmsqr[3][3]; 
    complex fU[3][3];     
    complex fUdagg[3][3];
    complex fUstar[3][3];
    complex fUtran[3][3];
    complex fA[3][3];     
  };
}
#endif


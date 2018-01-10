





#ifndef EARTHMODEL_H
#define EARTHMODEL_H
#include <list>
#include <vector>

namespace osc {
  class EarthModel 
  {
  public:
    EarthModel(const char* which, double tol);

    double Ne(double r);

    double Density(double r);

    double ZoverA(double r);
    
    void GetLayers(std::vector<double>& rlo, 
                   std::vector<double>& rhi,
                   std::vector<double>& ne);

    void LineProfile(double prodL, double cosQ, double rdet,
                     std::list<double>& Ls,
                     std::list<double>& Ns);
  private:
    void   InitPREM();
    void   InitStacey();
    double DensityPREM(double r);
    double DensityStacey(double r);
    double AveNe(double r1, double r2, int nsample);
    void   MakeLayers(double tol);
    int    IntersectLineAndCircle(double  x1, double  y1,
                                  double  x2, double  y2,
                                  double  r,
                                  double* xa, double* ya,
                                  double* xb, double* yb);
  private:
    int                 fModelID;    
    double              fRouterCore; 
    double              fRearth;     
    std::vector<double> fRregion;    

    std::vector<double> fRlo;    
    std::vector<double> fRhi;    
    std::vector<double> fNe;     
  };
}
#endif


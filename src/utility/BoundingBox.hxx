#ifndef BOUNDINGBOX_HXX_SEEN
#define BOUNDINGBOX_HXX_SEEN

#include "TVector3.h"

#include <array>
#include <string>

struct BoundingBox {
  std::array<double, 3> Max;
  std::array<double, 3> Min;

  //-x,+x,-y,+y,-z,+z
  std::array<std::array<double, 3>, 6> PlaneCenters;
  std::array<std::array<double, 3>, 6> PlaneNormals;
  std::array<std::array<double, 4>, 3> PlaneBBs;
  BoundingBox();
  BoundingBox(std::array<double,3> const &Max, std::array<double,3> const &Min);

  bool Contains(std::array<double,3> const &pos) const;

  std::string Print() const;

 private:
  void Sort();
  void CalcPlanes();
};

double CalculateToWall(BoundingBox const &BB, TVector3 const &TrStart,
                       TVector3 const &TrDir);

void TestBBIntersect();

#endif

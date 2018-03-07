#include "TVector3.h"

#include <array>
#include <iostream>

struct BoundingBox {
  std::array<Float_t, 3> Max;
  std::array<Float_t, 3> Min;

  //-x,+x,-y,+y,-z,+z
  std::array<std::array<Float_t, 3>, 6> PlaneCenters;
  std::array<std::array<Float_t, 3>, 6> PlaneNormals;
  std::array<std::array<Float_t, 4>, 3> PlaneBBs;
  BoundingBox() : Max({0, 0, 0}), Min({0, 0, 0}) {}
  BoundingBox(TVector3 const &Max, TVector3 const &Min) {
    this->Max[0] = Max[0];
    this->Max[1] = Max[1];
    this->Max[2] = Max[2];

    this->Min[0] = Min[0];
    this->Min[1] = Min[1];
    this->Min[2] = Min[2];

    Sort();
    CalcPlanes();
  }

 private:
  void Sort() {
    std::array<Float_t, 3> nMax;
    std::array<Float_t, 3> nMin;

    nMax[0] = std::max(Max[0], Min[0]);
    nMax[1] = std::max(Max[1], Min[1]);
    nMax[2] = std::max(Max[2], Min[2]);

    nMin[0] = std::min(Max[0], Min[0]);
    nMin[1] = std::min(Max[1], Min[1]);
    nMin[2] = std::min(Max[2], Min[2]);

    // Actually don't gain anything from this... but meh, #YOLO.
    Max = std::move(nMax);
    Min = std::move(nMin);
  }
  void CalcPlanes() {
    PlaneNormals[0] = std::array<Float_t, 3>{-1, 0, 0};
    PlaneNormals[1] = std::array<Float_t, 3>{1, 0, 0};

    PlaneNormals[2] = std::array<Float_t, 3>{0, -1, 0};
    PlaneNormals[3] = std::array<Float_t, 3>{0, 1, 0};

    PlaneNormals[4] = std::array<Float_t, 3>{0, 0, -1};
    PlaneNormals[5] = std::array<Float_t, 3>{0, 0, 1};

    std::array<Float_t, 3> Halfs{(Max[0] - Min[0]) / 2.0f,
                                 (Max[1] - Min[1]) / 2.0f,
                                 (Max[2] - Min[2]) / 2.0f};

    PlaneCenters[0] = std::array<Float_t, 3>{Min[0], Halfs[1], Halfs[2]};
    PlaneCenters[1] = std::array<Float_t, 3>{Max[0], Halfs[1], Halfs[2]};

    PlaneCenters[2] = std::array<Float_t, 3>{Halfs[0], Min[1], Halfs[2]};
    PlaneCenters[3] = std::array<Float_t, 3>{Halfs[0], Max[1], Halfs[2]};

    PlaneCenters[4] = std::array<Float_t, 3>{Halfs[0], Halfs[1], Min[2]};
    PlaneCenters[5] = std::array<Float_t, 3>{Halfs[0], Halfs[1], Max[2]};

    // The ordering of these may seem odd but it is so that the coordinate in
    // which
    // the BB doesn't care is 'surrounded' (in the %3 sense) by the coordinates
    // to
    // test against.
    PlaneBBs[0] = std::array<Float_t, 4>{Min[2], Min[1], Max[2], Max[1]};
    PlaneBBs[1] = std::array<Float_t, 4>{Min[0], Min[2], Max[0], Max[2]};
    PlaneBBs[2] = std::array<Float_t, 4>{Min[1], Min[0], Max[1], Max[0]};
  }
};

double CalculateToWall(BoundingBox const &BB, TVector3 const &TrStart,
              TVector3 const &TrDir) {
  double ToWall = 0;

  for (size_t i = 0; i < 6; ++i) {
    TVector3 PC(BB.PlaneCenters[i][0], BB.PlaneCenters[i][1],
                BB.PlaneCenters[i][2]);
    TVector3 PN(BB.PlaneNormals[i][0], BB.PlaneNormals[i][1],
                BB.PlaneNormals[i][2]);
    TVector3 TStrToP = TrStart - PC;

    float D = PN.Dot(TrDir);
    float N = -1.0 * PN.Dot(TStrToP);

    if (fabs(D) < 1E-8) {  // segment is parallel to plane
      if (N == 0) {        // segment lies in plane
        continue;
      }
    }
    // they are not parallel
    // compute intersect param
    float sI = N / D;

    TVector3 CrossingPoint = TrStart + sI * TrDir;

    int ci = i >> 1;  // divide by 2 and floor
    int cj = (ci + 2) % 3;
    int ck = (ci + 1) % 3;

    // These need to be strictly lt/gt to coincide with the PointIsInBox
    // definition (of not >= and not <=)
    bool CrossInPlaneBB = (BB.PlaneBBs[ci][0] < CrossingPoint[cj]) &&
                          (BB.PlaneBBs[ci][2] > CrossingPoint[cj]) &&
                          (BB.PlaneBBs[ci][1] < CrossingPoint[ck]) &&
                          (BB.PlaneBBs[ci][3] > CrossingPoint[ck]);

    if (CrossInPlaneBB) {
      if (sI > 0) {
        if (ToWall) {
          std::cout << "[ERROR]: Found ToWall of " << (sI * TrDir.Mag())
                    << ", but already found one of " << ToWall << std::endl;
          throw;
        }

        ToWall = sI * TrDir.Mag();
      }
    }
  }
  return ToWall;
}

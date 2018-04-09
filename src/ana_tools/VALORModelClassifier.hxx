#ifndef VALORMODELCLASSIFIER_HXX_SEEN
#define VALORMODELCLASSIFIER_HXX_SEEN

#include "InteractionModel.hxx"
#include "DepositsSummaryTreeReader.hxx"

#include <iostream>
#include <utility>
#include <vector>

namespace VALORModel {
enum class TrueClass;
}

#define VARLIST          \
  X(kNu_QE_1, 0.082)     \
  X(kNu_QE_2, 0.23)      \
  X(kNu_QE_3, 0.48)      \
  X(kNuBar_QE_1, 0.087)  \
  X(kNuBar_QE_2, 0.24)   \
  X(kNuBar_QE_3, 0.40)   \
  X(kNu_MEC, 1)          \
  X(kNuBar_MEC, 1)       \
  X(kNu_1Pi0_1, 0.13)    \
  X(kNu_1Pi0_2, 0.23)    \
  X(kNu_1Pi0_3, 0.35)    \
  X(kNu_1PiC_1, 0.13)    \
  X(kNu_1PiC_2, 0.24)    \
  X(kNu_1PiC_3, 0.40)    \
  X(kNuBar_1Pi0_1, 0.16) \
  X(kNuBar_1Pi0_2, 0.27) \
  X(kNuBar_1Pi0_3, 0.35) \
  X(kNuBar_1PiC_1, 0.16) \
  X(kNuBar_1PiC_2, 0.3)  \
  X(kNuBar_1PiC_3, 0.4)  \
  X(kNu_2Pi, 0.22)       \
  X(kNuBar_2Pi, 0.22)    \
  X(kNu_DIS_1, 0.035)    \
  X(kNu_DIS_2, 0.035)    \
  X(kNu_DIS_3, 0.027)    \
  X(kNuBar_DIS_1, 0.01)  \
  X(kNuBar_DIS_2, 0.017) \
  X(kNuBar_DIS_3, 0.017) \
  X(kNu_COH, 1.28)       \
  X(kNuBar_COH, 1.34)    \
  X(kNu_NC, 0.16)        \
  X(kNuBar_NC, 0.16)     \
  X(kNuMu_E_Ratio, 0.03) \
  X(kNVALORDials,0)

#define X(A, B) A,
enum class VALORModel::TrueClass { VARLIST };
#undef X
#define X(A, B)                    \
  case VALORModel::TrueClass::A: { \
    return os << #A;               \
  }

inline std::ostream &operator<<(std::ostream &os, VALORModel::TrueClass te) {
  switch (te) { VARLIST }
  return os;
}
#undef X
#define X(A, B)                    \
  case VALORModel::TrueClass::A: { \
    return B;                      \
  }

inline double GetDefaultDialTweak(VALORModel::TrueClass te) {
  switch (te) { VARLIST }
  return 0;
}
#undef X
#define X(A, B) dv.push_back(std::make_tuple(VALORModel::TrueClass::A, 1, B));
inline std::vector<std::tuple<VALORModel::TrueClass, double, double> >
GetDialValueVector(VALORModel::TrueClass te) {
  std::vector<std::tuple<VALORModel::TrueClass, double, double> > dv;
  VARLIST
  return dv;
}
#undef X
#undef VARLIST

std::vector<VALORModel::TrueClass> GetApplicableDials(DepositsSummary const &ed);

double GetVALORWeight(VALORModel::TrueClass te, double value, DepositsSummary const &ed);

#endif

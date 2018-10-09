#ifndef INTERACTIONMODEL_HXX_SEEN
#define INTERACTIONMODEL_HXX_SEEN

#include <iostream>

namespace InteractionModel {
  enum class TrueChannel;
}

#define VARLIST \
  Y(kQES, 1)    \
  X(k2p2h)       \
  X(kRES)       \
  X(kDIS)       \
  X(kCOH)       \
  X(kNuEEL)     \
  X(kIMD)       \
  X(kUnknown)

#define X(A) A,
#define Y(A, B) A = B,
enum class InteractionModel::TrueChannel { VARLIST };
#undef X
#undef Y
#define X(A)                         \
  case InteractionModel::TrueChannel::A: { \
    return os << #A;                 \
  }
#define Y(A, B)                      \
  case InteractionModel::TrueChannel::A: { \
    return os << #A;                 \
  }
inline std::ostream &operator<<(std::ostream &os, InteractionModel::TrueChannel te) {
  switch (te) { VARLIST }
  return os;
}
#undef X
#undef Y
#undef VARLIST

#endif

#pragma once

#include "TColor.h"

#include <array>

void defc(int &c, std::array<int, 3> cs, std::string const &n) {
  static int currC = 51;
  TColor *col = gROOT->GetColor(currC);
  if (!col) {
    col = new TColor(currC, Float_t(cs[0]) / 255.0, Float_t(cs[1]) / 255.0,
                     Float_t(cs[2]) / 255.0, n.c_str());
  } else {
    col->SetRGB(Float_t(cs[0]) / 255.0, Float_t(cs[1]) / 255.0,
                Float_t(cs[2]) / 255.0);
  }
  c = currC;
  currC++;
}

static int kDUNEOrange;
static int kDUNEBlue;
static int kMSUGreen;
static int kMSUPurple;
static int kMSUMud;
static int kLilac;

void DeclareColors(){
  defc(kDUNEOrange,{242, 159, 84},"kDUNEOrange");
  defc(kDUNEBlue,{125, 172, 213},"kDUNEBlue");
  defc(kMSUGreen,{13, 177, 75},"kMSUGreen");
  defc(kMSUPurple,{110, 0, 95},"kMSUPurple");
  defc(kMSUMud,{200, 154, 88},"kMSUMud");
  defc(kLilac,{191, 134, 191},"kLilac");
}

#pragma once

#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"

#include <functional>
#include <memory>
#include <array>

template <size_t dim> struct HistRank {};
template <> struct HistRank<1> {
  static constexpr size_t size = 1;
  using THType = TH1D;
  using THPtrType = THType *;
  using FillType = std::array<double, size + 1>;

  static void Fill(std::unique_ptr<THType> &h, FillType const &v) {
    h->Fill(v[0], v[1]);
  }
};
template <> struct HistRank<2> {
  static constexpr size_t size = 2;
  using THType = TH2D;
  using THPtrType = THType *;
  using FillType = std::array<double, size + 1>;

  static void Fill(std::unique_ptr<THType> &h, FillType const &v) {
    h->Fill(v[0], v[1], v[2]);
  }
};

template <size_t dim, class EvType> struct PlotProcessor {
  std::unique_ptr<typename HistRank<dim>::THType> hist;

  using PFuncType =
      std::function<typename HistRank<dim>::FillType(EvType const &)>;
  PFuncType processor;

  template <typename... Args>
  PlotProcessor(PFuncType &&proc, Args &&... hargs)
      : hist(new typename HistRank<dim>::THType(std::forward<Args>(hargs)...)),
        processor(std::move(proc)) {
    hist->SetDirectory(nullptr);
  }

  PlotProcessor() : hist(nullptr){};
  PlotProcessor(PlotProcessor &&other)
      : hist(std::move(other.hist)), processor(std::move(other.processor)) {}

  void Process(EvType const &evt) { HistRank<dim>::Fill(hist, processor(evt)); }

  void Write(TDirectory *d) { hist.release()->SetDirectory(d); }

  PlotProcessor Clone(std::string const &plotrename) const {
    PlotProcessor cln;
    cln.hist.reset(dynamic_cast<typename HistRank<dim>::THPtrType>(
        hist->Clone(plotrename.c_str())));
    cln.hist->SetDirectory(nullptr);
    cln.processor = processor;
    return cln;
  }

  std::string GetName() const { return hist->GetName(); }
};

using Plot1DCAF = PlotProcessor<1, CAFReader>;
using Plot2DCAF = PlotProcessor<2, CAFReader>;

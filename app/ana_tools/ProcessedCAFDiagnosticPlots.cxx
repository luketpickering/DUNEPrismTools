#include "CAFReader.hxx"
#include "ROOTUtility.hxx"

#include "DUNETDRNDHelper.hxx"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"

#include <functional>
#include <memory>

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

size_t const kAll = 0;
size_t const kFVSel = 1;
size_t const kNumuSel = 2;

// Expect a fhicl like
// InputFile: ND_FHC_CAF.root
int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cout << "[ERROR]: Expected to be passed a single FHiCL file "
                 "describing the inputs."
              << std::endl;
    return 1;
  }

  fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[1]);

  std::string inpf = ps.get<std::string>("InputFile");
  CAFReader rdr(inpf);

  size_t NMaxEvents =
      ps.get<size_t>("NMaxEvents", std::numeric_limits<size_t>::max());
  size_t fents = std::min(NMaxEvents, rdr.GetEntries());

  std::vector<Plot1DCAF> PlotDefinitions;

  PlotDefinitions.emplace_back(
      [](CAFReader const &ev) -> std::array<double, 2> {
        return {ev.det_x + ev.vtx_x * 1E-2, ev.POTWeight};
      },
      "OffAxisEvRate", "", 4500, -5, 40);

  std::vector<std::vector<Plot1DCAF>> AllPlots;
  AllPlots.resize(3);

  size_t idx = kAll;
  for (std::string const &suff : {"_All", "_FV", "_NumuSel"}) {
    for (Plot1DCAF const &p : PlotDefinitions) {
      AllPlots[idx].emplace_back(p.Clone(p.GetName() + suff));
    }
    idx++;
  }
  PlotDefinitions.clear();

  for (size_t ev_it = 0; ev_it < fents; ++ev_it) {
    rdr.GetEntry(ev_it);

    if (ev_it && (fents / 20) && !(ev_it % (fents / 20))) {
      std::cout << "\r[INFO]: Processing event " << ev_it << "/" << fents
                << " (" << std::round((double(ev_it) / double(fents)) * 100)
                << " %)." << std::flush;
    }

    if (!Ar_Select(rdr)) {
      continue;
    }

    for (Plot1DCAF &p : AllPlots[kAll]) {
      p.Process(rdr);
    }

    if (FV_Select(rdr)) {
      for (Plot1DCAF &p : AllPlots[kFVSel]) {
        p.Process(rdr);
      }
      if (ND_FHC_Select(rdr)) {
        for (Plot1DCAF &p : AllPlots[kNumuSel]) {
          p.Process(rdr);
        }
      }
    }
  }
  std::cout << "\n[INFO]: Processed " << fents << " events." << std::endl;

  std::string oupf = ps.get<std::string>("OutputFile");
  TFile *oupFile = CheckOpenFile(oupf, "RECREATE");

  for (auto &SelectionPlots : AllPlots) {
    for (Plot1DCAF &p : SelectionPlots) {
      p.Write(oupFile);
    }
  }

  oupFile->Write();
  oupFile->Close();
}

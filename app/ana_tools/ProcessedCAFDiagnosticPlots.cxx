#include "CAFReader.hxx"
#include "ROOTUtility.hxx"

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
  typedef TH1D THType;
  typedef std::array<double, size + 1> FillType;

  static void Fill(std::unique_ptr<THType> &h, FillType const &v) {
    h->Fill(v[0], v[1]);
  }
};
template <> struct HistRank<2> {
  static constexpr size_t size = 2;
  typedef TH2D THType;
  typedef std::array<double, size + 1> FillType;

  static void Fill(std::unique_ptr<THType> &h, FillType const &v) {
    h->Fill(v[0], v[1], v[2]);
  }
};

template <size_t dim, class EvType> struct PlotProcessor {
  std::unique_ptr<typename HistRank<dim>::THType> hist;

  typedef std::function<typename HistRank<dim>::FillType(EvType const &)>
      PFuncType;
  PFuncType processor;

  template <typename... Args>
  PlotProcessor(PFuncType &&proc, Args &&... hargs)
      : hist(new typename HistRank<dim>::THType(std::forward<Args>(hargs)...)),
        processor(std::move(proc)) {
    hist->SetDirectory(nullptr);
  }

  void Process(EvType const &evt) { HistRank<dim>::Fill(hist, processor(evt)); }

  void Write(TDirectory *d) { hist.release()->SetDirectory(d); }
};

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

  size_t NEvents = rdr.GetEntries();

  PlotProcessor<1, CAFReader> OffAxisEvRate(
      [](CAFReader const &ev) -> std::array<double, 2> {
        return {ev.det_x + ev.vtx_x * 1E-2, ev.POTWeight};
      },
      "OffAxisEvRate", "", 450, -5, 40);

  for (size_t ev_it = 0; ev_it < NEvents; ++ev_it) {
    rdr.GetEntry(ev_it);

    OffAxisEvRate.Process(rdr);
  }

  std::string oupf = ps.get<std::string>("OutputFile");
  TFile *oupFile = CheckOpenFile(oupf, "RECREATE");

  OffAxisEvRate.Write(oupFile);

  oupFile->Write();
  oupFile->Close();
}

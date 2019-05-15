
#include "GENIESplineReader.hxx"
#include "ROOTUtility.hxx"

#include "fhiclcpp/make_ParameterSet.h"

#include "TFile.h"
#include "TGraph.h"

int main(int argc, char const *argv[]) {

  if (argc < 2) {
    std::cout << "[INFO]: Expected to be run like: " << argv[0]
              << " <input_fhicl_file> <input_gxspl_file>" << std::endl;
    return 1;
  }

  fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[1]);

  std::vector<fhicl::ParameterSet> spline_groups;
  if (ps.has_key("spline_groups")) {
    spline_groups = ps.get<std::vector<fhicl::ParameterSet>>("spline_groups");
  } else {
    spline_groups.push_back(ps);
  }

  for (auto sg : spline_groups) {
    GENIEXSecReader rdr(
        sg.get<std::vector<std::pair<std::string, std::vector<std::string>>>>(
            "Splines"));

    rdr.Read(argv[2]);

    std::vector<std::pair<std::string, TGraph>> tg = rdr.GetTGraphs(0, 10);
    std::cout << "Got " << tg.size() << " Splines." << std::endl;

    TFile *oupF = CheckOpenFile(sg.get<std::string>("OutputFile"), "RECREATE");

    std::vector<double> EVals;
    std::vector<double> Y;

    for (size_t i = 0; i < tg.front().second.GetN(); ++i) {
      double x, y;
      tg.front().second.GetPoint(i, x, y);
      EVals.push_back(x);
      Y.push_back(0);
    }

    TGraph tot(EVals.size());

    for (auto g : tg) {
      for (size_t i = 0; i < EVals.size(); ++i) {
        Y[i] += g.second.Eval(EVals[i]);
      }

      std::string name = g.first;
      g.second.Write(name.c_str());
    }

    for (size_t i = 0; i < EVals.size(); ++i) {
      tot.SetPoint(i, EVals[i], Y[i]);
    }

    tot.Write("tot");

    oupF->Write();
    delete oupF;
  }
}

#include "HistogramBuilder.hxx"

#include "ROOTUtility.hxx"

std::vector<Double_t> ParseJaggedUniformAxis(fhicl::ParameterSet const &ps) {
  bool xisunform = ps.get<bool>("xisuniform", true);
  std::string const &uattr = ps.get<std::string>(xisunform ? "xattr" : "yattr");
  std::vector<fhicl::ParameterSet> const &binning =
      ps.get<std::vector<fhicl::ParameterSet>>("binning");

  std::vector<std::vector<Double_t>> UBinnings;
  for (fhicl::ParameterSet const &NUSlice : binning) {
    UBinnings.push_back(BuildBinEdges(NUSlice.get<std::string>(uattr)));
  }

  size_t NumUBinEdges =
      std::accumulate(UBinnings.begin(), UBinnings.end(), 0,
                      [](size_t sum, std::vector<Double_t> const &v) -> size_t {
                        return sum + v.size();
                      });

  std::vector<Double_t> UBinEdges;
  for (std::vector<Double_t> const &ub : UBinnings) {
    AppendVect(UBinEdges, ub);
  }
  UBinEdges = SanitizeBinEdges(UBinEdges);

  // Expect to have taken the last bin off every slice as a duplicate of the
  // first bin of the next slice, except the last slice
  if (UBinEdges.size() != (NumUBinEdges - (UBinnings.size() - 1))) {
    std::cout << "[ERROR]: Expected uniform binning to have "
              << (NumUBinEdges - (UBinnings.size() - 1)) << " edges, but found "
              << UBinEdges.size()
              << " make sure when specifying your jagged binning that the "
                 "uniform axis bin edges are contiuous and correctly ordered."
              << std::endl;
    throw;
  }

  return UBinEdges;
}

std::vector<TAxis> DetermineOptimalNonUniformBinning(
    TAxis const &UniformAxis,
    std::vector<std::tuple<float, float, double>> const &data, size_t NPerBin,
    double minNU, double maxNU, double minNuBinWidth) {
  std::vector<std::vector<std::pair<float, double>>> sorted_NUdata;
  sorted_NUdata.resize(UniformAxis.GetNbins());

  for (std::tuple<float, float, double> const &p : data) {
    size_t Bin = UniformAxis.FindFixBin(std::get<0>(p));

    if ((Bin == 0) || (Bin > sorted_NUdata.size())) {
      continue;
    }

    Bin--;
    sorted_NUdata[Bin].emplace_back(std::get<1>(p), std::get<2>(p));
  }

  std::vector<TAxis> axes;
  for (std::vector<std::pair<float, double>> &NUrow : sorted_NUdata) {
    std::stable_sort(
        NUrow.begin(), NUrow.end(),
        [](std::pair<float, double> const &l,
           std::pair<float, double> const &r) { return l.first < r.first; });

    std::vector<double> BinEdges;
    BinEdges.push_back((minNU == 0xd34db33f) ? NUrow.front().first : minNU);

    double bin_content = 0;
    double sumw2 = 0;
    // poisson error limit requested
    double fracerrlim = sqrt(double(NPerBin)) / double(NPerBin);
    for (std::pair<float, double> const &nu : NUrow) {
      if (nu.first > maxNU) {
        break;
      }
      if (nu.first < BinEdges.back()) {
        continue;
      }

      bin_content += nu.second;
      sumw2 += pow(nu.second, 2);

      if ((sqrt(sumw2) / bin_content) < fracerrlim) {
        double dist = ceil((nu.first - BinEdges.back()) / minNuBinWidth);

        BinEdges.push_back(BinEdges.back() + (dist * minNuBinWidth));

        bin_content = 0;
        sumw2 = 0;
      }
    }
    if (maxNU == 0xdeadbeef) {
      if (fabs(BinEdges.back() - NUrow.back().first) >
          (1E-8 * NUrow.back().first)) {
        double dist =
            ceil((NUrow.back().first - BinEdges.back()) / minNuBinWidth);
        BinEdges.push_back(BinEdges.back() + dist);
      }
    } else {
      if (fabs(BinEdges.back() - maxNU) > (1E-8 * maxNU)) {
        BinEdges.push_back(maxNU);
      }
    }

    SanitizeBinEdges(BinEdges);

    axes.emplace_back(BinEdges.size() - 1, BinEdges.data());
  }

  return axes;
}

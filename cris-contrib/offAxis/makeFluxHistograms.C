void makeFluxHistograms(std::string flux_directory,std::string pot_per_file, std::string user_position_x, std::string user_position_y,std::string user_position_z, std::string output_name) {

  //  gROOT->ProcessLine(".L eventRates.C++");
  gROOT->ProcessLine(".L /mnt/home/f0003917/DK2NU/v01_05_01/dk2nu/tree/dkmeta.cc++"); // Get DK2NU stuff
  gROOT->ProcessLine(".L /mnt/home/f0003917/DK2NU/v01_05_01/dk2nu/tree/dk2nu.cc++"); // Get DK2NU stuff
  gROOT->ProcessLine(".L eventRates_dk2nu.C++");

  std::string command = "eventRates_dk2nu t(\""+flux_directory+"\", "+pot_per_file+", "+user_position_x+", "+user_position_y+", "+user_position_z+",\""+output_name+"\")";

  std::cout<<"Executing: "<<command<<std::endl;

  gROOT->ProcessLine(command.c_str());

  gROOT->ProcessLine("t.Loop()");
}

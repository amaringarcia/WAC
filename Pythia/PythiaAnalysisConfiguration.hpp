#ifndef PYTHIAANALYSISCONFIGURATION_H
#define PYTHIAANALYSISCONFIGURATION_H

#include <TObject.h>
#include <TError.h>
#include "ParticleFilter.hpp"
#include "AnalysisConfiguration.hpp"

class PythiaAnalysisConfiguration : public TObject
{
 public:
  PythiaAnalysisConfiguration() {}
  ~PythiaAnalysisConfiguration() {}

 public:
  std::vector<float> abs_y = {10.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.8};
  float min_pt = 0.0;
  float max_pt = 1e6;
  std::vector<double> ptbins = {0.0,
                                0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                                0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4,
                                1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
                                3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 13.0, 20.0};
  int n_ptbins = ptbins.size() - 1;

  long nEventsRequested = 1000000;
  int nEventsReport = 10000;
  std::string logLevel = "info";
  std::vector<std::string> pythiaOptions = {
    "Init:showChangedSettings = on",
    "Init:showChangedParticleData = on",
    "Next:numberCount = 10000",
    "Next:numberShowInfo = 1",
    "Next:numberShowProcess = 0",
    "Next:numberShowEvent = 0",
    "SoftQCD:all = on",
    "Tune:pp = 14"};
  int projectileA = 2212;
  int projectileB = 2212;
  float energy = 7000.0;
  std::string geventfilter = "MB";
  std::string gparticlefilter = "All";
  std::string gchargefilter = "All";
  bool inrapidity = true;
  std::string outputfname = "PYTHIA8_Pairs_%%03d_";
  std::string taskname = "Rapidity%%03dAll";
  std::string teventfilter = "MB";
  std::vector<std::string> tpairs = {"AllP", "AllM"};
  std::vector<std::string> tsingles = {"AllA"};

  template <AnalysisConfiguration::RapidityPseudoRapidity r>
  static ParticleFilter<r>* particleFilter(std::string str, AnalysisConfiguration* ac)
  {
    if (str == "PiP") {
      return new ParticleFilter<r>(ParticleFilter<r>::Pion, ParticleFilter<r>::Positive, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
    } else if (str == "PiM") {
      return new ParticleFilter<r>(ParticleFilter<r>::Pion, ParticleFilter<r>::Negative, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
    } else if (str == "Pi0") {
      return new ParticleFilter<r>(ParticleFilter<r>::Pion, ParticleFilter<r>::Neutral, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
    } else if (str == "KaP") {
      return new ParticleFilter<r>(ParticleFilter<r>::Kaon, ParticleFilter<r>::Positive, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
    } else if (str == "KaM") {
      return new ParticleFilter<r>(ParticleFilter<r>::Kaon, ParticleFilter<r>::Negative, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
    } else if (str == "Ka0") {
      return new ParticleFilter<r>(ParticleFilter<r>::Kaon, ParticleFilter<r>::Neutral, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
    } else if (str == "PrP") {
      return new ParticleFilter<r>(ParticleFilter<r>::Proton, ParticleFilter<r>::Positive, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
    } else if (str == "PrM") {
      return new ParticleFilter<r>(ParticleFilter<r>::Proton, ParticleFilter<r>::Negative, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
    } else if (str == "AllP") {
      return new ParticleFilter<r>(ParticleFilter<r>::AllSpecies, ParticleFilter<r>::Positive, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
    } else if (str == "AllM") {
      return new ParticleFilter<r>(ParticleFilter<r>::AllSpecies, ParticleFilter<r>::Negative, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
    } else if (str == "AllA") {
      return new ParticleFilter<r>(ParticleFilter<r>::AllSpecies, ParticleFilter<r>::AllCharges, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
    } else {
      ::Error("main", "Paricle species %s still not supported for analysis. Please fix it!!", str.c_str());
      return (ParticleFilter<r>*)nullptr;
    }
  }

  ClassDef(PythiaAnalysisConfiguration, 1)
};

#endif // PYTHIAANALYSISCONFIGURATION_H

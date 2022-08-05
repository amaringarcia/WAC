//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include "Event.hpp"
#include "AnalysisConfiguration.hpp"
#include "TwoPartDiffCorrelationAnalyzer.hpp"
#include "TwoPartDiffCorrelationAnalyzerME.hpp"
#include "ParticleAnalyzer.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaConfiguration.hpp"
#include "PythiaEventGenerator.hpp"

int main(int argc, char* argv[])
{
  if (argc > 3 or argc < 2) {
    Fatal("main", "Wrong number of arguments. Use RunPythiaSimulationTwoParticleDiff jobix seed");
  }
  int jobix = stoi(argv[1]);
  long seed = stol(argv[2]);

  time_t begin, end; // time_t is a datatype to store time values.
  time(&begin);      // note time before execution

  // =========================
  // Short configuration
  // =========================
  float min_y = -2;
  float max_y = 2;
  float min_pt = 0.0;
  float max_pt = 1e6;
  int nBins_y = int((max_y - min_y) / 0.1);
  // int nBins_pt = int((max_pt - min_pt) / 0.1);

  long nEventsRequested = 1000000;
  int nEventsReport = 10000;
  MessageLogger::LogLevel logLevel = MessageLogger::Info;

  // ==========================
  // Event Section
  // ==========================
  Event* event = Event::getEvent();

  // ==========================
  // Generator Section
  // ==========================
  int nOptions = 0;
  TString** pythiaOptions = new TString*[50];
  /* no decay of strange barions */
  /* pythiaOptions[nOptions++] = new TString("130:mayDecay = off");
  pythiaOptions[nOptions++] = new TString("310:mayDecay = off");
  pythiaOptions[nOptions++] = new TString("311:mayDecay = off");
  pythiaOptions[nOptions++] = new TString("3112:mayDecay = off");
  pythiaOptions[nOptions++] = new TString("3222:mayDecay = off");
  pythiaOptions[nOptions++] = new TString("3212:mayDecay = off");
  pythiaOptions[nOptions++] = new TString("3322:mayDecay = off");
  pythiaOptions[nOptions++] = new TString("3312:mayDecay = off");
  pythiaOptions[nOptions++] = new TString("3334:mayDecay = off");
  pythiaOptions[nOptions++] = new TString("3122:mayDecay = off"); */

  pythiaOptions[nOptions++] = new TString("Init:showChangedSettings = on");      // list changed settings
  pythiaOptions[nOptions++] = new TString("Init:showChangedParticleData = on");  // list changed particle data
  pythiaOptions[nOptions++] = new TString("Next:numberCount = 10000");           // print message every n events
  pythiaOptions[nOptions++] = new TString("Next:numberShowInfo = 1");            // print event information n times
  pythiaOptions[nOptions++] = new TString("Next:numberShowProcess = 0");         // print process record n times
  pythiaOptions[nOptions++] = new TString("Next:numberShowEvent = 0");
  pythiaOptions[nOptions++] = new TString("Random:setSeed = on");
  pythiaOptions[nOptions++] = new TString(TString::Format("Random:seed = %ld", seed));
  pythiaOptions[nOptions++] = new TString("SoftQCD:all = on"); // Allow total sigma = elastic/SD/DD/ND
  pythiaOptions[nOptions++] = new TString("Tune:pp = 14");     // Monash 2013
  // pythiaOptions[nOptions++] = new TString("HardQCD:all = on");
  PythiaConfiguration* pc = new PythiaConfiguration(2212 /* p */,
                                                    2212 /* p */,
                                                    7000.0, /* energy in GeV */
                                                    nOptions,
                                                    pythiaOptions);
  EventFilter* eventFilterGen = new EventFilter(EventFilter::MinBias, 0.0, 0.0);
  ParticleFilter<AnalysisConfiguration::kRapidity>* particleFilterGen = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::AllSpecies,
                                                                                                                             ParticleFilter<AnalysisConfiguration::kRapidity>::AllCharges,
                                                                                                                             min_pt, max_pt,
                                                                                                                             min_y, max_y);
  PythiaEventGenerator<AnalysisConfiguration::kRapidity>* generator = new PythiaEventGenerator<AnalysisConfiguration::kRapidity>("PYTHIA", pc, event, eventFilterGen, particleFilterGen);
  generator->reportLevel = logLevel;

  // ==========================
  // Analysis Section
  // ==========================
  AnalysisConfiguration* ac = new AnalysisConfiguration("PYTHIA", "PYTHIA", "1.0");
  ac->loadHistograms = false;
  ac->createHistograms = true;
  ac->scaleHistograms = false;
  ac->calculateDerivedHistograms = false;
  ac->saveHistograms = true;
  ac->resetHistograms = false;
  ac->clearHistograms = false;
  ac->forceHistogramsRewrite = true;
  ac->inputPath = "Input/";
  ac->rootInputFileName = "";
  ac->outputPath = "Output/";
  ac->rootOuputFileName = TString::Format("PYTHIA8_Pairs_%03d_", jobix).Data();
  ac->histoBaseName = "TEST";

  ac->bin_edges_pt = {0.0,
                      0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                      0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4,
                      1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
                      3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 13.0, 20.0};
  ac->nBins_pt = ac->bin_edges_pt.size() - 1;
  ac->min_pt = min_pt;
  ac->max_pt = max_pt;
  ac->nBins_eta = nBins_y;
  ac->min_eta = min_y;
  ac->max_eta = max_y;
  ac->nBins_y = nBins_y;
  ac->min_y = min_y;
  ac->max_y = max_y;
  ac->nBins_phi = 72;
  ac->min_phi = 0.0;
  ac->max_phi = kTWOPI;

  ac->fillPairs = true;
  ac->fillPratt = false;
  ac->fill3D = false;
  ac->fillYorEta = AnalysisConfiguration::kRapidity;

  TString taskName;
  int nAnalysisTasks = 20;
  Task** analysisTasks = new Task*[nAnalysisTasks];

  EventFilter* eventFilter = new EventFilter(EventFilter::MinBias, 0.0, 0.0);
  std::vector<ParticleFilter<AnalysisConfiguration::kRapidity>*> particleFilters;
  particleFilters.push_back(new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Pion, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Pion, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Kaon, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Kaon, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Proton, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Proton, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));
  // particleFilters.push_back(new ParticleFilter(ParticleFilter::Lambda, ParticleFilter::Neutral, ac->min_pt, ac->max_pt, ac->min_eta, ac->max_eta, ac->min_y, ac->max_y));
  // particleFilters.push_back(new ParticleFilter(ParticleFilter::ALambda, ParticleFilter::Neutral, ac->min_pt, ac->max_pt, ac->min_eta, ac->max_eta, ac->min_y, ac->max_y));
  // particleFilters.push_back(new ParticleFilter(ParticleFilter::Photon, ParticleFilter::Neutral, ac->min_pt, ac->max_pt, ac->min_eta, ac->max_eta, ac->min_y, ac->max_y));

  int iTask = 0;
  /* the two-particle analyzer */
  analysisTasks[iTask++] = new TwoPartDiffCorrelationAnalyzerME<AnalysisConfiguration::kRapidity>("NarrowPiKaPrLa", ac, event, eventFilter, particleFilters);
  analysisTasks[iTask - 1]->reportLevel = logLevel;

  /* single particle analysis filters */
  /* int nParticleFilters = 1;
  ParticleFilter** singleParticleFilters = new ParticleFilter*[nParticleFilters];
  singleParticleFilters[0] = new ParticleFilter(ParticleFilter::AllSpecies, ParticleFilter::AllCharges, ac->min_pt, ac->max_pt, ac->min_eta, ac->max_eta, ac->min_y, ac->max_y);
  analysisTasks[iTask++] = new ParticleAnalyzer("SinglesPythia", ac, event, eventFilter, nParticleFilters, singleParticleFilters);
  analysisTasks[iTask - 1]->reportLevel = logLevel; */

  nAnalysisTasks = iTask;

  // ==========================
  // Event Loop
  // ==========================

  EventLoop* eventLoop = new EventLoop();
  eventLoop->addTask((Task*)generator);
  for (int iAnalysisTask = 0; iAnalysisTask < nAnalysisTasks; iAnalysisTask++) {
    eventLoop->addTask(analysisTasks[iAnalysisTask]);
  }
  eventLoop->run(nEventsRequested, nEventsReport);

  cout << "<INFO> PYTHIA Model Analysis - Pair Differential Correlations Histograms  - Completed" << endl;
  time(&end); // note time after execution
  double difference = difftime(end, begin);
  cout << "<INFO> in " << difference << " seconds" << endl;
}

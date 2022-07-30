//
//  ParticlePairCombinedDiffHistos.cpp
//  MyMC
//
//  Created by Claude Pruneau on 12/18/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//
#include "ParticlePairCombinedDiffHistos.hpp"

ParticlePairCombinedDiffHistos::ParticlePairCombinedDiffHistos(const TString& name,
                                                               AnalysisConfiguration* configuration,
                                                               LogLevel debugLevel)
  : Histograms(name, configuration, 100, debugLevel),
    h_R2_DetaDphi_shft(nullptr),
    h_P2_DetaDphi_shft(nullptr),
    h_G2_DetaDphi_shft(nullptr),
    h_R2_DyDphi_shft(nullptr),
    h_P2_DyDphi_shft(nullptr),
    h_G2_DyDphi_shft(nullptr)
{
  if (reportDebug())
    cout << "ParticlePairCombinedDiffHistos::CTOR() Started." << endl;
  initialize();
  if (reportDebug())
    cout << "ParticlePairCombinedDiffHistos::CTOR() Completed." << endl;
}

ParticlePairCombinedDiffHistos::ParticlePairCombinedDiffHistos(TFile* inputFile,
                                                               const TString& collectionName,
                                                               AnalysisConfiguration* analysisConfiguration,
                                                               LogLevel debugLevel)
  : Histograms(collectionName, analysisConfiguration, 100, debugLevel),
    h_R2_DetaDphi_shft(nullptr),
    h_P2_DetaDphi_shft(nullptr),
    h_G2_DetaDphi_shft(nullptr),
    h_R2_DyDphi_shft(nullptr),
    h_P2_DyDphi_shft(nullptr),
    h_G2_DyDphi_shft(nullptr)
{
  if (reportDebug())
    cout << "ParticlePairCombinedDiffHistos::CTOR(TFile * inputFile,...) Started." << endl;
  loadHistograms(inputFile);
  if (reportDebug())
    cout << "ParticlePairCombinedDiffHistos::CTOR(TFile * inputFile,...) Completed." << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////////////////////////////////////////////
ParticlePairCombinedDiffHistos::~ParticlePairCombinedDiffHistos()
{
  if (reportDebug())
    cout << "ParticlePairCombinedDiffHistos::DTOR() No ops." << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// load histograms
//////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticlePairCombinedDiffHistos::loadHistograms(TFile* inputFile)
{
  if (reportDebug())
    cout << "ParticlePairCombinedDiffHistos::loadHistograms(...) strted" << endl;

  TString bn = getHistoBaseName();

  if (configuration->fillYorEta == configuration->kPseudorapidity) {
    h_R2_DetaDphi_shft = loadH2(inputFile, bn + TString("R2_DetaDphi_shft"));
    h_P2_DetaDphi_shft = loadH2(inputFile, bn + TString("P2_DetaDphi_shft"));
    h_G2_DetaDphi_shft = loadH2(inputFile, bn + TString("G2_DetaDphi_shft"));
  }
  if (configuration->fillYorEta == configuration->kRapidity) {
    h_R2_DyDphi_shft = loadH2(inputFile, bn + TString("R2_DyDphi_shft"));
    h_P2_DyDphi_shft = loadH2(inputFile, bn + TString("P2_DyDphi_shft"));
    h_G2_DyDphi_shft = loadH2(inputFile, bn + TString("G2_DyDphi_shft"));
  }
  if (reportDebug())
    cout << "ParticlePairCombinedDiffHistos::loadHistograms(...) completed" << endl;
  /* the histograms are not owned */
  bOwnTheHistograms = false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// initialize histograms of this class
//////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticlePairCombinedDiffHistos::createHistograms()
{
  if (reportDebug())
    cout << "ParticlePairCombinedDiffHistos::initialize() started" << endl;
  TString name;

  AnalysisConfiguration& ac = *getConfiguration();
  TString bn = getHistoBaseName();

  if (ac.fillYorEta == ac.kPseudorapidity) {
    h_R2_DetaDphi_shft = createHistogram(bn + TString("R2_DetaDphi_shft"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "R_{2}", scaled, saved, plotted, printed);
    h_P2_DetaDphi_shft = createHistogram(bn + TString("P2_DetaDphi_shft"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "P_{2}", scaled, saved, plotted, printed);
    h_G2_DetaDphi_shft = createHistogram(bn + TString("G2_DetaDphi_shft"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "G_{2}", scaled, saved, plotted, printed);
  }
  if (ac.fillYorEta == ac.kRapidity) {
    h_R2_DyDphi_shft = createHistogram(bn + TString("R2_DyDphi_shft"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "R_{2}", scaled, saved, plotted, printed);
    h_P2_DyDphi_shft = createHistogram(bn + TString("P2_DyDphi_shft"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "P_{2}", scaled, saved, plotted, printed);
    h_G2_DyDphi_shft = createHistogram(bn + TString("G2_DyDphi_shft"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "G_{2}", scaled, saved, plotted, printed);
  }

  if (reportDebug())
    cout << "ParticlePairCombinedDiffHistos::initialize() completed" << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// calculate combinations of histograms with given coefficients
//////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticlePairCombinedDiffHistos::calculate(ParticlePairDerivedDiffHistos* pp,
                                               ParticlePairDerivedDiffHistos* mm,
                                               ParticlePairDerivedDiffHistos* pm,
                                               ParticlePairDerivedDiffHistos* mp,
                                               double app, double amm, double apm, double amp)
{
  if (reportDebug())
    cout << "ParticlePairCombinedDiffHistos::calculate() started" << endl;

  if (configuration->fillYorEta == configuration->kPseudorapidity) {
    h_R2_DetaDphi_shft->Reset();
    h_R2_DetaDphi_shft->Add(pp->h_R2_DetaDphi_shft, mm->h_R2_DetaDphi_shft, app, amm);
    h_R2_DetaDphi_shft->Add(pm->h_R2_DetaDphi_shft, apm);
    h_R2_DetaDphi_shft->Add(mp->h_R2_DetaDphi_shft, amp);

    h_P2_DetaDphi_shft->Reset();
    h_P2_DetaDphi_shft->Add(pp->h_P2_DetaDphi_shft, mm->h_P2_DetaDphi_shft, app, amm);
    h_P2_DetaDphi_shft->Add(pm->h_P2_DetaDphi_shft, apm);
    h_P2_DetaDphi_shft->Add(mp->h_P2_DetaDphi_shft, amp);

    h_G2_DetaDphi_shft->Reset();
    h_G2_DetaDphi_shft->Add(pp->h_G2_DetaDphi_shft, mm->h_G2_DetaDphi_shft, app, amm);
    h_G2_DetaDphi_shft->Add(pm->h_G2_DetaDphi_shft, apm);
    h_G2_DetaDphi_shft->Add(mp->h_G2_DetaDphi_shft, amp);
  }
  if (configuration->fillYorEta == configuration->kRapidity) {
    h_R2_DetaDphi_shft->Reset();
    h_R2_DyDphi_shft->Add(pp->h_R2_DyDphi_shft, mm->h_R2_DyDphi_shft, app, amm);
    h_R2_DyDphi_shft->Add(pm->h_R2_DyDphi_shft, apm);
    h_R2_DyDphi_shft->Add(mp->h_R2_DyDphi_shft, amp);

    h_P2_DetaDphi_shft->Reset();
    h_P2_DyDphi_shft->Add(pp->h_P2_DyDphi_shft, mm->h_P2_DyDphi_shft, app, amm);
    h_P2_DyDphi_shft->Add(pm->h_P2_DyDphi_shft, apm);
    h_P2_DyDphi_shft->Add(mp->h_P2_DyDphi_shft, amp);

    h_G2_DetaDphi_shft->Reset();
    h_G2_DyDphi_shft->Add(pp->h_G2_DyDphi_shft, mm->h_G2_DyDphi_shft, app, amm);
    h_G2_DyDphi_shft->Add(pm->h_G2_DyDphi_shft, apm);
    h_G2_DyDphi_shft->Add(mp->h_G2_DyDphi_shft, amp);
  }

  if (reportDebug())
    cout << "ParticlePairCombinedDiffHistos::calculate() completed" << endl;
}

ClassImp(ParticlePairCombinedDiffHistos)

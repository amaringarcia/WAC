// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining Task
 */

#include "ParticlePairDiffHistos.hpp"

ParticlePairDiffHistos::ParticlePairDiffHistos(const TString& name,
                                               AnalysisConfiguration* configuration,
                                               LogLevel debugLevel)
  : Histograms(name, configuration, 150, debugLevel),
    h_n2_ptPt(nullptr),
    h_n2_DetaDphi(nullptr),
    p_n2_DetaDphi(nullptr),
    h_ptpt_DetaDphi(nullptr),
    h_dptdpt_DetaDphi(nullptr),
    h_n2_DyDphi(nullptr),
    p_n2_DyDphi(nullptr),
    h_ptpt_DyDphi(nullptr),
    h_dptdpt_DyDphi(nullptr)
{
  if (reportDebug())
    cout << "ParticlePairDiffHistos::CTOR() Started." << endl;
  initialize();
  if (reportDebug())
    cout << "ParticlePairDiffHistos::CTOR() Completed." << endl;
}

ParticlePairDiffHistos::ParticlePairDiffHistos(TFile* inputFile,
                                               const TString& name,
                                               AnalysisConfiguration* configuration,
                                               LogLevel debugLevel)
  : Histograms(name, configuration, 150, debugLevel),
    h_n2_ptPt(nullptr),
    h_n2_DetaDphi(nullptr),
    p_n2_DetaDphi(nullptr),
    h_ptpt_DetaDphi(nullptr),
    h_dptdpt_DetaDphi(nullptr),
    h_n2_DyDphi(nullptr),
    p_n2_DyDphi(nullptr),
    h_ptpt_DyDphi(nullptr),
    h_dptdpt_DyDphi(nullptr)
{
  if (reportDebug())
    cout << "ParticlePairDiffHistos::CTOR() Started." << endl;
  loadHistograms(inputFile);
  if (reportDebug())
    cout << "ParticlePairDiffHistos::CTOR() Completed." << endl;
}

ParticlePairDiffHistos::~ParticlePairDiffHistos()
{
  // deleteHistograms();
}

void ParticlePairDiffHistos::initialize()
{
  if (reportDebug())
    cout << "ParticlePairDiffHistos::initialize() Started." << endl;
  AnalysisConfiguration& ac = *(AnalysisConfiguration*)getConfiguration();
  TString bn = getHistoBaseName();
  if (ac.bin_edges_pt.size() > 0) {
    ac.nBins_pt = ac.bin_edges_pt.size() - 1;
  }
  ac.range_pt = ac.max_pt - ac.min_pt;
  ac.range_phi = ac.max_phi - ac.min_phi;
  ac.range_eta = ac.max_eta - ac.min_eta;
  ac.range_y = ac.max_y - ac.min_y;
  ac.nBins_phiEta = ac.nBins_eta * ac.nBins_phi;
  ac.nBins_phiEtaPt = ac.nBins_eta * ac.nBins_phi * ac.nBins_pt;
  ac.nBins_phiY = ac.nBins_y * ac.nBins_phi;
  ac.nBins_phiYPt = ac.nBins_y * ac.nBins_phi * ac.nBins_pt;

  ac.nBins_Deta = ac.nBins_eta * 2 - 1;
  ac.min_Deta = ac.min_eta - ac.max_eta;
  ac.max_Deta = ac.max_eta - ac.min_eta;
  ac.width_Deta = (ac.max_eta - ac.min_eta) / ac.nBins_Deta;
  ac.nBins_Dphi = ac.nBins_phi;
  ac.nBins_Dphi_shft = ac.nBins_phi / 4;
  ac.width_Dphi = kTWOPI / ac.nBins_Dphi;
  ac.min_Dphi = -ac.width_Dphi / 2.;
  ac.max_Dphi = kTWOPI - ac.width_Dphi / 2.;
  ac.min_Dphi_shft = ac.min_Dphi - ac.width_Dphi * ac.nBins_Dphi_shft;
  ac.max_Dphi_shft = ac.max_Dphi - ac.width_Dphi * ac.nBins_Dphi_shft;
  ac.nBins_Dy = 2 * ac.nBins_y - 1;
  ac.width_Dy = (ac.max_y - ac.min_y) / ac.nBins_Dy;
  ac.min_Dy = ac.min_y - ac.max_y;
  ac.max_Dy = ac.max_y - ac.min_y;

  /* big histograms are forced to be created without sumw2 structure for it will not be used */
  bool defsumw2 = TH1::GetDefaultSumw2();
  TH1::SetDefaultSumw2(false);

  if (ac.bin_edges_pt.size() > 0) {
    h_n2_ptPt = createHistogram(bn + TString("n2_ptPt"), ac.bin_edges_pt, ac.bin_edges_pt, "p_{T,1}", "p_{T,2}", "N_{2}", scaled, saved, plotted, notPrinted);
  } else {
    h_n2_ptPt = createHistogram(bn + TString("n2_ptPt"), ac.nBins_pt, ac.min_pt, ac.max_pt, ac.nBins_pt, ac.min_pt, ac.max_pt, "p_{T,1}", "p_{T,2}", "N_{2}", scaled, saved, plotted, notPrinted);
  }
  if (ac.fillYorEta == ac.kPseudorapidity) {
    h_n2_DetaDphi = createHistogram(bn + TString("n2_DetaDphi"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta#eta", "#Delta#varphi", "n_{2}", scaled, saved, notPlotted, notPrinted, false);
    p_n2_DetaDphi = createProfile(bn + TString("p_n2_DetaDphi"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta#eta", "#Delta#varphi", "n_{2}", notScaled, saved, notPlotted, notPrinted, false);
    h_ptpt_DetaDphi = createHistogram(bn + TString("ptpt_DetaDphi"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta#eta", "#Delta#varphi", "p_{T}xp_{T}", scaled, saved, notPlotted, notPrinted, false);
    h_dptdpt_DetaDphi = createHistogram(bn + TString("dptdpt_DetaDphi"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta#eta", "#Delta#varphi", "(p_{T}-<p_{T}>)x(p_{T}-<p_{T}>)", scaled, saved, notPlotted, notPrinted, false);

    /* big histograms are forced to be created without sumw2 structure for it will not be used */
    h_n2_DetaDphi->SetBit(TH1::kIsNotW);
    h_n2_DetaDphi->Sumw2(false);
    p_n2_DetaDphi->SetBit(TH1::kIsNotW);
    p_n2_DetaDphi->Sumw2(false);
    h_ptpt_DetaDphi->SetBit(TH1::kIsNotW);
    h_ptpt_DetaDphi->Sumw2(false);
    h_dptdpt_DetaDphi->SetBit(TH1::kIsNotW);
    h_dptdpt_DetaDphi->Sumw2(false);
  }
  if (ac.fillYorEta == ac.kRapidity) {
    h_n2_DyDphi = createHistogram(bn + TString("n2_DyDphi"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta y", "#Delta#varphi", "<n_{2}>", scaled, saved, notPlotted, notPrinted, false);
    p_n2_DyDphi = createProfile(bn + TString("p_n2_DyDphi"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta y", "#Delta#varphi", "<n_{2}>", notScaled, saved, notPlotted, notPrinted, false);
    h_ptpt_DyDphi = createHistogram(bn + TString("ptpt_DyDphi"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta y", "#Delta#varphi", "p_{T}xp_{T}", scaled, saved, notPlotted, notPrinted, false);
    h_dptdpt_DyDphi = createHistogram(bn + TString("dptdpt_DyDphi"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta y", "#Delta#varphi", "p_{T}xp_{T}", scaled, saved, notPlotted, notPrinted, false);

    /* big histograms are forced to be created without sumw2 structure for it will not be used */
    h_n2_DyDphi->SetBit(TH1::kIsNotW);
    h_n2_DyDphi->Sumw2(false);
    p_n2_DyDphi->SetBit(TH1::kIsNotW);
    p_n2_DyDphi->Sumw2(false);
    h_ptpt_DyDphi->SetBit(TH1::kIsNotW);
    h_ptpt_DyDphi->Sumw2(false);
    h_dptdpt_DyDphi->SetBit(TH1::kIsNotW);
    h_dptdpt_DyDphi->Sumw2(false);
  }

  /* back to default behavior */
  TH1::SetDefaultSumw2(defsumw2);
  if (reportDebug())
    cout << "ParticlePairDiffHistos::initialize() Completed." << endl;
}

//________________________________________________________________________
void ParticlePairDiffHistos::loadHistograms(TFile* inputFile)
{
  if (!inputFile) {
    cout << "-Fatal- Attempting to load ParticleHistos from an invalid file pointer" << endl;
    return;
  }
  AnalysisConfiguration& ac = *(AnalysisConfiguration*)getConfiguration();
  TString bn = getHistoBaseName();
  if (ac.bin_edges_pt.data()) {
    ac.nBins_pt = ac.bin_edges_pt.size() - 1;
  }
  ac.range_pt = ac.max_pt - ac.min_pt;
  ac.range_phi = ac.max_phi - ac.min_phi;
  ac.range_eta = ac.max_eta - ac.min_eta;
  ac.range_y = ac.max_y - ac.min_y;
  ac.nBins_phiEta = ac.nBins_eta * ac.nBins_phi;
  ac.nBins_phiEtaPt = ac.nBins_eta * ac.nBins_phi * ac.nBins_pt;
  ac.nBins_phiY = ac.nBins_y * ac.nBins_phi;
  ac.nBins_phiYPt = ac.nBins_y * ac.nBins_phi * ac.nBins_pt;

  ac.nBins_Deta = ac.nBins_eta * 2 - 1;
  ac.width_Deta = (ac.max_eta - ac.min_eta) / ac.nBins_Deta;
  ac.min_Deta = ac.min_eta - ac.max_eta;
  ac.max_Deta = ac.max_eta - ac.min_eta;
  ac.nBins_Dphi = ac.nBins_phi;
  ac.nBins_Dphi_shft = ac.nBins_phi / 4;
  ac.width_Dphi = kTWOPI / ac.nBins_Dphi;
  ac.min_Dphi = -ac.width_Dphi / 2.;
  ac.max_Dphi = kTWOPI - ac.width_Dphi / 2.;
  ac.min_Dphi_shft = ac.min_Dphi - ac.width_Dphi * double(ac.nBins_Dphi / 4);
  ac.max_Dphi_shft = ac.max_Dphi - ac.width_Dphi * double(ac.nBins_Dphi / 4);
  ac.nBins_Dy = 2 * ac.nBins_y - 1;
  ac.width_Dy = (ac.max_y - ac.min_y) / ac.nBins_Dy;
  ac.min_Dy = ac.min_y - ac.max_y;
  ac.max_Dy = ac.max_y - ac.min_y;

  if (ac.fillYorEta == ac.kPseudorapidity) {
    h_n2_DetaDphi = loadH2(inputFile, bn + TString("n2_DetaDphi"), true);
    p_n2_DetaDphi = loadProfile2D(inputFile, bn + TString("p_n2_DetaDphi"), false); /* don't downscale the profiles */
    h_ptpt_DetaDphi = loadH2(inputFile, bn + TString("ptpt_DetaDphi"), true);
    h_dptdpt_DetaDphi = loadH2(inputFile, bn + TString("dptdpt_DetaDphi"), true);
  }

  if (ac.fillYorEta == ac.kRapidity) {
    h_n2_DyDphi = loadH2(inputFile, bn + TString("n2_DyDphi"), true);
    p_n2_DyDphi = loadProfile2D(inputFile, bn + TString("p_n2_DyDphi"), false); /* don't downscale the profiles */
    h_ptpt_DyDphi = loadH2(inputFile, bn + TString("ptpt_DyDphi"), true);
    h_dptdpt_DyDphi = loadH2(inputFile, bn + TString("dptdpt_DyDphi"), true);
  }
  /* the histograms are not owned */
  bOwnTheHistograms = false;
  return;
}

ClassImp(ParticlePairDiffHistos)

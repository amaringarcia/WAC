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

#include "ParticlePairHistos.hpp"

ClassImp(ParticlePairHistos)

  ParticlePairHistos::ParticlePairHistos(const TString& name,
                                         AnalysisConfiguration* configuration,
                                         LogLevel debugLevel)
  : Histograms(name, configuration, 150, debugLevel)
{
  if (reportDebug())
    cout << "ParticlePairHistos::CTOR() Started." << endl;
  initialize();
  if (reportDebug())
    cout << "ParticlePairHistos::CTOR() Completed." << endl;
}

ParticlePairHistos::ParticlePairHistos(TFile* inputFile,
                                       const TString& name,
                                       AnalysisConfiguration* configuration,
                                       LogLevel debugLevel)
  : Histograms(name, configuration, 150, debugLevel)
{
  if (reportDebug())
    cout << "ParticlePairHistos::CTOR() Started." << endl;
  loadHistograms(inputFile);
  if (reportDebug())
    cout << "ParticlePairHistos::CTOR() Completed." << endl;
}

ParticlePairHistos::~ParticlePairHistos()
{
  // deleteHistograms();
}

void ParticlePairHistos::initialize()
{
  if (reportDebug())
    cout << "ParticlePairHistos::initialize() Started." << endl;
  AnalysisConfiguration& ac = *(AnalysisConfiguration*)getConfiguration();
  TString bn = getHistoBaseName();
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
  ac.nBins_Dy = ac.nBins_y * 2 - 1;
  ac.min_Dy = ac.min_y - ac.max_y;
  ac.max_Dy = ac.max_y - ac.min_y;
  ac.nBins_Dphi = ac.nBins_phi;
  ac.nBins_Dphi_shft = ac.nBins_phi / 4;
  ac.width_Dphi = TMath::TwoPi() / ac.nBins_Dphi;
  ac.min_Dphi = -ac.width_Dphi / 2.;
  ac.max_Dphi = TMath::TwoPi() - ac.width_Dphi / 2.;
  ac.min_Dphi_shft = ac.min_Dphi - ac.width_Dphi * double(ac.nBins_Dphi / 4);
  ac.max_Dphi_shft = ac.max_Dphi - ac.width_Dphi * double(ac.nBins_Dphi / 4);

  h_n2_ptPt = createHistogram(bn + TString("n2_ptPt"), ac.nBins_pt, ac.min_pt, ac.max_pt, ac.nBins_pt, ac.min_pt, ac.max_pt, "p_{T,1}", "p_{T,2}", "N_{2}", scaled, saved, plotted, notPrinted);
  h_n2_phiPhi = createHistogram(bn + TString("n2_phiPhi"), ac.nBins_phi, ac.min_phi, ac.max_phi, ac.nBins_phi, ac.min_phi, ac.max_phi, "#varphi_{1}", "#varphi_{2}", "N_{2}", scaled, saved, notPlotted, notPrinted);
  h_npt_phiPhi = createHistogram(bn + TString("npt_phiPhi"), ac.nBins_phi, ac.min_phi, ac.max_phi, ac.nBins_phi, ac.min_phi, ac.max_phi, "#varphi_{1}", "#varphi_{2}", "n x p_{T}", scaled, saved, notPlotted, notPrinted);
  h_ptn_phiPhi = createHistogram(bn + TString("ptn_phiPhi"), ac.nBins_phi, ac.min_phi, ac.max_phi, ac.nBins_phi, ac.min_phi, ac.max_phi, "#varphi_{1}", "#varphi_{2}", "p_{T} x n", scaled, saved, notPlotted, notPrinted);
  h_ptpt_phiPhi = createHistogram(bn + TString("ptpt_phiPhi"), ac.nBins_phi, ac.min_phi, ac.max_phi, ac.nBins_phi, ac.min_phi, ac.max_phi, "#varphi_{1}", "#varphi_{2}", "p_{T}xp_{T}", scaled, saved, plotted, notPrinted);

  if (ac.fillYorEta == ac.kPseudorapidity) {
    h_n2_etaEta = createHistogram(bn + TString("n2_etaEta"), ac.nBins_eta, ac.min_eta, ac.max_eta, ac.nBins_eta, ac.min_eta, ac.max_eta, "#eta_{1}", "#eta_{2}", "N_{2}", scaled, saved, plotted, notPrinted);
#ifdef OPTIMIZEADDBINCONTENT
    /* big histograms are forced to be created without sumw2 structure for it will not be used */
    bool defsumw2 = TH1::GetDefaultSumw2();
    TH1::SetDefaultSumw2(false);
#endif // OPTIMIZEADDBINCONTENT
    h_n2_phiEtaPhiEta = createHistogram(bn + TString("n2_phiEtaPhiEta"), ac.nBins_phiEta, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiEta), ac.nBins_phiEta, 0.0, static_cast<double>(ac.nBins_phiEta), "#eta_{1}x#varphi_{1}", "#eta_{2}x#varphi_{2}", "N_{2}", scaled, saved, notPlotted, notPrinted, false);
    h_npt_phiEtaPhiEta = createHistogram(bn + TString("npt_phiEtaPhiEta"), ac.nBins_phiEta, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiEta), ac.nBins_phiEta, 0.0, static_cast<double>(ac.nBins_phiEta), "#eta_{1}x#varphi_{1}", "#eta_{2}x#varphi_{2}", "Nxp_{T}", scaled, saved, notPlotted, notPrinted, false);
    h_ptn_phiEtaPhiEta = createHistogram(bn + TString("ptn_phiEtaPhiEta"), ac.nBins_phiEta, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiEta), ac.nBins_phiEta, 0.0, static_cast<double>(ac.nBins_phiEta), "#eta_{1}x#varphi_{1}", "#eta_{2}x#varphi_{2}", "p_{T}xN", scaled, saved, notPlotted, notPrinted, false);
    h_ptpt_phiEtaPhiEta = createHistogram(bn + TString("ptpt_phiEtaPhiEta"), ac.nBins_phiEta, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiEta), ac.nBins_phiEta, 0.0, static_cast<double>(ac.nBins_phiEta), "#eta_{1}x#varphi_{1}", "#eta_{2}x#varphi_{2}", "p_{T}xp_{T}", scaled, saved, notPlotted, notPrinted, false);

    h_n2_detaDphi_o = createHistogram(bn + TString("n2_detaDphi_o"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "N_{2}", scaled, saved, notPlotted, notPrinted, false);
    h_npt_detaDphi_o = createHistogram(bn + TString("npt_detaDphi_o"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "Nxp_{T}", scaled, saved, notPlotted, notPrinted, false);
    h_ptn_detaDphi_o = createHistogram(bn + TString("ptn_detaDphi_o"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "p_{T}xN", scaled, saved, notPlotted, notPrinted, false);
    h_ptpt_detaDphi_o = createHistogram(bn + TString("ptpt_detaDphi_o"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "p_{T}xp_{T}", scaled, saved, notPlotted, notPrinted, false);
    h_n1n1_detaDphi_o = createHistogram(bn + TString("n1n1_detaDphi_o"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "n_{1}*n_{1}", scaled, saved, notPlotted, notPrinted, false);
    h_pt1pt1_detaDphi_o = createHistogram(bn + TString("pt1pt1_detaDphi_o"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "p_{T}*p_{T}", scaled, saved, notPlotted, notPrinted, false);
#ifdef OPTIMIZEADDBINCONTENT
    /* big histograms are forced to be created without sumw2 structure for it will not be used */
    h_n2_phiEtaPhiEta->SetBit(TH1::kIsNotW);
    h_n2_phiEtaPhiEta->Sumw2(false);
    h_npt_phiEtaPhiEta->SetBit(TH1::kIsNotW);
    h_npt_phiEtaPhiEta->Sumw2(false);
    h_ptn_phiEtaPhiEta->SetBit(TH1::kIsNotW);
    h_ptn_phiEtaPhiEta->Sumw2(false);
    h_ptpt_phiEtaPhiEta->SetBit(TH1::kIsNotW);
    h_ptpt_phiEtaPhiEta->Sumw2(false);
    h_n2_detaDphi_o->SetBit(TH1::kIsNotW);
    h_n2_detaDphi_o->Sumw2(false);
    h_npt_detaDphi_o->SetBit(TH1::kIsNotW);
    h_npt_detaDphi_o->Sumw2(false);
    h_ptn_detaDphi_o->SetBit(TH1::kIsNotW);
    h_ptn_detaDphi_o->Sumw2(false);
    h_ptpt_detaDphi_o->SetBit(TH1::kIsNotW);
    h_ptpt_detaDphi_o->Sumw2(false);
    h_n1n1_detaDphi_o->SetBit(TH1::kIsNotW);
    h_n1n1_detaDphi_o->Sumw2(false);
    h_pt1pt1_detaDphi_o->SetBit(TH1::kIsNotW);
    h_pt1pt1_detaDphi_o->Sumw2(false);
    TH1::SetDefaultSumw2(defsumw2);
#endif // OPTIMIZEADDBINCONTENT
    h_npt_etaEta = createHistogram(bn + TString("npt_etaEta"), ac.nBins_eta, ac.min_eta, ac.max_eta, ac.nBins_eta, ac.min_eta, ac.max_eta, "#eta_{1}", "#eta_{2}", "n x p_{T}", scaled, saved, notPlotted, notPrinted);
    h_ptn_etaEta = createHistogram(bn + TString("ptn_etaEta"), ac.nBins_eta, ac.min_eta, ac.max_eta, ac.nBins_eta, ac.min_eta, ac.max_eta, "#eta_{1}", "#eta_{2}", "p_{T} x n", scaled, saved, notPlotted, notPrinted);
    h_ptpt_etaEta = createHistogram(bn + TString("ptpt_etaEta"), ac.nBins_eta, ac.min_eta, ac.max_eta, ac.nBins_eta, ac.min_eta, ac.max_eta, "#eta_{1}", "#eta_{2}", "p_{T}xp_{T}", scaled, saved, plotted, notPrinted);
  }
  if (ac.fillYorEta == ac.kRapidity) {
    h_n2_yY = createHistogram(bn + TString("n2_yY"), ac.nBins_y, ac.min_y, ac.max_y, ac.nBins_y, ac.min_y, ac.max_y, "y_{1}", "y_{2}", "N_{2}", scaled, saved, notPlotted, notPrinted);
#ifdef OPTIMIZEADDBINCONTENT
    /* big histograms are forced to be created without sumw2 structure for it will not be used */
    bool defsumw2 = TH1::GetDefaultSumw2();
    TH1::SetDefaultSumw2(false);
#endif // OPTIMIZEADDBINCONTENT
    h_n2_phiYPhiY = createHistogram(bn + TString("n2_phiYPhiY"), ac.nBins_phiY, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiY), ac.nBins_phiY, 0.0, static_cast<double>(ac.nBins_phiY), "y_{1}x#varphi_{1}", "y_{2}x#varphi_{2}", "N_{2}", scaled, saved, notPlotted, notPrinted, false);
    h_npt_phiYPhiY = createHistogram(bn + TString("npt_phiYPhiY"), ac.nBins_phiY, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiY), ac.nBins_phiY, 0.0, static_cast<double>(ac.nBins_phiY), "y_{1}x#varphi_{1}", "y_{2}x#varphi_{2}", "Nxp_{T}", scaled, saved, notPlotted, notPrinted, false);
    h_ptn_phiYPhiY = createHistogram(bn + TString("ptn_phiYPhiY"), ac.nBins_phiY, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiY), ac.nBins_phiY, 0.0, static_cast<double>(ac.nBins_phiY), "y_{1}x#varphi_{1}", "y_{2}x#varphi_{2}", "p_{T}xN", scaled, saved, notPlotted, notPrinted, false);
    h_ptpt_phiYPhiY = createHistogram(bn + TString("ptpt_phiYPhiY"), ac.nBins_phiY, static_cast<double>(0.0), static_cast<double>(ac.nBins_phiY), ac.nBins_phiY, 0.0, static_cast<double>(ac.nBins_phiY), "y_{1}x#varphi_{1}", "y_{2}x#varphi_{2}", "p_{T}xp_{T}", scaled, saved, notPlotted, notPrinted, false);
#ifdef OPTIMIZEADDBINCONTENT
    /* big histograms are forced to be created without sumw2 structure for it will not be used */
    h_n2_phiYPhiY->SetBit(TH1::kIsNotW);
    h_n2_phiYPhiY->Sumw2(false);
    h_npt_phiYPhiY->SetBit(TH1::kIsNotW);
    h_npt_phiYPhiY->Sumw2(false);
    h_ptn_phiYPhiY->SetBit(TH1::kIsNotW);
    h_ptn_phiYPhiY->Sumw2(false);
    h_ptpt_phiYPhiY->SetBit(TH1::kIsNotW);
    h_ptpt_phiYPhiY->Sumw2(false);
    TH1::SetDefaultSumw2(defsumw2);
#endif // OPTIMIZEADDBINCONTENT
    h_npt_yY = createHistogram(bn + TString("npt_yY"), ac.nBins_y, ac.min_y, ac.max_y, ac.nBins_y, ac.min_y, ac.max_y, "y_{1}", "y_{2}", "n x p_{T}", scaled, saved, notPlotted, notPrinted);
    h_ptn_yY = createHistogram(bn + TString("ptn_yY"), ac.nBins_y, ac.min_y, ac.max_y, ac.nBins_y, ac.min_y, ac.max_y, "y_{1}", "y_{2}", "p_{T} x n", scaled, saved, notPlotted, notPrinted);
    h_ptpt_yY = createHistogram(bn + TString("ptpt_yY"), ac.nBins_y, ac.min_y, ac.max_y, ac.nBins_y, ac.min_y, ac.max_y, "y_{1}", "y_{2}", "p_{T}xp_{T}", scaled, saved, plotted, notPrinted);
  }
  if (reportDebug())
    cout << "ParticlePairHistos::initialize() Completed." << endl;
}

// complete filling the addicional histograms by projecting the
// higher dimensional ones
void ParticlePairHistos::completeFill()
{
  AnalysisConfiguration& ac = *(AnalysisConfiguration*)getConfiguration();

  if (ac.fillYorEta == ac.kPseudorapidity) {
    int nbinseta = h_n2_etaEta->GetNbinsX();
    int nbinsphi = h_n2_phiPhi->GetNbinsX();

    project_n2XYXY_n2XX(h_n2_phiEtaPhiEta, h_n2_etaEta, nbinseta, nbinsphi);
    project_n2XYXY_n2YY(h_n2_phiEtaPhiEta, h_n2_phiPhi, nbinseta, nbinsphi);

    project_n2XYXY_n2XX(h_ptn_phiEtaPhiEta, h_ptn_etaEta, nbinseta, nbinsphi);
    project_n2XYXY_n2XX(h_npt_phiEtaPhiEta, h_npt_etaEta, nbinseta, nbinsphi);
    project_n2XYXY_n2XX(h_ptpt_phiEtaPhiEta, h_ptpt_etaEta, nbinseta, nbinsphi);

    project_n2XYXY_n2YY(h_ptn_phiEtaPhiEta, h_ptn_phiPhi, nbinseta, nbinsphi);
    project_n2XYXY_n2YY(h_npt_phiEtaPhiEta, h_npt_phiPhi, nbinseta, nbinsphi);
    project_n2XYXY_n2YY(h_ptpt_phiEtaPhiEta, h_ptpt_phiPhi, nbinseta, nbinsphi);
  }

  if (ac.fillYorEta == ac.kRapidity) {
    int nbinsy = h_n2_yY->GetNbinsX();
    int nbinsphi = h_n2_phiPhi->GetNbinsX();

    project_n2XYXY_n2XX(h_n2_phiYPhiY, h_n2_yY, nbinsy, nbinsphi);
    project_n2XYXY_n2YY(h_n2_phiYPhiY, h_n2_phiPhi, nbinsy, nbinsphi);
    project_n2XYXY_n2XX(h_ptn_phiYPhiY, h_ptn_yY, nbinsy, nbinsphi);
    project_n2XYXY_n2XX(h_npt_phiYPhiY, h_npt_yY, nbinsy, nbinsphi);
    project_n2XYXY_n2XX(h_ptpt_phiYPhiY, h_ptpt_yY, nbinsy, nbinsphi);

    project_n2XYXY_n2YY(h_ptn_phiYPhiY, h_ptn_phiPhi, nbinsy, nbinsphi);
    project_n2XYXY_n2YY(h_npt_phiYPhiY, h_npt_phiPhi, nbinsy, nbinsphi);
    project_n2XYXY_n2YY(h_ptpt_phiYPhiY, h_ptpt_phiPhi, nbinsy, nbinsphi);
  }
}

//________________________________________________________________________
void ParticlePairHistos::loadHistograms(TFile* inputFile)
{
  if (!inputFile) {
    cout << "-Fatal- Attempting to load ParticleHistos from an invalid file pointer" << endl;
    return;
  }
  AnalysisConfiguration& ac = *(AnalysisConfiguration*)getConfiguration();
  TString bn = getHistoBaseName();
  h_n2_ptPt = loadH2(inputFile, bn + TString("n2_ptPt"), true);
  h_n2_phiPhi = loadH2(inputFile, bn + TString("n2_phiPhi"), true);
  h_ptn_phiPhi = loadH2(inputFile, bn + TString("ptn_phiPhi"), true);
  h_ptpt_phiPhi = loadH2(inputFile, bn + TString("ptpt_phiPhi"), true);
  h_npt_phiPhi = loadH2(inputFile, bn + TString("npt_phiPhi"), true);

  if (ac.fillYorEta == ac.kPseudorapidity) {
    h_n2_etaEta = loadH2(inputFile, bn + TString("n2_etaEta"), true);
    h_n2_phiEtaPhiEta = loadH2(inputFile, bn + TString("n2_phiEtaPhiEta"), true);
    h_npt_phiEtaPhiEta = loadH2(inputFile, bn + TString("npt_phiEtaPhiEta"), true);
    h_ptn_phiEtaPhiEta = loadH2(inputFile, bn + TString("ptn_phiEtaPhiEta"), true);
    h_ptpt_phiEtaPhiEta = loadH2(inputFile, bn + TString("ptpt_phiEtaPhiEta"), true);
    h_npt_etaEta = loadH2(inputFile, bn + TString("npt_etaEta"), true);
    h_ptn_etaEta = loadH2(inputFile, bn + TString("ptn_etaEta"), true);
    h_ptpt_etaEta = loadH2(inputFile, bn + TString("ptpt_etaEta"), true);
    h_n2_detaDphi_o = loadH2(inputFile, bn + TString("n2_detaDphi_o"), true);
    h_npt_detaDphi_o = loadH2(inputFile, bn + TString("npt_detaDphi_o"), true);
    h_ptn_detaDphi_o = loadH2(inputFile, bn + TString("ptn_detaDphi_o"), true);
    h_ptpt_detaDphi_o = loadH2(inputFile, bn + TString("ptpt_detaDphi_o"), true);
    h_n1n1_detaDphi_o = loadH2(inputFile, bn + TString("n1n1_detaDphi_o"), true);
    h_pt1pt1_detaDphi_o = loadH2(inputFile, bn + TString("pt1pt1_detaDphi_o"), true);
  }

  if (ac.fillYorEta == ac.kRapidity) {
    h_n2_yY = loadH2(inputFile, bn + TString("n2_yY"), true);
    h_n2_phiYPhiY = loadH2(inputFile, bn + TString("n2_phiYPhiY"), true);
    h_npt_phiYPhiY = loadH2(inputFile, bn + TString("npt_phiYPhiY"), true);
    h_ptn_phiYPhiY = loadH2(inputFile, bn + TString("ptn_phiYPhiY"), true);
    h_ptpt_phiYPhiY = loadH2(inputFile, bn + TString("ptpt_phiYPhiY"), true);
    h_npt_yY = loadH2(inputFile, bn + TString("npt_yY"), true);
    h_ptn_yY = loadH2(inputFile, bn + TString("ptn_yY"), true);
    h_ptpt_yY = loadH2(inputFile, bn + TString("ptpt_yY"), true);
  }
  /* the histograms are not owned */
  bOwnTheHistograms = false;
  return;
}

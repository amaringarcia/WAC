// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining ParticlePairHistos
 */

#ifndef WAC_ParticlePairHistos
#define WAC_ParticlePairHistos

#include "Histograms.hpp"

class ParticlePairHistos : public Histograms
{
 public:
  ParticlePairHistos(const TString& name,
                     AnalysisConfiguration* configuration,
                     LogLevel debugLevel);
  ParticlePairHistos(TFile* inputFile,
                     const TString& name,
                     AnalysisConfiguration* configuration,
                     LogLevel debugLevel);
  virtual ~ParticlePairHistos();
  void initialize();
  template <AnalysisConfiguration::RapidityPseudoRapidity r>
  void fill(Particle& particle1, Particle& particle2, double weight);
  void completeFill();
  void loadHistograms(TFile* inputFile);

  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////
  TH2* h_n2_ptPt;

  TH2* h_n2_phiEtaPhiEta;
  TH2* h_npt_phiEtaPhiEta;
  TH2* h_ptn_phiEtaPhiEta;
  TH2* h_ptpt_phiEtaPhiEta;
  TH2* h_n2_ptPhiEtaPtPhiEta;

  TH2* h_n2_detaDphi_o;
  TH2* h_npt_detaDphi_o;
  TH2* h_ptn_detaDphi_o;
  TH2* h_ptpt_detaDphi_o;
  TH2* h_n1n1_detaDphi_o;
  TH2* h_pt1pt1_detaDphi_o;

  TH2* h_n2_etaEta;
  TH2* h_npt_etaEta;
  TH2* h_ptn_etaEta;
  TH2* h_ptpt_etaEta;

  TH2* h_n2_phiPhi;
  TH2* h_ptn_phiPhi;
  TH2* h_npt_phiPhi;
  TH2* h_ptpt_phiPhi;

  TH2* h_n2_yY;
  TH2* h_n2_phiYPhiY;
  TH2* h_npt_phiYPhiY;
  TH2* h_ptn_phiYPhiY;
  TH2* h_ptpt_phiYPhiY;
  TH2* h_n2_ptPhiYPtPhiY;

  TH2* h_npt_yY;
  TH2* h_ptn_yY;
  TH2* h_ptpt_yY;

  TH3* h_n2_Q3D;   // number of pairs vs Q3D
  TH1* h_mInv_Lab; // lab fram calculation
  TH1* h_mInv_PF;  // pair frame
  TH1* h_beta;

  ClassDef(ParticlePairHistos, 1)
};

template <AnalysisConfiguration::RapidityPseudoRapidity r>
void ParticlePairHistos::fill(Particle& particle1, Particle& particle2, double weight)
{
  double pt1 = particle1.pt;
  double phi1 = particle1.phi;
  if (phi1 < 0)
    phi1 += TMath::TwoPi();

  double pt2 = particle2.pt;
  double phi2 = particle2.phi;
  if (phi2 < 0)
    phi2 += TMath::TwoPi();

  AnalysisConfiguration& ac = *(AnalysisConfiguration*)getConfiguration();
  h_n2_ptPt->Fill(pt1, pt2, weight);
  // delayed fill h_n2_etaEta  ->Fill(eta1, eta2, weight);
  // delayed fill h_n2_phiPhi  ->Fill(phi1, phi2, weight);

  // delayed fill h_ptn_etaEta ->Fill(eta1, eta2, weight*pt1);
  // delayed fill h_npt_etaEta ->Fill(eta1, eta2, weight*pt2);
  // delayed fill h_ptpt_etaEta->Fill(eta1, eta2, weight*pt1*pt2);
  // delayed fill h_ptn_phiPhi ->Fill(phi1, phi2, weight*pt1);
  // delayed fill h_npt_phiPhi ->Fill(phi1, phi2, weight*pt2);
  // delayed fill h_ptpt_phiPhi->Fill(phi1, phi2, weight*pt1*pt2);

  int iPhiYEta1 = particle1.ixYEtaPhi;
  int iPhiYEta2 = particle2.ixYEtaPhi;

  bool wrongix = iPhiYEta1 < 0 || iPhiYEta2 < 0;

  if (!wrongix) {
    if constexpr (r == AnalysisConfiguration::kPseudorapidity) {
      double eta1 = particle1.eta;
      double eta2 = particle2.eta;
      double deta = eta1 - eta2;
      double dphi = ac.getDphiShifted(phi1 - phi2);

      int binno = (iPhiYEta2 + 1) * (h_n2_phiEtaPhiEta->GetNbinsX() + 2) + (iPhiYEta1 + 1);
      double nentries = h_n2_phiEtaPhiEta->GetEntries() + 1;
      h_n2_phiEtaPhiEta->AddBinContent(binno, weight);
      h_ptn_phiEtaPhiEta->AddBinContent(binno, weight * pt1);
      h_npt_phiEtaPhiEta->AddBinContent(binno, weight * pt2);
      h_ptpt_phiEtaPhiEta->AddBinContent(binno, weight * pt1 * pt2);
      h_n2_phiEtaPhiEta->SetEntries(nentries);
      h_ptn_phiEtaPhiEta->SetEntries(nentries);
      h_npt_phiEtaPhiEta->SetEntries(nentries);
      h_ptpt_phiEtaPhiEta->SetEntries(nentries);

      binno = h_n2_detaDphi_o->Fill(deta, dphi, weight);
      nentries = h_n2_detaDphi_o->GetEntries();
      h_ptn_detaDphi_o->AddBinContent(binno, weight * pt1);
      h_npt_detaDphi_o->AddBinContent(binno, weight * pt2);
      h_ptpt_detaDphi_o->AddBinContent(binno, weight * pt1 * pt2);
      h_n1n1_detaDphi_o->AddBinContent(binno, weight);
      h_pt1pt1_detaDphi_o->AddBinContent(binno, weight * pt1 * pt2);
      h_ptn_detaDphi_o->SetEntries(nentries);
      h_npt_detaDphi_o->SetEntries(nentries);
      h_ptpt_detaDphi_o->SetEntries(nentries);
      h_n1n1_detaDphi_o->SetEntries(nentries);
      h_pt1pt1_detaDphi_o->SetEntries(nentries);
    }
  }

  if constexpr (r == AnalysisConfiguration::kRapidity) {

    // delayed fill h_n2_yY  ->Fill(y1, y2, weight);
    // delayed fill h_ptn_yY ->Fill(y1, y2, weight*pt1);
    // delayed fill h_npt_yY ->Fill(y1, y2, weight*pt2);
    // delayed fill h_ptpt_yY->Fill(y1, y2, weight*pt1*pt2);

    int binno = (iPhiYEta2 + 1) * (h_n2_phiYPhiY->GetNbinsX() + 2) + (iPhiYEta1 + 1);
    double nentries = h_n2_phiEtaPhiEta->GetEntries() + 1;
    h_n2_phiYPhiY->AddBinContent(binno, weight);
    h_ptn_phiYPhiY->AddBinContent(binno, weight * pt1);
    h_npt_phiYPhiY->AddBinContent(binno, weight * pt2);
    h_ptpt_phiYPhiY->AddBinContent(binno, weight * pt1 * pt2);
    h_n2_phiYPhiY->SetEntries(nentries);
    h_ptn_phiYPhiY->SetEntries(nentries);
    h_npt_phiYPhiY->SetEntries(nentries);
    h_ptpt_phiYPhiY->SetEntries(nentries);
  }
}

#endif /* WAC_ParticlePairHistos  */

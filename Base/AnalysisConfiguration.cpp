// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *************************************************************************/
/**
 \class Property
 \ingroup WAC

 Utility class used to defineAnalysisConfiguration
 */

#include "AnalysisConfiguration.hpp"

ClassImp(AnalysisConfiguration)

  AnalysisConfiguration::AnalysisConfiguration(const TString& name,
                                               const TString& type,
                                               const TString& version)
  : TaskConfiguration(name, type, version),
    histoBaseName(""),
    nBins_nPartTot(10),
    min_nPartTot(0.0),
    max_nPartTot(100.0),
    nBins_nCell(100),
    min_nCell(0.0),
    max_nCell(100.0),
    nBins_nPartPerCell(10),
    min_nPartPerCell(0.0),
    max_nPartPerCell(10.0),
    nBins_pt(40),
    min_pt(0.0),
    max_pt(2.0),
    range_pt(2.0),
    bin_edges_pt{},
    nBins_eta(20),
    min_eta(-2.0),
    max_eta(2.0),
    range_eta(4.0),
    nBins_y(20),
    min_y(-2.0),
    max_y(-2.0),
    range_y(4.0),
    nBins_phi(36),
    min_phi(0.0),
    max_phi(0.0),
    range_phi(0.0),
    nBins_phiEta(720),
    nBins_phiEtaPt(7200),
    nBins_phiY(720),
    nBins_phiYPt(7200),
    binCorrPM(1.0),
    binCorrMP(1.0),
    binCorrMM(1.0),
    binCorrPP(1.0),
    nBins_Qinv(40),
    min_Qinv(0.0),
    max_Qinv(2.0),
    nBins_Qlong(20),
    min_Qlong(-2.0),
    max_Qlong(2.0),
    nBins_Qside(20),
    min_Qside(-2.0),
    max_Qside(2.0),
    nBins_Qout(20),
    min_Qout(-2.0),
    max_Qout(2.0),
    nBins_kT(20),
    min_kT(0.0),
    max_kT(2.0),
    nBins_Dphi(36),
    min_Dphi(0.0),
    max_Dphi(0.0),
    width_Dphi(0.0),
    nBins_Dphi_shft(36),
    min_Dphi_shft(0.0),
    max_Dphi_shft(0.0),
    nBins_Deta(39),
    min_Deta(-2.0),
    max_Deta(2.0),
    width_Deta(4.0 / 39.0),
    nBins_Dy(39),
    min_Dy(-2.0),
    max_Dy(2.0),
    width_Dy(4.0 / 39.0),
    fillPairs(true),
    fillPratt(false),
    fill3D(false),
    fillYorEta(kPseudorapidity),
    nuDynVsMult(true),
    nuDynVsCent(false),
    nBins_mult(20),
    min_mult(0.0),
    max_mult(200.0),
    nBins_cent(20),
    min_cent(0.0),
    max_cent(100.0)
{
  max_phi = TMath::TwoPi();
  range_phi = TMath::TwoPi();
  max_Dphi = max_phi;
  width_Dphi = max_phi / double(nBins_Dphi);
  min_Dphi_shft = -TMath::Pi() / 4.0;
  max_Dphi_shft = 3.0 * TMath::Pi() / 4.0;
}

AnalysisConfiguration::AnalysisConfiguration(const AnalysisConfiguration& source)
  : TaskConfiguration(source),
    histoBaseName(source.histoBaseName),
    nBins_nPartTot(source.nBins_nPartTot),
    min_nPartTot(source.min_nPartTot),
    max_nPartTot(source.max_nPartTot),
    nBins_nCell(source.nBins_nCell),
    min_nCell(source.min_nCell),
    max_nCell(source.max_nCell),
    nBins_nPartPerCell(source.nBins_nPartPerCell),
    min_nPartPerCell(source.min_nPartPerCell),
    max_nPartPerCell(source.max_nPartPerCell),
    nBins_pt(source.nBins_pt),
    min_pt(source.min_pt),
    max_pt(source.max_pt),
    range_pt(source.range_pt),
    bin_edges_pt(source.bin_edges_pt),
    nBins_eta(source.nBins_eta),
    min_eta(source.min_eta),
    max_eta(source.max_eta),
    range_eta(source.range_eta),
    nBins_y(source.nBins_y),
    min_y(source.min_y),
    max_y(source.max_y),
    range_y(source.range_y),
    nBins_phi(source.nBins_phi),
    min_phi(source.min_phi),
    max_phi(source.max_phi),
    range_phi(source.range_phi),
    nBins_phiEta(source.nBins_phiEta),
    nBins_phiEtaPt(source.nBins_phiEtaPt),
    nBins_phiY(source.nBins_phiY),
    nBins_phiYPt(source.nBins_phiYPt),
    binCorrPM(source.binCorrPM),
    binCorrMP(source.binCorrMP),
    binCorrMM(source.binCorrMM),
    binCorrPP(source.binCorrPP),
    nBins_Qinv(source.nBins_Qinv),
    min_Qinv(source.min_Qinv),
    max_Qinv(source.max_Qinv),
    nBins_Qlong(source.nBins_Qlong),
    min_Qlong(source.min_Qlong),
    max_Qlong(source.max_Qlong),
    nBins_Qside(source.nBins_Qside),
    min_Qside(source.min_Qside),
    max_Qside(source.max_Qside),
    nBins_Qout(source.nBins_Qout),
    min_Qout(source.min_Qout),
    max_Qout(source.max_Qout),
    nBins_kT(source.nBins_kT),
    min_kT(source.min_kT),
    max_kT(source.max_kT),
    nBins_Dphi(source.nBins_Dphi),
    min_Dphi(source.min_Dphi),
    max_Dphi(source.max_Dphi),
    width_Dphi(source.width_Dphi),
    nBins_Dphi_shft(source.nBins_Dphi_shft),
    min_Dphi_shft(source.min_Dphi_shft),
    max_Dphi_shft(source.max_Dphi_shft),
    nBins_Deta(source.nBins_Deta),
    min_Deta(source.min_Deta),
    max_Deta(source.max_Deta),
    width_Deta(source.width_Deta),
    nBins_Dy(source.nBins_Dy),
    min_Dy(source.min_Dy),
    max_Dy(source.max_Dy),
    width_Dy(source.width_Dy),
    fillPairs(source.fillPairs),
    fillPratt(source.fillPratt),
    fill3D(source.fill3D),
    fillYorEta(source.fillYorEta),
    nuDynVsMult(source.nuDynVsMult),
    nuDynVsCent(source.nuDynVsCent),
    nBins_mult(source.nBins_mult),
    min_mult(source.min_mult),
    max_mult(source.max_mult),
    nBins_cent(source.nBins_cent),
    min_cent(source.min_cent),
    max_cent(source.max_cent)
{
}

AnalysisConfiguration& AnalysisConfiguration::operator=(const AnalysisConfiguration& source)
{
  TaskConfiguration::operator=(source);
  histoBaseName = source.histoBaseName;
  nBins_nPartTot = source.nBins_nPartTot;
  min_nPartTot = source.min_nPartTot;
  max_nPartTot = source.max_nPartTot;
  nBins_nCell = source.nBins_nCell;
  min_nCell = source.min_nCell;
  max_nCell = source.max_nCell;
  nBins_nPartPerCell = source.nBins_nPartPerCell;
  min_nPartPerCell = source.min_nPartPerCell;
  max_nPartPerCell = source.max_nPartPerCell;
  nBins_pt = source.nBins_pt;
  min_pt = source.min_pt;
  max_pt = source.max_pt;
  range_pt = source.range_pt;
  bin_edges_pt = source.bin_edges_pt;
  nBins_eta = source.nBins_eta;
  min_eta = source.min_eta;
  max_eta = source.max_eta;
  range_eta = source.range_eta;
  nBins_y = source.nBins_y;
  min_y = source.min_y;
  max_y = source.max_y;
  range_y = source.range_y;
  nBins_phi = source.nBins_phi;
  min_phi = source.min_phi;
  max_phi = source.max_phi;
  range_phi = source.range_phi;
  nBins_phiEta = source.nBins_phiEta;
  nBins_phiEtaPt = source.nBins_phiEtaPt;
  nBins_phiY = source.nBins_phiY;
  nBins_phiYPt = source.nBins_phiYPt;
  binCorrPM = source.binCorrPM;
  binCorrMM = source.binCorrMM;
  binCorrPP = source.binCorrPP;
  nBins_Qinv = source.nBins_Qinv;
  min_Qinv = source.min_Qinv;
  max_Qinv = source.max_Qinv;
  nBins_Qlong = source.nBins_Qlong;
  min_Qlong = source.min_Qlong;
  max_Qlong = source.max_Qlong;
  nBins_Qside = source.nBins_Qside;
  min_Qside = source.min_Qside;
  max_Qside = source.max_Qside;
  nBins_Qout = source.nBins_Qout;
  min_Qout = source.min_Qout;
  max_Qout = source.max_Qout;
  nBins_kT = source.nBins_kT;
  min_kT = source.min_kT;
  max_kT = source.max_kT;
  nBins_Dphi = source.nBins_Dphi;
  min_Dphi = source.min_Dphi;
  max_Dphi = source.max_Dphi;
  width_Dphi = source.width_Dphi;
  nBins_Dphi_shft = source.nBins_Dphi_shft;
  min_Dphi_shft = source.min_Dphi_shft;
  max_Dphi_shft = source.max_Dphi_shft;
  nBins_Deta = source.nBins_Deta;
  min_Deta = source.min_Deta;
  max_Deta = source.max_Deta;
  width_Deta = source.width_Deta;
  nBins_Dy = source.nBins_Dy;
  min_Dy = source.min_Dy;
  max_Dy = source.max_Dy;
  width_Dy = source.width_Dy;
  fillPairs = source.fillPairs;
  fillPratt = source.fillPratt;
  fill3D = source.fill3D;
  fillYorEta = source.fillYorEta;
  nuDynVsMult = source.nuDynVsMult;
  nuDynVsCent = source.nuDynVsCent;
  nBins_mult = source.nBins_mult;
  min_mult = source.min_mult;
  max_mult = source.max_mult;
  nBins_cent = source.nBins_cent;
  min_cent = source.min_cent;
  max_cent = source.max_cent;

  return *this;
}

////////////////////////////////////////////////////
// Print this configuration to the given stream
////////////////////////////////////////////////////
void AnalysisConfiguration::printConfiguration(ostream& os)
{
  printTaskConfiguration(os);
  os
    << "    Analysis   Parameters: " << endl
    << " ------------------------------------------------------------------------------------------" << endl
    << "      Histogram Base Name: " << histoBaseName << endl
    << "                  nBinsPt: " << nBins_pt << endl
    << "                    minPt: " << min_pt << endl
    << "                    maxPt: " << max_pt << endl
    << "                 nBinsEta: " << nBins_eta << endl
    << "                   minEta: " << min_eta << endl
    << "                   maxEta: " << max_eta << endl
    << "                   nBinsY: " << nBins_y << endl
    << "                     minY: " << min_y << endl
    << "                     maxY: " << max_y << endl
    << "                 nBinsPhi: " << nBins_phi << endl
    << "                   minPhi: " << min_phi << endl
    << "                   maxPhi: " << max_phi << endl
    << "                fillPairs: " << fillPairs << endl
    << "                fillPratt: " << fillPratt << endl
    << "                   fill3D: " << fill3D << endl
    << "             fillY or eta: " << ((fillYorEta == kRapidity) ? "Y" : "eta") << endl
    << "        +- bin correction: " << binCorrPM << endl
    << "        -+ bin correction: " << binCorrMP << endl
    << "        ++ bin correction: " << binCorrPP << endl
    << "        -- bin correction: " << binCorrMM << endl;
  if (bin_edges_pt.size() > 0) {
    os
      << "         -- Pt bins edges: ";
    int edgeno = 0;
    for (double edge : bin_edges_pt) {
      if ((edgeno + 1) % 8 == 0) {
        os << endl
           << "                           ";
      }
      os << edge << ", ";
      edgeno++;
    }
  }
}

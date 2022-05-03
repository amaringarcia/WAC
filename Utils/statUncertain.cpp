//
//  statUncertain.C
//  CodeV3
//
//  Created by Claude Pruneau on 8/22/17.
//  Copyright © 2017 Claude Pruneau. All rights reserved.
//

#include <stdio.h>
#include <cstring>

//
//  cm_figure9_6.c
//  MyMC
//
//  Created by Claude Pruneau on 4/10/17.
//  Copyright © 2017 Claude Pruneau. All rights reserved.
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <fstream>
#include <TError.h>
#include <TObjArray.h>
#include <TTimeStamp.h>
#include <TList.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TMath.h>
#include "TRandom.h"
#include "../Base/AnalysisConfiguration.hpp"
#include "../Base/TwoPartDiffCorrelationAnalyzer.hpp"

const int npart = 6;
const int ncomb = 2;
const char* partname[npart] = {"PiP", "PiM", "KaP", "KaM", "PrP", "PrM"};
const char* combname[ncomb] = {"CI", "CD"};
const int ncorrpart = 3;
const int ncorrfcomb = 4;
const char* corrfname[ncorrfcomb] = {"P2", "R2", "G2", "R2BF"};

std::vector<std::string> centfname = {
  "MB"};
int ncent = centfname.size();

std::vector<std::string> samplesfnames = {
  "BUNCH01/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH02/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH03/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH04/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH05/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH06/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH07/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH08/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH09/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH10/Output/PYTHIA8_PiKaPr_%s.root",
};
int nsamples = samplesfnames.size();

TFile *getSampleFile(int icent,int isample) {
  TFile *f = nullptr;

  std::string filename = TString::Format(samplesfnames[isample].c_str(),centfname[icent].c_str(),centfname[icent].c_str()).Data();
  f = new TFile(filename.c_str());
  if (f == nullptr or not f->IsOpen()) {
    Error("statUncertain::getSampleFile","File %s not found. ABORTING!!!", filename.c_str());
  }
  return f;
}


TH2 *extractHistoMeanAndStDevFromSubSets(const TObjArray &listsarray, Int_t ih, const TString &name) {

  TH2 *h2 = dynamic_cast<TH2*>(((TList*) listsarray[0])->At(ih));

  if (h2 != NULL) {
    /* extract the name of the result histogram */
    TString title = h2->GetTitle();
    TString xtitle = h2->GetXaxis()->GetTitle();
    TString ytitle = h2->GetYaxis()->GetTitle();
    TString ztitle = h2->GetZaxis()->GetTitle();

    /* let's first create the profile we use as support */
    TProfile2D *p = new TProfile2D("auxprofile","auxprofile",
        h2->GetNbinsX(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(),
        h2->GetNbinsY(), h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax());

    for (Int_t ix = 0; ix < h2->GetNbinsX(); ix++) {
      Double_t x = h2->GetXaxis()->GetBinCenter(ix+1);
      for (Int_t iy = 0; iy < h2->GetNbinsY(); iy++) {
        Double_t y = h2->GetYaxis()->GetBinCenter(iy+1);
        for (Int_t iset = 0; iset < listsarray.GetEntries(); iset++) {
          p->Fill(x,y,((TH2*)(((TList*) listsarray[iset])->At(ih)))->GetBinContent(ix+1,iy+1));
        }
      }
    }
    h2 = p->ProjectionXY(name.Data());
    h2->SetTitle(title.Data());
    h2->SetXTitle(xtitle.Data());
    h2->SetYTitle(ytitle.Data());
    h2->SetZTitle(ztitle.Data());

    /* set the proper number of entries */
    Double_t nentries = 0;
    for (Int_t iset = 0; iset < listsarray.GetEntries(); iset++) {
      nentries += ((TH2*)(((TList*) listsarray[iset])->At(ih)))->GetEntries();
    }
    h2->SetEntries(nentries);

    delete p;
    return h2;
  }
  else {
    Error("extractHistoMeanAndStDevFromSubSets", "Wrong histogram type");
    return NULL;
  }
}

TList* extractMeanAndStDevFromSubSets(const TObjArray& listsarray, const TString& pattern, char** name)
{
  /* basically we receive an array of histograms lists */
  /* each array item corresponds to a results subset   */
  /* we extract an equivalent list of histograms which */
  /* corresponds to the mean of the subsets histograms */
  /* and have as errors the standard deviation from    */
  /* the mean on each bin                              */

  /* first, some consistency checks */
  Int_t nhistos = ((TList*) listsarray[0])->GetEntries();
  if (nhistos != ncorrpart && nhistos != ncorrfcomb * ncomb)
    Error("extractMeanAndStDevFromSubSets", "Inconsistent number of histograms to average");
  for (Int_t iset = 0; iset < listsarray.GetEntries(); iset++) {
    if (nhistos != ((TList*) listsarray[iset])->GetEntries()) {
      Error("extractMeanAndStDevFromSubSets", "Inconsistent array of lists");
      return NULL;
    }
  }

  TList *list = new TList();
  list->SetOwner(kTRUE);
  for (Int_t ih = 0; ih < nhistos; ih++) {
    list->Add(extractHistoMeanAndStDevFromSubSets(listsarray,ih,TString::Format(pattern.Data(),name[ih])));
  }

  return list;
}

TList *extractSampleResults(AnalysisConfiguration *ac, int icent,int isample) {

  /* we will control what goes to the directory tree */
  Bool_t oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  std::vector<ParticleFilter*> particleFilters;
  particleFilters.push_back(new ParticleFilter(ParticleFilter::Pion, ParticleFilter::Positive, ac->min_pt, ac->max_pt, ac->min_eta, ac->max_eta, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter(ParticleFilter::Pion, ParticleFilter::Negative, ac->min_pt, ac->max_pt, ac->min_eta, ac->max_eta, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter(ParticleFilter::Kaon, ParticleFilter::Positive, ac->min_pt, ac->max_pt, ac->min_eta, ac->max_eta, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter(ParticleFilter::Kaon, ParticleFilter::Negative, ac->min_pt, ac->max_pt, ac->min_eta, ac->max_eta, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Positive, ac->min_pt, ac->max_pt, ac->min_eta, ac->max_eta, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Negative, ac->min_pt, ac->max_pt, ac->min_eta, ac->max_eta, ac->min_y, ac->max_y));

  EventFilter* eventFilter = new EventFilter(EventFilter::MinBias, 0, 0);

  Event *event = Event::getEvent();

  TwoPartDiffCorrelationAnalyzer* eventanalyzer = new TwoPartDiffCorrelationAnalyzer("NarrowPiKaP", ac, event, eventFilter, particleFilters);

  eventanalyzer->setReportLevel(MessageLogger::Error);

  TFile *myfile = getSampleFile(icent,isample);
  if (myfile == nullptr) return nullptr;

  eventanalyzer->loadBaseHistograms(myfile);
  eventanalyzer->finalize();

  TList *list = new TList();
  list->SetOwner(kFALSE);

  /* the pair single histos */
  for (int ipart = 0; ipart < npart; ++ipart) {   /* first component of the pair */
    for (int jpart = 0; jpart < npart; ++jpart) { /* second component of the pair */
      TList* plist = new TList();                 /* a list per pair */
      plist->SetOwner(kTRUE);
      for (int icf = 0; icf < ncorrpart; ++icf) { /* the individual pair cf */
        TH2* h2 = NULL;
        switch (icf) {
          case 0: /* P2 */
            h2 = eventanalyzer->pairs_Histos[ipart][jpart]->h_P2_DetaDphi_shft;
            break;
          case 1: /* R2 */
            h2 = eventanalyzer->pairs_Histos[ipart][jpart]->h_R2_DetaDphi_shft;
            break;
          case 2: /* G2 */
            h2 = eventanalyzer->pairs_Histos[ipart][jpart]->h_G2_DetaDphi_shft;
            break;
          default:
            Fatal("extractSampleResults", "Wrong correlator index");
        }
        plist->Add(h2->Clone(TString::Format("%s%s_Sub%02d", h2->GetName(), centfname[icent].c_str(), isample)));
      }
      list->Add(plist);
    }
  }

  /* the pair combined histos */
  for (int ipart = 0; ipart < npart; ++ipart) {                 /* first component of the pair */
    for (int jpart = 0; jpart < npart - (ipart + 1); ++jpart) { /* second component of the pair */
      TList* plist = new TList();                               /* a list per pair */
      plist->SetOwner(kTRUE);
      for (int icf = 0; icf < ncorrfcomb; ++icf) { /* the individual pair cf */
        TH2* h2ci = nullptr;
        TH2* h2cd = nullptr;
        switch (icf) {
          case 0: /* P2 */
            h2ci = eventanalyzer->pairs_CIHistos[ipart][jpart]->h_P2_DetaDphi_shft;
            h2cd = eventanalyzer->pairs_CDHistos[ipart][jpart]->h_P2_DetaDphi_shft;
            break;
          case 1: /* R2 */
            h2ci = eventanalyzer->pairs_CIHistos[ipart][jpart]->h_R2_DetaDphi_shft;
            h2cd = eventanalyzer->pairs_CDHistos[ipart][jpart]->h_R2_DetaDphi_shft;
            break;
          case 2: /* G2 */
            h2ci = eventanalyzer->pairs_CIHistos[ipart][jpart]->h_G2_DetaDphi_shft;
            h2cd = eventanalyzer->pairs_CDHistos[ipart][jpart]->h_G2_DetaDphi_shft;
            break;
          case 3: /* R2BF */
            h2ci = eventanalyzer->pairs_CIHistos[ipart][jpart]->h_R2BF_DetaDphi_shft;
            h2cd = eventanalyzer->pairs_CDHistos[ipart][jpart]->h_R2BF_DetaDphi_shft;
            break;
          default:
            Fatal("extractSampleResults", "Wrong correlator index");
        }
        plist->Add(h2ci->Clone(TString::Format("%s%s_Sub%02d", h2ci->GetName(), centfname[icent].c_str(), isample)));
        plist->Add(h2cd->Clone(TString::Format("%s%s_Sub%02d", h2cd->GetName(), centfname[icent].c_str(), isample)));
      }
      list->Add(plist);
    }
  }

  myfile->Close();

  delete myfile;
  delete eventanalyzer;
  delete eventFilter;
  for (auto pf : particleFilters) {
    delete pf;
  }

  /* back to normal */
  TH1::AddDirectory(oldstatus);

  return list;
}

/////////////////////////////////////////////////////////////////////////////////////////
// produce results with statistic uncertainties out of a set of sub-samples
/////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc > 2 or argc < 1) {
    Fatal("main", "Wrong number of arguments. Use statUncertain option");
  }

  Option_t* opt = argv[1];

  TTimeStamp now = TTimeStamp();
  if (!TString(opt).Contains("verb"))
    gErrorIgnoreLevel = kWarning;

  /* Cummulate errors by default */
  TH1::SetDefaultSumw2(kTRUE);

  TFile* outputfile = new TFile(TString::Format("Histograms%02dSub_%d_%d.root", nsamples, now.GetDate(), now.GetTime()), "RECREATE");

  for (int icent = 0; icent < ncent; icent++) {

    Info("statUncertain", "Setup configuration");

    AnalysisConfiguration* ac = new AnalysisConfiguration("PYTHIA8", "UCM-miniWAC", "1.0");

    ac->outputPath = "./";
    ac->configurationFileName = "Config_DUKE.txt";
    ac->rootOuputFileName = TString::Format("Histograms%02dSub", nsamples).Data();

    cout << "================================================================" << endl;
    cout << "" << endl;
    cout << "         running for = " << ac->getName() << endl;
    cout << "" << endl;
    cout << "================================================================" << endl;

    // =========================
    // Short configuration
    // =========================
    float min_eta = -1;
    float max_eta = 1;
    float min_pt = 0.2;
    float max_pt = 3.0;
    int nBins_eta = int((max_eta - min_eta) / 0.1);
    int nBins_pt = int((max_pt - min_pt) / 0.1);

    ac->nBins_pt = nBins_pt;
    ac->min_pt = min_pt;
    ac->max_pt = max_pt;
    ac->nBins_eta = nBins_eta;
    ac->min_eta = min_eta;
    ac->max_eta = max_eta;
    ac->nBins_y = nBins_eta;
    ac->min_y = min_eta;
    ac->max_y = max_eta;
    ac->nBins_phi = 72;
    ac->min_phi = 0.0;
    ac->max_phi = kTWOPI;

    ac->fillPairs = true;
    ac->fill3D = false;
    ac->fill6D = false;
    ac->fillQ3D = false;
    ac->fillY = false;

    ac->scaleHistograms = true;
    ac->createHistograms = false;
    ac->loadHistograms = true;
    ac->saveHistograms = false;
    ac->forceHistogramsRewrite = false;
    ac->calculateDerivedHistograms = true;

    /* the pair single histos plus the pair combined histos */
    int nmainlists = npart * npart + npart * (npart - 1) / 2;
    std::vector<TObjArray> pairslists(nmainlists, TObjArray(nsamples));
    for (int ilst = 0; ilst < nmainlists; ++ilst) {
      pairslists[ilst].SetOwner(kTRUE);
    }

    for (Int_t isamp = 0; isamp < nsamples; isamp++) {
      Warning("statUncertain", "Processing sample %d for centrality %s", isamp, centfname[icent].c_str());
      TList* list = extractSampleResults(ac, icent, isamp);

      for (Int_t ilst = 0; ilst < nmainlists; ++ilst) {
        pairslists[ilst][isamp] = list->At(ilst);
      }
      /* we write the individual results if required */
      if (TString(opt).Contains("savesub")) {
        outputfile->cd();
        for (Int_t ixh = 0; ixh < list->GetEntries(); ixh++) {
          list->At(ixh)->Write();
        }
      }
      delete list;
    }

    /* we will control what goes to the directory tree */
    Bool_t oldstatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    char* cfname[ncorrfcomb] = {nullptr};
    for (int i = 0; i < ncorrfcomb; ++i) {
      cfname[i] = new char[std::strlen(corrfname[i]) + 1];
      sprintf(cfname[i], "%s", corrfname[i]);
    }
    outputfile->cd();
    for (int ipart = 0; ipart < npart; ++ipart) {
      for (int jpart = 0; jpart < npart; ++jpart) {
        int ilst = ipart * npart + jpart;
        TString pattern = TString::Format("Pythia8_%s%s%%s_DetaDphi_shft_%s", partname[ipart], partname[jpart], centfname[icent].c_str());
        TList* meanhlist = extractMeanAndStDevFromSubSets(pairslists[ilst], pattern, cfname);
        for (int ixh = 0; ixh < meanhlist->GetEntries(); ixh++) {
          meanhlist->At(ixh)->Write();
        }
        delete meanhlist;
      }
    }
    char* cfnamecomb[ncorrfcomb * ncomb] = {nullptr};
    for (int i = 0; i < ncorrfcomb; ++i) {
      for (int j = 0; j < ncomb; ++j) {
        cfnamecomb[i * ncomb + j] = new char[std::strlen(corrfname[i]) + strlen(combname[j]) + 1];
        sprintf(cfnamecomb[i * ncomb + j], "%s%s", corrfname[i], combname[j]);
      }
    }
    int ilst = npart * npart;
    for (int ipart = 0; ipart < npart; ++ipart) {
      for (int jpart = 0; jpart < npart - (ipart + 1); ++jpart) {
        TString pattern = TString::Format("Pythia8_%s%s%%s_DetaDphi_shft_%s", partname[ipart], partname[ipart + 1 + jpart], centfname[icent].c_str());
        TList* meanhlist = extractMeanAndStDevFromSubSets(pairslists[ilst], pattern, cfnamecomb);
        for (int ixh = 0; ixh < meanhlist->GetEntries(); ixh++) {
          meanhlist->At(ixh)->Write();
        }
        delete meanhlist;
        ilst++;
      }
    }

    for (auto aname : cfname) {
      delete aname;
    }
    for (auto aname : cfnamecomb) {
      delete aname;
    }
    /* back to normal */
    TH1::AddDirectory(oldstatus);
  }
  outputfile->Close();
  delete outputfile;
}

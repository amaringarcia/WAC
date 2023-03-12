//
//  HistogramCollection.hpp
//  MyMC
//
//  Created by Claude Pruneau on 12/11/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//
#ifndef WAC_HistogramCollection
#define WAC_HistogramCollection
#include <stdio.h>
#include "TClass.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TList.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TPad.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom.h"
#include "TString.h"
#include "TLatex.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLegend.h"
#include "MessageLogger.hpp"

/// using namespace std;

class HistogramCollection : public MessageLogger
{
 public:
  ////////////////////////////////////////////////////////////////////////////
  // CTOR1
  ////////////////////////////////////////////////////////////////////////////
  HistogramCollection(const TString& name,
                      int histoCapacity = 100,
                      LogLevel debugStatus = Info);
  ////////////////////////////////////////////////////////////////////////////
  virtual ~HistogramCollection();

  ////////////////////////////////////////////////////////////////////////////
  // Add the given histogram to the list
  ////////////////////////////////////////////////////////////////////////////
  void addToList(TH1* h,
                 bool doScale = false,
                 bool doSave = true,
                 bool doPlot = true,
                 bool doPrint = true,
                 bool doSumw2 = false);
  TH1* getHisto(int i);
  void setDefaultOptions(bool color = 0);
  TH1* createHistogram(const TString& name,
                       int n, double min_x, double max_x,
                       const TString& title_x,
                       const TString& title_y,
                       bool scale = true,
                       bool save = true,
                       bool plot = true,
                       bool print = true,
                       bool sumw2 = true);
  TH1* createHistogram(const TString& name,
                       int n, double* bins,
                       const TString& title_x,
                       const TString& title_y,
                       bool scale = true,
                       bool save = true,
                       bool plot = true,
                       bool print = true,
                       bool sumw2 = true);
  TH1* createHistogram(const TString& name,
                       std::vector<double> edges,
                       const TString& title_x,
                       const TString& title_y,
                       bool scale = true,
                       bool save = true,
                       bool plot = true,
                       bool print = true,
                       bool sumw2 = true);
  TH2* createHistogram(const TString& name,
                       int n_x, double min_x, double max_x,
                       int n_y, double min_y, double max_y,
                       const TString& title_x,
                       const TString& title_y,
                       const TString& title_z,
                       bool scale = true,
                       bool save = true,
                       bool plot = true,
                       bool print = true,
                       bool sumw2 = true);
  TH2* createHistogram(const TString& name,
                       int n_x, double* xbins, int n_y, double min_y, double max_y,
                       const TString& title_x,
                       const TString& title_y,
                       const TString& title_z,
                       bool scale = true,
                       bool save = true,
                       bool plot = true,
                       bool print = true,
                       bool sumw2 = true);
  TH2* createHistogram(const TString& name,
                       std::vector<double> xedges, int n_y, double min_y, double max_y,
                       const TString& title_x,
                       const TString& title_y,
                       const TString& title_z,
                       bool scale = true,
                       bool save = true,
                       bool plot = true,
                       bool print = true,
                       bool sumw2 = true);
  TH2* createHistogram(const TString& name,
                       int n_x, double min_x, double max_x, std::vector<double> yedges,
                       const TString& title_x,
                       const TString& title_y,
                       const TString& title_z,
                       bool scale = true,
                       bool save = true,
                       bool plot = true,
                       bool print = true,
                       bool sumw2 = true);
  TH2* createHistogram(const TString& name,
                       std::vector<double> xedges,
                       std::vector<double> yedges,
                       const TString& title_x,
                       const TString& title_y,
                       const TString& title_z,
                       bool scale = true,
                       bool save = true,
                       bool plot = true,
                       bool print = true,
                       bool sumw2 = true);
  TH3* createHistogram(const TString& name,
                       int n_x, double min_x, double max_x,
                       int n_y, double min_y, double max_y,
                       int n_z, double min_z, double max_z,
                       const TString& title_x,
                       const TString& title_y,
                       const TString& title_z,
                       const TString& title_w,
                       bool scale = true,
                       bool save = true,
                       bool plot = true,
                       bool print = true,
                       bool sumw2 = true);
  TProfile* createProfile(const TString& name,
                          int n_x, double min_x, double max_x,
                          const TString& title_x,
                          const TString& title_y,
                          bool save = true,
                          bool plot = true,
                          bool print = true);

  TProfile* createProfile(const TString& name,
                          int n_x, double* bins,
                          const TString& title_x,
                          const TString& title_y,
                          bool scale = true,
                          bool save = true,
                          bool plot = true,
                          bool print = true);

  TProfile2D* createProfile(const TString& title,
                            int n_x, double min_x, double max_x,
                            int n_y, double min_y, double max_y,
                            const TString& title_x,
                            const TString& title_y,
                            const TString& title_z,
                            bool scale = true,
                            bool save = true,
                            bool plot = true,
                            bool print = true,
                            bool sumw2 = true);

  void addHistogramsToExtList(TList* list, bool all = false);
  void saveHistograms(TDirectory* dir, bool saveAll = false);
  void scale(double factor, bool scaleAll = false);
  void addHistos(HistogramCollection* c1, double a1);
  void addHistos(HistogramCollection* c1, HistogramCollection* c2, double a1, double a2);
  void addHistos(HistogramCollection* c1, HistogramCollection* c2, HistogramCollection* c3, double a1, double a2, double a3);
  void addHistos(HistogramCollection* c1, HistogramCollection* c2, HistogramCollection* c3, HistogramCollection* c4,
                 double a1, double a2, double a3, double a4);
  double calculateN1N1(const TH1* h_1, const TH1* h_2, TH1* h_12, double a1, double a2);
  double calculateN1N1_H1H1H2(const TH1* h_1, const TH1* h_2, TH2* h_12, double a1, double a2);
  double calculateN1N1_H2H2H2(const TH2* h_1, const TH2* h_2, TH2* h_12, double a1, double a2);
  double calculateN1N1N1(const TH1* h_1, const TH1* h_2, const TH1* h_3, TH1* h_123);
  double calculateN1N1N1(const TH1* h_1, const TH1* h_2, const TH1* h_3, TH3* h_123);
  void calculateN2N1(const TH2* s2, const TH1* s1, TH2* target, int single);
  void calculateN2N1x(const TH2* s2, const TH1* s1, TH3* target, int single);
  bool sameDimensions1D(const TH1* h1, const TH1* h2);
  bool sameDimensions2D(const TH2* h1, const TH2* h2);
  bool sameDimensions3D(const TH3* h1, const TH3* h2);
  void calculateDptDpt(const TH2* spp, const TH2* spn, const TH2* snp, const TH2* snn,
                       const TH2* avgp1, const TH2* avgp2, TH2* s2dptdpt, TH2* dptdpt,
                       bool ijNormalization, int nEta, int nPhi);
  void calculateDptDpt(const TH2* spp, const TH2* spn, const TH2* snp, const TH2* snn,
                       const TH1* avgp1, const TH1* avgp2,
                       TH2* dptdpt,
                       bool ijNormalization, int nBins);
  void calculateSc(const TH1* spp, const TH1* n1n1, const TH1* pt1pt1, TH1* sean, bool ijNormalization);
  void calculateG2_H2H2H2H2(const TH2* spp, const TH2* n1n1, const TH2* pt1pt1, TH2* sean, double avgpt1_1, double avgpt1_2, bool ijNormalization, double a1, double a2);
  void calculateBf(const TH2* n2, const TH2* n1_1, const TH2* n1_2, TH2* bf_12, TH2* bf_21);
  void calculateBfR2(const TH2* r2, const TH2* n1_1, const TH2* n1_2, TH2* bf_12, TH2* bf_21);
  void calculateSean_H1H2H2H2(const TH1* spp, const TH2* n1n1, const TH2* pt1pt1, TH2* sean, bool ijNormalization, double a1, double a2);
  int calculateQ3DwPtPhiEta(double pt1, double phi1, double eta1,
                            double pt2, double phi2, double eta2,
                            double& Qlong, double& Qout, double& Qside);
  int calculateQ3DwPtPhiY(double pt1, double phi1, double y1,
                          double pt2, double phi2, double y2,
                          double& Qlong, double& Qout, double& Qside);
  void calculateN1N1H2H2_Q3D_MCY(TH2* n1_1, TH2* n1_2, TH3* n1n1_Q3D, double a1, double a2);
  void calculateN1N1H2H2_Q3D_MCEta(TH2* n1_1, TH2* n1_2, TH3* n1n1_Q3D, double a1, double a2);
  void calculateN1N1H2H2_Q3D(const TH2* n1_1, const TH2* n1_2, TH3* n1n1_Q3D, double a1, double a2);
  void calculateN1N1H3H3_Q3D(const TH3* n1_1, const TH3* n1_2, TH3* n1n1_Q3D, double a1, double a2);
  void calculateR2_Q3D(const TH3* n2_Q3D, const TH3* n1n1_Q3D, TH3* R2_Q3D, double a1, double a2);
  double avgValue(TH1* h);

  void setHistogram(TH1* h, double v, double ev);
  void setHistogram(TH2* h, double v, double ev);
  TH1* loadH1(TDirectory* dir,
              const TString& histoName,
              bool scale = false,
              bool save = false,
              bool plot = true,
              bool print = true,
              bool sumw2 = false);
  TH2* loadH2(TDirectory* dir,
              const TString& histoName,
              bool scale = false,
              bool save = false,
              bool plot = true,
              bool print = true,
              bool sumw2 = false);
  TH3* loadH3(TDirectory* dir,
              const TString& histoName,
              bool scale = false,
              bool save = false,
              bool plot = true,
              bool print = true,
              bool sumw2 = false);
  TProfile* loadProfile(TDirectory* dir,
                        const TString& histoName,
                        bool scale = false,
                        bool save = false,
                        bool plot = true,
                        bool print = true,
                        bool sumw2 = false);
  TProfile2D* loadProfile2D(TDirectory* dir,
                            const TString& histoName,
                            bool scale = false,
                            bool save = false,
                            bool plot = true,
                            bool print = true,
                            bool sumw2 = false);
  TH1* cloneH1(const TH1* h1, const TString& histoName);
  TH2* cloneH2(const TH2* h2, const TString& histoName);
  TH3* cloneH3(const TH3* h3, const TString& histoName);
  void findMaximum(TH1* h, int xFirstBin, int xLastBin, int& xMaxValueBin, double& xMaxValue);
  void findMinimum(TH1* h, int xFirstBin, int xLastBin, int& xMinValueBin, double& xMinValue);
  void scaleByBinWidth1D(TH1* h, double scale);
  void scaleByBinWidth2D(TH2* h, double scale);
  void scaleByBinWidth3D(TH3* h, double scale);
  void scaleByBinWidth(TH1* h, double scale);
  void scaleAllHistoByBinWidth(double scale);
  void sumw2All();
  void unpack_vsXY_to_vsXVsY(const TH1* source, TH2* target);
  void correctMerging(TH1* h, int nEta, int nPhi, bool reverse);
  void calculateR2_H1H1H1(const TH1* n2_12, const TH1* n1n1_12, TH1* r2_12, bool ijNormalization, double a1, double a2);
  void calculateR2_H2H2H2(const TH2* n2_12, const TH2* n1n1_12, TH2* r2_12, bool ijNormalization, double a1, double a2);
  void calculateR2_H1H2H2(const TH1* n2_12, const TH2* n1n1_12, TH2* r2_12, bool ijNormalization, double a1, double a2);
  void calculateR2VsM(const TProfile* h1, const TProfile* h2, const TProfile* h12, TH1* r2VsM, TH1* intR2, bool sameFilter);
  void calculateBinCorr(const TProfile* h1, const TProfile* h2, TH1* intBinCorrVsM1, bool sameFilter);
  void calculateAveragePt(const TH1* h1, const TH1* h2, TH1* h3);
  void calculateAveragePtH2(const TH2* h1, const TH2* h2, TH2* h3);
  void calculateAveragePtH1H2(const TH1* h1, const TH1* h2, TH2* h3);
  void calculateAveragePt(const TProfile* h1, const TProfile* h2, TH1* h3);
  void symmetrize3D(TH3* h);
  void symmetrizeDeltaEtaDeltaPhi(TH2* h, bool ijNormalization);
  void symmetrizeXX(TH2* h, bool ijNormalization);
  void reduce_n2xEtaPhi_n2DetaDphi(const TH2* source, TH2* target, int nEtaBins, int nPhiBins, Option_t* opt = "");
  void reduce_n2xEtaPhi_n2EtaEta(const TH1* source, TH2* target, int nEtaBins, int nPhiBins);
  void convolve_n2xEtaPhi_n2DetaDphi(const TH2* source1, const TH2* source2, TH2* target, double a1, double a2);
  void project_n2XYXY_n2XX(const TH2* source, TH2* target, int nXBins, int nYBins);
  void project_n2XYXY_n2YY(const TH2* source, TH2* target, int nXBins, int nYBins);
  TH2* symmetrize(TH2* h);
  TH2* transpose(const TH2* h);
  void shiftY(const TH2& source, TH2& target, int nbins);
  float* getFloatArray(int size, float v);
  double* getDoubleArray(int size, double v);
  void resetDoubleArray(int n, double* array, double value);

  int index2(int i1, int i2);
  int index3(int i1, int i2, int i3);
  int index4(int i1, int i2, int i3, int i4);
  void calculateF2R2(TH1* h_f1_1, TH1* h_f1_2, TH1* h_f2_12, TH1* h_F2_12, TH1* h_R2_12);
  void calculateNudyn(TH1* h_R2_11, TH1* h_R2_12, TH1* h_R2_22, TH1* h_nudyn_12);
  void calculateF3R3(TH1* h_f1_1, TH1* h_f1_2, TH1* h_f1_3,
                     TH1* h_f2_12, TH1* h_f2_13, TH1* h_f2_23,
                     TH1* h_f3_123,
                     TH1* h_F3_123, TH1* h_R2_123);
  void calculateF4R4(TH1* h_f1_1, TH1* h_f1_2, TH1* h_f1_3, TH1* h_f1_4,
                     TH1* h_f2_12, TH1* h_f2_13, TH1* h_f2_14, TH1* h_f2_23, TH1* h_f2_24, TH1* h_f2_34,
                     TH1* h_f3_123, TH1* h_f3_124, TH1* h_f3_134, TH1* h_f3_234,
                     TH1* h_f4_1234,
                     TH1* h_F4_1234, TH1* h_R4_1234);
  void calculateF2R2(double f1_1, double ef1_1, double f1_2, double ef1_2, double f2_12, double ef2_12, double& F2_12, double& eF2_12, double& R2_12, double& eR2_12);
  void calculateF3R3(double f1_1, double ef1_1, double f1_2, double ef1_2, double f1_3, double ef1_3,
                     double f2_12, double ef2_12, double f2_13, double ef2_13, double f2_23, double ef2_23,
                     double f3_123, double ef3_123,
                     double& F3_123, double& eF3_123, double& R3_123, double& eR3_123);
  void calculateF4R4(double f1_1, double ef1_1, double f1_2, double ef1_2, double f1_3, double ef1_3, double f1_4, double ef1_4,
                     double f2_12, double ef2_12, double f2_13, double ef2_13, double f2_14, double ef2_14, double f2_23, double ef2_23,
                     double f2_24, double ef2_24, double f2_34, double ef2_34,
                     double f3_123, double ef3_123, double f3_124, double ef3_124, double f3_134, double ef3_134, double f3_234, double ef3_234,
                     double f4_1234, double ef4_1234,
                     double& F4_1234, double& eF4_1234, double& R4_1234, double& eR4_1234);
  void calculateNudyn(double r2_11, double er2_11, double r2_12, double er2_12, double r2_22, double er2_22, double& nudyn, double& enudyn);

  TString getName() const
  {
    return collectionName;
  }
  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Inputs
  ////////////////////////////////////////////////////////////////////////////
  int nHistoCapacity;
  int nHistograms;
  TH1** histograms;
  bool* isScaled;  // saved to root file at the end
  bool* isSaved;   // saved to root file at the end
  bool* isPlotted; // plot on screen  at the end
  bool* isPrinted; // save histogram content to file at end
  TString collectionName;
  TRandom* randomGenerator;
  bool bOwnTheHistograms; // if the instance owns its own histograms
                          // use to be the case unless the histograms are incorporated to an external list
                          // when true at instance destruction time the histograms are deleted

  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Constants
  ////////////////////////////////////////////////////////////////////////////

  bool scaled;     // = true;
  bool saved;      // = true;
  bool plotted;    // = true;
  bool printed;    //= true;
  bool notScaled;  //= false;
  bool notSaved;   //= false;
  bool notPlotted; //= false;
  bool notPrinted; //= false;

  ClassDef(HistogramCollection, 0);

}; // HistogramCollection

#endif /* WAC_HistogramCollection */

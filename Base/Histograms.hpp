// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *************************************************************************/
#ifndef WAC_Histograms
#define WAC_Histograms
#include "HistogramCollection.hpp"
#include "AnalysisConfiguration.hpp"

class Histograms : public HistogramCollection
{
 public:
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // CTOR
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  Histograms(const TString& name,
             AnalysisConfiguration* configuration,
             int nHistogramsMax,
             LogLevel debugLevel);
  Histograms(TFile* inputFile,
             const TString& name,
             AnalysisConfiguration* configuration,
             int nHistogramsMax,
             LogLevel debugLevel);
  ~Histograms();
  virtual void initialize();
  virtual void finalize();
  virtual void reset();
  virtual void clear();
  virtual void createHistograms();
  virtual void loadHistograms(TDirectory* dir);

  AnalysisConfiguration* getConfiguration() const;
  void setConfiguration(AnalysisConfiguration* config);

  virtual TString getHistoBaseName() const;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  AnalysisConfiguration* configuration;

  ClassDef(Histograms, 0)
};

#endif /* WAC_Histograms */

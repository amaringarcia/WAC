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

 Class defining Event Histograms
 */

#include "Histograms.hpp"
ClassImp(Histograms)

  Histograms::Histograms(const TString& name, AnalysisConfiguration* config,
                         int nHistogramsMax, LogLevel debugLevel)
  : HistogramCollection(name, nHistogramsMax, debugLevel),
    configuration(config)
{
  if (reportDebug())
    cout << "Histograms::CTOR() No OPS." << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// CTOR2
////////////////////////////////////////////////////////////////////////////////////////////////////////
Histograms::Histograms(TFile*, const TString& name,
                       AnalysisConfiguration* configuration, int nHistogramsMax,
                       LogLevel debugLevel)
  : HistogramCollection(name, nHistogramsMax, debugLevel),
    configuration(configuration)
{
  if (reportDebug())
    cout << "Histograms::CTOR() No OPS." << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// DTOR
////////////////////////////////////////////////////////////////////////////////////////////////////////
Histograms::~Histograms()
{
  if (reportDebug())
    cout << "Histograms::DTOR() No ops." << endl;
}

// overload this method only if you need to do more than
// create histogram - overload createHistograms to create histograms.
void Histograms::initialize()
{
  if (reportDebug())
    cout << "Histograms::initialize() Started." << endl;
  createHistograms();
  if (reportDebug())
    cout << "Histograms::initialize() Completed." << endl;
}

void Histograms::finalize()
{
  if (reportDebug())
    cout << "Histograms::finalize() No ops." << endl;
}

void Histograms::reset()
{
  if (reportDebug())
    cout << "Histograms::reset()" << endl;
  for (int iHisto = 0; iHisto < nHistograms; iHisto++) {
    histograms[iHisto]->Reset();
  }
}

void Histograms::clear()
{
  if (reportDebug())
    cout << "Histograms::clear()" << endl;
  for (int iHisto = 0; iHisto < nHistograms; iHisto++) {
    histograms[iHisto]->Reset();
  }
}

// overload this class to create histograms...
void Histograms::createHistograms()
{
  if (reportDebug())
    cout << "Histograms::clear()" << endl;
}

// ==============================================================
// load the cluster histograms from the given file and base name
// ==============================================================
void Histograms::loadHistograms(TDirectory*)
{
  if (reportDebug())
    cout << "Histograms::loadHistograms(...) No ops." << endl;
}

AnalysisConfiguration* Histograms::getConfiguration() const
{
  return configuration;
}

void Histograms::setConfiguration(AnalysisConfiguration* config)
{
  configuration = config;
}

TString Histograms::getHistoBaseName() const
{
  //  AnalysisConfiguration & ac = *(AnalysisConfiguration*) getConfiguration();
  TString bn; // ac.histoBaseName;
  //  bn += "_";
  bn = getName();
  bn += "_";
  return bn;
}

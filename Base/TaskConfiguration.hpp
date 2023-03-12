// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                 *
 * All rights reserved.                                                *
 * Based on the ROOT package and environment                           *
 *                                                                     *
 * For the licensing terms see LICENSE.                                *
 **********************************************************************/

#ifndef WAC_TaskConfiguration
#define WAC_TaskConfiguration
#include <ostream>
#include "TString.h"
using std::endl;
using std::ostream;

////////////////////////////////////////////////////////////////////////////////////////////////////////
// TaskConfiguration of a given object or task
////////////////////////////////////////////////////////////////////////////////////////////////////////
class TaskConfiguration
{
 public:
  TaskConfiguration();
  TaskConfiguration(const TString& theName, const TString& theType, const TString& theVersion);
  TaskConfiguration(const TaskConfiguration& source);
  virtual ~TaskConfiguration();
  TaskConfiguration& operator=(const TaskConfiguration& source);

  TString getName() const;
  TString getType() const;
  TString getVersion() const;
  virtual void printTaskConfiguration(ostream& os);

  TString getInputRootFileName();
  TString getOutputRootFileName();
  TString getOutputDirectory();

  ////////////////////////////////////////////////////
  // Data members
  ////////////////////////////////////////////////////
  TString name;
  TString type;
  TString version;
  bool loadHistograms;
  bool createHistograms;
  bool scaleHistograms;
  bool calculateDerivedHistograms;
  bool saveHistograms;
  bool resetHistograms;
  bool clearHistograms;
  bool forceHistogramsRewrite;

  TString inputPath;
  TString configurationFileName;
  TString rootInputFileName;
  TString outputPath;
  TString rootOuputFileName;
  TString outputDirectory;

  TString dataInputPath;
  TString dataInputFileName;
  TString dataInputTreeName;
  int dataInputFileMinIndex;
  int dataInputFileMaxIndex;
};

#endif /* TaskConfiguration_hpp */

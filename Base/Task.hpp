// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_Task
#define WAC_Task
#include <iostream>
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
#include "TFile.h"
#include "TaskConfiguration.hpp"
#include "MessageLogger.hpp"
#include "Event.hpp"

using namespace std;

class Task : public MessageLogger
{
 public:
  ////////////////////////////////////////////
  // Data members
  ////////////////////////////////////////////
  TString taskName;
  TaskConfiguration* taskConfiguration;
  TRandom* taskRandomGenerator;
  Event* event;

  long eventsProcessed;

  ////////////////////////////////////////////
  // Member functions
  ////////////////////////////////////////////

  Task(const TString& name,
       TaskConfiguration* configuration,
       Event* event);
  virtual ~Task();
  virtual void initialize();
  virtual void execute();
  virtual void finalize();
  virtual void reset();
  virtual void clear();

  virtual void createHistograms();
  virtual void loadHistograms();
  virtual void loadHistograms(TFile* inputFile);
  virtual void resetHistograms();
  virtual void clearHistograms();
  virtual void calculateDerivedHistograms();
  virtual void scaleHistograms();
  virtual void scaleHistograms(double factor);
  virtual void saveHistogramsAsText();
  virtual void saveHistograms();
  virtual void saveHistograms(TDirectory* dir);
  virtual void addHistogramsToExtList(TList* list, bool all = false);
  Event* getEvent();

  TaskConfiguration* getTaskConfiguration();
  void setTaskConfiguration(TaskConfiguration* config);
  virtual void printConfiguration(ostream& output);
  TString getName() const;
  void setName(const TString& name);

  TRandom* getRandomGenerator()
  {
    return taskRandomGenerator;
  }

  void setRandomGenerator(TRandom* randomGenerator);

  enum TaskStatus { Unknown,
                    TaskOk,
                    TaskEof,
                    TaskEod,
                    TaskWarning,
                    TaskError,
                    TaskFatal };

 private:
  static TaskStatus taskStatus;

 public:
  static TaskStatus getTaskStatus();
  static void setTaskStatus(TaskStatus taskStatus);
  static void postTaskOk();
  static void postTaskEof(); // end of file
  static void postTaskEod(); // end of data
  static void postTaskWarning();
  static void postTaskError();
  static void postTaskFatal();
  static bool isTaskOk();
  static TString getTaskStatusName();

  ClassDef(Task, 0)
};

#endif /* WAC_Task */

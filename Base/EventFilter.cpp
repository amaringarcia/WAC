//
//  EventFilter.cpp
//  MyMC
//
//  Created by Claude Pruneau on 12/12/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//

#include "EventFilter.hpp"

ClassImp(EventFilter)

  //////////////////////////////////////////////////////////////////////////////////////////
  // CTOR
  //////////////////////////////////////////////////////////////////////////////////////////
  EventFilter::EventFilter(EventSelection selection, double min, double max)
  : eventSelectionRequested(selection), minimum(min), maximum(max)
{
  // no ops
}

//////////////////////////////////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////////////////////////////////
EventFilter::~EventFilter()
{
  // no ops
}

//////////////////////////////////////////////////////////////////////////////////////////
// accept/reject the given particle based on filter parameter
//////////////////////////////////////////////////////////////////////////////////////////
bool EventFilter::accept(Event& event)
{
  if (eventSelectionRequested == MinBias)
    return true;
  double condition = 0.0;
  switch (eventSelectionRequested) {
    case MinBias:
      break;
    case Centrality:
      condition = event.getMultiplicityClass();
      break;
    case Multiplicity:
      condition = event.getMultiplicity();
      break;
    case ImpactParameter:
      condition = event.getImpactParameter();
      break;
    case Other:
      condition = event.getOther();
      break;
  }

  return (condition >= minimum) && (condition < maximum);
}

// enum EventSelection   { MinBias, Centrality, Multiplicity, ImpactParameter, Other };

TString EventFilter::getName()
{
  TString name;
  switch (eventSelectionRequested) {
    case MinBias:
      name = "MB";
      break;
    case Centrality:
      name = "C";
      break;
    case Multiplicity:
      name = "M";
      break;
    case ImpactParameter:
      name = "B";
      break;
    case Other:
      name = "O";
      break;
  }
  if (eventSelectionRequested == MinBias)
    return name;

  TString minName;
  int minLim = int(minimum);
  int maxLim = int(maximum);
  name += minLim;
  name += "to";
  name += maxLim;
  return name;
}

TString EventFilter::getTitle()
{
  TString name;
  switch (eventSelectionRequested) {
    case MinBias:
      name = "MinBias";
      break;
    case Centrality:
      name = "Centrality";
      break;
    case Multiplicity:
      name = "Mult";
      break;
    case ImpactParameter:
      name = "b";
      break;
    case Other:
      name = "O";
      break;
  }
  if (eventSelectionRequested == MinBias)
    return name;

  TString minName;
  double minLim = double(int(1000.0 * minimum)) / 1000.0;
  double maxLim = double(int(1000.0 * maximum)) / 1000.0;
  name += minLim;
  name += " to ";
  name += maxLim;
  return name;
}

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PythiaConfiguration + ;
#pragma link C++ class PythiaEventReader + ;
#pragma link C++ class PythiaEventGenerator < AnalysisConfiguration::kRapidity> + ;
#pragma link C++ class PythiaEventGenerator < AnalysisConfiguration::kPseudorapidity> + ;

#endif

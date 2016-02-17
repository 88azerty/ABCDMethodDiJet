#include "LLGAnalysis.h"
void LLGAnalysis::SetupTriggerCheck() {

    // setup the _cutFlow
    _cutFlow.insert( pair<string,double>( "trig", 0. ) );
    // and the histograms

    return;
}

void LLGAnalysis::TriggerCheckSelection() {


    for( unsigned int iTrig = 0; iTrig < triggerNames->size(); ++iTrig ) {
      if( triggerNames->at(iTrig).find( "HLT_PFMET170_NoiseCleaned" ) != std::string::npos ) {
        if( triggerBits->at(iTrig) == 1 ) _cutFlow.at("trig") += evtWeight;
      }
    }
    return;
}

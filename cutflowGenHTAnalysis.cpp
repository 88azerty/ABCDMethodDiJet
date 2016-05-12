#include "LLGAnalysis.h"
void LLGAnalysis::SetupGenHTAnalysis() {

    // setup the cutflow

    // and the histograms
    makeHist( "GenLevel_HT", 48, 390, 610, "Gen Level HT", "Number of Events" );
    return;
}

void LLGAnalysis::GenHTAnalysisSelection() {

    double HT = 0.;
    for( unsigned int ip = 0; ip < mct_px -> size(); ++ip ) {
      if( mct_status->at(ip) != 23 ) continue;
      if( mct_id->at(ip) == 21 || abs(mct_id->at(ip)) <= 6 ) {
        HT += sqrt( mct_px->at(ip)*mct_px->at(ip) + mct_py->at(ip)*mct_py->at(ip) );
      }
    }
    cout << "GOT HT " << HT << endl;
    _histograms1D.at("GenLevel_HT").Fill( HT );
    return;
}

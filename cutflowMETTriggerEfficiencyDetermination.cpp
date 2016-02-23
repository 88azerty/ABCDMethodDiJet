#include "LLGAnalysis.h"
void LLGAnalysis::SetupMETTriggerEfficiencyDetermination() {

    // setup the cutflow
    _cutFlow.insert(pair<string,int>("01_MuonTrigger", 0) );
    _cutFlow.insert(pair<string,int>("02_MuonMultiplicity", 0 ) );
    _cutFlow.insert(pair<string,int>("03_MuonTriggerMatching", 0) );
  
    // and the histograms
    makeHist("allEvents_MET_MHT", 1000, 0., 1000., 30, 0., 1000., "MET [GeV]", "MHT [GeV]", "Number of Events", "COLZ", 1.5, 1.5, 1.0 );
    makeHist("passedEvents_MET_MHT", 100, 0., 1000., 30, 0., 1000., "MET [GeV]", "MHT [GeV]", "Number of Events", "COLZ", 1.5, 1.5, 1.0 );
    makeHist("allEvents_MET", 100, 0., 1000., "MET [GeV]", "Number of Events");
    makeHist("passedEvents_MET", 100, 0., 1000., "MET [GeV]", "Number of Events");
    makeHist("allEvents_MHT", 100, 0., 1000., "MHT [GeV]", "Number of Events");
    makeHist("passedEvents_MHT", 100, 0., 1000., "MHT [GeV]", "Number of Events");

    return;
}

void LLGAnalysis::METTriggerEfficiencyDeterminationSelection() {

    bool passTrigger = false;
    for( unsigned int iTrig = 0; iTrig < triggerNames->size(); ++iTrig ) {
      if( triggerNames->at(iTrig).find("HLT_Mu50") != std::string::npos && triggerBits->at(iTrig) == 1 ) passTrigger = true;
    }
    if( !passTrigger ) return;
    _cutFlow.at("01_MuonTrigger") += 1;

    if( muon_px->size() != 1 ) return; 
    _cutFlow.at("02_MuonMultiplicity") += 1;

    //double m_pt = sqrt( muon_px->at(0)*muon_px->at(0) + muon_py->at(0)*muon_py->at(0) );
    double m_eta = muon_eta->at(0);
    double m_phi = muon_phi->at(0);
    bool thisMuonFiredTrigger = false;

    for( unsigned int itrig = 0; itrig < to_TriggerNames->size(); ++itrig ) {
      if( to_TriggerNames->at(itrig).find("HLT_Mu50") == std::string::npos ) continue;
      for( unsigned int ito = 0; ito < to_pt->at(itrig).size(); ++ito ) {
        double deta = to_eta->at(itrig).at(ito) - m_eta;
        double dphi = to_phi->at(itrig).at(ito) - m_phi;
        if( dphi > M_PI ) dphi = 2*M_PI - dphi;
        double dr = sqrt(deta*deta + dphi*dphi);
        if( dr < 0.1 ) thisMuonFiredTrigger = true;
      }
    }

    if( !thisMuonFiredTrigger) return; 
    _cutFlow.at("03_MuonTriggerMatching") += 1;

    double mht = 0;
    double mht_x = 0.;
    double mht_y = 0.;
    for( unsigned int iJet = 0; iJet < recoJet_pt->size(); ++iJet ) {
      if( recoJet_pt->at(iJet).at(SYSJET) < 20. ) continue;
      if( fabs(recoJet_eta->at(iJet)) > 5 ) continue;
      mht_x -= recoJet_pt->at(iJet).at(SYSJET)*TMath::Cos(recoJet_phi->at(iJet));
      mht_y -= recoJet_pt->at(iJet).at(SYSJET)*TMath::Sin(recoJet_phi->at(iJet));
    }
    mht = sqrt( mht_x*mht_x + mht_y*mht_y );
    double rec_met_x = met_x->at(0) + muon_px->at(0);
    double rec_met_y = met_y->at(0) + muon_py->at(0);
    double rec_met = sqrt( rec_met_x*rec_met_x + rec_met_y*rec_met_y );

    _histograms2D.at("allEvents_MET_MHT").Fill(rec_met, mht);
    if( rec_met > 200. ) _histograms1D.at("allEvents_MHT").Fill(mht);
    if( mht > 200. ) _histograms1D.at("allEvents_MET").Fill(rec_met);
   
    bool passTargetTrigger = false;
    for( unsigned int iTrig = 0; iTrig < triggerNames->size(); ++iTrig ) {
      if(     (triggerNames->at(iTrig).find("HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight") != std::string::npos
           ||  triggerNames->at(iTrig).find("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight") != std::string::npos
           ||  triggerNames->at(iTrig).find("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight") != std::string::npos )
           && triggerBits->at(iTrig) == 1 ) passTargetTrigger = true;
    }

    if( !passTargetTrigger ) return; 
    _histograms2D.at("passedEvents_MET_MHT").Fill(rec_met, mht);
    if( rec_met > 200. ) _histograms1D.at("passedEvents_MHT").Fill(mht);
    if( mht > 200. ) _histograms1D.at("passedEvents_MET").Fill(rec_met);


    return;
}

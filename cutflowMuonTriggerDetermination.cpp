#include "LLGAnalysis.h"
#include "TLorentzVector.h"
void LLGAnalysis::SetupMuonTriggerDetermination() {

    // setup the cutflow
    _cutFlow.insert(pair<string,int>("0_NoCut", 0) );
    _cutFlow.insert(pair<string,int>("1_Trigger", 0) );
    _cutFlow.insert(pair<string,int>("2_ElectronVeto", 0) );
    _cutFlow.insert(pair<string,int>("3_2TightMuons", 0) );
    _cutFlow.insert(pair<string,int>("4_ZWindow", 0) );

    // and the histograms
    makeHist("allMuonsPt", 50, 0., 100., "Muon p_{T} [GeV]", "Number of Muons" );
    makeHist("passedMuonsPt", 50, 0., 100., "Muon p_{T} [GeV]", "Number of Muons" );
    makeHist("allMuonsEta", 50, -2.6, 2.6, "Muon #eta", "Number of Muons" );
    makeHist("passedMuonsEta", 50, -2.6, 2.6, "Muon #eta", "Number of Muons" );
    makeHist("allMuonsPhi", 50, -M_PI, M_PI, "Muon #eta", "Number of Muons" );
    makeHist("passedMuonsPhi", 50, -M_PI, M_PI, "Muon #eta", "Number of Muons" );
    return;
}

void LLGAnalysis::MuonTriggerDeterminationSelection() {
    _cutFlow.at("0_NoCut") += 1;



    bool passTrigger = false;
    for( unsigned int iTrig = 0; iTrig < triggerNames->size(); ++iTrig ) {
        if( (triggerNames->at(iTrig).find("HLT_Mu50") != string::npos ) && triggerBits->at(iTrig) == 1 ) passTrigger = true;
    }

    if( !passTrigger ) return;
    _cutFlow.at("1_Trigger") += 1;

    if( vetoElectrons.size() > 0 ) return;
    _cutFlow.at("2_ElectronVeto") += 1;

    if( tightMuons.size() != 2 ) return;
    _cutFlow.at("3_2TightMuons") += 1;

    TLorentzVector LVM1, LVM2;
    double pt1 = sqrt( muon_px->at(tightMuons.at(0))*muon_px->at(tightMuons.at(0)) + muon_py->at(tightMuons.at(0))*muon_py->at(tightMuons.at(0)) );
    double pt2 = sqrt( muon_px->at(tightMuons.at(1))*muon_px->at(tightMuons.at(1)) + muon_py->at(tightMuons.at(1))*muon_py->at(tightMuons.at(1)) );
    LVM1.SetPtEtaPhiM( pt1, muon_eta->at(tightMuons.at(0)), muon_phi->at(tightMuons.at(0)), 0.105 );
    LVM2.SetPtEtaPhiM( pt2, muon_eta->at(tightMuons.at(1)), muon_phi->at(tightMuons.at(1)), 0.105 );

    TLorentzVector LVZC = LVM1 + LVM2;
    if( fabs( LVZC.M() - 91.19 ) > 10. ) return; 
    _cutFlow.at("4_ZWindow") += 1;

    //find matches for the muons;
    vector<int> foundMatch(2, 0);
    for( unsigned int k = 0; k < to_TriggerNames->size(); ++k ) {
      if( to_TriggerNames->at(k) == "HLT_Mu50" ) {
        for( unsigned int to = 0; to < to_pt->at(k).size(); ++to ) {
          for( unsigned int iMuon = 0; iMuon < 2; ++iMuon ) {
            double deta = fabs(muon_eta->at(tightMuons.at(iMuon)) - to_eta->at(k).at(to) );
            double dphi = fabs(muon_phi->at(tightMuons.at(iMuon)) - to_phi->at(k).at(to) );
            if( dphi > M_PI ) dphi = 2*M_PI - dphi;
            double dr = sqrt( dphi*dphi + deta*deta );
            if( dr < 0.4 ) foundMatch.at(iMuon) += 1;
          }
        }
      }
    }
    
    for( unsigned int iMuon2 = 0; iMuon2 < tightMuons.size(); ++iMuon2 ) {
      if( !foundMatch.at(iMuon2) ) continue;
      
      for( unsigned int iMuon = 0; iMuon < tightMuons.size(); ++iMuon ) {
        if( iMuon == iMuon2 ) continue;
        double pt = sqrt( muon_px->at(tightMuons.at(iMuon))*muon_px->at(tightMuons.at(iMuon)) + muon_py->at(tightMuons.at(iMuon))*muon_py->at(tightMuons.at(iMuon)) );
        _histograms1D.at("allMuonsPt").Fill( pt );
        if( pt > 55. ) {
          _histograms1D.at("allMuonsEta").Fill( muon_eta->at(tightMuons.at(iMuon)));
          _histograms1D.at("allMuonsPhi").Fill( muon_phi->at(tightMuons.at(iMuon)));
        }
        if( foundMatch.at(iMuon) > 0 ) _histograms1D.at("passedMuonsPt").Fill(pt);
        if( pt > 55. ) {
          if( foundMatch.at(iMuon) > 0 ) _histograms1D.at("passedMuonsEta").Fill(muon_eta->at(tightMuons.at(iMuon)));
          if( foundMatch.at(iMuon) > 0 ) _histograms1D.at("passedMuonsPhi").Fill(muon_phi->at(tightMuons.at(iMuon)));
        }
      }

    }


    return;
}

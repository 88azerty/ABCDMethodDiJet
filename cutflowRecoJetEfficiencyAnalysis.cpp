#include "LLGAnalysis.h"
#include "TLorentzVector.h"
void LLGAnalysis::SetupRecoJetEfficiencyAnalysis() {

    // setup the cutflow
    _cutFlow.insert( pair<string, int>( "total", 0 ) );
    _cutFlow.insert( pair<string, int>( "matchedCHS", 0 ) );
    _cutFlow.insert( pair<string, int>( "matchedNoCHS", 0 ) );
    // and the histograms
    
    makeHist( "allQOI_pt",              100, 0., 250., "QOI p_{T} [GeV]", "Number of QOIs" );
    makeHist( "chsMatchedQOI_pt",       100, 0., 250., "QOI p_{T} [GeV]", "Number of QOIs" );
    makeHist( "nochsMatchedQOI_pt",     100, 0., 250., "QOI p_{T} [GeV]", "Number of QOIs" );
    makeHist( "allQOI_eta",             60, -3., 3., "QOI #eta", "Number of QOIs" );
    makeHist( "chsMatchedQOI_eta",      60, -3., 3., "QOI #eta", "Number of QOIs" );
    makeHist( "nochsMatchedQOI_eta",    60, -3., 3., "QOI #eta", "Number of QOIs" );
    makeHist( "allQOI_drpvsv",          100, 0., 650., "QOI dr PVSV [mm]", "Number of QOIs" );
    makeHist( "chsMatchedQOI_drpvsv",   100, 0., 650., "QOI dr PVSV [mm]", "Number of QOIs" );
    makeHist( "nochsMatchedQOI_drpvsv", 100, 0., 650., "QOI dr PVSV [mm]", "Number of QOIs" );
    
    makeHist( "chsMatchedSJQOI_pt",       100, 0., 250., "QOI p_{T} [GeV]", "Number of QOIs" );
    makeHist( "nochsMatchedSJQOI_pt",     100, 0., 250., "QOI p_{T} [GeV]", "Number of QOIs" );
    makeHist( "chsMatchedSJQOI_eta",      60, -3., 3., "QOI #eta", "Number of QOIs" );
    makeHist( "nochsMatchedSJQOI_eta",    60, -3., 3., "QOI #eta", "Number of QOIs" );
    makeHist( "chsMatchedSJQOI_drpvsv",   100, 0., 650., "QOI dr PVSV [mm]", "Number of QOIs" );
    makeHist( "nochsMatchedSJQOI_drpvsv", 100, 0., 650., "QOI dr PVSV [mm]", "Number of QOIs" );
    
    
    makeHist( "chsMatchedQOI_dpt",      100, -10., 1.1, "(QOI p_{T} - MatchedJet p_{T})/QOI p_{T}", "Number of QOIs/RecoJets" );
    makeHist( "nochsMatchedQOI_dpt",    100, -10., 1.1, "(QOI p_{T} - MatchedJet p_{T})/QOI p_{T}", "Number of QOIs/RecoJets" );
    makeHist( "bothMatched_chsMatchedQOI_dpt",      100, -10., 1.1, "(QOI p_{T} - MatchedJet p_{T})/QOI p_{T}", "Number of QOIs/RecoJets" );
    makeHist( "bothMatched_nochsMatchedQOI_dpt",    100, -10., 1.1, "(QOI p_{T} - MatchedJet p_{T})/QOI p_{T}", "Number of QOIs/RecoJets" );
  
    makeHist( "NoCHSJets_FromPV_ds_PV", 100, 0., 10., "Distances of NoCHS Jet Constitutents (fromPV) to matched PV [mm]", "Number of Jet Constituents" );
    makeHist( "NoCHSJets_FromPV_ds_trueSV", 100, 0., 10., "Distances of NoCHS Jet Constitutents (fromPV) to true SV [mm]", "Number of Jet Constituents" );
    makeHist( "NoCHSJets_NotFromPV_ds_PV", 100, 0., 10., "Distances of NoCHS Jet Constitutents (not fromPV) to matched PV [mm]", "Number of Jet Constituents" );
    makeHist( "NoCHSJets_NotFromPV_ds_trueSV", 100, 0., 10., "Distances of NoCHS Jet Constitutents (not fromPV) to true SV [mm]", "Number of Jet Constituents" );
    makeHist( "NoCHSJets_FromPV_ds_PV_vs_ds_trueSV", 100, 0., 10., 100, 0., 10., "Distances of NoCHS Jet Constitutents (fromPV) to matched PV [mm]", 
                                                                                    "Distances of NoCHS Jet Constitutents (fromPV) to true SV [mm]", "Number of Jet Constituents", "COLZ" );
    makeHist( "NoCHSJets_NotFromPV_ds_PV_vs_ds_trueSV", 100, 0., 10., 100, 0., 10., "Distances of NoCHS Jet Constitutents (not fromPV) to matched PV [mm]", 
                                                                                    "Distances of NoCHS Jet Constitutents (not fromPV) to true SV [mm]", "Number of Jet Constituents", "COLZ" );

    makeHist( "chsMatched_DeltaDTruth", 100, 0., 20., "QOI ds SV-SV [mm]", "Number of QOIs" );
    makeHist( "nochsMatched_DeltaDTruth", 100, 0., 20., "QOI ds SV-SV [mm]", "Number of QOIs" );
    makeHist( "chsMatched_fromPV", 4, -0.5, 3.5, "Matched CHS Jets Constitutents from PV", "Number of Constitutents" );
    makeHist( "nochsMatched_fromPV", 4, -0.5, 3.5, "Matched NoCHS Jets Constitutents from PV", "Number of Constitutents" );
    makeHist( "chsMatched_fromPV_drpvsv", 4, -0.5, 3.5, 100, 0., 10., "Matched CHS Jets Constitutents from PV", "QOI dr PVSV [mm]", "Number of Constituents", "COLZ" );
    makeHist( "nochsMatched_fromPV_drpvsv", 4, -0.5, 3.5, 100, 0., 10., "Matched NoCHS Jets Constitutents from PV", "QOI dr PVSV [mm]", "Number of Constituents", "COLZ" );
    makeHist( "allMatchedCHSJets_PVSV", 3, -0.5, 2.5, "Matched CHS Jet associated to Vertex?", "Number of CHS Jets" );    
    makeHist( "allMatchedNoCHSJets_PVSV", 3, -0.5, 2.5, "Matched No-CHS Jet associated to Vertex?", "Number of No-CHS Jets" );    
    makeHist( "matchedCHSJets_PVAss_dxAssCalc", 200, -10., 10., "#Delta x - Associated PV, Matched PV [mm]", "Number of Constituents" );
    makeHist( "matchedCHSJets_PVAss_dyAssCalc", 200, -10., 10., "#Delta y - Associated PV, Matched PV [mm]", "Number of Constituents" );
    makeHist( "matchedCHSJets_PVAss_dzAssCalc", 200, -20., 20., "#Delta z - Associated PV, Matched PV [mm]", "Number of Constituents" );
    makeHist( "matchedCHSJets_PVAss_drAssCalc", 200, 0., 10., "#Delta r - Associated PV, Matched PV [mm]", "Number of Constituents" );
    makeHist( "matchedCHSJets_PVAss_dsAssCalc", 200, 0., 10., "#Delta s - Associated PV, Matched PV [mm]", "Number of Constituents" );
    makeHist( "matchedCHSJets_SVAss_dxAssCalc", 200, -10., 10., "#Delta x - Associated PV, Matched SV [mm]", "Number of Constituents" );
    makeHist( "matchedCHSJets_SVAss_dyAssCalc", 200, -10., 10., "#Delta y - Associated PV, Matched SV [mm]", "Number of Constituents" );
    makeHist( "matchedCHSJets_SVAss_dzAssCalc", 200, -10., 10., "#Delta z - Associated PV, Matched SV [mm]", "Number of Constituents" );
    makeHist( "matchedCHSJets_SVAss_drAssCalc", 200, 0., 50., "#Delta r - Associated PV, Matched SV [mm]", "Number of Constituents" );
    makeHist( "matchedCHSJets_SVAss_dsAssCalc", 200, 0., 50., "#Delta s - Associated PV, Matched SV [mm]", "Number of Constituents" );
    makeHist( "matchedNoCHSJets_PVAss_dxAssCalc", 200, -10., 10., "#Delta x - Associated PV, Matched PV [mm]", "Number of Constituents" );
    makeHist( "matchedNoCHSJets_PVAss_dyAssCalc", 200, -10., 10., "#Delta y - Associated PV, Matched PV [mm]", "Number of Constituents" );
    makeHist( "matchedNoCHSJets_PVAss_dzAssCalc", 200, -20., 20., "#Delta z - Associated PV, Matched PV [mm]", "Number of Constituents" );
    makeHist( "matchedNoCHSJets_PVAss_drAssCalc", 200, 0., 10., "#Delta r - Associated PV, Matched PV [mm]", "Number of Constituents" );
    makeHist( "matchedNoCHSJets_PVAss_dsAssCalc", 200, 0., 10., "#Delta s - Associated PV, Matched PV [mm]", "Number of Constituents" );
    makeHist( "matchedNoCHSJets_SVAss_dxAssCalc", 200, -10., 10., "#Delta x - Associated PV, Matched SV [mm]", "Number of Constituents" );
    makeHist( "matchedNoCHSJets_SVAss_dyAssCalc", 200, -10., 10., "#Delta y - Associated PV, Matched SV [mm]", "Number of Constituents" );
    makeHist( "matchedNoCHSJets_SVAss_dzAssCalc", 200, -10., 10., "#Delta z - Associated PV, Matched SV [mm]", "Number of Constituents" );
    makeHist( "matchedNoCHSJets_SVAss_drAssCalc", 200, 0., 50., "#Delta r - Associated PV, Matched SV [mm]", "Number of Constituents" );
    makeHist( "matchedNoCHSJets_SVAss_dsAssCalc", 200, 0., 50., "#Delta s - Associated PV, Matched SV [mm]", "Number of Constituents" );
    makeHist( "matchedCHSJets_DeltaS_TrueCalc", 200, 0., 100., "#DeltaS - CHS Jets, True vs Calculated [mm]", "Number of Jets" );
    makeHist( "matchedCHSJets_DeltaR_TrueCalc", 200, 0., 100., "#DeltaR - CHS Jets, True vs Calculated [mm]", "Number of Jets" );
    makeHist( "matchedCHSJets_DeltaX_TrueCalc", 200, -25., 25., "#DeltaX - CHS Jets, True vs Calculated [mm]", "Number of Jets" );
    makeHist( "matchedCHSJets_DeltaY_TrueCalc", 200, -25., 25., "#DeltaY - CHS Jets, True vs Calculated [mm]", "Number of Jets" );
    makeHist( "matchedCHSJets_DeltaZ_TrueCalc", 200, -25., 25., "#DeltaZ - CHS Jets, True vs Calculated [mm]", "Number of Jets" );
    makeHist( "matchedNoCHSJets_DeltaS_TrueCalc", 200, 0, 100., "#DeltaS - NoCHS Jets, True vs Calculated [mm]", "Number of Jets" );
    makeHist( "matchedNoCHSJets_DeltaR_TrueCalc", 200, 0., 100., "#DeltaR - NoCHS Jets, True vs Calculated [mm]", "Number of Jets" );
    makeHist( "matchedNoCHSJets_DeltaX_TrueCalc", 200, -25., 25., "#DeltaX - NoCHS Jets, True vs Calculated [mm]", "Number of Jets" );
    makeHist( "matchedNoCHSJets_DeltaY_TrueCalc", 200, -25., 25., "#DeltaY - NoCHS Jets, True vs Calculated [mm]", "Number of Jets" );
    makeHist( "matchedNoCHSJets_DeltaZ_TrueCalc", 200, -25., 25., "#DeltaZ - NoCHS Jets, True vs Calculated [mm]", "Number of Jets" );
    makeHist( "GluinoDecayLength", 2000, 0., 200., "Gluino Decay Length", "Number of Events" );
    makeHist( "matchedCHSJet_ds_SVAss_drmin5_PCAPV", 200, 0., 50, "Distance at PCA wrt Associated PV [mm]", "Number of CHS Jets" );
    makeHist( "matchedCHSJet_ds_vs_constPt_SVAss_drmin5_PCAPV", 200, 0., 50, 200, 0., 20., "Distance at PCA wrt PV", "Constituent pT [mm]", "Number of CHS Jets", "COLZ" );
    makeHist( "matchedCHSJet_dr_SVAss_drmin5_PCAPV", 200, 0., 50, "Transverse Distance at PCA wrt Associated PV [mm]", "Number of CHS Jets" );
    makeHist( "matchedCHSJet_dx_SVAss_drmin5_PCAPV", 200, -50., 50, "Distance at PCA in x wrt Associated PV [mm]", "Number of CHS Jets" );
    makeHist( "matchedCHSJet_dy_SVAss_drmin5_PCAPV", 200, -50., 50, "Distance at PCA in y wrt Associated PV [mm]", "Number of CHS Jets" );
    makeHist( "matchedCHSJet_dz_SVAss_drmin5_PCAPV", 200, -50., 50, "Distance at PCA in z wrt Associated PV [mm]", "Number of CHS Jets" );


    makeHist( "PVisGluinoProductionVertex", 2, -0.5, 1.5, "PV is Gluino Production Vertex", "Number of Events" );
    makeHist( "PVGluinoProudctionVertex_dx", 100, -1., 1., "Distance in x between PV and GluinoProductionVertex [mm]", "Number of Events" );
    makeHist( "PVGluinoProudctionVertex_dy", 100, -1., 1., "Distance in y between PV and GluinoProductionVertex [mm]", "Number of Events" );
    makeHist( "PVGluinoProudctionVertex_dz", 100, -1., 1., "Distance in z between PV and GluinoProductionVertex [mm]", "Number of Events" );
    makeHist( "PVGluinoProudctionVertex_dr", 100, 0., 1., "Distance in r between PV and GluinoProductionVertex [mm]", "Number of Events" );
    makeHist( "PVGluinoProudctionVertex_ds", 100, 0., 1., "Distance in s between PV and GluinoProductionVertex [mm]", "Number of Events" );

    // and the truthtree
    lastTruthEntry = -1; 
    fTruth = new TFile( GenFileName.c_str(), "OPEN");
    tTruth = (TTree*)fTruth->Get("outputSorted");
    tmct_px = new vector<double>;
    tmct_py = new vector<double>;
    tmct_pz = new vector<double>;
    tmct_vx = new vector<double>;
    tmct_vy = new vector<double>;
    tmct_vz = new vector<double>;
    tmct_e = new vector<double>;
    tmct_id = new vector<int>;
    tmct_status = new vector<int>;
    tmct_parent = new vector<int>;
    tmct_daughters = new vector<vector<int> >;
    tTruth->SetBranchAddress("EventNumber", &GenEventNumber );
    tTruth->SetBranchAddress("LumiBlock", &GenLumiBlock );
    tTruth->SetBranchAddress("RunNumber", &GenRunNumber );
    tTruth->SetBranchAddress("Particle_px", &tmct_px );
    tTruth->SetBranchAddress("Particle_py", &tmct_py );
    tTruth->SetBranchAddress("Particle_pz", &tmct_pz );
    tTruth->SetBranchAddress("Particle_vx", &tmct_vx );
    tTruth->SetBranchAddress("Particle_vy", &tmct_vy );
    tTruth->SetBranchAddress("Particle_vz", &tmct_vz );
    tTruth->SetBranchAddress("Particle_E", &tmct_e );
    tTruth->SetBranchAddress("Particle_ID", &tmct_id );
    tTruth->SetBranchAddress("Particle_Status", &tmct_status );
    tTruth->SetBranchAddress("Particle_Parent", &tmct_parent );
    tTruth->SetBranchAddress("Particle_Daughter", &tmct_daughters );
    return;
}

void LLGAnalysis::RecoJetEfficiencyAnalysisSelection() {
   
   
    vector<int> recoCHSJet_assignedPV(   recoCHSJet_pt->size(), -1 );
    vector<int> recoNoCHSJet_assignedPV( recoNoCHSJet_pt->size(), -1 );
    vector<int> recoCHSJet_assignedSV(   recoCHSJet_pt->size(), -1 );
    vector<int> recoNoCHSJet_assignedSV( recoNoCHSJet_pt->size(), -1 );

    for( unsigned int iJet = 0; iJet < recoCHSJet_pt->size(); ++iJet ) {
        //calculate jet vertex position:
        unsigned int nCons;
        double weightednCons;
        vector<double> error(3,0.);
        vector<double> position = CalculateVertex( recoCHSJet_constVertex_x->at(iJet), recoCHSJet_constVertex_y->at(iJet), recoCHSJet_constVertex_z->at(iJet), recoCHSJet_const_pt->at(iJet), recoCHSJet_const_charge->at(iJet), recoCHSJet_const_closestVertex_d->at(iJet), nCons, weightednCons, error );
        int nMatch = 0;
        
        //std::cout << "this is jet # " << iJet << " at position " << position.at(0) << " " << position.at(1) << " " << position.at(2) << std::endl; 
        for( unsigned int iVtx = 0; iVtx < vertex_x->size(); ++iVtx ) {
            if( fabs(position.at(0) - vertex_x->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(1) - vertex_y->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(2) - vertex_z->at(iVtx) ) < 1.e-10 ) {
                recoCHSJet_assignedPV.at(iJet) = iVtx;
                nMatch += 1;
            }
        }
        for( unsigned int iVtx = 0; iVtx < secVertex_x->size(); ++iVtx ) {
            //std::cout << "checking SV with " << secVertex_x->at(iVtx) << " " << secVertex_y->at(iVtx) << " " << secVertex_z->at(iVtx) << std::endl;
            if( fabs(position.at(0) - secVertex_x->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(1) - secVertex_y->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(2) - secVertex_z->at(iVtx) ) < 1.e-10 ) {
                recoCHSJet_assignedSV.at(iJet) = iVtx;    
                nMatch += 1;
            }
        }
        if( nMatch > 1 ) {
            cout << "WARNING! ASSOCIATED JET TO MORE THAN 1 VERTEX ?!" << endl;
        }
    }
    
    for( unsigned int iJet = 0; iJet < recoNoCHSJet_pt->size(); ++iJet ) {
        //calculate jet vertex position:
        unsigned int nCons;
        double weightednCons;
        vector<double> error(3,0.);
        vector<double> position = CalculateVertex( recoNoCHSJet_constVertex_x->at(iJet), recoNoCHSJet_constVertex_y->at(iJet), recoNoCHSJet_constVertex_z->at(iJet), recoNoCHSJet_const_pt->at(iJet), recoNoCHSJet_const_charge->at(iJet), recoNoCHSJet_const_closestVertex_d->at(iJet), nCons, weightednCons, error );
        int nMatch = 0;
        
        //std::cout << "this is jet # " << iJet << " at position " << position.at(0) << " " << position.at(1) << " " << position.at(2) << std::endl; 
        for( unsigned int iVtx = 0; iVtx < vertex_x->size(); ++iVtx ) {
            if( fabs(position.at(0) - vertex_x->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(1) - vertex_y->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(2) - vertex_z->at(iVtx) ) < 1.e-10 ) {
                recoNoCHSJet_assignedPV.at(iJet) = iVtx;
                nMatch += 1;
            }
        }
        for( unsigned int iVtx = 0; iVtx < secVertex_x->size(); ++iVtx ) {
            //std::cout << "checking SV with " << secVertex_x->at(iVtx) << " " << secVertex_y->at(iVtx) << " " << secVertex_z->at(iVtx) << std::endl;
            if( fabs(position.at(0) - secVertex_x->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(1) - secVertex_y->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(2) - secVertex_z->at(iVtx) ) < 1.e-10 ) {
                recoNoCHSJet_assignedSV.at(iJet) = iVtx;    
                nMatch += 1;
            }
        }
        if( nMatch > 1 ) {
            cout << "WARNING! ASSOCIATED JET TO MORE THAN 1 VERTEX ?!" << endl;
        }
    }


    bool FoundTruthEvent = false;
    int idxTruthEvent = -1;
    if( lastTruthEntry < tTruth->GetEntries() ) {
    tTruth->GetEntry(lastTruthEntry+1 ); 
    if( RunNumber == GenRunNumber && LumiBlock == GenLumiBlock && EventNumber == GenEventNumber ) {
      lastTruthEntry += 1;
      FoundTruthEvent = true;
    }
    }
    if( !FoundTruthEvent ) {
      for( unsigned int itruth = 0; itruth < tTruth->GetEntries(); ++itruth ) {
        tTruth->GetEntry(itruth);
        if( RunNumber == GenRunNumber && LumiBlock == GenLumiBlock && EventNumber == GenEventNumber ) {
          FoundTruthEvent = true;
          idxTruthEvent = itruth;
          break;
        }
      }
      lastTruthEntry = idxTruthEvent;
    }
    if( !FoundTruthEvent ) {
      cout << "DID NOT FIND TRUTH INFO " << endl;
      getchar();
    }
    // now find the quarks from the gluino decay
    vector<double> GluinoProductionVertex(3, -10000. );
    vector<vector<double> > GluinoDecayVertex;
    vector<double> decV1(3, -10000.);
    GluinoDecayVertex.push_back( decV1 );
    GluinoDecayVertex.push_back( decV1 );
    vector<int> qoi;
    int nGluino = 0;
    /*
    for( unsigned int ipart = 0; ipart < tmct_px->size(); ++ipart ) {
      if( tmct_id->at(ipart) == 1000021 ) {
        if( abs( tmct_id->at( tmct_parent->at(ipart) )) < 100  ) {
          double gpv_x = tmct_vx->at(ipart); double gpv_y = tmct_vy->at(ipart); double gpv_z = tmct_vz->at(ipart);
          TLorentzVector GluinoP4;
          GluinoP4.SetPxPyPzE( tmct_px->at(ipart), tmct_py->at(ipart), tmct_pz->at(ipart), tmct_e->at(ipart) );
          int gluinoID = ipart;
          int SUSYDaughter = -1;
          for( unsigned int ic = 0; ic < tmct_daughters->at(gluinoID).size(); ++ic ) {
            if( abs(tmct_id->at(tmct_daughters->at(gluinoID).at(ic))) > 1000000 ) SUSYDaughter = ic;
          }
          if( SUSYDaughter < 0 ) {
            cout << "RPV!" << endl;
            getchar();
          }
          while( tmct_id->at( tmct_daughters->at(gluinoID).at(SUSYDaughter)) != 1000022 ) {
            gluinoID = tmct_daughters->at(gluinoID).at(SUSYDaughter);
            for( unsigned int ic = 0; ic < tmct_daughters->at(gluinoID).size(); ++ic ) {
              if( abs(tmct_id->at(tmct_daughters->at(gluinoID).at(ic))) > 1000000 ) SUSYDaughter = ic;
            }
            std::cout << "\tsubsequent gluino (ID = " << tmct_id->at(gluinoID) << " found with " << tmct_daughters->at(gluinoID).size() << " daughters : ";
            if( tmct_daughters->at(gluinoID).size() == 0 ) {
              break;
            }
            for( unsigned int ic = 0; ic < tmct_daughters->at(gluinoID).size(); ++ic ) {
              std::cout << tmct_id->at(tmct_daughters->at(gluinoID).at(ic)) << " ";
            }
            std::cout << endl;
            
          }
          if( tmct_daughters->at(gluinoID).size() == 0 ) continue;
          int neutralinoID = tmct_daughters->at(gluinoID).at(SUSYDaughter);
          if( tmct_id->at(neutralinoID) != 1000022 ) continue;
          double gdv_x = tmct_vx->at(neutralinoID); double gdv_y = tmct_vy->at(neutralinoID); double gdv_z = tmct_vz->at(neutralinoID);
          double dx = gpv_x - gdv_x; double dy = gpv_y - gdv_y; double dz = gpv_z - gdv_z;
          _histograms1D.at("GluinoDecayLength").Fill( sqrt( dx*dx + dy*dy + dz*dz )*GluinoP4.M()/GluinoP4.P() );
        }
      }
    }
    */
    for( unsigned int ipart = 0; ipart < tmct_px->size(); ++ipart ) {
      if( tmct_id->at(ipart) == 1000021 ) {
        int gluinoParent = tmct_parent->at(ipart);
        if( abs( tmct_id->at(gluinoParent)) < 100  ) {
          if( GluinoProductionVertex.at(0) < -9999. ) {
            GluinoProductionVertex.at(0) = tmct_vx->at(ipart);
            GluinoProductionVertex.at(1) = tmct_vy->at(ipart);
            GluinoProductionVertex.at(2) = tmct_vz->at(ipart);
            double distMin = 100000.;
            int bestMatch = -1;
            for( unsigned int iVtx = 0; iVtx < vertex_x->size(); ++iVtx ) {
              double ldx = GluinoProductionVertex.at(0) - vertex_x->at(iVtx);
              double ldy = GluinoProductionVertex.at(1) - vertex_y->at(iVtx);
              double ldz = GluinoProductionVertex.at(2) - vertex_z->at(iVtx);
              double dist = sqrt( ldx*ldx + ldy*ldy + ldz*ldz );
              if( dist < distMin ) {
                bestMatch = iVtx;
                distMin = dist;
              }
            }
            double ldxmin = GluinoProductionVertex.at(0) - vertex_x->at(bestMatch);
            double ldymin = GluinoProductionVertex.at(1) - vertex_y->at(bestMatch);
            double ldzmin = GluinoProductionVertex.at(2) - vertex_z->at(bestMatch);
            bool hasPVJet150 = false;
            for( unsigned int iJet = 0; iJet < recoCHSJet_pt->size(); ++iJet ) {
              if( recoCHSJet_pt->at(iJet) > 150. && recoCHSJet_assignedPV.at(iJet) >= 0. ) hasPVJet150 = true;
            }
            if( hasPVJet150 ) {
              _histograms1D.at( "PVisGluinoProductionVertex" ).Fill( ((bestMatch == 0 ) ? 1 : 0 ));
              _histograms1D.at( "PVGluinoProudctionVertex_dx" ).Fill( 10.*ldxmin );
              _histograms1D.at( "PVGluinoProudctionVertex_dy" ).Fill( 10.*ldymin );
              _histograms1D.at( "PVGluinoProudctionVertex_dz" ).Fill( 10.*ldzmin );
              _histograms1D.at( "PVGluinoProudctionVertex_dr" ).Fill( 10.*sqrt( ldxmin*ldxmin + ldymin*ldymin ) );
              _histograms1D.at( "PVGluinoProudctionVertex_ds" ).Fill( 10.*sqrt( ldxmin*ldxmin + ldymin*ldymin + ldzmin*ldzmin) );
            }
          }
          nGluino += 1;
        }
      }
      if( tmct_id->at(ipart) == 1000021 && tmct_daughters->at(ipart).size() == 3 && tmct_id->at(tmct_daughters->at(ipart).at(0)) == 1000022 ) {
        //cout << endl << "got a gluino decaying into " << tmct_id->at(tmct_daughters->at(ipart).at(0)) << " " << tmct_id->at(tmct_daughters->at(ipart).at(1)) << " " << tmct_id->at(tmct_daughters->at(ipart).at(2)) << endl;       
        int dau1 = tmct_daughters->at(ipart).at(1);
        int dau2 = tmct_daughters->at(ipart).at(2);
        if( qoi.size() == 0 ) {
          GluinoDecayVertex.at(0).at(0) = tmct_vx->at(dau1);
          GluinoDecayVertex.at(0).at(1) = tmct_vy->at(dau1);
          GluinoDecayVertex.at(0).at(2) = tmct_vz->at(dau1);
        
          double dec_dx = GluinoProductionVertex.at(0) - GluinoDecayVertex.at(0).at(0);
          double dec_dy = GluinoProductionVertex.at(1) - GluinoDecayVertex.at(0).at(1);
          double dec_dz = GluinoProductionVertex.at(2) - GluinoDecayVertex.at(0).at(2);
          TLorentzVector GluinoP4;
          GluinoP4.SetPxPyPzE( tmct_px->at(ipart), tmct_py->at(ipart), tmct_pz->at(ipart), tmct_e->at(ipart) );
          double d0 = sqrt( dec_dx*dec_dx + dec_dy*dec_dy + dec_dz*dec_dz )*GluinoP4.M()/GluinoP4.P();
          _histograms1D.at("GluinoDecayLength").Fill( d0 );
        }
        else {
          GluinoDecayVertex.at(1).at(0) = tmct_vx->at(dau1);
          GluinoDecayVertex.at(1).at(1) = tmct_vy->at(dau1);
          GluinoDecayVertex.at(1).at(2) = tmct_vz->at(dau1);
          
          double dec_dx = GluinoProductionVertex.at(0) - GluinoDecayVertex.at(1).at(0);
          double dec_dy = GluinoProductionVertex.at(1) - GluinoDecayVertex.at(1).at(1);
          double dec_dz = GluinoProductionVertex.at(2) - GluinoDecayVertex.at(1).at(2);
          TLorentzVector GluinoP4;
          GluinoP4.SetPxPyPzE( tmct_px->at(ipart), tmct_py->at(ipart), tmct_pz->at(ipart), tmct_e->at(ipart) );
          double d0 = sqrt( dec_dx*dec_dx + dec_dy*dec_dy + dec_dz*dec_dz )*GluinoP4.M()/GluinoP4.P();
          _histograms1D.at("GluinoDecayLength").Fill( d0 );
          
        }
        bool dau1Final = true;
        bool dau2Final = true;
         
        do {
          dau1Final = true;
          if( tmct_daughters->at(dau1).size() == 0 || tmct_daughters->at(dau1).at(0) < 0 || tmct_daughters->at(dau1).at(0) == dau1 ) break;
          for( unsigned int icc = 0; icc < tmct_daughters->at(dau1).size(); ++icc ) {
            if( tmct_daughters->at(dau1).at(icc) < 0) continue;
            if( tmct_id->at(tmct_daughters->at(dau1).at(icc)) == tmct_id->at(dau1) ) {
              dau1Final = false;
              dau1 = tmct_daughters->at(dau1).at(icc);
              continue;
            }
          }
        } while( !dau1Final );
        do {
          dau2Final = true;
          if( tmct_daughters->at(dau2).size() == 0 || tmct_daughters->at(dau2).at(0) < 0 || tmct_daughters->at(dau2).at(0) == dau2 ) break;
          for( unsigned int icc = 0; icc < tmct_daughters->at(dau2).size(); ++icc ) {
            if( tmct_daughters->at(dau2).at(icc) < 0) continue;
            if( tmct_id->at(tmct_daughters->at(dau2).at(icc)) == tmct_id->at(dau2) ) {
              dau2Final = false;
              dau2 = tmct_daughters->at(dau2).at(icc);
              continue;
            }
          }
        } while( !dau2Final );
        qoi.push_back( dau1 );
        qoi.push_back( dau2 );
      }
    }
    //got the quarks of interest. now try to find a jet matching them;
    for( unsigned int ipoi = 0; ipoi < qoi.size(); ++ipoi ) {
      _cutFlow.at("total") += 1;
      int part = qoi.at(ipoi);
      double px = tmct_px->at(part);
      double py = tmct_py->at(part);
      double pz = tmct_pz->at(part);
      double energy = sqrt( px*px + py*py + pz*pz );
      TLorentzVector qVec;
      qVec.SetPxPyPzE( px, py, pz, energy );
      double quark_phi = qVec.Phi();
      double quark_eta = qVec.Eta();
      
      double ds_x = GluinoProductionVertex.at(0) - tmct_vx->at(part);
      double ds_y = GluinoProductionVertex.at(1) - tmct_vy->at(part);
      double ds_z = GluinoProductionVertex.at(2) - tmct_vz->at(part);
      double ds_r = 10.*sqrt( ds_x*ds_x + ds_y*ds_y );
      _histograms1D.at("allQOI_pt").Fill( qVec.Pt() );
      if( qVec.Pt() > 50 ) {
        _histograms1D.at("allQOI_eta").Fill( qVec.Eta() );
        _histograms1D.at("allQOI_drpvsv").Fill( ds_r ); 
      }
      double drmin = 10000.;
      int matchidxjet = -1;
      int matchCHSID = -1;
      int matchNoCHSID = -1;
      for( unsigned int iJet = 0; iJet < recoCHSJet_pt->size(); ++iJet ) {
        //if( recoCHSJet_pt->at(iJet) < 30. ) continue;
        double deta = quark_eta - recoCHSJet_eta->at(iJet);
        double dphi = quark_phi - recoCHSJet_phi->at(iJet);
        if( dphi > M_PI ) dphi = 2*M_PI-dphi;
        double dr = sqrt( deta*deta + dphi*dphi );
        if( dr < drmin ) { drmin = dr; matchidxjet = iJet; matchCHSID = iJet; }
      }
      if( drmin < 0.2 ) {
        _cutFlow.at("matchedCHS") += 1;
        if( recoCHSJet_assignedPV.at(matchidxjet) >= 0 && recoCHSJet_assignedSV.at(matchidxjet) >= 0 ) {
          cout << "Jet assigned to PV and SV - this should not happen" << endl;
          getchar();
        }
        else if( recoCHSJet_assignedPV.at(matchidxjet) >= 0 ) _histograms1D.at("allMatchedCHSJets_PVSV").Fill( 1 );
        else if( recoCHSJet_assignedSV.at(matchidxjet) >= 0 ) _histograms1D.at("allMatchedCHSJets_PVSV").Fill( 2 );
        else _histograms1D.at("allMatchedCHSJets_PVSV").Fill( 0 );

        _histograms1D.at("chsMatchedQOI_pt").Fill( qVec.Pt() );
        if( qVec.Pt() > 50. ) {
          _histograms1D.at("chsMatchedQOI_eta").Fill( qVec.Eta() );
          _histograms1D.at("chsMatchedQOI_drpvsv").Fill( ds_r ); 
          _histograms1D.at("chsMatchedQOI_dpt").Fill( (qVec.Pt() - recoCHSJet_pt->at(matchidxjet))/qVec.Pt() );
        }
        for( unsigned int iConst = 0; iConst <  recoCHSJet_const_fromPV->at(matchidxjet).size(); ++iConst ) {
          if( recoCHSJet_const_charge->at(matchidxjet).at(iConst) == 0 ) continue; 
          _histograms1D.at("chsMatched_fromPV").Fill( recoCHSJet_const_fromPV->at(matchidxjet).at(iConst));
          _histograms2D.at("chsMatched_fromPV_drpvsv").Fill( recoCHSJet_const_fromPV->at(matchidxjet).at(iConst), ds_r );
        }
        

        if( recoCHSJet_pt->at(matchidxjet) > 30. && recoCHSJet_assignedSV.at(matchidxjet) >= 0 ) {
          int thisJetMatchedIdx = recoCHSJet_assignedSV.at(matchidxjet);
          for( unsigned int iJet2 = 0; iJet2 < recoCHSJet_pt->size(); ++iJet2 ) {
            if( iJet2 == matchidxjet ) continue;
            if( recoCHSJet_assignedSV.at(iJet2) == thisJetMatchedIdx && recoCHSJet_pt->at(iJet2) > 30. ) {
              _histograms1D.at("chsMatchedSJQOI_pt").Fill( qVec.Pt() );
              if( qVec.Pt() > 50 ) {
                _histograms1D.at("chsMatchedSJQOI_eta").Fill( qVec.Eta() );
                _histograms1D.at("chsMatchedSJQOI_drpvsv").Fill( ds_r ); 
              }
            }
          }
        }


        unsigned int nCons;
        double weightednCons;
        vector<double> error(3,0.);
        vector<double> jetVertex = CalculateVertex( recoCHSJet_constVertex_x->at(matchidxjet), recoCHSJet_constVertex_y->at(matchidxjet), recoCHSJet_constVertex_z->at(matchidxjet), recoCHSJet_const_pt->at(matchidxjet), recoCHSJet_const_charge->at(matchidxjet), recoCHSJet_const_closestVertex_d->at(matchidxjet), nCons, weightednCons, error );
       



        double tr_dx = jetVertex.at(0) - GluinoDecayVertex.at( (ipoi < 2 ) ? 0 : 1 ).at(0);
        double tr_dy = jetVertex.at(1) - GluinoDecayVertex.at( (ipoi < 2 ) ? 0 : 1 ).at(1);
        double tr_dz = jetVertex.at(2) - GluinoDecayVertex.at( (ipoi < 2 ) ? 0 : 1 ).at(2);
        /* 
        std::cout << "True Gluino Decay Vertices : " << endl;
        cout << GluinoDecayVertex.at(0).at(0) << " " << GluinoDecayVertex.at(0).at(1) << " " << GluinoDecayVertex.at(0).at(2) << endl;
        cout << GluinoDecayVertex.at(1).at(0) << " " << GluinoDecayVertex.at(1).at(1) << " " << GluinoDecayVertex.at(1).at(2) << endl;
        cout << "this matched jet vertex : " << endl;
        cout << jetVertex.at(0) << " " << jetVertex.at(1) << " " << jetVertex.at(2) << endl;
        */
        if( recoCHSJet_assignedSV.at(matchidxjet) >= 0 ) {

          double ds = sqrt( tr_dx*tr_dx + tr_dy*tr_dy + tr_dz*tr_dz );
          double dr = sqrt( tr_dx*tr_dx + tr_dy*tr_dy );
          if( ds < dr ) {
            cout << "THIS CAN't BE: " << tr_dx << " " << tr_dy << " " << tr_dz << " and " << ds << " vs " << dr << "??" << endl;
            getchar();
          } 
          _histograms1D.at("matchedCHSJets_DeltaS_TrueCalc").Fill( 10.*sqrt( tr_dx*tr_dx + tr_dy*tr_dy + tr_dz*tr_dz ) );   
          _histograms1D.at("matchedCHSJets_DeltaR_TrueCalc").Fill( 10.*sqrt( tr_dx*tr_dx + tr_dy*tr_dy ) );   
          _histograms1D.at("matchedCHSJets_DeltaX_TrueCalc").Fill( 10.*tr_dx );
          _histograms1D.at("matchedCHSJets_DeltaY_TrueCalc").Fill( 10.*tr_dy );
          _histograms1D.at("matchedCHSJets_DeltaZ_TrueCalc").Fill( 10.*tr_dz );
        }




        for( unsigned int iConst = 0; iConst < recoCHSJet_const_fromPV->at(matchidxjet).size(); ++iConst ) {
          if( recoCHSJet_const_charge->at(matchidxjet).at(iConst) == 0 ) continue; 
          if( recoCHSJet_const_fromPV->at(matchidxjet).at(iConst) == 0 ) continue;
          double dx = jetVertex.at(0) - recoCHSJet_constVertexRef_x->at(matchidxjet).at(iConst);
          double dy = jetVertex.at(1) - recoCHSJet_constVertexRef_y->at(matchidxjet).at(iConst);
          double dz = jetVertex.at(2) - recoCHSJet_constVertexRef_z->at(matchidxjet).at(iConst);
          double dr = sqrt( dx*dx + dy*dy );
          double ds = sqrt( dr*dr + dz*dz );

          if( recoCHSJet_assignedPV.at(matchidxjet) >= 0 ) {
            _histograms1D.at("matchedCHSJets_PVAss_dxAssCalc").Fill( 10.*dx );
            _histograms1D.at("matchedCHSJets_PVAss_dyAssCalc").Fill( 10.*dy );
            _histograms1D.at("matchedCHSJets_PVAss_dzAssCalc").Fill( 10.*dz );
            _histograms1D.at("matchedCHSJets_PVAss_drAssCalc").Fill( 10.*dr );
            _histograms1D.at("matchedCHSJets_PVAss_dsAssCalc").Fill( 10.*ds );
          }
          if( recoCHSJet_assignedSV.at(matchidxjet) >= 0 ) {
            // jet is associated to a SV
            // jet constituent is still from PV
            // what's the distance between the PCA to the PV and the PV?
            // what's the distance between the PCA to the SV and the SV?
            if( recoCHSJet_const_fromPV->at(matchidxjet).at(iConst) != 0 ) {
              TLorentzVector constP4; constP4.SetPtEtaPhiM( recoCHSJet_const_pt->at(matchidxjet).at(iConst), 
                                                          recoCHSJet_const_eta->at(matchidxjet).at(iConst),
                                                          recoCHSJet_const_phi->at(matchidxjet).at(iConst),
                                                          0. );
              double px = constP4.Px();
              double py = constP4.Py();
              double pz = constP4.Pz();
              double p = sqrt( px*px + py*py + pz*pz);
              double px_e = px/p;
              double py_e = py/p;
              double pz_e = pz/p;
              double or_x = recoCHSJet_constVertex_x->at(matchidxjet).at(iConst);
              double or_y = recoCHSJet_constVertex_y->at(matchidxjet).at(iConst);
              double or_z = recoCHSJet_constVertex_z->at(matchidxjet).at(iConst);
            
              double upv_x = vertex_x->at(0);
              double upv_y = vertex_y->at(0);
              double upv_z = vertex_z->at(0);

              double tmin = - ( px_e*(or_x-upv_x) + py_e*(or_y-upv_y) + pz_e*(or_z-upv_z) );
              double dx_min = upv_x - or_x - tmin*px_e;
              double dy_min = upv_y - or_y - tmin*py_e;
              double dz_min = upv_z - or_z - tmin*pz_e;
              if( sqrt( (upv_x-jetVertex.at(0))*(upv_x-jetVertex.at(0)) + (upv_y-jetVertex.at(1))*(upv_y-jetVertex.at(1)) + (upv_z-jetVertex.at(2))*(upv_z-jetVertex.at(2))) > 0.5 ) {
                _histograms1D.at("matchedCHSJet_ds_SVAss_drmin5_PCAPV").Fill( 10.* sqrt(dx_min*dx_min + dy_min*dy_min + dz_min*dz_min) );
                _histograms2D.at("matchedCHSJet_ds_vs_constPt_SVAss_drmin5_PCAPV").Fill( 10.*sqrt(dx_min*dx_min + dy_min*dy_min + dz_min*dz_min), constP4.P() );
                _histograms1D.at("matchedCHSJet_dr_SVAss_drmin5_PCAPV").Fill( 10.* sqrt(dx_min*dx_min + dy_min*dy_min) );
                _histograms1D.at("matchedCHSJet_dx_SVAss_drmin5_PCAPV").Fill( 10.* dx_min );
                _histograms1D.at("matchedCHSJet_dy_SVAss_drmin5_PCAPV").Fill( 10.* dy_min );
                _histograms1D.at("matchedCHSJet_dz_SVAss_drmin5_PCAPV").Fill( 10.* dz_min );
              }
            }
            _histograms1D.at("matchedCHSJets_SVAss_dxAssCalc").Fill( 10.*dx );
            _histograms1D.at("matchedCHSJets_SVAss_dyAssCalc").Fill( 10.*dy );
            _histograms1D.at("matchedCHSJets_SVAss_dzAssCalc").Fill( 10.*dz );
            _histograms1D.at("matchedCHSJets_SVAss_drAssCalc").Fill( 10.*dr );
            _histograms1D.at("matchedCHSJets_SVAss_dsAssCalc").Fill( 10.*ds );
          }
        }
      }
      
      drmin = 10000.;
      matchidxjet = -1;
      for( unsigned int iJet = 0; iJet < recoNoCHSJet_pt->size(); ++iJet ) {
        //if( recoNoCHSJet_pt->at(iJet) < 30. ) continue;
        double deta = quark_eta - recoNoCHSJet_eta->at(iJet);
        double dphi = quark_phi - recoNoCHSJet_phi->at(iJet);
        if( dphi > M_PI ) dphi = 2*M_PI-dphi;
        double dr = sqrt( deta*deta + dphi*dphi );
        if( dr < drmin ) { drmin = dr; matchidxjet = iJet; matchNoCHSID = iJet;}
      }
      if( drmin < 0.2 ) {
        _cutFlow.at("matchedNoCHS") += 1;
        if( recoNoCHSJet_assignedPV.at(matchidxjet) >= 0 && recoNoCHSJet_assignedSV.at(matchidxjet) >= 0 ) {
          cout << "Jet assigned to PV and SV - this should not happen" << endl;
          getchar();
        }
        else if( recoNoCHSJet_assignedPV.at(matchidxjet) >= 0 ) _histograms1D.at("allMatchedNoCHSJets_PVSV").Fill( 1 );
        else if( recoNoCHSJet_assignedSV.at(matchidxjet) >= 0 ) _histograms1D.at("allMatchedNoCHSJets_PVSV").Fill( 2 );
        else _histograms1D.at("allMatchedNoCHSJets_PVSV").Fill( 0 );
        _histograms1D.at("nochsMatchedQOI_pt").Fill( qVec.Pt() );
        if( qVec.Pt() > 50 ) {
          _histograms1D.at("nochsMatchedQOI_eta").Fill( qVec.Eta() );
          _histograms1D.at("nochsMatchedQOI_drpvsv").Fill( ds_r ); 
          _histograms1D.at("nochsMatchedQOI_dpt").Fill( (qVec.Pt() - recoNoCHSJet_pt->at(matchidxjet))/qVec.Pt() );
        }
        for( unsigned int iConst = 0; iConst <  recoNoCHSJet_const_fromPV->at(matchidxjet).size(); ++iConst ) {
          if( recoNoCHSJet_const_charge->at(matchidxjet).at(iConst) == 0 ) continue; 
          _histograms1D.at("nochsMatched_fromPV").Fill( recoNoCHSJet_const_fromPV->at(matchidxjet).at(iConst));
          _histograms2D.at("nochsMatched_fromPV_drpvsv").Fill( recoNoCHSJet_const_fromPV->at(matchidxjet).at(iConst), ds_r );
        }
        
        if( recoNoCHSJet_pt->at(matchidxjet) > 30. && recoNoCHSJet_assignedSV.at(matchidxjet) >= 0 ) {
          int thisJetMatchedIdx = recoNoCHSJet_assignedSV.at(matchidxjet);
          for( unsigned int iJet2 = 0; iJet2 < recoNoCHSJet_pt->size(); ++iJet2 ) {
            if( iJet2 == matchidxjet ) continue;
            if( recoNoCHSJet_assignedSV.at(iJet2) == thisJetMatchedIdx && recoNoCHSJet_pt->at(iJet2) > 30. ) {
              _histograms1D.at("nochsMatchedSJQOI_pt").Fill( qVec.Pt() );
              if( qVec.Pt() > 50 ) {
                _histograms1D.at("nochsMatchedSJQOI_eta").Fill( qVec.Eta() );
                _histograms1D.at("nochsMatchedSJQOI_drpvsv").Fill( ds_r ); 
              }
            }
          }
        }
      
        unsigned int nCons;
        double weightednCons;
        vector<double> error(3,0.);
        vector<double> jetVertex = CalculateVertex( recoNoCHSJet_constVertex_x->at(matchidxjet), recoNoCHSJet_constVertex_y->at(matchidxjet), recoNoCHSJet_constVertex_z->at(matchidxjet), recoNoCHSJet_const_pt->at(matchidxjet), recoNoCHSJet_const_charge->at(matchidxjet), recoNoCHSJet_const_closestVertex_d->at(matchidxjet), nCons, weightednCons, error );
        
        
        double tr_dx = jetVertex.at(0) - GluinoDecayVertex.at( (ipoi < 2 ) ? 0 : 1 ).at(0);
        double tr_dy = jetVertex.at(1) - GluinoDecayVertex.at( (ipoi < 2 ) ? 0 : 1 ).at(1);
        double tr_dz = jetVertex.at(2) - GluinoDecayVertex.at( (ipoi < 2 ) ? 0 : 1 ).at(2);

        if( recoNoCHSJet_assignedSV.at(matchidxjet) >= 0 ) {
          _histograms1D.at("matchedNoCHSJets_DeltaS_TrueCalc").Fill( 10.*sqrt( tr_dx*tr_dx + tr_dy*tr_dy + tr_dz*tr_dz ) );   
          _histograms1D.at("matchedNoCHSJets_DeltaR_TrueCalc").Fill( 10.*sqrt( tr_dx*tr_dx + tr_dy*tr_dy ) );   
          _histograms1D.at("matchedNoCHSJets_DeltaX_TrueCalc").Fill( 10.*tr_dx );
          _histograms1D.at("matchedNoCHSJets_DeltaY_TrueCalc").Fill( 10.*tr_dy );
          _histograms1D.at("matchedNoCHSJets_DeltaZ_TrueCalc").Fill( 10.*tr_dz );
        }
        
        
        
        for( unsigned int iConst = 0; iConst < recoNoCHSJet_const_fromPV->at(matchidxjet).size(); ++iConst ) {
          if( recoNoCHSJet_const_charge->at(matchidxjet).at(iConst) == 0 ) continue; 
          
          // consider only those matched to a PV
          if( recoNoCHSJet_assignedPV.at(matchidxjet) >= 0 ) {
            if(  fabs(jetVertex.at(0) - vertex_x->at(recoNoCHSJet_assignedPV.at(matchidxjet))) > 1.e-10
              || fabs(jetVertex.at(1) - vertex_y->at(recoNoCHSJet_assignedPV.at(matchidxjet))) > 1.e-10 
              || fabs(jetVertex.at(2) - vertex_z->at(recoNoCHSJet_assignedPV.at(matchidxjet))) > 1.e-10 ) {
              cout << "VERTEX POSITION MATCH MISMATCH!" << endl;
              getchar();
            }

            if( 10.*sqrt( tr_dx*tr_dx + tr_dy*tr_dy ) > 5. ) {
              // take only jets associated to a clearly off PV
              // calculate the distance at the point of closest approach to this PV 
              // and compare it to the distance at the point of closest approach to the actual Vertex (if there is a SV that more or less matches the actual vertex).
              vector<double> refPoint;
              refPoint.push_back( recoNoCHSJet_constVertex_x->at(matchidxjet).at(iConst) );
              refPoint.push_back( recoNoCHSJet_constVertex_y->at(matchidxjet).at(iConst) );
              refPoint.push_back( recoNoCHSJet_constVertex_z->at(matchidxjet).at(iConst) );
              TLorentzVector constP4; constP4.SetPtEtaPhiM( recoNoCHSJet_const_pt->at(matchidxjet).at(iConst),
                                                            recoNoCHSJet_const_eta->at(matchidxjet).at(iConst),
                                                            recoNoCHSJet_const_phi->at(matchidxjet).at(iConst),
                                                            0. );

              vector<double> momentum;
              momentum.push_back( constP4.Px() ); 
              momentum.push_back( constP4.Py() );
              momentum.push_back( constP4.Pz() );
              vector<double> genVertex = GluinoDecayVertex.at( (ipoi < 2 ) ? 0 : 1 );

              vector<double> pcaPV = CalculatePCA( &refPoint, &momentum, &jetVertex );
              vector<double> pcaSV = CalculatePCA( &refPoint, &momentum, &genVertex );

              if( recoNoCHSJet_const_fromPV->at(matchidxjet).at(iConst) != 0 ) {
                _histograms1D.at( "NoCHSJets_FromPV_ds_PV" ).Fill( 10.*sqrt( pcaPV.at(3)*pcaPV.at(3) + pcaPV.at(4)*pcaPV.at(4) + pcaPV.at(5)*pcaPV.at(5) ) );
                _histograms1D.at( "NoCHSJets_FromPV_ds_trueSV" ).Fill( 10.*sqrt( pcaSV.at(3)*pcaSV.at(3) + pcaSV.at(4)*pcaSV.at(4) + pcaSV.at(5)*pcaSV.at(5) ) );
                _histograms2D.at( "NoCHSJets_FromPV_ds_PV_vs_ds_trueSV").Fill( 10.*sqrt( pcaPV.at(3)*pcaPV.at(3) + pcaPV.at(4)*pcaPV.at(4) + pcaPV.at(5)*pcaPV.at(5) ),
                                                                               10.*sqrt( pcaSV.at(3)*pcaSV.at(3) + pcaSV.at(4)*pcaSV.at(4) + pcaSV.at(5)*pcaSV.at(5) ) );
              }
              else {
                _histograms1D.at( "NoCHSJets_NotFromPV_ds_PV" ).Fill( 10.*sqrt( pcaPV.at(3)*pcaPV.at(3) + pcaPV.at(4)*pcaPV.at(4) + pcaPV.at(5)*pcaPV.at(5) ) );
                _histograms1D.at( "NoCHSJets_NotFromPV_ds_trueSV" ).Fill( 10.*sqrt( pcaSV.at(3)*pcaSV.at(3) + pcaSV.at(4)*pcaSV.at(4) + pcaSV.at(5)*pcaSV.at(5) ) );
                _histograms2D.at( "NoCHSJets_NotFromPV_ds_PV_vs_ds_trueSV").Fill( 10.*sqrt( pcaPV.at(3)*pcaPV.at(3) + pcaPV.at(4)*pcaPV.at(4) + pcaPV.at(5)*pcaPV.at(5) ),
                                                                               10.*sqrt( pcaSV.at(3)*pcaSV.at(3) + pcaSV.at(4)*pcaSV.at(4) + pcaSV.at(5)*pcaSV.at(5) ) );
              }
              
            }
           
          }

          if( recoNoCHSJet_const_fromPV->at(matchidxjet).at(iConst) == 0 ) continue;
          double dx = jetVertex.at(0) - recoNoCHSJet_constVertexRef_x->at(matchidxjet).at(iConst);
          double dy = jetVertex.at(1) - recoNoCHSJet_constVertexRef_y->at(matchidxjet).at(iConst);
          double dz = jetVertex.at(2) - recoNoCHSJet_constVertexRef_z->at(matchidxjet).at(iConst);
          double dr = sqrt( dx*dx + dy*dy );
          double ds = sqrt( dr*dr + dz*dz );

          if( recoNoCHSJet_assignedPV.at(matchidxjet) >= 0 ) {
            _histograms1D.at("matchedNoCHSJets_PVAss_dxAssCalc").Fill( 10.*dx );
            _histograms1D.at("matchedNoCHSJets_PVAss_dyAssCalc").Fill( 10.*dy );
            _histograms1D.at("matchedNoCHSJets_PVAss_dzAssCalc").Fill( 10.*dz );
            _histograms1D.at("matchedNoCHSJets_PVAss_drAssCalc").Fill( 10.*dr );
            _histograms1D.at("matchedNoCHSJets_PVAss_dsAssCalc").Fill( 10.*ds );
          }
          if( recoNoCHSJet_assignedSV.at(matchidxjet) >= 0 ) {
            _histograms1D.at("matchedNoCHSJets_SVAss_dxAssCalc").Fill( 10.*dx );
            _histograms1D.at("matchedNoCHSJets_SVAss_dyAssCalc").Fill( 10.*dy );
            _histograms1D.at("matchedNoCHSJets_SVAss_dzAssCalc").Fill( 10.*dz );
            _histograms1D.at("matchedNoCHSJets_SVAss_drAssCalc").Fill( 10.*dr );
            _histograms1D.at("matchedNoCHSJets_SVAss_dsAssCalc").Fill( 10.*ds );
          }
        }
      
      
      
      }
      if( matchCHSID != -1 && matchNoCHSID != -1 ) {
        double dptCHS   = (qVec.Pt() - recoCHSJet_pt->at(matchCHSID))/qVec.Pt();
        double dptNoCHS = (qVec.Pt() - recoNoCHSJet_pt->at(matchNoCHSID))/qVec.Pt();
        _histograms1D.at("bothMatched_chsMatchedQOI_dpt").Fill( dptCHS );
        _histograms1D.at("bothMatched_nochsMatchedQOI_dpt").Fill( dptNoCHS );
      }
    
    }

    return;
}

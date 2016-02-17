#include "LLGAnalysis.h"
#include "TLorentzVector.h"

void LLGAnalysis::SetupSignalRegionTruthAnalysis() {

    // setup the cutflow
    _cutFlow.insert(pair<string,int>("0_NoCut", 0) );
    _cutFlow.insert(pair<string,int>("1_Trigger", 0) );
    _cutFlow.insert(pair<string,int>("2_MuonVeto", 0) );
    _cutFlow.insert(pair<string,int>("3_ElectronVeto", 0) );
    _cutFlow.insert(pair<string,int>("4_HasPVJet", 0) );
    _cutFlow.insert(pair<string,int>("5_HasSV20", 0) );
    _cutFlow.insert(pair<string,int>("6_MET", 0) );
    _cutFlow.insert(pair<string,int>("6a_SVJetPT", 0 ) );
    _cutFlow.insert(pair<string,int>("7_DiJetMass", 0 ) );
    _cutFlow.insert(pair<string,int>("8_SVPVDistance", 0) );
    _cutFlow.insert(pair<string,int>("9_SVResolution", 0 ) );

    // and the histograms 
    makeHist( "nBjetAtSV", 5, -0.5, 4.5, "Number of b-jets associated to SV", "Number of SV" );
    makeHist( "mJJSV", 100, 0., 500., "DiJet mass at SV", "Number of Jet Pairs" );
    makeHist( "nJetsSV", 7, -0.5, 6.5, "Number of Jets associated to SV", "Number of SV" );
    makeHist( "SVJet1Pt", 50, 0., 500., "SV Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "SVJet2Pt", 50, 0., 500., "SV 2^{nd} Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "SVJet3Pt", 50, 0., 500., "SV 3^{rd} Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "SVJet4Pt", 50, 0., 500., "SV 4^{th} Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "PVJet1Pt", 50, 0., 1000., "PV Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "PVJet2Pt", 50, 0., 500., "PV 2^{nd} Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "PVJet3Pt", 50, 0., 500., "PV 3^{rd} Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "PVJet4Pt", 50, 0., 500., "PV 4^{th} Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "nJetsTotal", 26, -0.5, 25.5, "Number of Jets with p_{T} > 30 GeV", "Number of Events" );
    makeHist( "BJet1Pt", 50, 0., 500., "Leading B-tagged Jet p_{T} [GeV]", "Number of Events" );
    makeHist( "BJet2Pt", 50, 0., 500., "Subleading B-tagged Jet p_{T} [GeV]", "Number of Events" );
    makeHist( "BJet3Pt", 50, 0., 500., "3^{rd} leading B-tagged Jet p_{T} [GeV]", "Number of Events" );
    makeHist( "BJet4Pt", 50, 0., 500., "4^{th} leading B-tagged Jet p_{T} [GeV]", "Number of Events" );
    makeHist( "JetLeptonDr", 50, 0., 6., "#DeltaR(jet, electron)", "# Jet-Electron Pairs" ); 
    makeHist( "nPVWithJet75", 10, -0.5, 9.5, "Number of PV with >= 1 Jet > 75 GeV", "# events" );
    makeHist( "nSVWith2Jets30", 10, -0.5, 9.5, "Number of SV with >= 2 Jets > 30 GeV", "# events" );
    makeHist( "mjjvsleadingjetpt", 100, 0., 500., 100, 0., 500., "DiJet Mass at SV [GeV]", "Leading Jet p_{T} at SV [GeV]", "Number of Events", "COLZ" );
    makeHist( "distancePVSV", 40, 0., 40., "Distance between leading PV and SV [mm]", "Number of PV-SV pairs" );
    makeHist( "jetOrigin", 10, -0.5, 9.5, "Gen Object", "Number of Events" );
    makeHist( "jetOrigin2D", 10, -0.5, 9.5, 10, -0.5, 9.5, "Gen Object", "Gen Object", "Number of SV", "COLZ" );
    makeHist( "jetOriginDistance", 102, -1.5, 100.5, "Distance between Reco and Gen Vertex", "Number of Events" );
    makeHist( "jetOriginResolution", 10, -0.5, 9.5, 20, 0., 1.e-3, "Gen Object", "SV #sigma_{r}", "Number of SV", "COLZ" );
    makeHist( "jetOriginNconsidered", 10, -0.5, 9.5, 71, -0.5, 70.5, "GenObject", "Considered Jet Constitutents", "Number of SV jets", "COLZ" );
    makeHist( "jetOriginNormScore", 10, -0.5, 9.5, 20, 0., 10000, "GenObject", "Normalised Vertex Score", "Number of SV jets", "COLZ" );

    // set some bin labels
    _histograms1D.at("jetOrigin").GetXaxis()->SetBinLabel( 1, "unmachted" );
    _histograms1D.at("jetOrigin").GetXaxis()->SetBinLabel( 2, "match unknown" );
    _histograms1D.at("jetOrigin").GetXaxis()->SetBinLabel( 3, "light flavour" ); 
    _histograms1D.at("jetOrigin").GetXaxis()->SetBinLabel( 4, "charm" ); 
    _histograms1D.at("jetOrigin").GetXaxis()->SetBinLabel( 5, "bottom" ); 
    _histograms1D.at("jetOrigin").GetXaxis()->SetBinLabel( 6, "tau" ); 
    _histograms1D.at("jetOrigin").GetXaxis()->SetBinLabel( 7, "muon" ); 
    _histograms1D.at("jetOrigin").GetXaxis()->SetBinLabel( 8, "electron" ); 
    _histograms1D.at("jetOrigin").GetXaxis()->SetBinLabel( 9, "photon" ); 
    _histograms1D.at("jetOrigin").GetXaxis()->SetBinLabel(10, "SUSY" ); 
    
    _histograms2D.at("jetOrigin2D").GetXaxis()->SetBinLabel( 1, "unmachted" );
    _histograms2D.at("jetOrigin2D").GetXaxis()->SetBinLabel( 2, "match unknown" );
    _histograms2D.at("jetOrigin2D").GetXaxis()->SetBinLabel( 3, "light flavour" ); 
    _histograms2D.at("jetOrigin2D").GetXaxis()->SetBinLabel( 4, "charm" ); 
    _histograms2D.at("jetOrigin2D").GetXaxis()->SetBinLabel( 5, "bottom" ); 
    _histograms2D.at("jetOrigin2D").GetXaxis()->SetBinLabel( 6, "tau" ); 
    _histograms2D.at("jetOrigin2D").GetXaxis()->SetBinLabel( 7, "muon" ); 
    _histograms2D.at("jetOrigin2D").GetXaxis()->SetBinLabel( 8, "electron" ); 
    _histograms2D.at("jetOrigin2D").GetXaxis()->SetBinLabel( 9, "photon" ); 
    _histograms2D.at("jetOrigin2D").GetXaxis()->SetBinLabel(10, "SUSY" ); 
    _histograms2D.at("jetOrigin2D").GetYaxis()->SetBinLabel( 1, "unmachted" );
    _histograms2D.at("jetOrigin2D").GetYaxis()->SetBinLabel( 2, "match unknown" );
    _histograms2D.at("jetOrigin2D").GetYaxis()->SetBinLabel( 3, "light flavour" ); 
    _histograms2D.at("jetOrigin2D").GetYaxis()->SetBinLabel( 4, "charm" ); 
    _histograms2D.at("jetOrigin2D").GetYaxis()->SetBinLabel( 5, "bottom" ); 
    _histograms2D.at("jetOrigin2D").GetYaxis()->SetBinLabel( 6, "tau" ); 
    _histograms2D.at("jetOrigin2D").GetYaxis()->SetBinLabel( 7, "muon" ); 
    _histograms2D.at("jetOrigin2D").GetYaxis()->SetBinLabel( 8, "electron" ); 
    _histograms2D.at("jetOrigin2D").GetYaxis()->SetBinLabel( 9, "photon" ); 
    _histograms2D.at("jetOrigin2D").GetYaxis()->SetBinLabel(10, "SUSY" ); 

    for( double jptpv = 30.; jptpv <= 401.; jptpv += 5. ) {
      std::vector<double> yield;
      //for( double jptsv = 40.; jptsv <= 401.; jptsv += 10. ) {
      //for( double jptsv = 15.; jptsv <= 201.; jptsv += 5. ) {
      for( double jptsv = 0.; jptsv <= 50.; jptsv += 2. ) {
        yield.push_back( 0. );
      }
      _yields2DOptimisation.push_back( yield );
    }
    return;
}

void LLGAnalysis::SignalRegionTruthAnalysisSelection() {


    int IDXFIRSTCUT = -1;
    int IDXSECONDCUT = -1;


    int leadingPV = -1;
    double leadingVertexPt = 0.;

    _cutFlow.at("0_NoCut") += 1;

    
    
    bool passTrigger = false;
    for( unsigned int iTrig = 0; iTrig < triggerNames->size(); ++iTrig ) {
        if( triggerNames->at(iTrig) == "HLT_PFMET170_NoiseCleaned_v1" && triggerBits->at(iTrig) == 1 ) passTrigger = true;
    }

    if( !passTrigger ) return; 
    _cutFlow.at("1_Trigger") += 1;

    // lepton veto:
    if( vetoMuons.size() > 0 ) return; 
    _cutFlow.at("2_MuonVeto") += 1;
        
    
    if( vetoElectrons.size() > 0 ) return;
    _cutFlow.at("3_ElectronVeto") += 1;


    // now assign jets to the vertices:
    vector<int> nJetsToPV( vertex_x->size(), 0 );
    vector<int> nJetsToSV( secVertex_x->size(), 0 );
    vector<vector<int> > idJetsToPV;
    vector<vector<int> > idJetsToSV;
    
    for( unsigned int iVtx = 0; iVtx < vertex_x->size(); ++iVtx ) {
        vector<int> idx;
        idJetsToPV.push_back( idx );
    }
    for( unsigned int iVtx = 0; iVtx < secVertex_x->size(); ++iVtx ) {
        vector<int> idx;
        idJetsToSV.push_back( idx );
    }
    

    for( unsigned int iselJet = 0; iselJet < selectedJets.size(); ++iselJet ) {
        int iJet = selectedJets.at(iselJet);
        //calculate jet vertex position:
        //vector<double> error(3,0.);
        vector<double> position(3,0.);
        position.at(0) = recoJet_vertex_x->at(iJet);
        position.at(1) = recoJet_vertex_y->at(iJet);
        position.at(2) = recoJet_vertex_z->at(iJet);
        //vector<double> position = CalculateVertex( recoJet_constVertex_x->at(iJet), recoJet_constVertex_y->at(iJet), recoJet_constVertex_z->at(iJet), recoJet_const_pt->at(iJet), recoJet_const_charge->at(iJet), recoJet_const_closestVertex_d->at(iJet), nCons, weightednCons, error );
        int nMatch = 0;
        
        //std::cout << "this is jet # " << iJet << " at position " << position.at(0) << " " << position.at(1) << " " << position.at(2) << std::endl; 
        for( unsigned int iVtx = 0; iVtx < vertex_x->size(); ++iVtx ) {
            if( fabs(position.at(0) - vertex_x->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(1) - vertex_y->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(2) - vertex_z->at(iVtx) ) < 1.e-10 ) {
                nJetsToPV.at(iVtx) += 1;
                idJetsToPV.at(iVtx).push_back( iJet );
                if( recoJet_pt->at(iJet).at(SYSJET) > leadingVertexPt ) {
                  leadingVertexPt = recoJet_pt->at(iJet).at(SYSJET);
                  leadingPV = iVtx;
                }
                nMatch += 1;
            }
        }
        for( unsigned int iVtx = 0; iVtx < secVertex_x->size(); ++iVtx ) {
            //std::cout << "checking SV with " << secVertex_x->at(iVtx) << " " << secVertex_y->at(iVtx) << " " << secVertex_z->at(iVtx) << std::endl;
            if( fabs(position.at(0) - secVertex_x->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(1) - secVertex_y->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(2) - secVertex_z->at(iVtx) ) < 1.e-10 ) {
                
                nJetsToSV.at(iVtx) += 1;
                if( recoJet_pt->at(iJet).at(SYSJET) > JET_PT_CUT_SV ) { 
                    idJetsToSV.at(iVtx).push_back( iJet );
                }
                nMatch += 1;

            }
        }
        if( nMatch > 1 ) {
            cout << "WARNING! ASSOCIATED JET TO MORE THAN 1 VERTEX ?!" << endl;
        }
    }

    // now count the number of vertices with jets:
    vector<int> PVWithJet;
    vector<int> SVWithJets;
    vector<int> SVWith2Jets;

 
    double allPVLeadingJetPt = -1.;
    double allSVLeadingJetPt = -1.;
    double allSVLeadingmJJ = -1.;
    for( unsigned int iPV = 0; iPV < vertex_x -> size(); ++iPV ) {
      bool hasJetPV = false;
      double leadingJetPt = 0.;
      for( unsigned int iiJet = 0; iiJet < idJetsToPV.at(iPV).size(); ++iiJet ) {
          int iJet = idJetsToPV.at(iPV).at(iiJet);
          if( recoJet_pt->at(iJet).at(SYSJET) > JET_PT_CUT_PV ) hasJetPV = true;
          if( recoJet_pt->at(iJet).at(SYSJET) > leadingJetPt ) leadingJetPt = recoJet_pt->at(iJet).at(SYSJET);
          if( leadingJetPt > allPVLeadingJetPt ) allPVLeadingJetPt = leadingJetPt;
      }
      if( hasJetPV ) {
        PVWithJet.push_back( iPV );
      }
    }
    
    if( PVWithJet.size() > 0 ) {
    int counter = -1;
    for( double jptpv = 30.; jptpv <= 401.; jptpv += 5. ) {
      counter++;
      if( allPVLeadingJetPt > jptpv && allPVLeadingJetPt <= jptpv + 5. ) break;
    }
    IDXFIRSTCUT = counter;
    
    for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
        if( idJetsToSV.at(iSV).size() > 0 ) SVWithJets.push_back( iSV );
        if( idJetsToSV.at(iSV).size() >= 2 ) SVWith2Jets.push_back( iSV );
    }
    }

    if( SVWith2Jets.size() > 0 ) _outputTree->Fill();


    // do the n-1 plots here:
    // 1st: Leading Lepton pT from PV:
    if( SVWith2Jets.size() > 0 && met->at(SYSMET) > MET_CUT ) {
      for( unsigned int iPV = 0; iPV < vertex_x -> size(); ++iPV ) {
        int idxLeadingJetPV = -1;
        double ptLeadingJetPV = -1.;
        int idxSubLeadingJetPV = -1;
        double ptSubLeadingJetPV = -1;
        int idxThirdLeadingJetPV = -1;
        double ptThirdLeadingJetPV = -1;
        int idxFourthLeadingJetPV = -1;
        double ptFourthLeadingJetPV = -1;
        bool hasJetPV = false;
        for( unsigned int iiJet = 0; iiJet < idJetsToPV.at(iPV).size(); ++iiJet ) {
            int iJet = idJetsToPV.at(iPV).at(iiJet);
            if( recoJet_pt->at(iJet).at(SYSJET) > JET_PT_CUT_PV ) hasJetPV = true;
            
            if( recoJet_pt->at(iJet).at(SYSJET) > ptLeadingJetPV ) {
                idxFourthLeadingJetPV = idxThirdLeadingJetPV;
                ptFourthLeadingJetPV = ptThirdLeadingJetPV;
                idxThirdLeadingJetPV = idxSubLeadingJetPV;
                ptThirdLeadingJetPV = ptSubLeadingJetPV;
                idxSubLeadingJetPV = idxLeadingJetPV;
                ptSubLeadingJetPV = ptLeadingJetPV;
                idxLeadingJetPV = iJet;
                ptLeadingJetPV = recoJet_pt->at(iJet).at(SYSJET);
            }
            else if( recoJet_pt->at(iJet).at(SYSJET) > ptSubLeadingJetPV ) {
                idxFourthLeadingJetPV = idxThirdLeadingJetPV;
                ptFourthLeadingJetPV = ptThirdLeadingJetPV;
                idxThirdLeadingJetPV = idxSubLeadingJetPV;
                ptThirdLeadingJetPV = ptSubLeadingJetPV;
                idxSubLeadingJetPV = iJet;
                ptSubLeadingJetPV = recoJet_pt->at(iJet).at(SYSJET);
            }
            else if( recoJet_pt->at(iJet).at(SYSJET) > ptThirdLeadingJetPV ) {
                idxFourthLeadingJetPV = idxThirdLeadingJetPV;
                ptFourthLeadingJetPV = ptThirdLeadingJetPV;
                idxThirdLeadingJetPV = iJet;
                ptThirdLeadingJetPV = recoJet_pt->at(iJet).at(SYSJET);
            }
            else if( recoJet_pt->at(iJet).at(SYSJET) > ptThirdLeadingJetPV ) {
                idxFourthLeadingJetPV = iJet;
                ptFourthLeadingJetPV = recoJet_pt->at(iJet).at(SYSJET);
            }
        }
        _histograms1D.at("PVJet1Pt").Fill( ptLeadingJetPV, evtWeight );
        _histograms1D.at("PVJet2Pt").Fill( ptSubLeadingJetPV, evtWeight );
        _histograms1D.at("PVJet3Pt").Fill( ptThirdLeadingJetPV, evtWeight );
        _histograms1D.at("PVJet4Pt").Fill( ptFourthLeadingJetPV, evtWeight );
      }
      _histograms1D.at("nPVWithJet75").Fill( PVWithJet.size(), evtWeight ); 
    }
    
    if( PVWithJet.size() >= 1 && met->at(SYSMET) > MET_CUT ) {
      
      for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
        _histograms1D.at("nJetsSV").Fill( idJetsToSV.at(iSV).size(), evtWeight );
      }
      _histograms1D.at("nSVWith2Jets30").Fill( SVWith2Jets.size(), evtWeight );
    }
    

    // and run the selection:
    if( PVWithJet.size() >= 1 ) {
        _cutFlow.at("4_HasPVJet") += 1;
        
        if( SVWith2Jets.size() > 0 ) {
            _cutFlow.at("5_HasSV20") += 1;
                
            vector<double> allDistances;
            double maxDist = 0.;
            for( unsigned int iPV = 0; iPV < PVWithJet.size(); ++iPV ) {
                double thispv_x = vertex_x->at(PVWithJet.at(iPV));
                double thispv_y = vertex_y->at(PVWithJet.at(iPV));
                double thispv_z = vertex_z->at(PVWithJet.at(iPV));
                for( unsigned int iSV = 0; iSV < SVWithJets.size(); ++iSV ) {
                    double thissv_x = secVertex_x->at(SVWithJets.at(iSV));
                    double thissv_y = secVertex_y->at(SVWithJets.at(iSV));
                    double thissv_z = secVertex_z->at(SVWithJets.at(iSV));
                    double dx = thissv_x - thispv_x;
                    double dy = thissv_y - thispv_y;
                    double dz = thissv_z - thispv_z;
                    double dist = 10.*sqrt(dx*dx + dy*dy + dz*dz );
                    allDistances.push_back( dist );
                    if( dist > maxDist ) maxDist = dist;
                }
            }
            if( met->at(SYSMET) > MET_CUT ) {
                _cutFlow.at("6_MET") += 1;
                
                for( unsigned int iDist = 0; iDist < allDistances.size(); ++iDist ) {
                  _histograms1D.at("distancePVSV").Fill(allDistances.at(iDist), evtWeight );
                }

                bool hasDiJetPair100 = false;
                bool hasSVJetCut = false;
                bool hasSVFarAway = false;
                bool hasSVPassDR = false;

                vector<double> svDistanceToPV;
                for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
                  double dx = secVertex_x->at(iSV) - vertex_x->at(leadingPV);
                  double dy = secVertex_y->at(iSV) - vertex_y->at(leadingPV);
                  double dz = secVertex_z->at(iSV) - vertex_z->at(leadingPV);
                  double dr = 10.*sqrt( dx*dx + dy*dy );
                  svDistanceToPV.push_back( dr );
                  if( idJetsToSV.at(iSV).size() <= 1 ) continue;
                  _histograms1D.at("nJetsSV").Fill( idJetsToSV.at(iSV).size(), evtWeight );
                  
                  int idxLeadingJet = -1;
                  double ptLeadingJet = -1.;
                  int idxSubLeadingJet = -1;
                  double ptSubLeadingJet = -1;
                  int idxThirdLeadingJet = -1;
                  double ptThirdLeadingJet = -1;
                  int idxFourthLeadingJet = -1;
                  double ptFourthLeadingJet = -1;

                  for( unsigned int iJToSV = 0; iJToSV < idJetsToSV.at(iSV).size(); ++iJToSV ) {
                    int jIdx = idJetsToSV.at(iSV).at(iJToSV);
                    if( recoJet_pt->at(jIdx).at(SYSJET) > ptLeadingJet ) {
                      idxFourthLeadingJet = idxThirdLeadingJet;
                      ptFourthLeadingJet = ptThirdLeadingJet;
                      idxThirdLeadingJet = idxSubLeadingJet;
                      ptThirdLeadingJet = ptSubLeadingJet;
                      idxSubLeadingJet = idxLeadingJet;
                      ptSubLeadingJet = ptLeadingJet;
                      idxLeadingJet = jIdx;
                      ptLeadingJet = recoJet_pt->at(jIdx).at(SYSJET);
                    }
                    else if ( recoJet_pt->at(jIdx).at(SYSJET) > ptSubLeadingJet ) {
                      idxFourthLeadingJet = idxThirdLeadingJet;
                      ptFourthLeadingJet = ptThirdLeadingJet;
                      idxThirdLeadingJet = idxSubLeadingJet;
                      ptThirdLeadingJet = ptSubLeadingJet;
                      idxSubLeadingJet = jIdx;
                      ptSubLeadingJet = recoJet_pt->at(jIdx).at(SYSJET);
                    }
                    else if ( recoJet_pt->at(jIdx).at(SYSJET) > ptThirdLeadingJet ) {
                      idxFourthLeadingJet = idxThirdLeadingJet;
                      ptFourthLeadingJet = ptThirdLeadingJet;
                      idxThirdLeadingJet = jIdx;
                      ptThirdLeadingJet = recoJet_pt->at(jIdx).at(SYSJET);
                    }
                    else if ( recoJet_pt->at(jIdx).at(SYSJET) > ptFourthLeadingJet ) {
                      idxFourthLeadingJet = jIdx;
                      ptFourthLeadingJet = recoJet_pt->at(jIdx).at(SYSJET);
                    }
                  }

                  _histograms1D.at("SVJet1Pt").Fill( ptLeadingJet, evtWeight );
                  _histograms1D.at("SVJet2Pt").Fill( ptSubLeadingJet, evtWeight );
                  _histograms1D.at("SVJet3Pt").Fill( ptThirdLeadingJet, evtWeight );
                  _histograms1D.at("SVJet4Pt").Fill( ptFourthLeadingJet, evtWeight );
                  TLorentzVector p4Jet1, p4Jet2;
                  p4Jet1.SetPtEtaPhiM( recoJet_pt->at(idxLeadingJet).at(SYSJET), recoJet_eta->at(idxLeadingJet), recoJet_phi->at(idxLeadingJet), 0. );
                  p4Jet2.SetPtEtaPhiM( recoJet_pt->at(idxSubLeadingJet).at(SYSJET), recoJet_eta->at(idxSubLeadingJet), recoJet_phi->at(idxSubLeadingJet), 0. );
                  TLorentzVector p4DiJet = p4Jet1 + p4Jet2;
                  _histograms1D.at("mJJSV").Fill( p4DiJet.M(), evtWeight );
                  
                  if( p4DiJet.M() > MJJ_CUT ) hasDiJetPair100 = true;
                  if( ptLeadingJet > LEADING_SV_JET_CUT ) hasSVJetCut = true;
                  if( p4DiJet.M() <= MJJ_CUT && svDistanceToPV.at(iSV) > 5. ) hasSVFarAway = true;
                  if( sqrt( secVertex_dx->at(iSV)*secVertex_dx->at(iSV) + secVertex_dy->at(iSV)*secVertex_dy->at(iSV) ) < 1.e-4 && p4DiJet.M() < MJJ_CUT && svDistanceToPV.at(iSV) > 5 ) hasSVPassDR = true; 
                  _histograms2D.at("mjjvsleadingjetpt").Fill( p4DiJet.M(), ptLeadingJet, evtWeight ); 
                  if( ptLeadingJet > allSVLeadingJetPt ) allSVLeadingJetPt = ptLeadingJet;
                  if( p4DiJet.M() > allSVLeadingmJJ ) allSVLeadingmJJ = p4DiJet.M();
                }
                int counter = -1;
                /*
                for( double jptsv = 40.; jptsv <= 401.; jptsv += 10. ) {
                  counter++;
                  if( (allSVLeadingJetPt < jptsv && allSVLeadingJetPt >= jptsv - 10.) || allSVLeadingJetPt < 40. ) break;
                }*/
                /*for( double jptsv = 15.; jptsv <= 201.; jptsv += 5. ) {
                  counter ++;
                  if( (allSVLeadingmJJ < jptsv && allSVLeadingmJJ >= jptsv - 5. ) || allSVLeadingmJJ < 15. ) break;
                }*/
                for( double jptsv = 0.; jptsv <= 50.; jptsv += 2. ) {
                  counter++;
                  if( maxDist >= jptsv && maxDist < jptsv+2. ) break;
                }
                IDXSECONDCUT = counter;
                for( unsigned int firstCut = IDXFIRSTCUT; firstCut >= 0; --firstCut ) {
                  //for( unsigned int secondCut = IDXSECONDCUT; secondCut < _yields2DOptimisation.at(firstCut).size(); ++secondCut ) {
                  for( unsigned int secondCut = IDXSECONDCUT; secondCut >=0; --secondCut ) {
                    //std::cout << "attempting to attac at " << firstCut << " and " << secondCut << endl;
                    _yields2DOptimisation.at(firstCut).at(secondCut) += evtWeight;
                    if( secondCut == 0 ) break;
                  }
                  if( firstCut == 0  ) break;
                }
                if( hasSVJetCut )  return;
                _cutFlow.at("6a_SVJetPT") += 1;
                if( hasDiJetPair100 ) return;
                _cutFlow.at("7_DiJetMass") += 1;
                //if( !hasSVFarAway ) return;
                //_cutFlow.at("8_SVPVDistance" )+= 1;
                //if( !hasSVPassDR ) return;
                //_cutFlow.at("9_SVResolution") += 1;
                TFile *fTruth = new TFile( GenFileName.c_str(), "OPEN");
                TTree *tTruth = (TTree*)fTruth->Get("output");
                int GenRunNumber, GenEventNumber, GenLumiBlock;
                vector<double> *mct_px = new vector<double>;
                vector<double> *mct_py = new vector<double>;
                vector<double> *mct_pz = new vector<double>;
                vector<double> *mct_vx = new vector<double>;
                vector<double> *mct_vy = new vector<double>;
                vector<double> *mct_vz = new vector<double>;
                vector<double> *mct_e = new vector<double>;
                vector<int> *mct_id = new vector<int>;
                vector<int> *mct_status = new vector<int>;
                vector<int> *mct_parent = new vector<int>;
                vector<vector<int> > *mct_daughters = new vector<vector<int> >;
                tTruth->SetBranchAddress("EventNumber", &GenEventNumber );
                tTruth->SetBranchAddress("LumiBlock", &GenLumiBlock );
                tTruth->SetBranchAddress("RunNumber", &GenRunNumber );
                tTruth->SetBranchAddress("Particle_px", &mct_px );
                tTruth->SetBranchAddress("Particle_py", &mct_py );
                tTruth->SetBranchAddress("Particle_pz", &mct_pz );
                tTruth->SetBranchAddress("Particle_vx", &mct_vx );
                tTruth->SetBranchAddress("Particle_vy", &mct_vy );
                tTruth->SetBranchAddress("Particle_vz", &mct_vz );
                tTruth->SetBranchAddress("Particle_E", &mct_e );
                tTruth->SetBranchAddress("Particle_ID", &mct_id );
                tTruth->SetBranchAddress("Particle_Status", &mct_status );
                tTruth->SetBranchAddress("Particle_Parent", &mct_parent );
                tTruth->SetBranchAddress("Particle_Daughter", &mct_daughters );
               
                bool FoundTruthEvent = false;
                for( unsigned int itruth = 0; itruth < tTruth->GetEntries(); ++itruth ) {
                  tTruth->GetEntry(itruth);
                  if( RunNumber == GenRunNumber && LumiBlock == GenLumiBlock && EventNumber == GenEventNumber ) {
                    FoundTruthEvent = true;
                    break;
                  }
                }
                if( !FoundTruthEvent ) {
                  cout << "DID NOT FIND TRUTH INFO " << endl;
                  getchar();
                }
                /*
                for( unsigned int ipart = 0; ipart < mct_px->size(); ++ipart ) {
                  if( mct_id->at(ipart) == 1000021 && mct_daughters->at(ipart).size() == 3 ) {
                    TLorentzVector lvg;
                    lvg.SetPxPyPzE( mct_px->at(ipart), mct_py->at(ipart), mct_pz->at(ipart), mct_e->at(ipart) );
                    TLorentzVector lvc1, lvc2, lvc3;
                    lvc1.SetPxPyPzE( mct_px->at( mct_daughters->at(ipart).at(0)),
                                     mct_py->at( mct_daughters->at(ipart).at(0)),
                                     mct_pz->at( mct_daughters->at(ipart).at(0)),
                                     mct_e->at( mct_daughters->at(ipart).at(0)) );
                    lvc2.SetPxPyPzE( mct_px->at( mct_daughters->at(ipart).at(1)),
                                     mct_py->at( mct_daughters->at(ipart).at(1)),
                                     mct_pz->at( mct_daughters->at(ipart).at(1)),
                                     mct_e->at( mct_daughters->at(ipart).at(1)) );
                    lvc3.SetPxPyPzE( mct_px->at( mct_daughters->at(ipart).at(2)),
                                     mct_py->at( mct_daughters->at(ipart).at(2)),
                                     mct_pz->at( mct_daughters->at(ipart).at(2)),
                                     mct_e->at( mct_daughters->at(ipart).at(2)) );
                    cout << "Gluino (" << lvg.Pt() << " " << lvg.Eta() << " " << lvg.Phi() << endl;
                    cout << "\t\t" << mct_id->at(mct_daughters->at(ipart).at(0)) << " ( " << lvc1.Pt() << " " << lvc1.Eta() << " " << lvc1.Phi() << endl;
                    cout << "\t\t" << mct_id->at(mct_daughters->at(ipart).at(1)) << " ( " << lvc2.Pt() << " " << lvc2.Eta() << " " << lvc2.Phi() << endl;
                    cout << "\t\t" << mct_id->at(mct_daughters->at(ipart).at(2)) << " ( " << lvc3.Pt() << " " << lvc3.Eta() << " " << lvc3.Phi() << endl;
                  }
                }
                */
                for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
                  std::vector<int> jetMatch(idJetsToSV.at(iSV).size(), -1 );
                  double SVDeltaR = sqrt( secVertex_dx->at(iSV)*secVertex_dx->at(iSV) + secVertex_dy->at(iSV)*secVertex_dz->at(iSV) ); 
                  if( jetMatch.size() != 2 ) continue;
                  vector<int> jetOrigins(2,-1);
                  for( unsigned int iijet = 0; iijet < idJetsToSV.at(iSV).size(); ++iijet ) {
                    int ijet = idJetsToSV.at(iSV).at(iijet);
                    cout << endl << "now running matching for jet " << ijet << " to sv " << iSV << endl;
                    double jetEta = recoJet_eta->at(ijet);
                    double jetPt = recoJet_pt->at(ijet).at(SYSJET);
                    double jetPhi = recoJet_phi->at(ijet);
                    double drmin = 10000.;
                    double dptmin = 10000.;
                    int idxmin = -1;
                    for( unsigned int ipart = 0; ipart < mct_px->size(); ++ipart ) {
                      if( fabs(mct_px->at(ipart)) < 1.e-10 && fabs(mct_py->at(ipart)) < 1.e-10 ) continue; 
                      //if( mct_parent->at(ipart) < 0 ) continue; 
                      if( abs(mct_id->at(ipart)) > 6 && abs(mct_id->at(ipart)) != 11 && abs(mct_id->at(ipart)) != 13 && abs(mct_id->at(ipart)) != 15 && mct_id->at(ipart) != 22 && mct_id->at(ipart) != 21 ) continue; 
                     
                      /*
                      if( mct_daughters->at(ipart).size() == 1 ) {
                        if( mct_daughters->at(ipart).at(0) > 0 ) {
                          if( mct_id->at(ipart) == mct_id->at(mct_daughters->at(ipart).at(0) ) ) continue;
                        }
                      }*/
                      
                      TLorentzVector vec;
                      vec.SetPxPyPzE( mct_px->at(ipart), mct_py->at(ipart), mct_pz->at(ipart), mct_e->at(ipart) );
                      
                      double tEta = vec.Eta();
                      double tPhi = vec.Phi();
                      double tPt = vec.Pt();
                      double deta = fabs( tEta - jetEta );
                      double dphi = fabs( tPhi - jetPhi );
                      double dpt = fabs( tPt - jetPt );
                      if( dphi > M_PI ) dphi = 2*M_PI - dphi;
                      double dr = sqrt( deta*deta + dphi*dphi );
                      if( dr < drmin && dpt < 0.8*jetPt ) {
                        drmin = dr;
                        dptmin = dpt;
                        idxmin = ipart;
                      }
                    }
                    double reco_mct_distance = -1.;
                    double reco_x = secVertex_x->at(iSV);
                    double reco_y = secVertex_y->at(iSV);
                    double reco_z = secVertex_z->at(iSV);
                    if( drmin < 0.4) { //&& dptmin < 0.5*jetPt ) {
                      jetMatch.at(iijet) = idxmin;
                      double gen_x = mct_vx->at(idxmin);
                      double gen_y = mct_vy->at(idxmin);
                      double gen_z = mct_vz->at(idxmin);
                      double dx = gen_x - reco_x;
                      double dy = gen_y - reco_y;
                      double dz = gen_z - reco_z;
                      reco_mct_distance = sqrt( dx*dx + dy*dy + dz*dz ); 
                      cout << endl;
                      cout << "matching summary for jet " << ijet << " matched to SV " << iSV << endl;
                      cout << "reco jet: " << jetPt << " " << jetEta << " " << jetPhi << "               ORIGIN: " << reco_x << " " << reco_y << " " << reco_z << endl;
                      TLorentzVector t4vec;
                      t4vec.SetPxPyPzE( mct_px->at(idxmin), mct_py->at(idxmin), mct_pz->at(idxmin), mct_e->at(idxmin) );
                       cout << "match  " << idxmin << "   : " << t4vec.Pt() << " " << t4vec.Eta() << " " << t4vec.Phi() << "    WITH ID: " << mct_id->at(idxmin) << "    and STATUS " << mct_status->at(idxmin) << "         ORIGIN: " << gen_x << " " << gen_y << " " << gen_z <<  endl;
                      cout << "\t origin summary: "<< endl;
                      
                      int mother = mct_parent->at(idxmin);
                      
                      int mid = mct_id->at(idxmin); //abs(mct_id->at(mother));
                      
                      bool isSUSY = false;
                      bool isCharm = false;
                      bool isBottom = false;
                      bool isTau = false;
                      while( mother >= 0 ) {
                        TLorentzVector motherVec;
                        motherVec.SetPxPyPzE( mct_px->at(mother), mct_py->at(mother), mct_pz->at(mother), mct_e->at(mother) );
                        int nDaugh = mct_daughters->at(mother).size();
                        cout << "\t\t idx " << mother << " " << mct_id->at(mother) << " " << "with mother " << mct_parent->at(mother) << " with " << mct_daughters->at(mother).size() << " children and status " << mct_status->at(mother) << " (" << motherVec.Pt() << " " << motherVec.Eta() << " " << motherVec.Phi() << ")  -> ";
                        for( unsigned int child = 0; child < mct_daughters->at(mother).size(); ++child ) {
                          if( mct_daughters->at(mother).at(child) < 0 ) continue;
                          TLorentzVector this4Vec;
                          this4Vec.SetPxPyPzE( mct_px->at( mct_daughters->at(mother).at(child) ), mct_py->at( mct_daughters->at(mother).at(child) ), mct_pz->at( mct_daughters->at(mother).at(child) ), mct_e->at( mct_daughters->at(mother).at(child) ) );
                          cout << " " << mct_id->at( mct_daughters->at(mother).at(child)  ) << "(" << this4Vec.Pt() << "," << this4Vec.Eta() << "," << this4Vec.Phi() << ")" << mct_daughters->at(mother).at(child) << " ";// << "(" << this4Vec.Pt() << " " << this4Vec.Eta() << " " << this4Vec.Phi() << ")";
                        }
                        cout << endl;
                        if( mct_id->at(mother) > 1000000 ) isSUSY = true;
                        if( abs(mct_id->at(mother)) == 5 ) isBottom = true;
                        if( abs(mct_id->at(mother)) == 4 ) isCharm = true;
                        if( abs(mct_id->at(mother)) == 15 ) isTau = true;
                        if( mother == mct_parent->at(mother) ) break;
                        //if( mct_parent->at(mother) >= 0 ) {
                          //if( mother == mct_parent->at(mct_parent->at(mother)) ) break;
                        //}
                        mother = mct_parent->at(mother); 
                      }
                      int origin = 1;
                      // 0 - unmatched
                      // 1 - matched unknown
                      // 2 - light flavour
                      // 3 - charm
                      // 4 - bottom
                      // 5 - tau 
                      // 6 - muon  
                      // 7 - electron
                      // 8 - photon
                      // 9 - SUSY
                      if( isSUSY ) origin = 9;
                      else if( abs(mid) == 15 || isTau ) origin = 5;
                      else if( abs(mid) == 5 || isBottom ) origin = 4;
                      else if( abs(mid) == 4 || isCharm ) origin = 3;
                      else {
                        if( (abs(mid) > 0 && abs(mid) < 4) || mid == 21 ) origin = 2;
                        else if( abs(mid) == 13 ) origin = 6;
                        else if( abs(mid) == 11 ) origin = 7;
                        else if( mid == 22 ) origin = 8;
                      }
                      cout << "THIS ONE WAS MATCHED TO " << origin << endl;
                      jetOrigins.at(iijet) = origin;
                      _histograms1D.at("jetOrigin").Fill( origin, evtWeight );
                      _histograms1D.at("jetOriginDistance").Fill( 10.*reco_mct_distance, evtWeight );         
                      //_histograms2D.at("jetOriginNormScore").Fill( origin, recoJet_vertex_score->at(ijet)/(float)recoJet_nConsidered->at(ijet), evtWeight );
                      _histograms2D.at("jetOriginNconsidered").Fill( origin, recoJet_nConsidered->at(ijet), evtWeight );
                      _histograms2D.at("jetOriginResolution").Fill( origin, SVDeltaR, evtWeight );
                      cout << endl << endl;
                    }
                    else {
                      std::cout << "NO MATCH FOUND" << std::endl;
                      jetOrigins.at(iijet) = 0;
                      _histograms1D.at("jetOrigin").Fill( 0., evtWeight );
                      _histograms1D.at("jetOriginDistance").Fill( -1., evtWeight );
                      _histograms2D.at("jetOriginNconsidered").Fill( 0., recoJet_nConsidered->at(ijet), evtWeight );
                      _histograms2D.at("jetOriginResolution").Fill( 0., SVDeltaR, evtWeight );
                      //_histograms2D.at("jetOriginNormScore").Fill( 0., recoJet_vertex_score->at(ijet)/(float)recoJet_nConsidered->at(ijet), evtWeight );
                    }
                  }
                  std::sort( jetOrigins.begin(), jetOrigins.end() );
                  _histograms2D.at("jetOrigin2D").Fill( jetOrigins.at(0), jetOrigins.at(1), evtWeight );
                } 

            }
        }
    }   
    return;
}

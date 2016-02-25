#include "LLGAnalysis.h"
#include "TLorentzVector.h"

void LLGAnalysis::SetupMakeROOTTrees() {

    // setup the cutflow

    // and the histograms

    return;
}

void LLGAnalysis::MakeROOTTreesSelection() {

    int triggerBitsCounted_PFMET170 = 0;
    int triggerBitsCounted_PFMET90_MHT = 0;
    int triggerBitsCounted_PFMET90_MHT_NoiseCleaned = 0;
    int triggerBitsCounted_PFMET90_MHT_JetCleaned = 0;
    int triggerBitsCounted_Mu50 = 0;
    int triggerBitsCounted_Mu45_eta2p1 = 0;
    int triggerBitsCounted_Ele27_WP85_Gsf = 0;
    for( unsigned int i = 0; i < triggerNames->size(); ++i ) {
      if( triggerNames->at(i).find("HLT_PFMET170_NoiseCleaned") != std::string::npos ) {
        _RT_HLT_PFMET170_NoiseCleaned = triggerBits->at(i);
        triggerBitsCounted_PFMET170 += 1;
      }
      if( triggerNames->at(i).find("HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight") != std::string::npos ) {
        _RT_HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight = triggerBits->at(i);
        triggerBitsCounted_PFMET90_MHT_JetCleaned += 1;
      }
      if( triggerNames->at(i).find("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight") != std::string::npos ) {
        _RT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight = triggerBits->at(i);
        triggerBitsCounted_PFMET90_MHT_NoiseCleaned += 1;
      }
      if( triggerNames->at(i).find("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight") != std::string::npos ) {
        _RT_HLT_PFMETNoMu90_PFMHTNoMu90_IDTight = triggerBits->at(i);
        triggerBitsCounted_PFMET90_MHT += 1;
      }
      if( triggerNames->at(i).find("HLT_Mu50") != std::string::npos ) {
        _RT_HLT_Mu50 = triggerBits->at(i);
        triggerBitsCounted_Mu50 += 1;
      }
      if( triggerNames->at(i).find("HLT_Mu45_eta2p1") != std::string::npos ) {
        _RT_HLT_Mu45_eta2p1 = triggerBits->at(i);
        triggerBitsCounted_Mu45_eta2p1 += 1;
      }
      if( triggerNames->at(i).find("HLT_Ele27_WP85_Gsf") != std::string::npos ) {
        _RT_HLT_Ele27_WP85_Gsf = triggerBits->at(i);
        triggerBitsCounted_Ele27_WP85_Gsf += 1;
      }

    }
    _RT_LeadingMuonPt = 0.;
    _RT_LeadingMuonIso = -1.;
    _RT_LeadingElectronPt = 0.;
    _RT_LeadingElectronIso = -1.;
    
    _RT_RunNumber = RunNumber;
    _RT_EventNumber = EventNumber;
    _RT_LumiBlock = LumiBlock;
    _RT_evtWeight = evtWeight;
    _RT_lumiWeight = lumiWeight;
    _RT_generatorWeight = generatorWeight;
    _RT_pileupWeight = pileupWeight;


    _RT_met->clear();
    _RT_met_y->clear();
    _RT_met_x->clear();

    for( unsigned int i = 0; i < met->size(); ++i ) {
      _RT_met->push_back( met->at(i) );
      _RT_met_x->push_back( met_x->at(i) );
      _RT_met_y->push_back( met_y->at(i) );

    }

    _RT_nJets10 = 0;
    _RT_nJets20 = 0;
    _RT_nJets30 = 0;
    
    _RT_nLooseBJets10 = 0;
    _RT_nLooseBJets20 = 0;
    _RT_nLooseBJets30 = 0;
    _RT_nMediumBJets10 = 0;
    _RT_nMediumBJets20 = 0;
    _RT_nMediumBJets30 = 0;
    _RT_nTightBJets10 = 0;
    _RT_nTightBJets20 = 0;
    _RT_nTightBJets30 = 0;
    _RT_nSVWith2Jets = 0;
    _RT_nPVWithJet150 = 0;
    _RT_PV_LeadingJetPt = 0.;
    _RT_SV_LeadingJetPt->clear();
    _RT_SV_SubLeadingJetPt->clear();
    _RT_SV_LeadingJetEta->clear();
    _RT_SV_SubLeadingJetEta->clear();
    _RT_SV_LeadingJetPhi->clear();
    _RT_SV_SubLeadingJetPhi->clear();
    
    _RT_SV_LeadingDiJetMass = -10.;
    _RT_SV_MaxDistance = -10.;
    _RT_SV_MaxDistanceR = -10.;
    _RT_SV_MaxDistanceZ = -10.;
    _RT_SV_MaxDistance_Uncert = -10.;
    _RT_SV_MaxDistanceR_Uncert = -10.;
    _RT_SV_MaxDistanceZ_Uncert = -10.;

    _RT_nVetoElectrons = vetoElectrons.size();
    _RT_nLooseElectrons = looseElectrons.size();
    _RT_nMediumElectrons = mediumElectrons.size();
    _RT_nTightElectrons = tightElectrons.size();
    _RT_nVetoMuons = vetoMuons.size();
    _RT_nTightMuons = tightMuons.size();
    for( unsigned int im = 0; im < tightMuons.size(); ++im ) {
      int itm = tightMuons.at(im);
      double pt = sqrt( muon_px->at(itm)*muon_px->at(itm) + muon_py->at(itm)*muon_py->at(itm) );
      if( pt > _RT_LeadingMuonPt ) {
        _RT_LeadingMuonPt = pt;
        _RT_LeadingMuonIso = muon_iso->at(itm);
      }
    }
    
    for( unsigned int ie = 0; ie < tightElectrons.size(); ++ie ) {
      int ite = tightElectrons.at(ie);
      double pt = sqrt( electron_px->at(ite)*electron_px->at(ite) + electron_py->at(ite)*electron_py->at(ite) );
      if( pt > _RT_LeadingElectronPt ) {
        _RT_LeadingElectronPt = pt;
        _RT_LeadingElectronIso = electron_iso->at(ite);
      }
    }
    double leadingPV_x = -10000.;
    double leadingPV_y = -10000.;
    double leadingPV_z = -10000.;
    double leadingPV_dx = -10000.;
    double leadingPV_dy = -10000.;
    double leadingPV_dz = -10000.;
    double SVmaxDistR = 0.;

    for( unsigned int iiJet = 0; iiJet < selectedJets.size(); ++iiJet ) {
      int iJet = selectedJets.at(iiJet);
      
      if( recoJet_pt->at(iJet).at(SYSJET) > 10. ) {
        _RT_nJets10 += 1;
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.605 ) _RT_nLooseBJets10 += 1;     
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.890 ) _RT_nMediumBJets10 += 1;     
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.970 ) _RT_nTightBJets10 += 1;     
      }

      if( recoJet_pt->at(iJet).at(SYSJET) > 20. ) {
        _RT_nJets20 += 1;
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.605 ) _RT_nLooseBJets20 += 1;     
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.890 ) _RT_nMediumBJets20 += 1;     
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.970 ) _RT_nTightBJets20 += 1;     
      }

      if( recoJet_pt->at(iJet).at(SYSJET) > 30. ) {
        _RT_nJets30 += 1;
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.605 ) _RT_nLooseBJets30 += 1;     
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.890 ) _RT_nMediumBJets30 += 1;     
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.970 ) _RT_nTightBJets30 += 1;     
      }
    }

    std::vector<std::vector<int> > idxAssocPV(vertex_x->size());
    std::vector<std::vector<int> > idxAssocSV(secVertex_x->size());
    std::vector<std::vector<int> > idxAssocSV30(secVertex_x->size());
    for( unsigned int iiJet = 0; iiJet < selectedJets.size(); ++iiJet ) {
      int iJet = selectedJets.at(iiJet);
      int nMatch = 0;
      for( unsigned int iPV = 0; iPV < vertex_x->size(); ++iPV ) {
        if( fabs( vertex_x->at(iPV) - recoJet_vertex_x->at(iJet) ) < 1.e-10 &&
            fabs( vertex_y->at(iPV) - recoJet_vertex_y->at(iJet) ) < 1.e-10 && 
            fabs( vertex_z->at(iPV) - recoJet_vertex_z->at(iJet) ) < 1.e-10 ) {
            
            idxAssocPV.at(iPV).push_back(iJet);
            nMatch += 1;
        }
      }
      for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
        if( fabs( secVertex_x->at(iSV) - recoJet_vertex_x->at(iJet) ) < 1.e-10 &&
            fabs( secVertex_y->at(iSV) - recoJet_vertex_y->at(iJet) ) < 1.e-10 && 
            fabs( secVertex_z->at(iSV) - recoJet_vertex_z->at(iJet) ) < 1.e-10 ) {
            
            idxAssocSV.at(iSV).push_back(iJet);
            if( recoJet_pt->at(iJet).at(SYSJET) > JET_PT_CUT_SV ) idxAssocSV30.at(iSV).push_back(iJet);
            nMatch += 1;
        }
      }
      
      if( nMatch > 1 ) {
        std::cout << "matched a jet to MORE than 1 vertex!!! " << std::endl;
      }
    }


    for( unsigned int iPV = 0; iPV < vertex_x->size(); ++iPV ) {
      int nJets150 = 0;
      for( unsigned int iiJet = 0; iiJet < idxAssocPV.at(iPV).size(); ++iiJet ) {
        int iJet = idxAssocPV.at(iPV).at(iiJet);
        if( recoJet_pt->at(iJet).at(SYSJET) > _RT_PV_LeadingJetPt ) {
          _RT_PV_LeadingJetPt = recoJet_pt->at(iJet).at(SYSJET);
          leadingPV_x = vertex_x->at(iPV);
          leadingPV_y = vertex_y->at(iPV);
          leadingPV_z = vertex_z->at(iPV);
          leadingPV_dx = vertex_dx->at(iPV);
          leadingPV_dy = vertex_dy->at(iPV);
          leadingPV_dz = vertex_dz->at(iPV);
        }
        if( recoJet_pt->at(iJet).at(SYSJET) > JET_PT_CUT_PV ) nJets150 += 1;
      }
      if( nJets150 >= 1 ) _RT_nPVWithJet150 += 1;
    }
    
    for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
      if( idxAssocSV30.at(iSV).size() >= 2 ) {
        _RT_nSVWith2Jets += 1;
        double leadingJetPt = 0.;
        double subLeadingJetPt = 0.;
        int idxLeadingJet = -1;
        int idxSubLeadingJet = -1;
        for( unsigned int iiJet = 0; iiJet < idxAssocSV30.at(iSV).size(); ++iiJet ) {
          int iJet = idxAssocSV30.at(iSV).at(iiJet);
          if( recoJet_pt->at(iJet).at(SYSJET) > leadingJetPt ) {
            subLeadingJetPt = leadingJetPt;
            idxSubLeadingJet = idxLeadingJet;
            leadingJetPt = recoJet_pt->at(iJet).at(SYSJET);
            idxLeadingJet = iJet;
          }
          else if( recoJet_pt->at(iJet).at(SYSJET) > subLeadingJetPt ) {
            subLeadingJetPt = recoJet_pt->at(iJet).at(SYSJET);
            idxSubLeadingJet = iJet;
          }
        }
        _RT_SV_LeadingJetPt->push_back(recoJet_pt->at(idxLeadingJet).at(SYSJET) );
        _RT_SV_SubLeadingJetPt->push_back(recoJet_pt->at(idxSubLeadingJet).at(SYSJET) );
        _RT_SV_LeadingJetEta->push_back(recoJet_eta->at(idxLeadingJet) );
        _RT_SV_SubLeadingJetEta->push_back(recoJet_eta->at(idxSubLeadingJet) );
        _RT_SV_LeadingJetPhi->push_back(recoJet_phi->at(idxLeadingJet) );
        _RT_SV_SubLeadingJetPhi->push_back(recoJet_phi->at(idxSubLeadingJet ) );
        
        TLorentzVector p4Jet1, p4Jet2;
        p4Jet1.SetPtEtaPhiM( recoJet_pt->at(idxLeadingJet).at(SYSJET), recoJet_eta->at(idxLeadingJet), recoJet_phi->at(idxLeadingJet), 0. );
        p4Jet2.SetPtEtaPhiM( recoJet_pt->at(idxSubLeadingJet).at(SYSJET), recoJet_eta->at(idxSubLeadingJet), recoJet_phi->at(idxSubLeadingJet), 0. );
        TLorentzVector p4DiJet = p4Jet1 + p4Jet2;
        if( p4DiJet.M() > _RT_SV_LeadingDiJetMass ) _RT_SV_LeadingDiJetMass =  p4DiJet.M();
        
        if( p4DiJet.M() < MJJ_CUT ) {
          double dx = leadingPV_x - secVertex_x->at(iSV);
          double dy = leadingPV_y - secVertex_y->at(iSV);
          double dz = fabs( leadingPV_z - secVertex_z->at(iSV) );
          double dr = sqrt( dx*dx + dy*dy );
          double d = sqrt( dx*dx + dy*dy + dz*dz );
          double edx = sqrt( leadingPV_dx*leadingPV_dx + secVertex_dx->at(iSV)*secVertex_dx->at(iSV) );          
          double edy = sqrt( leadingPV_dy*leadingPV_dy + secVertex_dy->at(iSV)*secVertex_dy->at(iSV) );          
          double edz = sqrt( leadingPV_dz*leadingPV_dz + secVertex_dz->at(iSV)*secVertex_dz->at(iSV) );          
          if( dr > SVmaxDistR ) {
            _RT_SV_MaxDistance = d;
            _RT_SV_MaxDistanceR = dr;
            _RT_SV_MaxDistanceZ = dz;
            _RT_SV_MaxDistance_Uncert = sqrt( dx/d*dx/d*edx*edx + dy/d*dy/d*edy*edy + dz/d*dz/d*edz*edz );
            _RT_SV_MaxDistanceR_Uncert = sqrt( dx/dr*dx/dr*edx*edx + dy/dr*dy/dr*edy*edy );
            _RT_SV_MaxDistanceZ_Uncert = edz; 
          }
        }
      }
    }
    _RT_outputTree->Fill();
    return;
}

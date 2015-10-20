#include "LLGAnalysis.h"
#include "TLorentzVector.h"

void LLGAnalysis::SetupMakeROOTTrees() {

    // setup the cutflow

    // and the histograms

    return;
}

void LLGAnalysis::MakeROOTTreesSelection() {

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
    _RT_PV_LeadingJetPt = 0.;
    _RT_SV_LeadingDiJetMass = 10000.;

    _RT_nVetoElectrons = vetoElectrons.size();
    _RT_nLooseElectrons = looseElectrons.size();
    _RT_nMediumElectrons = mediumElectrons.size();
    _RT_nTightElectrons = tightElectrons.size();
    _RT_nVetoMuons = vetoMuons.size();
    _RT_nTightMuons = tightMuons.size();
    
    for( unsigned int iiJet = 0; iiJet < selectedJets.size(); ++iiJet ) {
      int iJet = selectedJets.at(iiJet);
      
      if( recoJet_pt->at(iJet) > 10. ) {
        _RT_nJets10 += 1;
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.605 ) _RT_nLooseBJets10 += 1;     
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.890 ) _RT_nMediumBJets10 += 1;     
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.970 ) _RT_nTightBJets10 += 1;     
      }

      if( recoJet_pt->at(iJet) > 20. ) {
        _RT_nJets20 += 1;
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.605 ) _RT_nLooseBJets20 += 1;     
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.890 ) _RT_nMediumBJets20 += 1;     
        if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.970 ) _RT_nTightBJets20 += 1;     
      }

      if( recoJet_pt->at(iJet) > 30. ) {
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
            if( recoJet_pt->at(iJet) > JET_PT_CUT_SV ) idxAssocSV30.at(iSV).push_back(iJet);
            nMatch += 1;
        }
      }
      
      if( nMatch > 1 ) {
        std::cout << "matched a jet to MORE than 1 vertex!!! " << std::endl;
      }
    }

    
    for( unsigned int iPV = 0; iPV < vertex_x->size(); ++iPV ) {
      for( unsigned int iiJet = 0; iiJet < idxAssocPV.at(iPV).size(); ++iiJet ) {
        int iJet = idxAssocPV.at(iPV).at(iiJet);
        if( recoJet_pt->at(iJet) > _RT_PV_LeadingJetPt ) _RT_PV_LeadingJetPt = recoJet_pt->at(iJet);
      }
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
          if( recoJet_pt->at(iJet) > leadingJetPt ) {
            subLeadingJetPt = leadingJetPt;
            idxSubLeadingJet = idxLeadingJet;
            leadingJetPt = recoJet_pt->at(iJet);
            idxLeadingJet = iJet;
          }
          else if( recoJet_pt->at(iJet) > subLeadingJetPt ) {
            subLeadingJetPt = recoJet_pt->at(iJet);
            idxSubLeadingJet = iJet;
          }
        }
        TLorentzVector p4Jet1, p4Jet2;
        p4Jet1.SetPtEtaPhiM( recoJet_pt->at(idxLeadingJet), recoJet_eta->at(idxLeadingJet), recoJet_phi->at(idxLeadingJet), 0. );
        p4Jet2.SetPtEtaPhiM( recoJet_pt->at(idxSubLeadingJet), recoJet_eta->at(idxSubLeadingJet), recoJet_phi->at(idxSubLeadingJet), 0. );
        TLorentzVector p4DiJet = p4Jet1 + p4Jet2;
        if( p4DiJet.M() > _RT_SV_LeadingDiJetMass ) _RT_SV_LeadingDiJetMass =  p4DiJet.M();
      }
    }


    _RT_outputTree->Fill();
    return;
}

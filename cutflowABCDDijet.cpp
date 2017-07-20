#include "LLGAnalysis.h"
#include "TLorentzVector.h"
void LLGAnalysis::SetupABCDDijet() {

	RegionA=0;
	RegionAWeighted=0;
	RegionB=0;
	RegionBWeighted=0;
	RegionC=0;
	RegionCWeighted=0;
	RegionD=0;
	RegionDWeighted=0;

	makeHist( "mJJSV", 100, 0., 500., "DiJet mass at SV", "Number of Jet Pairs" );

	return;
}

void LLGAnalysis::ABCDDijetSelection() {
	//assign jets to vertices
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
		vector<double> position(3,0.);
		position.at(0) = recoJet_vertex_x->at(iJet);
		position.at(1) = recoJet_vertex_y->at(iJet);
		position.at(2) = recoJet_vertex_z->at(iJet);
		int nMatch = 0;

		for( unsigned int iVtx = 0; iVtx < vertex_x->size(); ++iVtx ) {
			if( fabs(position.at(0) - vertex_x->at(iVtx) ) < 1.e-10 &&
			fabs(position.at(1) - vertex_y->at(iVtx) ) < 1.e-10 &&
			fabs(position.at(2) - vertex_z->at(iVtx) ) < 1.e-10 ) {
				nJetsToPV.at(iVtx) += 1;
				idJetsToPV.at(iVtx).push_back( iJet );
				nMatch += 1;
			}
		}

		for( unsigned int iVtx = 0; iVtx < secVertex_x->size(); ++iVtx ) {
			if( fabs(position.at(0) - secVertex_x->at(iVtx) ) < 1.e-10 &&
			fabs(position.at(1) - secVertex_y->at(iVtx) ) < 1.e-10 &&
			fabs(position.at(2) - secVertex_z->at(iVtx) ) < 1.e-10 ) {
				nJetsToSV.at(iVtx) += 1;
				idJetsToSV.at(iVtx).push_back( iJet );
				nMatch += 1;
			}
		}

		if( nMatch > 1 ) {
			cout << "WARNING! ASSOCIATED JET TO MORE THAN 1 VERTEX ?!" << endl;
		}

	}

	//sort Secondary Vertex jets
	int idxLeadingJet = -1;
	double ptLeadingJet = -1.;
	int idxSubLeadingJet = -1;
	double ptSubLeadingJet = -1;
	int idxThirdLeadingJet = -1;
	double ptThirdLeadingJet = -1;
	int idxFourthLeadingJet = -1;
	double ptFourthLeadingJet = -1;

	for (unsigned int iSV = 0; iSV < secVertex_x->size(); iSV++) {cout<<"Firstfor";
		for( unsigned int iJToSV = 0; iJToSV < idJetsToSV.at(iSV).size(); ++iJToSV ) {cout<<"Secondfor";
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

		//Calculate DiJet Mass
		TLorentzVector p4Jet1, p4Jet2;
		p4Jet1.SetPtEtaPhiM( recoJet_pt->at(idxLeadingJet).at(SYSJET), recoJet_eta->at(idxLeadingJet), recoJet_phi->at(idxLeadingJet), 0. );
		p4Jet2.SetPtEtaPhiM( recoJet_pt->at(idxSubLeadingJet).at(SYSJET), recoJet_eta->at(idxSubLeadingJet), recoJet_phi->at(idxSubLeadingJet), 0. );
		TLorentzVector p4DiJet = p4Jet1 + p4Jet2;
		_histograms1D.at("mJJSV").Fill( p4DiJet.M(), evtWeight );

	}

	return;

}

#include "LLGAnalysis.h"
#include "TLorentzVector.h"
#include <cmath>
void LLGAnalysis::SetupABCDDijet() {

	_cutFlow.insert(pair<string,int>("0_NoCut", 0 ) );
	_cutFlow.insert(pair<string,int>("1_MET", 0 ) );
	RegionA=0;
	RegionAWeighted=0;
	RegionAError=0;
	RegionB=0;
	RegionBWeighted=0;
	RegionBError=0;
	RegionC=0;
	RegionCWeighted=0;
	RegionCError=0;
	RegionD=0;
	RegionDWeighted=0;
	RegionDError=0;

	makeHist( "mJJSV", 100, 0., 500., "DiJet mass at SV", "Number of Jet Pairs" );

	return;
}

void LLGAnalysis::ABCDDijetSelection() {
	_cutFlow.at("0_NoCut") += 1;
	if (sqrt(met_x->at(SYSMET) * met_x->at(SYSMET) + met_y->at(SYSMET) * met_y->at(SYSMET)) > MET_CUT ) {
		return;
	}
	_cutFlow.at("1_MET") += 1;
	int leadingPV = -1;
	double leadingVertexPt = 0.;
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
				if( recoJet_pt->at(iJet).at(SYSJET) > leadingVertexPt ) {
					leadingVertexPt = recoJet_pt->at(iJet).at(SYSJET);
					leadingPV = iVtx;
				}
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

	for (unsigned int iSV = 0; iSV < secVertex_x->size(); iSV++) {
		if ( idJetsToSV.at(iSV).size() <= 1 ){
			if ( leadingVertexPt<150 ) {
				RegionA++;
				RegionAWeighted+=evtWeight;
			} else {
				RegionD++;
				RegionDWeighted+=evtWeight;
			}
			continue; //skip if vertex has fewer than 2 jets
		}
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

		//Calculate DiJet Mass
		TLorentzVector p4Jet1, p4Jet2;
		p4Jet1.SetPtEtaPhiM( recoJet_pt->at(idxLeadingJet).at(SYSJET), recoJet_eta->at(idxLeadingJet), recoJet_phi->at(idxLeadingJet), 0. );
		p4Jet2.SetPtEtaPhiM( recoJet_pt->at(idxSubLeadingJet).at(SYSJET), recoJet_eta->at(idxSubLeadingJet), recoJet_phi->at(idxSubLeadingJet), 0. );
		TLorentzVector p4DiJet = p4Jet1 + p4Jet2;
		_histograms1D.at("mJJSV").Fill( p4DiJet.M(), evtWeight );

		if ( p4DiJet.M() <60 ) {
			if ( leadingVertexPt <150 ) {
				RegionA++;
				RegionAWeighted+=evtWeight;
				RegionAError+=pow(evtWeight,2);
			} else {
				RegionD++;
				RegionDWeighted+=evtWeight;
				RegionDError+=pow(evtWeight,2);
			}
		} else {
			if ( leadingVertexPt<150 ) {
				RegionB++;
				RegionBWeighted+=evtWeight;
				RegionBError+=pow(evtWeight,2);
			} else {
				RegionC++;
				RegionCWeighted+=evtWeight;
				RegionCError+=pow(evtWeight,2);
			}
		}
	}

	return;

}

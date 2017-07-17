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

	//check conditions
	

	return;
}

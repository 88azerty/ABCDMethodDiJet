#include "LLGAnalysis.h"
void LLGAnalysis::SetupABCDMethod() {

	// setup the cutflow

	// and the histograms
	makeHist( "ww", 10, 0, 10, "Region", "Count");
	makeHist( "nw", 10, 0, 10, "Region", "Count");
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

void LLGAnalysis::ABCDMethodSelection() {

	int leadingPV = -1;
	double leadingVertexPt = 0.;

	// assign jets to the vertices:
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

	//Fill the regions
	unsigned int numberVtx1Jet75 = 0;
	unsigned int numberVtx2Jet30 = 0;
	for ( unsigned int p = 0 ; p < idJetsToPV.size() ; p++ ){ //se pasea por los vértices
		unsigned int numberJet75 = 0;
		for ( unsigned int q = 0; q < idJetsToPV.at(p).size() ; q++ ){ //se pasea por los jets
			unsigned int idcurrentjet = idJetsToPV.at(p).at(q);
			if ( recoJet_pt->at(idcurrentjet).at(SYSJET) > 75 ){
				numberJet75++;					//numero de jets con pt>75 en este vertice
			}
		}
		if ( numberJet75 == 1 ){
			numberVtx1Jet75++;				//numero de vertices con exactamente 1 jet con pt>75
		}
	}
	for ( unsigned int p = 0 ; p < idJetsToSV.size() ; p++ ){
		unsigned int numberJet30 = 0;
		for ( unsigned int q = 0; q<idJetsToSV.at(p).size(); q++ ){
			unsigned idcurrentjet = idJetsToSV.at(p).at(q);
			if ( recoJet_pt->at(idcurrentjet).at(SYSJET) > 30 ){
				numberJet30++;
			}
		}
		if ( numberJet30 == 2){
			numberVtx2Jet30++;
		}
	}

	if ( numberVtx1Jet75 == 0 && numberVtx2Jet30 == 0){
		RegionA++;
		RegionAWeighted+=evtWeight;
		_histograms1D.at("ww").Fill("RegionA",evtWeight);
		_histograms1D.at("nw").Fill("RegionA",1);
	}
	if ( numberVtx1Jet75 == 0 && numberVtx2Jet30 == 1){
		RegionC++;
		RegionCWeighted+=evtWeight;
		_histograms1D.at("ww").Fill("RegionC",evtWeight);
		_histograms1D.at("nw").Fill("RegionC",1);
	}
	if ( numberVtx1Jet75 == 1 && numberVtx2Jet30 == 0){
		RegionB++;
		RegionBWeighted+=evtWeight;
		_histograms1D.at("ww").Fill("RegionB",evtWeight);
		_histograms1D.at("nw").Fill("RegionB",1);
	}
	if ( numberVtx1Jet75 == 1 && numberVtx2Jet30 == 1){
		RegionD++;
		RegionDWeighted+=evtWeight;
		_histograms1D.at("ww").Fill("RegionD",evtWeight);
		_histograms1D.at("nw").Fill("RegionD",1);
	}
	return;
}

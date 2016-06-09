#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <stdlib.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TColor.h>
#include <math.h>
#include <THnSparse.h>
#include <TChain.h>
#include <map>
#include <string>
#include <vector>
#include <TRandom3.h>
#include "THnSparse.h"
#include "TF1.h"
#include "TSystem.h"

using namespace std;

class LLGAnalysis {
    public:
        static LLGAnalysis* GetInstance( char* configFileName );
        ~LLGAnalysis() {}
        
        vector<double> CalculateVertex( vector<double> x, vector<double> y, vector<double> z, vector<double> weight, vector<int> charge, vector<double> distance, unsigned int &nConsidered, double &weightednConsidered, vector<double> &error ); 
        
        bool Init();
        void RunEventLoop( int nEventsMax = -1);
        void FinishRun();
    
    private:
        static LLGAnalysis* _instance;

        void RunObjectID();

        void makeHist( string nametitle, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, string xtitle, string ytitle, string ztitle, string drawOption = "", double xAxisOffset = 1., double yAxisOffset = 1.2, double zAxisOffset = 1. ); 
        void makeHist( string nametitle, int nbins, double xmin, double xmax, string xtitle, string ytitle, string drawOption = "", double xAxisOffset = 1., double yAxisOffset = 1.2 );
        void setStyle(double ytoff = 1.0, bool marker = true, double left_margin = 0.15); 
        void MakeEfficiencyPlot( TH1D hpass, TH1D htotal, TCanvas *c, string triggerName = "");
        void FillEfficiencyHistograms();

        std::vector<double> CalculatePCA( std::vector<double> *refPoint, std::vector<double> *momentum, std::vector<double> *vertex );

        void SetupSignalRegion();
        void SignalRegionSelection();
        void SetupWJetsCR();
        void WJetsCRSelection();
        void SetupMakeROOTTrees();
        void MakeROOTTreesSelection();
        void SetupSignalRegionTruthAnalysis();
        void SignalRegionTruthAnalysisSelection();
        void SetupWJetsTestRegion();
        void WJetsTestRegionSelection();
        void SetupZJetsTestRegion();
        void ZJetsTestRegionSelection();
        void SetupMuonTriggerDetermination();
        void MuonTriggerDeterminationSelection();
        void SetupVertexStudy();
        void VertexStudySelection();
        void SetupTriggerCheck();
        void TriggerCheckSelection();
        void SetupMETTriggerEfficiencyDetermination();
        void METTriggerEfficiencyDeterminationSelection();
        void SetupGenHTAnalysis();
        void GenHTAnalysisSelection();
        void SetupRecoJetEfficiencyAnalysis();
        void RecoJetEfficiencyAnalysisSelection();
        // INSERT YOUR SELECTION HERE


    private:
        LLGAnalysis() {}
        LLGAnalysis( char* configFileName );
        
        map<string, int>        _cutFlow;
        map<string, TH1D>       _histograms1D;
        map<string, TH2D>       _histograms2D;
        map<string, string>     _histograms1DDrawOptions;
        map<string, string>     _histograms2DDrawOptions;
        vector<string>          _inputFileNames;
        string                  _inputTreeName;
        TChain                  *_inputTree;
        string                  _outputFileName;
        string                  _RT_outputFileName;
        TTree                   *_outputTree;
        TFile                   *_outputFile;
        TTree                   *_RT_outputTree;
        TFile                   *_RT_outputFile;
        bool                    _writeOutputTree;

        // variables for the 'normal' output trees
        vector<double> *recoCHSJet_pt;
        vector<double> *recoCHSJet_eta;
        vector<double> *recoCHSJet_phi;
        vector<double> *recoNoCHSJet_pt;
        vector<double> *recoNoCHSJet_eta;
        vector<double> *recoNoCHSJet_phi;
        vector<vector<double> > *recoJet_pt; 
        vector<double> *recoJet_phi; 
        vector<double> *recoJet_eta; 
        vector<double> *recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags;
        vector<double> *recoJet_btag_jetBProbabilityBJetTags;
        vector<double> *recoJet_btag_jetProbabilityBJetTags;
        vector<double> *recoJet_btag_trackCountingHighPurBJetTags;
        vector<double> *recoJet_btag_trackCountingHighEffBJetTags;
        vector<double> *recoJet_vertex_x;
        vector<double> *recoJet_vertex_y;
        vector<double> *recoJet_vertex_z;
        vector<double> *recoJet_vertex_score;
        vector<int> *recoJet_nConsidered;
        vector<double> *recoJet_averageDistance;;
        vector<double> *recoJet_rmsDistance;
        vector<double> *muon_px; 
        vector<double> *muon_py; 
        vector<double> *muon_pz; 
        vector<double> *muon_phi; 
        vector<double> *muon_eta; 
        vector<double> *muon_iso;
        vector<bool>   *muon_isLooseMuon;
        vector<bool>   *muon_isTightMuon;
        vector<double> *electron_px; 
        vector<double> *electron_py; 
        vector<double> *electron_pz; 
        vector<double> *electron_phi; 
        vector<double> *electron_eta; 
        vector<double> *electron_iso;
        vector<bool> *electron_isVeto;
        vector<bool> *electron_isLoose;
        vector<bool> *electron_isMedium;
        vector<bool> *electron_isTight;
        vector<bool> *electron_isHEEP;
        vector<int> *triggerBits; 
        vector<string> *triggerNames; 
        /*
        vector<vector<double> > *recoJet_constVertex_x; 
        vector<vector<double> > *recoJet_constVertex_y; 
        vector<vector<double> > *recoJet_constVertex_z; 
        vector<vector<double> > *recoJet_const_pt; 
        vector<vector<double> > *recoJet_const_closestVertex_dxy; 
        vector<vector<double> > *recoJet_const_closestVertex_dz; 
        vector<vector<double> > *recoJet_const_closestVertex_d;
        vector<vector<int> > *recoJet_const_charge; 
        */
        std::vector<std::vector<double> >* recoCHSJet_constVertex_x; 
        std::vector<std::vector<double> >* recoCHSJet_constVertex_y; 
        std::vector<std::vector<double> >* recoCHSJet_constVertex_z; 
        std::vector<std::vector<double> >* recoCHSJet_constVertexRef_x; 
        std::vector<std::vector<double> >* recoCHSJet_constVertexRef_y;
        std::vector<std::vector<double> >* recoCHSJet_constVertexRef_z; 
        std::vector<std::vector<double> >* recoCHSJet_const_pt;     
        std::vector<std::vector<double> >* recoCHSJet_const_eta;     
        std::vector<std::vector<double> >* recoCHSJet_const_phi;     
        std::vector<std::vector<int> >*    recoCHSJet_const_charge;  
        std::vector<std::vector<int> >*    recoCHSJet_const_fromPV;  
        std::vector<std::vector<double> >* recoCHSJet_const_pca0_x; 
        std::vector<std::vector<double> >* recoCHSJet_const_pca0_y; 
        std::vector<std::vector<double> >* recoCHSJet_const_pca0_z; 
        std::vector<std::vector<double> >* recoCHSJet_const_closestVertex_dxy; 
        std::vector<std::vector<double> >* recoCHSJet_const_closestVertex_dz; 
        std::vector<std::vector<double> >* recoCHSJet_const_closestVertex_d; 
        
        std::vector<std::vector<double> >* recoNoCHSJet_constVertex_x; 
        std::vector<std::vector<double> >* recoNoCHSJet_constVertex_y; 
        std::vector<std::vector<double> >* recoNoCHSJet_constVertex_z; 
        std::vector<std::vector<double> >* recoNoCHSJet_constVertexRef_x; 
        std::vector<std::vector<double> >* recoNoCHSJet_constVertexRef_y; 
        std::vector<std::vector<double> >* recoNoCHSJet_constVertexRef_z;
        std::vector<std::vector<double> >* recoNoCHSJet_const_pt;      
        std::vector<std::vector<double> >* recoNoCHSJet_const_eta;     
        std::vector<std::vector<double> >* recoNoCHSJet_const_phi;     
        std::vector<std::vector<int> >*    recoNoCHSJet_const_charge;  
        std::vector<std::vector<int> >*    recoNoCHSJet_const_fromPV;  
        std::vector<std::vector<double> >* recoNoCHSJet_const_pca0_x; 
        std::vector<std::vector<double> >* recoNoCHSJet_const_pca0_y; 
        std::vector<std::vector<double> >* recoNoCHSJet_const_pca0_z; 
        std::vector<std::vector<double> >* recoNoCHSJet_const_closestVertex_dxy; 
        std::vector<std::vector<double> >* recoNoCHSJet_const_closestVertex_dz; 
        std::vector<std::vector<double> >* recoNoCHSJet_const_closestVertex_d; 
   
        std::vector<double>* tightJet_pt;
        std::vector<double>* tightJet_eta;
        std::vector<double>* tightJet_phi;

        std::vector<double>* signalJets_pt;
        std::vector<double>* signalJets_eta; 
        std::vector<double>* signalJets_phi; 
        std::vector<std::vector<double> >* signalJets_constVertex_x; 
        std::vector<std::vector<double> >* signalJets_constVertex_y; 
        std::vector<std::vector<double> >* signalJets_constVertex_z; 
        std::vector<std::vector<double> >* signalJets_constVertexRef_x; 
        std::vector<std::vector<double> >* signalJets_constVertexRef_y; 
        std::vector<std::vector<double> >* signalJets_constVertexRef_z; 
        std::vector<std::vector<double> >* signalJets_const_pt; 
        std::vector<std::vector<int> >* signalJets_const_charge;
        std::vector<std::vector<int> >* signalJets_const_fromPV; 
        std::vector<std::vector<double> >* signalJets_const_pca0_x; 
        std::vector<std::vector<double> >* signalJets_const_pca0_y;
        std::vector<std::vector<double> >* signalJets_const_pca0_z; 
        std::vector<std::vector<double> >* signalJets_const_eta;
        std::vector<std::vector<double> >* signalJets_const_phi; 
        double gluinoProdVertex_x, gluinoProdVertex_y, gluinoProdVertex_z; 
        std::vector<double> *gluinoDecVertex_x;
        std::vector<double> *gluinoDecVertex_y;
        std::vector<double> *gluinoDecVertex_z;
 
        vector<bool> *recoJet_isLeptonLike;

        vector<double> *vertex_x;
        vector<double> *vertex_y;
        vector<double> *vertex_z; 
        vector<double> *vertex_nTracks; 
        vector<double> *vertex_pt; 
        vector<double> *vertex_ndof;
        vector<double> *vertex_dx;
        vector<double> *vertex_dy;
        vector<double> *vertex_dz;

        vector<double> *secVertex_x;
        vector<double> *secVertex_y; 
        vector<double> *secVertex_z;
        vector<double> *secVertex_ndof;
        vector<double> *secVertex_chi2;
        vector<double> *secVertex_pt;
        vector<double> *secVertex_dx;
        vector<double> *secVertex_dy; 
        vector<double> *secVertex_dz; 
         
        vector<string> *to_TriggerNames;
        vector<vector<double> > *to_pt;
        vector<vector<double> > *to_eta;
        vector<vector<double> > *to_phi;


        vector<double> *mct_px;
        vector<double> *mct_py;
        vector<double> *mct_pz;
        vector<double> *mct_e;
        vector<int> *mct_id;
        vector<int> *mct_status;
        vector<vector<int> > *mct_parentId;
        vector<vector<int> > *mct_parentStatus;
        vector<vector<double> > *mct_parentPx;
        vector<vector<double> > *mct_parentPy;
        vector<vector<double> > *mct_parentPz;
        vector<vector<double> > *mct_parentE;

        
        vector<int> vetoElectrons;
        vector<int> looseElectrons;
        vector<int> mediumElectrons;
        vector<int> tightElectrons;
        vector<int> heepElectrons;
        vector<int> vetoMuons;
        vector<int> tightMuons;
        vector<int> selectedJets;

        vector<double> *met;
        vector<double> *met_x;
        vector<double> *met_y;
        int RunNumber;
        int EventNumber;
        int LumiBlock;
        int NumberOfObservedInteractions;
        float NumberOfTrueInteractions;
        double generatorWeight;
        double genLevel_HT;
        bool applyGenLevelHTCut;
        double GENLEVEL_HT_CUT;
        double lumiWeight;
        double pileupWeight;
        string GenFileName;


        TFile *fTruth; 
        TTree *tTruth;
        int GenRunNumber, GenEventNumber, GenLumiBlock;
        vector<double> *tmct_px; 
        vector<double> *tmct_py; 
        vector<double> *tmct_pz; 
        vector<double> *tmct_vx; 
        vector<double> *tmct_vy; 
        vector<double> *tmct_vz; 
        vector<double> *tmct_e; 
        vector<int> *tmct_id; 
        vector<int> *tmct_status; 
        vector<int> *tmct_parent; 
        vector<vector<int> > *tmct_daughters; 
        int lastTruthEntry;


        string PUTYPE;
        // variables for systematics
        int SYSJET;
        int SYSMET;
        string SYSPILEUP;
        string RUNSYS;

        ofstream passedLogFile;

        // variables for ROOT trees
        double            _RT_generatorWeight;
        double            _RT_pileupWeight;
        double            _RT_lumiWeight;
        vector<double> *  _RT_met;
        vector<double> *  _RT_met_x;
        vector<double> *  _RT_met_y;
        double            _RT_evtWeight;
        double            _RT_LeadingMuonPt;
        double            _RT_LeadingMuonIso;
        double            _RT_LeadingElectronPt;
        double            _RT_LeadingElectronIso;
        int               _RT_RunNumber;
        int               _RT_EventNumber;
        int               _RT_LumiBlock;
        int               _RT_HLT_PFMET170_NoiseCleaned;
        int               _RT_HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight;
        int               _RT_HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight;
        int               _RT_HLT_PFMETNoMu90_PFMHTNoMu90_IDTight;
        int               _RT_HLT_Mu50;
        int               _RT_HLT_Mu45_eta2p1;
        int               _RT_HLT_Ele27_WP85_Gsf;
        int               _RT_HLT_Ele105_CaloIdVT_GsfTrkIdT;
        int               _RT_nVetoElectrons;
        int               _RT_nLooseElectrons;
        int               _RT_nMediumElectrons;
        int               _RT_nTightElectrons;
        int               _RT_nVetoMuons;
        int               _RT_nTightMuons;
        int               _RT_nJets10;
        int               _RT_nJets20;
        int               _RT_nJets30;
        int               _RT_nLooseBJets10;
        int               _RT_nLooseBJets20;
        int               _RT_nLooseBJets30;
        int               _RT_nMediumBJets10;
        int               _RT_nMediumBJets20;
        int               _RT_nMediumBJets30;
        int               _RT_nTightBJets10;
        int               _RT_nTightBJets20;
        int               _RT_nTightBJets30;
        int               _RT_nSVWith2Jets;
        int               _RT_nPVWithJet150;
        double            _RT_PV_LeadingJetPt;
        std::vector<double>*          _RT_SV_LeadingJetPt;
        std::vector<double>*          _RT_SV_SubLeadingJetPt;
        std::vector<double>*          _RT_SV_LeadingJetEta;
        std::vector<double>*          _RT_SV_SubLeadingJetEta;
        std::vector<double>*          _RT_SV_LeadingJetPhi;
        std::vector<double>*          _RT_SV_SubLeadingJetPhi;
        double            _RT_SV_LeadingDiJetMass;
        double            _RT_SV_MaxDistance;
        double            _RT_SV_MaxDistanceR;
        double            _RT_SV_MaxDistanceZ;
        double            _RT_SV_MaxDistance_Uncert;
        double            _RT_SV_MaxDistanceR_Uncert;
        double            _RT_SV_MaxDistanceZ_Uncert;
        std::vector<std::string>*     _RT_AllTriggerNames;
        std::vector<int>*             _RT_AllTriggerBits;

        double evtWeight;
        double JET_PT_CUT_SV;
        double JET_PT_CUT_PV;
        double JET_ETA_CUT;
        double MUON_PT_CUT;
        double ELECTRON_PT_CUT;
        double MET_CUT;
        double LEADING_SV_JET_CUT;
        double MJJ_CUT;
        std::string SELECTION;
        std::string metadataFileName;
        std::string datasetName;
        std::vector<std::vector<double> > _yields2DOptimisation;

        // the total number of events per sample (take it from DAS!)
        double PROC_NTOT;
        // the total cross section for the process in pb
        double PROC_XSEC;
        // the target luminosity in fb-1
        double TARGET_LUMI;

        bool applyEventWeights;
        bool applyPileupWeights;
        

        //for pileup reweighting
        std::string PUFILE;
        TH1D *hPU_weights;

        bool requireGenBranches; 
        vector<string> _plotFormats;

};


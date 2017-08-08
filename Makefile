ROOTGLIBS =`root-config --glibs` -lMinuit
LIBS = `root-config --glibs` -L/usr/X11R6/lib -lm -ldl -lstdc++ -Wl,--no-as-needed $(ROOTGLIBS)

CXXFLAGS = -O3 --exceptions `root-config --ldflags --cflags` -I./ -w

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

RunAnalysis: LLGAnalysis.o RunAnalysis.o ObjectID.o cutflowSignalRegion.o cutflowABCDDijet.o cutflowABCDMethod.o cutflowRecoJetEfficiencyAnalysis.o cutflowGenHTAnalysis.o cutflowMETTriggerEfficiencyDetermination.o cutflowTriggerCheck.o cutflowVertexStudy.o cutflowMuonTriggerDetermination.o cutflowZJetsTestRegion.o cutflowWJetsTestRegion.o cutflowSignalRegionTruthAnalysis.o cutflowMakeROOTTrees.o cutflowWJetsCR.o
	$(CXX) -o RunAnalysis RunAnalysis.o LLGAnalysis.o ObjectID.o cutflowSignalRegion.o cutflowABCDDijet.o cutflowABCDMethod.o cutflowRecoJetEfficiencyAnalysis.o cutflowGenHTAnalysis.o cutflowMETTriggerEfficiencyDetermination.o cutflowTriggerCheck.o cutflowVertexStudy.o cutflowMuonTriggerDetermination.o cutflowZJetsTestRegion.o cutflowWJetsTestRegion.o cutflowSignalRegionTruthAnalysis.o cutflowMakeROOTTrees.o cutflowWJetsCR.o $(LDFLAGS) $(LIBS)
clean:
	rm -f RunAnalysis
	rm -f *.o
	rm -f *~
	rm -f *.so
	rm -f Loader_C.d
	rm -f *.root
	rm -f *.pdf
	rm -f *.log
	rm -f *.eps
	rm -f *.png

void runAll(){
  gROOT->ProcessLine(".L AliAnalysisTaskHyperTriton3VtxPerf.cxx+");
  gROOT->ProcessLine(".L AliAnalysisTaskHyperTriton3VtxPerf.h+");
  gROOT->ProcessLine(".L runTaskGrid.cpp+");
  gROOT->ProcessLine("runTaskGrid()");
}
#include <TFile.h>
#include <TTree.h>
#include "AliSelectorTaskHyperTriton3VtxPerf.h"

void runSelectorTaskHyperTriton3VtxPerf(TString inputFile = "../HyperTritonTree.root", TString outputName = "selector_results.root", int vertexer = 0, TString outputPath = "./") {
    TFile lFile(inputFile.Data());
    TTree* lTree = (TTree*)lFile.Get("Hyp3");
    AliSelectorTaskHyperTriton3VtxPerf lSelector(outputName.Data(), outputPath.Data(), vertexer);
    lTree->Process(&lSelector);
}
void runSelectorTaskHyperTriton3VtxPerfAllVertexer() {
    runSelectorTaskHyperTriton3VtxPerf("Output.root","selector_resultsStd.root",2,"./");
    runSelectorTaskHyperTriton3VtxPerf("Output.root","selector_resultsKFChi2_50.root",0,"./");
    runSelectorTaskHyperTriton3VtxPerf("Output.root","selector_resultsO2.root",1,"./");
}
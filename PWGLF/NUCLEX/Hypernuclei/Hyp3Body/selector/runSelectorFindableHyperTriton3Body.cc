#include <TFile.h>
#include <TTree.h>
#include <AliSelectorFindableHyperTriton3Body.h>

void runSelectorFindableHyperTriton3Body(TString inputFile = "Output.root", TString outputName = "selector_results.root", int vertexer = 0, TString outputPath = "./") {
    TFile lFile(inputFile.Data());
    TTree* lTree = (TTree*)lFile.Get("FindableTree/fTreeHyperTriton3Body");
    AliSelectorFindableHyperTriton3Body lSelector(outputName.Data(), outputPath.Data(), vertexer);
    lTree->Process(&lSelector);
}
void runSelectorFindableHyperTriton3BodyAllVertexer() {
    runSelectorFindableHyperTriton3Body("Output.root","selector_resultsStd.root",2,"./");
    runSelectorFindableHyperTriton3Body("Output.root","selector_resultsKFChi2_50.root",0,"./");
    runSelectorFindableHyperTriton3Body("Output.root","selector_resultsO2.root",1,"./");
}
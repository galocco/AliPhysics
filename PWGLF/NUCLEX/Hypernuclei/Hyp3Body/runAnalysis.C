#include "AliAnalysisTaskHyperTriton3VtxPerf.h"

void runAnalysis(){
  // header location
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  // create the analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("AnalysisMyTask");
  AliESDInputHandler* aodH = new AliESDInputHandler();
  mgr->SetInputEventHandler(aodH);
  // compile the class ( locally )
  gROOT->LoadMacro("AliAnalysisTaskHyperTriton3VtxPerf.h ++g") ;
  gROOT->LoadMacro("AliAnalysisTaskHyperTriton3VtxPerf.cxx ++g") ;
  // load the addtask macro
  gROOT->LoadMacro("AddMyTaskHyperTritonVtxPerf.C");
  // create an instance of your analysis task
  AliAnalysisTaskHyperTriton3VtxPerf * task;
  task->AddTask(true,"");
  // if you want to run locally , we need to define some input
  TChain * chain = new TChain("aodTree");
  chain->Add("001/AliESDs.root");
  // start the analysis locally
  mgr->StartAnalysis("local", chain );
}
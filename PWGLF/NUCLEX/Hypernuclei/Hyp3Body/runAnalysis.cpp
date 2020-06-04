//#include "AliAnalysisTaskHyperTriton3VtxPerf.h"
//#include "AddTaskHyperTritonVtxPerf.h"

//#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"

void runAnalysis(){
  // header location
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  // create the analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("AnalysisMyTask");
  AliESDInputHandler* esdH = new AliESDInputHandler();
  
  mgr->SetInputEventHandler(esdH);
  //???
  esdH->SetNeedField();
  // compile the class ( locally );
  gROOT->ProcessLine(".L AliAnalysisTaskHyperTriton3VtxPerf.cxx+");
  gROOT->ProcessLine(".L AliAnalysisTaskHyperTriton3VtxPerf.h+");
  gROOT->ProcessLine(".L AddTaskHyperTritonVtxPerf.C+");
  gROOT->LoadMacro("AliAnalysisTaskHyperTriton3VtxPerf.cxx ++g") ;
  gROOT->LoadMacro("AliAnalysisTaskHyperTriton3VtxPerf.h") ;
  // load the addtask macro and create the task
  AliAnalysisTaskHyperTriton3VtxPerf *task = reinterpret_cast<AliAnalysisTaskHyperTriton3VtxPerf*>(gInterpreter->ExecuteMacro("AddTaskHyperTritonVtxPerf.C"));
  // if you want to run locally , we need to define some input
  TChain * chain = new TChain("esdTree");
  chain->Add("001/AliESDs.root");
  // start the analysis locally
  mgr->StartAnalysis("local", chain );
  
}
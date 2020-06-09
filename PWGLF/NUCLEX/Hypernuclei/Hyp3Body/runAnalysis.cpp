#include "TSystem.h"
#include "TROOT.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisAlien.h"
#include "AliPhysicsSelectionTask.h"
#include "AliMultSelectionTask.h"
#include "AliAnalysisGrid.h"
#include "AliMultSelectionTask.h"
#include "AliVEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskHyperTriton3VtxPerf.h"
#include "TChain.h"

void runAnalysis(){
  // header location
  gROOT->ProcessLine(".include $ROOTSYS/include");
  // add aliroot include path
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/STEER");
  gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS");
  gROOT->ProcessLine(".include $ALICE_ROOT/ESD");
  gROOT->ProcessLine(".include $ALICE_ROOT/STEER/ESD");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  //gROOT->ProcessLine(".include $ALICE_ROOT/include");
  // create the analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("AnalysisMyTask");
  // Input
  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);
  // compile the class ( locally );
  gROOT->ProcessLine(".L AliAnalysisTaskHyperTriton3VtxPerf.cxx++g");
  //gROOT->ProcessLine(".L AliAnalysisTaskHyperTriton3VtxPerf.h+");
  gROOT->ProcessLine(".L AddTaskHyperTritonVtxPerf.C");
  //gROOT->LoadMacro("AliAnalysisTaskHyperTriton3VtxPerf.cxx ++g") ;
  //gROOT->LoadMacro("AliAnalysisTaskHyperTriton3VtxPerf.h") ;
  // load the addtask macro and create the task
  AliAnalysisTaskHyperTriton3VtxPerf *task = reinterpret_cast<AliAnalysisTaskHyperTriton3VtxPerf*>(gInterpreter->ExecuteMacro("AddTaskHyperTritonVtxPerf.C"));
  // if you want to run locally , we need to define some input
  TChain * chain = new TChain("esdTree");
  chain->Add("AliESDs.root");
  // start the analysis locally
  mgr->StartAnalysis("local", chain );
  
}
//#include "Common.h"

#include <iostream> 
using namespace std;

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

//______________________________________________________________________________
void runTaskGrid() {
  bool kIsMC=true;
  if (kIsMC) {
    printf("MC analysis chosen\n");
  }
  else printf("Data analysis chosen\n");

  // load libraries
  gSystem->Load("libCore");        
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGPP");
  gSystem->Load("libANALYSISaliceBase");
  gSystem->Load("libCORRFW");
  gSystem->Load("libOADB");
  
  // add aliroot include path
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/STEER");
  gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS");
  gROOT->ProcessLine(".include $ALICE_ROOT/ESD");
  gROOT->ProcessLine(".include $ALICE_ROOT/STEER/ESD");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  gROOT->SetStyle("Plain");
  
  //gROOT->ProcessLine(".L AliAnalysisTaskHyperTriton3VtxPerf.h+");
  // analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("AnalysisMyTask");
  
  // create the alien handler and attach it to the manager
  // AliAnalysisGrid *plugin = CreateAlienHandler(); 
  // mgr->SetGridHandler(plugin);

  // Input
  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);
  //AliAODInputHandler* iH = new AliAODInputHandler("handler","handler for my analisys");
  //mgr->SetInputEventHandler(iH);
  
  // mc event handler
  
  if(kIsMC) {
    AliMCEventHandler* mchandler = new AliMCEventHandler();
    // Not reading track references
    mchandler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mchandler);
  }  
  
  //==========================================
  // MULTIPLICITY SELECTION Task 
  //==========================================
  AliMultSelectionTask *mult_task = reinterpret_cast<AliMultSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"))));

  //Select the desired tasks
  //=========================================
  // PHYSICS SELECTION Task
  //=========================================
  AliPhysicsSelectionTask *physSelTask = reinterpret_cast<AliPhysicsSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s(%d)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"), kIsMC)));
  
  //==========================================
  // PID Task
  //==========================================
  AliAnalysisTask *PIDTask = reinterpret_cast<AliAnalysisTask *>(gInterpreter->ProcessLine(Form(".x %s(%d)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"), kIsMC)));

  //gROOT->ProcessLine(".L AliAnalysisTaskHyperTriton3VtxPerf.cxx++g");
  AliAnalysisTaskHyperTriton3VtxPerf *task = reinterpret_cast<AliAnalysisTaskHyperTriton3VtxPerf*>(gInterpreter->ExecuteMacro("AddTaskHyperTritonVtxPerf.C"));

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  mgr->SetDebugLevel(10);
  TChain* chain = new TChain("esdTree");
  chain->Add("AliESDs.root");
  mgr->StartAnalysis("LOCAL", chain); 
}

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////
//
// Author:
// P. Fecchio, pfecchio@cern.ch
///////////////////////////////////////////////////////////////////////////

#include <array>
#include <climits>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

/// ROOT includes
#include <Riostream.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TList.h>
#include <TObjArray.h>
#include <TVector3.h>

/// AliRoot icludes
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultVariable.h"
#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliPDG.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "TLorentzVector.h"

#include "AliAnalysisTaskFindableHypertriton3.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskFindableHypertriton3);

namespace {

const AliPID::EParticleType kSpecies[3] = {AliPID::kDeuteron, AliPID::kProton, AliPID::kPion};

bool IsHyperTriton3(const AliVParticle *vPart, AliMCEvent *mcEvent) {
  int nDaughters = 0;

  int vPartPDG   = vPart->PdgCode();
  int vPartLabel = vPart->GetLabel();

  if (!mcEvent->IsPhysicalPrimary(vPartLabel) || (std::abs(vPartPDG) != 1010010030)) return false;

  for (int iD = vPart->GetDaughterFirst(); iD <= vPart->GetDaughterLast(); iD++) {
    AliVParticle *dPart = mcEvent->GetTrack(iD);

    int dPartPDG = dPart->PdgCode();
    if (std::abs(dPartPDG) != 11) nDaughters++;
  }
  if (nDaughters == 3) return true;
  return false;
}
} // namespace

//________________________________________________________________________
AliAnalysisTaskFindableHypertriton3::AliAnalysisTaskFindableHypertriton3(TString taskname)
    : AliAnalysisTaskSE(taskname.Data()),
      // support objects
      fEventCuts{},
      fPIDResponse{nullptr},
      fESDtrackCuts{nullptr},
      fPrimaryVertex{nullptr},
      // setting parameters
      fCosPoiningAngleLimit{0},
      // output objects
      fOutputList{nullptr},
      fFindableTree{nullptr},
      fTreeHyp3BodyVarTracks{nullptr},
      fTreeHyp3BodyVarNsigmaTPC{0},
      fTreeHyp3BodyVarNsigmaTOF{0},
      fTreeHyp3BodyVarTruePx{0},
      fTreeHyp3BodyVarTruePy{0},
      fTreeHyp3BodyVarTruePz{0},
      fTreeHyp3BodyVarDecayVx{0},
      fTreeHyp3BodyVarDecayVy{0},
      fTreeHyp3BodyVarDecayVz{0},
      fTreeHyp3BodyVarDecayT{0},
      fTreeHyp3BodyVarPVx{0},
      fTreeHyp3BodyVarPVy{0},
      fTreeHyp3BodyVarPVz{0},
      fTreeHyp3BodyVarPVt{0},
      fTreeHyp3BodyVarMagneticField{0},
      fTreeHyp3BodyVarCentrality{0},
      fHistEventCounter{nullptr},
      fHistCentrality{nullptr},
      fHistGeneratedPtVsYVsCentralityHypTrit{nullptr},
      fHistGeneratedPtVsYVsCentralityAntiHypTrit{nullptr} {

  // Standard Output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskFindableHypertriton3::~AliAnalysisTaskFindableHypertriton3() {
  // destructor
  if (fOutputList) {
    delete fOutputList;
    fOutputList = nullptr;
  }
  if (fFindableTree) {
    delete fFindableTree;
    fFindableTree = nullptr;
  }
}

//________________________________________________________________________
void AliAnalysisTaskFindableHypertriton3::UserCreateOutputObjects() {

  AliAnalysisManager *fMgr = AliAnalysisManager::GetAnalysisManager();
  if (!fMgr) AliFatal("Could not find analysis manager.");
  AliInputEventHandler *fHandl = (AliInputEventHandler *)fMgr->GetInputEventHandler();
  if (!fHandl) AliFatal("No input event handler.");
  fPIDResponse = fHandl->GetPIDResponse();
  fHandl->SetNeedField();

  // Multiplicity
  if (!fESDtrackCuts) {
    fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, kFALSE);
    fESDtrackCuts->SetPtRange(0.15); // adding pt cut
    fESDtrackCuts->SetEtaRange(-1.0, 1.0);
  }

  // Create a TList with Histograms
  fOutputList = new TList();
  fOutputList->SetOwner();
  // fEventCuts.AddQAplotsToList(fOutputList);

  // Histogram Output: Event-by-Event
  fHistCentrality =
      new TH1D("fHistCentrality", "WARNING: no pileup rejection applied!;Centrality;Event Count", 100, 0, 100);
  fOutputList->Add(fHistCentrality);

  // Histogram Output: Efficiency Denominator
  fHistGeneratedPtVsYVsCentralityHypTrit = new TH3D("fHistGeneratedPtVsYVsCentralityHypTrit", ";#it{p}_{T} (GeV/#it{c});y;centrality", 500, 0, 25, 40, -1.0, 1.0, 100, 0, 100);
  fOutputList->Add(fHistGeneratedPtVsYVsCentralityHypTrit);
  fHistGeneratedPtVsYVsCentralityAntiHypTrit = new TH3D("fHistGeneratedPtVsYVsCentralityAntiHypTrit", ";#it{p}_{T} (GeV/#it{c});y;centrality", 500, 0, 25, 40, -1.0, 1.0, 100, 0, 100);
  fOutputList->Add(fHistGeneratedPtVsYVsCentralityAntiHypTrit);

  // Histogram Output: Event-by-Event
  fHistEventCounter = new TH1D("fHistEventCounter", ";Evt. Sel. Step;Count", 2, 0, 2);
  fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
  fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected");
  fOutputList->Add(fHistEventCounter);

  OpenFile(2);
  fFindableTree = new TTree("fTreeHyperTriton3Body", "HyperTriton3BodyCandidates");

  fFindableTree->Branch("fTreeHyp3BodyVarTrack0", &fTreeHyp3BodyVarTracks[0]);
  fFindableTree->Branch("fTreeHyp3BodyVarTrack1", &fTreeHyp3BodyVarTracks[1]);
  fFindableTree->Branch("fTreeHyp3BodyVarTrack2", &fTreeHyp3BodyVarTracks[2]);

  fFindableTree->Branch("fPrimaryVertex", &fPrimaryVertex);

  fFindableTree->Branch("fTreeHyp3BodyVarNsigmaTPC0", &fTreeHyp3BodyVarNsigmaTPC[0], "fTreeHyp3BodyVarNsigmaTPC0/F");
  fFindableTree->Branch("fTreeHyp3BodyVarNsigmaTPC1", &fTreeHyp3BodyVarNsigmaTPC[1], "fTreeHyp3BodyVarNsigmaTPC1/F");
  fFindableTree->Branch("fTreeHyp3BodyVarNsigmaTPC2", &fTreeHyp3BodyVarNsigmaTPC[2], "fTreeHyp3BodyVarNsigmaTPC2/F");

  fFindableTree->Branch("fTreeHyp3BodyVarNsigmaTOF0", &fTreeHyp3BodyVarNsigmaTOF[0], "fTreeHyp3BodyVarNsigmaTOF0/F");
  fFindableTree->Branch("fTreeHyp3BodyVarNsigmaTOF1", &fTreeHyp3BodyVarNsigmaTOF[1], "fTreeHyp3BodyVarNsigmaTOF1/F");
  fFindableTree->Branch("fTreeHyp3BodyVarNsigmaTOF2", &fTreeHyp3BodyVarNsigmaTOF[2], "fTreeHyp3BodyVarNsigmaTOF2/F");

  fFindableTree->Branch("fTreeHyp3BodyVarTruePx", &fTreeHyp3BodyVarTruePx, "fTreeHyp3BodyVarTruePx/F");
  fFindableTree->Branch("fTreeHyp3BodyVarTruePy", &fTreeHyp3BodyVarTruePy, "fTreeHyp3BodyVarTruePy/F");
  fFindableTree->Branch("fTreeHyp3BodyVarTruePz", &fTreeHyp3BodyVarTruePz, "fTreeHyp3BodyVarTruePz/F");

  fFindableTree->Branch("fTreeHyp3BodyVarDecayVx", &fTreeHyp3BodyVarDecayVx, "fTreeHyp3BodyVarDecayVx/F");
  fFindableTree->Branch("fTreeHyp3BodyVarDecayVy", &fTreeHyp3BodyVarDecayVy, "fTreeHyp3BodyVarDecayVy/F");
  fFindableTree->Branch("fTreeHyp3BodyVarDecayVz", &fTreeHyp3BodyVarDecayVz, "fTreeHyp3BodyVarDecayVz/F");
  fFindableTree->Branch("fTreeHyp3BodyVarDecayT", &fTreeHyp3BodyVarDecayT, "fTreeHyp3BodyVarDecayT/F");

  fFindableTree->Branch("fTreeHyp3BodyVarPVx", &fTreeHyp3BodyVarPVx, "fTreeHyp3BodyVarPVx/F");
  fFindableTree->Branch("fTreeHyp3BodyVarPVy", &fTreeHyp3BodyVarPVy, "fTreeHyp3BodyVarPVy/F");
  fFindableTree->Branch("fTreeHyp3BodyVarPVz", &fTreeHyp3BodyVarPVz, "fTreeHyp3BodyVarPVz/F");
  fFindableTree->Branch("fTreeHyp3BodyVarPVt", &fTreeHyp3BodyVarPVt, "fTreeHyp3BodyVarPVt/F");

  fFindableTree->Branch("fTreeHyp3BodyVarMagneticField", &fTreeHyp3BodyVarMagneticField,"fTreeHyp3BodyVarMagneticField/F");
  fFindableTree->Branch("fTreeHyp3BodyVarCentrality", &fTreeHyp3BodyVarCentrality, "fTreeHyp3BodyVarCentrality/F");



  PostData(1, fOutputList);
  PostData(2, fFindableTree);
}

//________________________________________________________________________
void AliAnalysisTaskFindableHypertriton3::UserExec(Option_t *) {
  // main loop called for each analized event

  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent) {
    ::Fatal("AliAnalysisTaskStrangenessLifetimes::UserExec", "AliESDEvent not found.");
    return;
  }

  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent) {
    ::Fatal("AliAnalysisTaskStrangenessLifetimes::UserExec", "Could not retrieve MC event");
    return;
  }

  auto ComputeRapidity = [](double rE, double rPz) {
    double rValue = -100;
    if ((rE - rPz + 1.e-13) != 0 && (rE + rPz) != 0) {
      rValue = 0.5 * TMath::Log((rE + rPz) / (rE - rPz + 1.e-13));
    }
    return rValue;
  };

  const long vNTracks           = esdEvent->GetNumberOfTracks();
  fTreeHyp3BodyVarMagneticField = esdEvent->GetMagneticField();

  // total number of analyzed events
  fHistEventCounter->Fill(0.5);

  if (!fEventCuts.AcceptEvent(esdEvent)) {
    PostData(1, fOutputList);
    PostData(2, fFindableTree);
    return;
  }
  fTreeHyp3BodyVarCentrality = fEventCuts.GetCentrality();
  fPrimaryVertex = (AliESDVertex*)fEventCuts.GetPrimaryVertex();
  fHistCentrality->Fill(fTreeHyp3BodyVarCentrality);

  // number of selected events
  fHistEventCounter->Fill(1.5); // selected events

  //--------------------------------------------------------------------------------
  // Part 1: fill the histograms of the MC hypertritons
  //--------------------------------------------------------------------------------
  for (Int_t iPart = 0; iPart < mcEvent->GetNumberOfTracks(); iPart++) {
    AliVParticle *vPart = mcEvent->GetTrack(iPart);
    if (!vPart) {
      ::Warning("AliAnalysisTaskHyperTriton2He3piML::UserExec",
                "Generated loop %i - MC TParticle pointer to current stack particle = 0x0 ! Skipping.", iPart);
      continue;
    }
    if (!mcEvent->IsPhysicalPrimary(iPart)) continue;

    // fill the histos of generated particles for efficiency denominator
    int vPartPDG    = vPart->PdgCode();
    double vPartPt  = vPart->Pt();
    double vPartRap = ComputeRapidity(vPart->E(), vPart->Pz());
    if (vPartPDG == 1010010030)
      fHistGeneratedPtVsYVsCentralityHypTrit->Fill(vPartPt, vPartRap, fTreeHyp3BodyVarCentrality);
    if (vPartPDG == -1010010030)
      fHistGeneratedPtVsYVsCentralityAntiHypTrit->Fill(vPartPt, vPartRap, fTreeHyp3BodyVarCentrality);
  }

  //--------------------------------------------------------------------------------
  // Part 2: establish list of hypertritons in the 3 body channel
  //--------------------------------------------------------------------------------

  std::vector<AliVParticle*> vGenHyp3;
  for (Long_t iTrack = 0; iTrack < vNTracks; iTrack++) {
    AliESDtrack *esdTrack = esdEvent->GetTrack(iTrack);
    if (!esdTrack) continue;
    int lLabel          = (int)TMath::Abs(esdTrack->GetLabel());
    AliVParticle *vPart = mcEvent->GetTrack(lLabel);
    if(IsHyperTriton3(vPart, mcEvent)){
      vGenHyp3.push_back(vPart);
    }
  }

  //--------------------------------------------------------------------------------
  // Part 3: find the triplets of reconstructed daughters and fill the tree
  //--------------------------------------------------------------------------------
  

  if (!vGenHyp3.empty()) {
    fTreeHyp3BodyVarPVt = vGenHyp3.back()->Tv();
    fTreeHyp3BodyVarPVx = vGenHyp3.back()->Xv();
    fTreeHyp3BodyVarPVy = vGenHyp3.back()->Yv();
    fTreeHyp3BodyVarPVz = vGenHyp3.back()->Zv();

    for (size_t iGenTrack = 0; iGenTrack < vGenHyp3.size(); iGenTrack++) {
      //loop
      for (int iD = vGenHyp3[iGenTrack]->GetDaughterFirst(); iD <= vGenHyp3[iGenTrack]->GetDaughterLast(); iD++) {
        AliESDtrack *dPart = esdEvent->GetTrack(iD);
        int dPartPDG        = dPart->PdgCode();
        int sTrack = 0;
        if (std::abs(dPartPDG) == 11){ continue;}
        else if(std::abs(dPartPDG) == 211){ sTrack = 2;} //pion
        else if(std::abs(dPartPDG) == 2212){ sTrack = 1;} //proton
        else{ //deuteron
          sTrack = 0;
          int lLabel          = (int)TMath::Abs(dPart->GetLabel());
          AliVParticle *vPart = mcEvent->GetTrack(lLabel);
          fTreeHyp3BodyVarDecayVx = vPart->Xv();
          fTreeHyp3BodyVarDecayVy = vPart->Yv();
          fTreeHyp3BodyVarDecayVz = vPart->Zv();
          fTreeHyp3BodyVarDecayT  = vPart->Tv();
        }

        fTreeHyp3BodyVarTracks[sTrack] = dPart;
        fTreeHyp3BodyVarNsigmaTPC[sTrack] = fPIDResponse->NumberOfSigmasTPC(dPart,kSpecies[sTrack]);
        fTreeHyp3BodyVarNsigmaTOF[sTrack] = (HasTOF(dPart)) ? fPIDResponse->NumberOfSigmasTOF(dPart,kSpecies[sTrack]): -999.;
      }

      fTreeHyp3BodyVarTruePx     = vGenHyp3[iGenTrack]->Px();
      fTreeHyp3BodyVarTruePy     = vGenHyp3[iGenTrack]->Py();
      fTreeHyp3BodyVarTruePz     = vGenHyp3[iGenTrack]->Pz();
      fFindableTree->Fill();
    }
  }

  PostData(1, fOutputList);
  PostData(2, fFindableTree);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskFindableHypertriton3::Terminate(Option_t *) {
  // Merge output
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList *>(GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: fOutputList not available\n");
    return;
  }

  printf("end of Terminate");
  return;

} // end of Terminate


bool AliAnalysisTaskFindableHypertriton3::HasTOF(AliESDtrack *track) {
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len = track->GetIntegratedLength();
  return hasTOFout && hasTOFtime && (len > 350.);
}

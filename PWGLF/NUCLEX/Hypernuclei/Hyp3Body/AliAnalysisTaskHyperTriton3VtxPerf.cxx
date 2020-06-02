#include "AliAnalysisTaskHyperTriton3VtxPerf.h"

#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliPDG.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVVertex.h"

#include <Riostream.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <vector>
#include <unordered_map>

#include "Math/Vector3Dfwd.h"
#include "Math/Vector3D.h"

#include "AliDataFile.h"
#include <TFile.h>
#include <TSpline.h>

#define HomogeneousField
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPVertex.h"
#include "KFPTrack.h"

ClassImp(AliAnalysisTaskHyperTriton3VtxPerf);

namespace {

struct HelperParticle {
  AliESDtrack* track = nullptr;
  float        nSigmaTPC = -1.f;
  float        nSigmaTOF = -1.f;
};

constexpr float kHyperTritonMass{2.99131};

constexpr float kDeuMass{1.87561};
constexpr float kPMass{0.938272};
constexpr float kPiMass{0.13957};
constexpr float kMasses[3]{kDeuMass,kPMass,kPiMass};
constexpr AliPID::EParticleType kAliPID[3]{AliPID::kDeuteron,AliPID::kProton,AliPID::kPion};
const int kPDGs[3]{AliPID::ParticleCode(kAliPID[0]),AliPID::ParticleCode(kAliPID[1]),AliPID::ParticleCode(kAliPID[2])};

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

int IsTrueHyperTriton3Candidate(AliESDtrack *t1, AliESDtrack *t2, AliESDtrack *t3, AliMCEvent *mcEvent) {
  if (!mcEvent) return 0;

  int lab1 = std::abs(t1->GetLabel());
  int lab2 = std::abs(t2->GetLabel());
  int lab3 = std::abs(t3->GetLabel());

  if (mcEvent->IsPhysicalPrimary(lab1)) return -1;
  if (mcEvent->IsPhysicalPrimary(lab2)) return -1;
  if (mcEvent->IsPhysicalPrimary(lab3)) return -1;

  AliVParticle *part1 = mcEvent->GetTrack(lab1);
  AliVParticle *part2 = mcEvent->GetTrack(lab2);
  AliVParticle *part3 = mcEvent->GetTrack(lab3);

  if (!part1 || !part2 || !part3) 
    return -1;

  int mom1 = part1->GetMother();
  int mom2 = part2->GetMother();
  int mom3 = part3->GetMother();

  if (mom1 != mom2 || mom1 != mom3 || mom2 != mom3) return -1;

  AliVParticle *mom = mcEvent->GetTrack(mom1);
  if (!mom) return -1;

  return (IsHyperTriton3(mom, mcEvent)) ? mom1 : -1;
}

bool HasTOF(AliVTrack *track) {
  const bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  const bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  return hasTOFout && hasTOFtime;
}

int GetIndex(int PdgCode){
  else if(std::abs(dPartPDG) == 211){ return 2; } //pion
  else if(std::abs(dPartPDG) == 2212){ return 1; } //proton
  return 0; //deuteron
}

/// helper functions
template <typename T> double Sq(T a) { return a * a; }
template <typename F> double Hypot(F a, F b, F c) { return std::sqrt(Sq(a) + Sq(b) + Sq(c)); }
template <typename F> double Hypot(F a, F b, F c, F d) { return std::sqrt(Sq(a) + Sq(b) + Sq(c) + Sq(d)); }

}    // namespace

AliAnalysisTaskHyperTriton3VtxPerf::AliAnalysisTaskHyperTriton3VtxPerf(bool mc, std::string name)
    : AliAnalysisTaskSE(name.data()), fEventCuts{}, fMC{mc}, fREvent{}, fGenHyp{} {
  fTrackCuts.SetMinNClustersTPC(0);
  fTrackCuts.SetEtaRange(-0.9,0.9);
  /// Settings for the custom vertexer

  /// Standard output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());    // Basic Histograms
  DefineOutput(2, TTree::Class());    // Hypertriton Candidates Tree output
}

AliAnalysisTaskHyperTriton3VtxPerf::~AliAnalysisTaskHyperTriton3VtxPerf() {
  if (fListHist) {
    delete fListHist;
    fListHist = nullptr;
  }

  if (fTreeHyp3) {
    delete fTreeHyp3;
    fTreeHyp3 = nullptr;
  }

  if (fCosPAsplineFile)
    delete fCosPAsplineFile;

}

void AliAnalysisTaskHyperTriton3VtxPerf::UserCreateOutputObjects() {
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  fInputHandler           = (AliInputEventHandler *)(man->GetInputEventHandler());
  fPIDResponse            = fInputHandler->GetPIDResponse();

  fInputHandler->SetNeedField();

  fListHist = new TList();
  fListHist->SetOwner(true);
  fEventCuts.AddQAplotsToList(fListHist);

  OpenFile(2);
  fTreeHyp3 = new TTree("Hyp3KF","Hypetriton 3 Body with the KFParticle");
  fTreeHyp3->Branch("RCollision", &fREvent);

  if (man->GetMCtruthEventHandler()) {
    fTreeHyp3->Branch("SHyperTriton", &fGenHyp);
    fTreeHyp3->Branch("SGenRecMap", &fGenRecMap);
    fTreeHyp3->Branch("RecDe", &fRecDe);
    fTreeHyp3->Branch("RecPr", &fRecPr);
    fTreeHyp3->Branch("RecPi", &fRecPi);
  }

  fCosPAsplineFile = TFile::Open(AliDataFile::GetFileName(fCosPAsplineName).data());
  if (fCosPAsplineFile) {
    fCosPAspline = (TSpline3*)fCosPAsplineFile->Get("cutSpline");
  }

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);

  AliPDG::AddParticlesToPdgDataBase();

}    /// end UserCreateOutputObjects

void AliAnalysisTaskHyperTriton3VtxPerf::UserExec(Option_t *) {
    // set Magnetic field for KF
  KFParticle::SetField(fInputEvent->GetMagneticField());
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent) {
    ::Fatal("AliAnalysisTaskHyperTriton3VtxPerf::UserExec", "AliESDEvent not found.");
    return;
  }

  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent && fMC) {
    ::Fatal("AliAnalysisTaskHyperTriton3VtxPerf::UserExec", "Could not retrieve MC event");
    return;
  }

  if (!fEventCuts.AcceptEvent(esdEvent)) {
    PostData(1, fListHist);
    PostData(2, fTreeHyp3);
    return;
  }

  fGenHyp.clear();
  fGenRecMap.clear();
  fRecDe.clear();
  fRecPr.clear();
  fRecPi.clear();

  double pvPos[3], pvCov[6];
  fEventCuts.GetPrimaryVertex()->GetXYZ(pvPos);
  fEventCuts.GetPrimaryVertex()->GetCovarianceMatrix(pvCov);
  fREvent.fX = pvPos[0];
  fREvent.fY = pvPos[1];
  fREvent.fZ = pvPos[2];
  fREvent.fCent = fEventCuts.GetCentrality();

  fREvent.fTrigger = 0u;
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) fREvent.fTrigger |= kINT7;
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral) fREvent.fTrigger |= kCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral) fREvent.fTrigger |= kSemiCentral;
  fREvent.fTrigger |= esdEvent->GetMagneticField() > 0 ? kPositiveB : 0;

  std::unordered_map<int, int> mcMap;
  if (fMC) {
    double mcVtx[4];
    mcEvent->GetPrimaryVertex()->GetXYZ(mcVtx);
    for (int iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
      AliVParticle *part = mcEvent->GetTrack(iTrack);

      if (!part) {
        ::Warning("AliAnalysisTaskHyperTriton3VtxPerf::UserExec",
                  "Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skipping.", iTrack);
        continue;
      }

      if (std::abs(part->Y()) > 1.) continue;
      if (!IsHyperTriton3(part, mcEvent)) continue;

      double decayVtx[4]{0.0, 0.0, 0.0, 0.0};

      for (int iD = part->GetDaughterFirst(); iD <= part->GetDaughterLast(); ++iD) {
        AliVParticle *daughter = mcEvent->GetTrack(iD);

        if (mcEvent->IsSecondaryFromWeakDecay(iD) && daughter && std::abs(daughter->PdgCode()) != 11) {
          decayVtx[0] = daughter->Xv();
          decayVtx[1] = daughter->Yv();
          decayVtx[2] = daughter->Zv();
          decayVtx[3] = daughter->Tv();
          mcVtx[3] = part->Tv();
          break;
        }
      }
      SHyperTriton3 genHyp;
      genHyp.pt = part->Pt();
      genHyp.phi = std::atan2(part->Py(),part->Px());
      genHyp.pz = part->Pz();
      for(int iCoord=0; iCoord<4; iCoord++){
        genHyp.prim_vert = mcVtx[iCoord];
        genHyp.dec_vert = decayVtx[iCoord];
      }
      genHyp.positive = part->PdgCode() > 0;
      mcMap[iTrack] = fGenHyp.size();
      fGenHyp.emplace_back(genHyp);
    }
  }

  std::vector<HelperParticle> helpers[3];

  for (int iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) {
    AliESDtrack *track = esdEvent->GetTrack(iTrack);
    if (!track) continue;

    if (!fTrackCuts.AcceptTrack(track)) continue;

    int iPart = 0;
    if (fMC && fOnlyTrueCandidates) {
      int lab = std::abs(track->GetLabel());
      if (!mcEvent->IsSecondaryFromWeakDecay(lab)) continue;
      AliVParticle *part = mcEvent->GetTrack(lab);
      AliVParticle *moth = mcEvent->GetTrack(part->GetMother());
      if (std::abs(moth->PdgCode()) != 1010010030) continue;
      if (std::abs(part->PdgCode()) == 11) continue;
      iPart = GetIndex(part->PdgCode());
    }

    HelperParticle helper;
    helper.track = track;
    helper.nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track, kAliPID[iPart]);
    helper.nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, kAliPID[iPart]);
    helpers[iPart].push_back(helper);
    
  }

  /// downscaling the output tree saving only a fDownscalingFactorByEvent fraction of the events
  if (!fMC && fDownscaling && ((int)helpers[kDeuteron].size() > 0)) {
    if (gRandom->Rndm() > fDownscalingFactorByEvent) return;
  }

  /// if event mixing is enabled takes deuteron from the event mixing pool
  // if (fEnableEventMixing && fApplyML) {
  //   deuterons = GetEventMixingTracks(fREvent.fCent, fREvent.fZ);
  // }

  for (const auto &de : helpers[kDeuteron]) {
    for (const auto &pr : helpers[kProton]) {
      if (de.track == pr.track || pr.particle.GetQ() * de.particle.GetQ() < 0)
        continue;
      for (const auto &pi : helpers[kPion]) {
        if (pr.track == pi.track || de.track == pi.track || pi.particle.GetQ() * pr.particle.GetQ() > 0) 
          continue;

        bool record{!fMC || !fOnlyTrueCandidates};
        if (fMC) {
          int momId = IsTrueHyperTriton3Candidate(deu.track, pr.track, pi.track, mcEvent);
          record = record || momId >=0;
          if (record) {
            fGenRecMap.push_back(mcMap[momId]);
            fRecDe.push_back(de.track);
            fRecPr.push_back(pr.track);
            fRecPi.push_back(pi.track);
          }
        }
      }
    }
  }

  /// if event mixing is enabled fill the event mixing pool with deuterons
  // if (fEnableEventMixing && fApplyML) {
  //   FillEventMixingPool(fREvent.fCent, fREvent.fZ, fDeuVector);
  // }

  fTreeHyp3->Fill();

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);
}

void AliAnalysisTaskHyperTriton3VtxPerf::Terminate(Option_t *) {}

int AliAnalysisTaskHyperTriton3VtxPerf::FindEventMixingCentBin(const float centrality) {
  if (centrality > 90) return -999;
  return static_cast<int>(centrality / 10);
}

int AliAnalysisTaskHyperTriton3VtxPerf::FindEventMixingZBin(const float zvtx) {
  if (zvtx > 10. || zvtx < -10.) return -999.;
  return static_cast<int>((zvtx + 10.) / 2);
}

void AliAnalysisTaskHyperTriton3VtxPerf::FillEventMixingPool(const float centrality, const float zvtx,
                                                        std::vector<AliESDtrack *> tracks) {
  int centBin = FindEventMixingCentBin(centrality);
  int zBin    = FindEventMixingZBin(zvtx);

  auto &trackVector = fEventMixingPool[centBin][zBin];

  for (auto &t : tracks) {
    trackVector.emplace_back(AliESDtrack{*t});
  }

  if (trackVector.size() - fEventMixingPoolDepth > 0) trackVector.pop_front();

  return;
}

std::vector<AliESDtrack *> AliAnalysisTaskHyperTriton3VtxPerf::GetEventMixingTracks(const float centrality,
                                                                               const float zvtx) {
  int centBin = FindEventMixingCentBin(centrality);
  int zBin    = FindEventMixingZBin(zvtx);

  std::vector<AliESDtrack *> tmpVector;

  for (auto &v : fEventMixingPool[centBin][zBin]) {
    tmpVector.emplace_back(&v);
  }

  return tmpVector;
}

AliAnalysisTaskHyperTriton3VtxPerf *AliAnalysisTaskHyperTriton3VtxPerf::AddTask(bool isMC, TString suffix) {
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHyperTriton2BodyML", "No analysis manager found.");
    return nullptr;
  }
  mgr->SetDebugLevel(2);

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHyperTriton3VtxPerf", "This task requires an input event handler");
    return nullptr;
  }

  TString tskname = "AliAnalysisTaskHyperTriton3VtxPerf";
  tskname.Append(suffix.Data());
  AliAnalysisTaskHyperTriton3VtxPerf *task = new AliAnalysisTaskHyperTriton3VtxPerf(isMC, tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("HyperTritonTree%s", suffix.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, Form("HyperTritonTree.root:%s", suffix.Data()));
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  return task;
}

#ifndef AliAnalysisTaskHyperTriton3VtxPerf_H
#define AliAnalysisTaskHyperTriton3VtxPerf_H

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "Math/Vector4D.h"

#include <TString.h>

#include <list>
#include <map>
#include <string>
#include <vector>

class TH1D;
class TH2D;
class TList;
class TTree;
class TFile;
class TSpline3;

class AliPIDResponse;
class AliESDtrack;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LVector_t;

struct REvent3KF
{
  float fX = -999.f;
  float fY = -999.f;
  float fZ = -999.f;
  float fCent = -999.f;
  float fMag_field = -999.f;
  unsigned char fTrigger;
};

struct SHyperTriton3 {
  float pt = -999.f;
  float phi = -999.f;
  float pz = -999.f;
  float dec_vert[4] = {-999.f, -999.f, -999.f, -999.f};
  bool positive = false;
};

class AliAnalysisTaskHyperTriton3VtxPerf : public AliAnalysisTaskSE {

public:
  enum kReducedTrigger { kINT7 = BIT(0), kCentral = BIT(1), kSemiCentral = BIT(2), kPositiveB = BIT(3) };

  AliAnalysisTaskHyperTriton3VtxPerf(bool mc = false, std::string name = "HyperTriton3");
  virtual ~AliAnalysisTaskHyperTriton3VtxPerf();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  static AliAnalysisTaskHyperTriton3VtxPerf *AddTask(bool isMC = false, TString suffix = "");

  void SetDownscaling(bool down) { fDownscaling = down; }

  void SetDownscalingFactorByEvent(float fraction) { fDownscalingFactorByEvent = fraction; }
  void SetDownscalingFactorByCandidate(float fraction) { fDownscalingFactorByCandidate = fraction; }

  void SetEnableEventMixing(bool enableEM) { fEnableEventMixing = enableEM; }
  void SetEventMixingPoolDepth(int maxDepth) { fEventMixingPoolDepth = maxDepth; }

  AliEventCuts fEventCuts;                  /// Event cuts class
  AliESDtrackCuts fTrackCuts = *AliESDtrackCuts::GetStandardV0DaughterCuts(); /// Track cuts Object

  enum kProng { kDeuteron = 0, kProton = 1, kPion = 2 };
  bool  fSwapSign = false;
  bool  fUseAbsCosPAcut = true;
  bool  fOnlyTrueCandidates = false;
  std::string fCosPAsplineName = "PWGLF/NUCLEX/HypertritonAnalysis/Cuts/spline3.root";


private:

  int FindEventMixingCentBin(const float centrality);
  int FindEventMixingZBin(const float zVtx);
  void FillEventMixingPool(const float centrality, const float xVtx, std::vector<AliESDtrack *> tracks);
  std::vector<AliESDtrack *> GetEventMixingTracks(const float centrality, const float zvtx);

  AliInputEventHandler* fInputHandler = nullptr; //!
  AliPIDResponse*       fPIDResponse = nullptr;  //!

  TList *fListHist = nullptr;    //! List of Cascade histograms
  TTree *fTreeHyp3 = nullptr;    //! Output Tree, V0s

  TFile *fCosPAsplineFile = nullptr;  //! Pointer to the spline file
  TSpline3 *fCosPAspline = nullptr;   //! Pointer to the cosPA cut spline

  bool fMC = false;
  bool fDownscaling = false;
  bool fEnableEventMixing = false;

  /// Control histograms to monitor the filtering
  TH2D *fHistNSigmaDeu = nullptr;    //! # sigma TPC for the deuteron
  TH2D *fHistNSigmaP = nullptr;      //! # sigma TPC proton for the positive prong
  TH2D *fHistNSigmaPi = nullptr;     //! # sigma TPC pion for the negative prong
  TH2D *fHistInvMass = nullptr;      //! # Invariant mass histogram

  float fDownscalingFactorByEvent = 1.;        // fraction of the events saved in the tree
  float fDownscalingFactorByCandidate = 1.;    // fraction of the candidates saved in the tree

  std::list<AliESDtrack> fEventMixingPool[10][10];    /// container for the ESD used fot event mixing
  int fEventMixingPoolDepth = 0;                      /// max depth of the event mixing pool

  REvent3KF                    fREvent;
  std::vector<SHyperTriton3> fGenHyp;
  std::vector<int>             fGenRecMap;
  std::vector<AliESDtrack*>    fRecDe;
  std::vector<AliESDtrack*>    fRecPr;
  std::vector<AliESDtrack*>    fRecPi;

  AliAnalysisTaskHyperTriton3VtxPerf(const AliAnalysisTaskHyperTriton3VtxPerf &);               // not implemented
  AliAnalysisTaskHyperTriton3VtxPerf &operator=(const AliAnalysisTaskHyperTriton3VtxPerf &);    // not implemented

  ClassDef(AliAnalysisTaskHyperTriton3VtxPerf, 1);
};

#endif
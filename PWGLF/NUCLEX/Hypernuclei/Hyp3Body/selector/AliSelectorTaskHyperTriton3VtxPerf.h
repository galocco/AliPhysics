#ifndef AliSelectorTaskHyperTriton3VtxPerf_h
#define AliSelectorTaskHyperTriton3VtxPerf_h

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliVertexerHyperTriton3Body.h"
#include "AliVertexerTracks.h"
#include <TVector3.h>
#include <TStopwatch.h>
#include <vector>
#include "../AliAnalysisTaskHyperTriton3VtxPerf.h"
#include <KFParticle.h>
#include "DCAFitterN.h"

class TH1D;
class TH2D;
class TH3D;

struct RParticles {
  KFParticle pion;
  KFParticle proton;
  KFParticle deuteron;
  KFParticle hypertriton;
};

struct RHyperTritonO2 {
  float fCent = -1.;
  float pt = -999.f;
  float phi = -999.f;
  float pz = -999.f;
  float ct = -1.f;
  float r = -1.f;
  float cosPA = -2.f;
  float m = -1;
  float cosPA_Lambda = -2.; 
  Double32_t mppi_vert = -1.; //[1.077,1.203,8]
  Double32_t mppi = -1.;      //[1.077,1.203,8]
  Double32_t dca_lambda_hyper = -1.0; //[0.0,8.0,8]
  Double32_t dca_de = -1.0; //[0.0,8.0,8]
  Double32_t dca_pr = -1.0; //[0.0,8.0,8]
  Double32_t dca_pi = -1.0; //[0.0,8.0,8]
  Double32_t tpcNsig_de = -4.0; //[-4.0,4.0,8]
  Double32_t tpcNsig_pr = -4.0; //[-4.0,4.0,8]
  Double32_t tpcNsig_pi = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_de = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_pr = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_pi = -4.0; //[-4.0,4.0,8]
  Double32_t dca_de_pr = -4.0; //[0.0,8.0,8]
  Double32_t dca_de_pi = -4.0; //[0.0,8.0,8]
  Double32_t dca_pr_pi = -4.0; //[0.0,8.0,8]
  Double32_t dca_de_sv = -4.0; //[0.0,8.0,8]
  Double32_t dca_pr_sv = -4.0; //[0.0,8.0,8]
  Double32_t dca_pi_sv = -4.0; //[0.0,8.0,8]
  Double32_t chi2 = -1.f;      //[0.0,16.,16]
  UChar_t tpcClus_de = 0u;
  UChar_t tpcClus_pr = 0u;
  UChar_t tpcClus_pi = 0u;
  UChar_t candidates = 0u;
  UChar_t fTrigger = 0u;
  bool hasTOF_de = false;
  bool hasTOF_pr = false;
  bool hasTOF_pi = false;
  bool positive = false;
};

class AliSelectorTaskHyperTriton3VtxPerf : public TSelector {
public:
  TTreeReader fReader; //! the tree reader
  TTree *fChain = 0;   //! pointer to the analyzed TTree or TChain

  // Readers to access the data (delete the ones you do not need).

  TTreeReaderValue<REvent> fREvent = {fReader, "REvent"};

  TTreeReaderArray<int> fHypGenMap = {fReader, "SHypGenMap"};

  TTreeReaderArray<RCandidate> fCandidate = {fReader, "RCandidate"};

  TTreeReaderArray<SHyperTriton3> fTrueHyp = {fReader, "SHypeTriton"};

  AliSelectorTaskHyperTriton3VtxPerf(TString outputName = "output.root", TString outputPath = ".", int vertexer = 0,
                                      TTree * /*tree*/ = 0);
  virtual ~AliSelectorTaskHyperTriton3VtxPerf() {}
  AliSelectorTaskHyperTriton3VtxPerf(const AliSelectorTaskHyperTriton3VtxPerf &) = delete;
  AliSelectorTaskHyperTriton3VtxPerf &operator=(const AliSelectorTaskHyperTriton3VtxPerf &other) = delete;
  virtual Int_t Version() const { return 2; }
  virtual void Begin(TTree *tree);
  virtual void SlaveBegin(TTree *tree);
  virtual void Init(TTree *tree);
  virtual Bool_t Notify();
  virtual Bool_t Process(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
    return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
  }
  virtual void SetOption(const char *option) { fOption = option; }
  virtual void SetObject(TObject *obj) { fObject = obj; }
  virtual void SetInputList(TList *input) { fInput = input; }
  virtual TList *GetOutputList() const { return fOutput; }
  virtual void SlaveTerminate();
  virtual void Terminate();

  TString fOutputFileName;
  TString fOutputFilePath;
  int fAlg;
  //fAlg select the vertexer:
  //0 ->kf
  //1 ->02
  //2 ->std
  //
  AliVertexerHyperTriton3Body fHypertritonVertexer;
  AliESDtrackCuts *fESDtrackCuts = nullptr;
  float fNsigmaTPC[3];
  float fNsigmaTOF[3];
  int fNclones = -1;
  int fNcycles = 0;
  int fNrec = 0;
  TH1D* fOperations = nullptr;
  TH1D* fTotTime = nullptr;
  TStopwatch* fTimer = nullptr;
  REvent fEvent;
  // Histogram for efficiencies
  TH1D *fHistCtRec[2] = {nullptr};
  TH3D *fHistClonesVsCuts[3][2] = {{nullptr}};

  TH3D *fHistResolutionVsCuts[3][2] = {{nullptr}};
  TH3D *fHistSingleRecVsCuts[3][2] = {{nullptr}};
  TH2D *fHistGen[2] = {nullptr};

  TH1D *fHistGenPt = nullptr;
  TH1D *fHistRecPt = nullptr;
  // Histograms for selection
  TH2D *fHistInvMassPt[2][3] = {{nullptr}};
  TH2D *fHistInvMassPtSel[2][3] = {{nullptr}};
  TH2D *fHistMassResPt[2] = {nullptr};//[Matter]
  TH2D *fHistMassResCt[2] = {nullptr};//[Matter]
  TH2D *fHistCtResCt[2] = {nullptr};//[Matter]
  TH2D *fHistCtResCtTrueP[2] = {nullptr};//[Matter]
  TH2D *fHistPtResPt[2] = {nullptr};//[Matter]
  TH2D *fHistPResP[2] = {nullptr};//[Matter]
  TH2D *fHistXResX[2] = {nullptr};//[Matter]
  TH2D *fHistYResY[2] = {nullptr};//[Matter]
  TH2D *fHistZResZ[2] = {nullptr};//[Matter]
  TH1D *fHistPDG[3] = {nullptr};//[Matter]

  TH1D *fHistDaughterPt[3][3] = {{nullptr}};
  TH1D *fHistDaughterTPCchi2[3][3] = {{nullptr}};
  TH1D *fHistDaughterITSchi2[3][3] = {{nullptr}};
  TH1D *fHistNclsITS[3][3] = {{nullptr}};
  TH1D *fHistNclsTPC[3][3] = {{nullptr}};
  TH1D *fHistNSigmaTPC[3][3] = {{nullptr}};
  TH1D *fHistNSigmaTOF[3][3] = {{nullptr}};

  // Histograms after vertexer
  TH1D *fHistVertexChi2[2] = {nullptr};
  TH1D *fHistCosPAngle[2] = {nullptr};
  TH1D *fHist2ProngChi2[2] = {nullptr};
  TH1D *fHist3ProngChi2[2] = {nullptr};
  TH1D *fHistVertChi2[2] = {nullptr};
  TH1D *fHistDCA2pvXY[2][3] = {{nullptr}};
  TH1D *fHistDCA2pvZ[2][3] = {{nullptr}};
  TH1D *fHistDCA2pv[2][3] = {{nullptr}};
  TH1D *fHistDCA2dvXY[2][3] = {{nullptr}};
  TH1D *fHistDCA2dvZ[2][3] = {{nullptr}};
  TH1D *fHistDCA2dv[2][3] = {{nullptr}};
  TH1D *fHistTrackDistance[2][3] = {{nullptr}};

  bool AcceptCandidate(RCandidate, int, int);
  bool KFVertexer(AliESDtrack* [], RParticles &, double [], float []);
  bool O2Vertexer(AliESDtrack* [], RHyperTritonO2 &, double [],float ,float []);

  ClassDef(AliSelectorTaskHyperTriton3VtxPerf, 0);
};

#endif

//______________________________________________________________________________
#ifdef AliSelectorTaskHyperTriton3VtxPerf_cxx

void AliSelectorTaskHyperTriton3VtxPerf::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // AliSelectorTaskHyperTriton3VtxPerf, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);
}

//______________________________________________________________________________
Bool_t AliSelectorTaskHyperTriton3VtxPerf::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated AliSelectorTaskHyperTriton3VtxPerf, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef AliSelectorTaskHyperTriton3VtxPerf_cxx

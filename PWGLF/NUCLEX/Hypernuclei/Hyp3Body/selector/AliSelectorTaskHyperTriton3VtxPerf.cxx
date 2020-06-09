#define AliSelectorTaskHyperTriton3VtxPerf_cxx

#include "Hyp3FindConfig.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliPID.h"
#include "AliVTrack.h"
#include <TCanvas.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>
#include <cmath>
#include <KFVertex.h>
#include <KFParticle.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include <stdio.h>
#include "DCAFitterN.h"
#include "AliSelectorTaskHyperTriton3VtxPerf.h"

const float kHypMass{2.99131};
const float kDeuMass{1.87561};
const float kPMass{0.938272};
const float kPiMass{0.13957};
/// usefull functions

template <typename T> double Sq(T a) { return a * a; }

template <typename T> double Point2PointDistance(T *p0, T *p1) {
  double d2 = 0.;
  for (int iDim = 0; iDim < 3; ++iDim) {
    d2 += Sq(p0[iDim] - p1[iDim]);
  }
  return std::sqrt(d2);
}

template <typename T> double Norm(T x, T y) { return std::sqrt(Sq(x) + Sq(y)); }

template <typename T> double Norm(T x, T y, T z) { return std::sqrt(Sq(x) + Sq(y) + Sq(z)); }

AliSelectorTaskHyperTriton3VtxPerf::AliSelectorTaskHyperTriton3VtxPerf(TString outputName, TString outputPath, int vertexer ,TTree *):
fOutputFileName{outputName},
fOutputFilePath{outputPath},
fAlg{vertexer},
fHypertritonVertexer(){
  fESDtrackCuts = AliESDtrackCuts::GetStandardV0DaughterCuts();
  fESDtrackCuts->SetMinNClustersTPC(0);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);
  //test to compare with the other selector

}

//______________________________________________________________________________
void AliSelectorTaskHyperTriton3VtxPerf::Begin(TTree * /*tree*/){
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  TString option = GetOption();
}

//______________________________________________________________________________
void AliSelectorTaskHyperTriton3VtxPerf::SlaveBegin(TTree * /*tree*/){
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  std::cout<<"ciao";
  const char lAM[3]{"AM"};
  const char lRL[3]{"RL"};
  const char *lSpecies[3]{"d", "p", "pi"};
  const char *lCuts[3]{"tree", "selection", "vertexer"};
  const char *lType[3]{"true", "fake", "clones"};
  const char *lProj[2]={"pT","ct"};
  const char *lVarName[3]={"NsigmaTPC", "NclusterTPC","NclusterITS"};
  
  fTimer = new TStopwatch();
  fTimer->Reset();
  fHistGenPt= new TH1D("fHistGenPt",";p_{T} [GeV/c];",10,0,10);
  GetOutputList()->Add(fHistGenPt);
  fHistRecPt= new TH1D("fHistRecPt",";p_{T} [GeV/c];",10,0,10);
  GetOutputList()->Add(fHistRecPt);
  fTotTime= new TH1D("fTotTime",";;time [s]",1,0,1);
  GetOutputList()->Add(fTotTime);
  fOperations= new TH1D("fOperations",";;cycles",2,0,2);
  GetOutputList()->Add(fOperations);
  /// histograms for efficiencies
  for(int iMatter=0; iMatter<2; iMatter++){
    fHistCtRec[iMatter] = new TH1D(Form("fHistCtRec%c",lAM[iMatter]),";;",10,0,50);

    fHistMassResPt[iMatter] = new TH2D(Form("fHistMassRes_%s_%c",lProj[0], lAM[iMatter]), ";#it{m_{rec}}-#it{m_{gen}} (dp#pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", 100,-0.05, 0.05, 10, 0, 10);
    fHistMassResCt[iMatter] = new TH2D(Form("fHistMassRes_%s_%c",lProj[1], lAM[iMatter]), ";#it{m_{rec}}-#it{m_{gen}} (dp#pi) (GeV/#it{c}^{2}); #it{c}t (cm);", 100, -0.05, 0.05, 10, 0, 50);
    fHistPtResPt[iMatter] = new TH2D(Form("fHistPtRes_pT_%c", lAM[iMatter]), ";#it{p}_{Trec}-#it{p}_{Tgen} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});", 200, -0.5, 0.5, 10, 0, 10);
    fHistCtResCt[iMatter] = new TH2D(Form("fHistCtRes_ct_%c", lAM[iMatter]), ";#it{c}t_{rec}-#it{c}t_{gen} (cm); #it{c}t (cm);", 400, -20., 20., 10, 0, 50);
    fHistCtResCtTrueP[iMatter] = new TH2D(Form("fHistCtRes_ct_TrueP_%c", lAM[iMatter]), ";#it{c}t_{rec}-#it{c}t_{gen} (cm); #it{c}t (cm);", 400, -20., 20., 10, 0, 50);
    fHistPResP[iMatter] = new TH2D(Form("fHistPRes_p_%c", lAM[iMatter]), ";#it{p}_{rec}-#it{p}_{gen} (GeV/#it{c});#it{p} (GeV/#it{c});", 200, -0.5, 0.5, 10, 0, 10);
    fHistXResX[iMatter] = new TH2D(Form("fHistXRes_X_%c", lAM[iMatter]), ";X_{rec}-X_{gen} (cm); X (cm);", 400, -20., 20., 10, 0, 50);
    fHistYResY[iMatter] = new TH2D(Form("fHistYRes_Y_%c", lAM[iMatter]), ";Y_{rec}-Y_{gen} (cm); Y (cm);", 400, -20., 20., 10, 0, 50);
    fHistZResZ[iMatter] = new TH2D(Form("fHistZRes_Z_%c", lAM[iMatter]), ";Z_{rec}-Z_{gen} (cm); Z (cm);", 400, -20., 20., 10, 0, 50);

    fHistGen[iMatter] = new TH2D(Form("fHistGen_%c",lAM[iMatter]),";#it{p}_{T} (GeV/#it{c});#it{ct} (cm)",kNBinPt,0.,10.,kNBinCt,0.,50.);

    GetOutputList()->Add(fHistMassResPt[iMatter]);
    GetOutputList()->Add(fHistMassResCt[iMatter]);
    GetOutputList()->Add(fHistPtResPt[iMatter]);
    GetOutputList()->Add(fHistCtResCt[iMatter]);
    GetOutputList()->Add(fHistCtResCtTrueP[iMatter]);
    GetOutputList()->Add(fHistPResP[iMatter]);
    GetOutputList()->Add(fHistXResX[iMatter]);
    GetOutputList()->Add(fHistYResY[iMatter]);
    GetOutputList()->Add(fHistZResZ[iMatter]);
    GetOutputList()->Add(fHistGen[iMatter]);     
    GetOutputList()->Add(fHistCtRec[iMatter]);

    for(int iCut=0; iCut<kNcuts; iCut++){
      fHistClonesVsCuts[iCut][iMatter] = new TH3D(Form("fHistClones_%s_%c",lVarName[iCut],lAM[iMatter]),";#it{p}_{T} (GeV/#it{c});#it{ct} (cm);",kNBinPt,0.,10.,kNBinCt,0.,50.,kNvariations,0,kNvariations);
      fHistSingleRecVsCuts[iCut][iMatter] = new TH3D(Form("fHistRec_%s_%c",lVarName[iCut],lAM[iMatter]),";#it{p}_{T} (GeV/#it{c});#it{ct} (cm);",kNBinPt,0.,10.,kNBinCt,0.,50.,kNvariations,0,kNvariations);
      fHistResolutionVsCuts[iCut][iMatter] = new TH3D(Form("fHistRes_%s_%c",lVarName[iCut],lAM[iMatter]),";#Delta#it{p}_{T} (GeV/#it{c});#Delta#it{ct} (cm);",100,-2.,2.,100,-15.,15.,kNvariations,0,kNvariations);

      GetOutputList()->Add(fHistSingleRecVsCuts[iCut][iMatter]);
      GetOutputList()->Add(fHistResolutionVsCuts[iCut][iMatter]); 
      GetOutputList()->Add(fHistClonesVsCuts[iCut][iMatter]);
    }
  }

  /// Histograms for selection

  for(int iCuts = 0; iCuts < 3; iCuts++){
    for(int iMatter = 0; iMatter < 2; iMatter++){
      fHistInvMassPt[iMatter][iCuts] = new TH2D(Form("fHistInvMassPt_%c_%s", lAM[iMatter], lType[iCuts]), ";#it{m} (dp#pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", 200, 2.95, 3.35, 20, 0, 10);
      GetOutputList()->Add(fHistInvMassPt[iMatter][iCuts]);
      fHistInvMassPtSel[iMatter][iCuts] = new TH2D(Form("fHistInvMassPt_%c_%s", lAM[iMatter], lCuts[iCuts]), ";#it{m} (dp#pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", 200, 2.95, 3.35, 20, 0, 10);
      GetOutputList()->Add(fHistInvMassPtSel[iMatter][iCuts]);
    }
    for(int iSpecies = 0; iSpecies < 3; iSpecies++){
      fHistDaughterPt[iSpecies][iCuts] = new TH1D(Form("fHistPt_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";#it{p}_{T} (GeV/#it{c});Counts", 100, 0, 10);
      fHistDaughterTPCchi2[iSpecies][iCuts] = new TH1D(Form("fHistTPCchi2_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";#chi^{2};Counts", 100, 0, 10);
      fHistDaughterITSchi2[iSpecies][iCuts] = new TH1D(Form("fHistITSchi2_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";#chi^{2};Counts", 100, 0, 10);
      fHistNclsITS[iSpecies][iCuts] = new TH1D(Form("fHistNclsITS_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";N_{cluster ITS};Counts", 7, -0.5, 6.5);
      fHistNclsTPC[iSpecies][iCuts] = new TH1D(Form("fHistNclsTPC_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";N_{cluster TPC};Counts", 101, -0.5, 200.5);
      fHistNSigmaTPC[iSpecies][iCuts] = new TH1D(Form("fHistNSigmaTPC_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), "n#sigma_{TPC}", 10, 0, 10);
      fHistNSigmaTOF[iSpecies][iCuts] = new TH1D(Form("fHistNSigmaTOF_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), "n#sigma_{TOF}", 10, 0, 10);
      GetOutputList()->Add(fHistDaughterPt[iSpecies][iCuts]);
      GetOutputList()->Add(fHistDaughterTPCchi2[iSpecies][iCuts]);
      GetOutputList()->Add(fHistDaughterITSchi2[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNclsITS[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNclsTPC[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNSigmaTPC[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNSigmaTOF[iSpecies][iCuts]);
    } 
  }
  /// Histograms after vertexer


  for(int iSide = 0; iSide < 2; iSide++){
    fHistVertexChi2[iSide] = new TH1D(Form("fHistVertexChi2_%c",lRL[iSide]), "", 100, 0, 200);
    fHistCosPAngle[iSide]  = new TH1D(Form("fCosPointingAngle_%c",lRL[iSide]), ";#it{cos#theta_{pointing}};Counts", 5000, 0.5, 1.);
    if(fAlg==0){
      fHist2ProngChi2[iSide]  = new TH1D(Form("fHist2ProngChi2_%c",lRL[iSide]), ";#chi^{2};Counts", 5000, 0., 50.);
      fHist3ProngChi2[iSide]  = new TH1D(Form("fHist3ProngChi2_%c",lRL[iSide]), ";#chi^{2};Counts", 5000, 0., 50.);
      fHistVertChi2[iSide]  = new TH1D(Form("fHistVertChi2_%c",lRL[iSide]), ";#chi^{2};Counts", 5000, 0., 50.);
      GetOutputList()->Add(fHist2ProngChi2[iSide]);
      GetOutputList()->Add(fHist3ProngChi2[iSide]);
      GetOutputList()->Add(fHistVertChi2[iSide]);
    }
    GetOutputList()->Add(fHistVertexChi2[iSide]);
    GetOutputList()->Add(fHistCosPAngle[iSide]);
    /// DCA to primary vertex
    for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
      fHistDCA2pvXY[iSide][iSpecies] = new TH1D(Form("fDCA2PrimaryvtxXY_%s_%c", lSpecies[iSpecies],lRL[iSide]), ";DCA_{xy} (mm);Counts", 600, 0, 30);
      fHistDCA2pvZ[iSide][iSpecies]  = new TH1D(Form("fDCA2PrimaryvtxZ_%s_%c", lSpecies[iSpecies],lRL[iSide]), ";DCA_{z} (mm);Counts", 600, 0, 30);
      fHistDCA2pv[iSide][iSpecies]   = new TH1D(Form("fDCA2Primaryvtx_%s_%c", lSpecies[iSpecies],lRL[iSide]), ";DCA (mm);Counts", 600, 0, 30);
      fHistDCA2dvXY[iSide][iSpecies] = new TH1D(Form("fDCA2DecayvtxXY_%s_%c", lSpecies[iSpecies],lRL[iSide]), ";DCA_{xy} (mm);Counts", 600, 0, 30);
      fHistDCA2dvZ[iSide][iSpecies]  = new TH1D(Form("fDCA2DecayvtxZ_%s_%c", lSpecies[iSpecies],lRL[iSide]), ";DCA_{z} (mm);Counts", 600, 0, 30);
      fHistDCA2dv[iSide][iSpecies]   = new TH1D(Form("fDCA2Decayvtx_%s_%c", lSpecies[iSpecies],lRL[iSide]), ";DCA (mm);Counts", 600, 0, 30);
      fHistTrackDistance[iSide][iSpecies] = new TH1D(Form("lTrackDistance_%s-%s_%c",lSpecies[iSpecies], lSpecies[(iSpecies + 1) % 3],lRL[iSide]), ";distance (mm);counts", 200, 0, 200);
      GetOutputList()->Add(fHistDCA2pvXY[iSide][iSpecies]);
      GetOutputList()->Add(fHistDCA2pvZ[iSide][iSpecies]);
      GetOutputList()->Add(fHistDCA2pv[iSide][iSpecies]);
      GetOutputList()->Add(fHistDCA2dvXY[iSide][iSpecies]);
      GetOutputList()->Add(fHistDCA2dvZ[iSide][iSpecies]);
      GetOutputList()->Add(fHistDCA2dv[iSide][iSpecies]);
      GetOutputList()->Add(fHistTrackDistance[iSide][iSpecies]);
    }
  }
  fTimer->Continue();
}

//______________________________________________________________________________
Bool_t AliSelectorTaskHyperTriton3VtxPerf::Process(Long64_t entry){
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetEntry(entry);
  //if(fNcycles>1000000)
  //  return true;
  // Support variables

  


  const float lMasses[3]{AliPID::ParticleMass(AliPID::kDeuteron), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion)};
  constexpr AliPID::EParticleType kAliPID[3]{AliPID::kDeuteron,AliPID::kProton,AliPID::kPion};

  //------------------------------------------------------------
  // Get main observables from tree
  //------------------------------------------------------------
  fEvent = *fREvent;
  float lMagField = fEvent.fMag_field;
  double lTruePrimaryVtx[3]{fEvent.fX,fEvent.fY,fEvent.fZ};
  double lPrimaryVtxCov[6] = {0.};
  double lRecPrimaryVtx[3] = {0.};
  double lRecDecayVtx[3] = {0.};
  double lRecDecayLenght[3] = {0.};

  //------------------------------------------------------------
  //                  Generated hypertritons                   
  //------------------------------------------------------------

  //TODO: put right name
  for (auto &SHyper : fTrueHyp)
  {  
    double lHypPtGen = std::hypot(SHyper.px,SHyper.py); 
    double lHypPGen = std::hypot(lHypPtGen,SHyper.pz);
    double lHypCtGen = Norm(SHyper.dec_vert[0]-lTruePrimaryVtx[0],SHyper.dec_vert[1]-lTruePrimaryVtx[1],SHyper.dec_vert[2]-lTruePrimaryVtx[2])*kHypMass/lHypPGen;
    int lCharge = (SHyper.positive) ? 0 : 1;
    fHistGen[lCharge]->Fill(lHypPtGen,lHypCtGen);
    fHistGenPt->Fill(lHypPtGen);
  }

  //loop on all the triples
  int index = 0;
  std::vector<int> OldMap;
  for (auto &RCand : fCandidate)
  { 
    int clone = 0;
    index++;
    for(auto rec_index : OldMap){
      if(fHypGenMap[index] == rec_index)
        clone++;
    }

    double lTrueDecayVtx[3]{fTrueHyp[fHypGenMap[index]].dec_vert[0],fTrueHyp[fHypGenMap[index]].dec_vert[1],fTrueHyp[fHypGenMap[index]].dec_vert[2]};
    double lHypPtGen = std::hypot(fTrueHyp[fHypGenMap[index]].px,fTrueHyp[fHypGenMap[index]].px); 
    double lHypPGen = std::hypot(lHypPtGen,fTrueHyp[fHypGenMap[index]].pz);  
    double lHypCtGen = Norm(fTrueHyp[fHypGenMap[index]].dec_vert[0]-lTruePrimaryVtx[0],fTrueHyp[fHypGenMap[index]].dec_vert[1]-lTruePrimaryVtx[1],fTrueHyp[fHypGenMap[index]].dec_vert[2]-lTruePrimaryVtx[2])*kHypMass/lHypPGen;

    TLorentzVector lLVhyp = {0., 0., 0., 0.};
    TLorentzVector lLVdaughter[3];

    AliESDtrack *lTrack[3]{&*RCand.track[0],&*RCand.track[1],&*RCand.track[2]};
    int lDaughterNclsITS[3]{lTrack[0]->GetITSNcls(), lTrack[1]->GetITSNcls(), lTrack[2]->GetITSNcls()};
    int lDaughterNclsTPC[3]{lTrack[0]->GetTPCNcls(), lTrack[1]->GetTPCNcls(), lTrack[2]->GetTPCNcls()};
    double lDaughterChi2TPC[3]{lTrack[0]->GetTPCchi2(), lTrack[1]->GetTPCchi2(), lTrack[2]->GetTPCchi2()};
    double lDaughterChi2ITS[3]{lTrack[0]->GetITSchi2(), lTrack[1]->GetITSchi2(), lTrack[2]->GetITSchi2()};
    float lDaughterNsigmaTPC[3] = {RCand.NsigmaTPC[0],RCand.NsigmaTPC[1],RCand.NsigmaTPC[2]};
    float lDaughterNsigmaTOF[3] = {RCand.NsigmaTOF[0],RCand.NsigmaTOF[1],RCand.NsigmaTOF[2]};
    int lCharge = (lTrack[0]->GetSign() < 0)? -1 : 1;


    float lDaughterPt[3] = {0.};
    for (int iTrack = 0; iTrack < 3; iTrack++) {
      lLVdaughter[iTrack].SetXYZM(lTrack[iTrack]->Px(), lTrack[iTrack]->Py(), lTrack[iTrack]->Pz(), lMasses[iTrack]);
      lDaughterPt[iTrack] = lLVdaughter[iTrack].Pt();
      lLVhyp += lLVdaughter[iTrack];
    }


    double lHypMassRec = lLVhyp.M();
    double lHypPtRec   = lLVhyp.Pt();
    //------------------------------------------------------------
    //                    Before selection
    //------------------------------------------------------------

    fHistInvMassPtSel[lCharge][0]->Fill(lHypMassRec,lHypPtRec);
    for(int iTrack = 0; iTrack < 3; iTrack++){
      fHistDaughterPt[iTrack][0]->Fill(lDaughterPt[iTrack]);
      fHistNclsITS[iTrack][0]->Fill(lDaughterNclsITS[iTrack]);
      fHistNclsTPC[iTrack][0]->Fill(lDaughterNclsTPC[iTrack]);
      fHistNSigmaTPC[iTrack][0]->Fill(lDaughterNsigmaTPC[iTrack]);
      fHistNSigmaTOF[iTrack][0]->Fill(lDaughterNsigmaTOF[iTrack]);
      fHistDaughterTPCchi2[iTrack][0]->Fill(lDaughterChi2TPC[iTrack]);
      fHistDaughterITSchi2[iTrack][0]->Fill(lDaughterChi2ITS[iTrack]);
    }

    //------------------------------------------------------------
    //                      After selection
    //------------------------------------------------------------
    
    if(AcceptCandidate(RCand,0,0)){
      fHistInvMassPtSel[lCharge][1]->Fill(lHypMassRec,lHypPtRec);
      for(int iTrack = 0; iTrack < 3; iTrack++){
        fHistDaughterPt[iTrack][1]->Fill(lDaughterPt[iTrack]);
        fHistNclsITS[iTrack][1]->Fill(lDaughterNclsITS[iTrack]);
        fHistNclsTPC[iTrack][1]->Fill(lDaughterNclsTPC[iTrack]);
        fHistNSigmaTPC[iTrack][1]->Fill(lDaughterNsigmaTPC[iTrack]);
        fHistNSigmaTOF[iTrack][1]->Fill(lDaughterNsigmaTOF[iTrack]);
        fHistDaughterTPCchi2[iTrack][1]->Fill(lDaughterChi2TPC[iTrack]);
        fHistDaughterITSchi2[iTrack][1]->Fill(lDaughterChi2ITS[iTrack]);
      }
    }

    //------------------------------------------------------------
    // Secondary vertex reconstruction
    //------------------------------------------------------------

    // reconstruct the decay vertex with the dedicated vertexer
    RParticles lCollection;
    RHyperTritonO2 recHyp;
    fNcycles++;
    bool vertexer_result = false;
    float decVert[3],kfchi2[3];
    if(fAlg==0){//kalman filter vertexer
      
      fTimer->Continue();
      vertexer_result=KFVertexer(lTrack, lCollection,lRecPrimaryVtx,kfchi2);
      fTimer->Stop();  
    }
    else if(fAlg==1){// 02 vertexer
      fTimer->Continue();
      vertexer_result=O2Vertexer(lTrack, recHyp, lRecPrimaryVtx, lMagField, decVert);
      fTimer->Stop();
    }
    else if(fAlg==2){//standard vertexer
      fTimer->Continue();
      vertexer_result=fHypertritonVertexer.FindDecayVertex(lTrack[0], lTrack[1], lTrack[2], lMagField);
      fTimer->Stop();
    }
    if(!vertexer_result) continue;
    fNrec++;
    
    AliESDVertex *lDecayVertex = nullptr;
    TVector3 vRecDecayLenght;
    float cospa=0;
    float lHypPRec=0,lHypCtRec=0,lHypCtRecTrueP=0,lHypXRec=0,lHypYRec=0,lHypZRec=0;
    float lHypXGen = lTrueDecayVtx[0];
    float lHypYGen = lTrueDecayVtx[1];
    float lHypZGen = lTrueDecayVtx[2];
    
    //the cut on the cosPA is applied only to the histograms to compute efficiency and resolution
      if(fAlg==2){//std vertexer
      lDecayVertex = static_cast<AliESDVertex *>(fHypertritonVertexer.GetCurrentVertex());   
      lDecayVertex->GetXYZ(lRecDecayVtx);
      for (int iCoord = 0; iCoord < 3; iCoord++) {
        lRecDecayLenght[iCoord] = lRecDecayVtx[iCoord] - lRecPrimaryVtx[iCoord];
      }
      vRecDecayLenght = TVector3(lRecDecayLenght[0], lRecDecayLenght[1], lRecDecayLenght[2]);
      lHypCtRec = vRecDecayLenght.Mag()*kHypMass/lLVhyp.P();
      lHypCtRecTrueP = vRecDecayLenght.Mag()*kHypMass/lHypPGen;
      fHistCtRec[lCharge]->Fill(lHypCtRec);

      lHypPRec = lLVhyp.P();
      lHypXRec = lRecDecayVtx[0];
      lHypYRec = lRecDecayVtx[1];
      lHypZRec = lRecDecayVtx[2];
    }
    else if(fAlg==0){//kf vertexer
      ROOT::Math::XYZVectorF decayl{lCollection.hypertriton.X()-(float)lRecPrimaryVtx[0], lCollection.hypertriton.Y()-(float)lRecPrimaryVtx[1], lCollection.hypertriton.Z()-(float)lRecPrimaryVtx[2]};
      ROOT::Math::XYZVectorF mom{lCollection.hypertriton.Px(), lCollection.hypertriton.Py(), lCollection.hypertriton.Pz()};
      lHypPRec = std::sqrt(mom.Mag2());
      cospa = mom.Dot(decayl) / std::sqrt(decayl.Mag2() * mom.Mag2());
      lHypPtRec = TMath::Sqrt(lCollection.hypertriton.Px()*lCollection.hypertriton.Px()+lCollection.hypertriton.Py()*lCollection.hypertriton.Py());
      lHypCtRec = std::sqrt(decayl.Mag2())*kHypMass/lHypPRec;
      lHypCtRecTrueP = std::sqrt(decayl.Mag2())*kHypMass/lHypPGen;
      lHypMassRec = lCollection.hypertriton.GetMass();
      fHistCtRec[lCharge]->Fill(lHypCtRec);
      
      lHypXRec = lCollection.hypertriton.X();
      lHypYRec = lCollection.hypertriton.Y();
      lHypZRec = lCollection.hypertriton.Z();
    }
    else {//O2 vertexer
      fHistCtRec[lCharge]->Fill(recHyp.ct);
      lHypMassRec=recHyp.m;
      lHypCtRec=recHyp.ct;
      lHypCtRecTrueP = recHyp.ct/lHypPGen*recHyp.pz;
      lHypPtRec=recHyp.pt;
      lHypPRec = recHyp.pz;//not a bug but a way to save time
      lHypXRec = decVert[0];
      lHypYRec = decVert[1];
      lHypZRec = decVert[2];
    }


    int ctside = ((lHypCtRec-lHypCtGen)>=0) ? 0 : 1;

    if(AcceptCandidate(RCand,0,0)){
    
      //this part is not implemented for the o2
      if(fAlg==0){

        if(lHypCtGen<20 && lHypCtGen>10){

        fHist2ProngChi2[ctside]->Fill(kfchi2[0]);
        fHist3ProngChi2[ctside]->Fill(kfchi2[1]);
        fHistVertChi2[ctside]->Fill(kfchi2[2]);
        //float dec_vertex[3] = {lCollection.hypertriton.X(), lCollection.hypertriton.Y(), lCollection.hypertriton.Z()};
        float dec_vertex[3] = {(float)lTrueDecayVtx[0], (float)lTrueDecayVtx[1], (float)lTrueDecayVtx[2]};
        //float pri_vertex[3] = {lRecPrimaryVtx[0], lRecPrimaryVtx[1], lRecPrimaryVtx[2]};
        float pri_vertex[3] = {(float)lTruePrimaryVtx[0], (float)lTruePrimaryVtx[1], (float)lTruePrimaryVtx[2]};
        double dcaXY[3],dca[3],dcaXYpv[3],dcapv[3],dcaPP[3];
        dcaXY[0] = 10.*lCollection.deuteron.GetDistanceFromVertexXY(dec_vertex);
        dca[0] = 10.*lCollection.deuteron.GetDistanceFromVertex(dec_vertex);
        dcaXY[1] = 10.*lCollection.proton.GetDistanceFromVertexXY(dec_vertex);
        dca[1] = 10.*lCollection.proton.GetDistanceFromVertex(dec_vertex);
        dcaXY[2] = 10.*lCollection.pion.GetDistanceFromVertexXY(dec_vertex);
        dca[2] = 10.*lCollection.pion.GetDistanceFromVertex(dec_vertex);

        dcaXY[0] = 10.*lCollection.deuteron.GetDistanceFromVertexXY(dec_vertex);
        dca[0] = 10.*lCollection.deuteron.GetDistanceFromVertex(dec_vertex);
        dcaXY[1] = 10.*lCollection.proton.GetDistanceFromVertexXY(dec_vertex);
        dca[1] = 10.*lCollection.proton.GetDistanceFromVertex(dec_vertex);
        dcaXY[2] = 10.*lCollection.pion.GetDistanceFromVertexXY(dec_vertex);
        dca[2] = 10.*lCollection.pion.GetDistanceFromVertex(dec_vertex);

        dcaPP[0] = 10.*lCollection.deuteron.GetDistanceFromParticle(lCollection.proton);
        dcaPP[1] = 10.*lCollection.pion.GetDistanceFromParticle(lCollection.deuteron);
        dcaPP[2] = 10.*lCollection.pion.GetDistanceFromParticle(lCollection.proton);
        
        dcaXYpv[0] = 10.*lCollection.deuteron.GetDistanceFromVertexXY(pri_vertex);
        dcapv[0] = 10.*lCollection.deuteron.GetDistanceFromVertex(pri_vertex);
        dcaXYpv[1] = 10.*lCollection.proton.GetDistanceFromVertexXY(pri_vertex);
        dcapv[1] = 10.*lCollection.proton.GetDistanceFromVertex(pri_vertex);
        dcaXYpv[2] = 10.*lCollection.pion.GetDistanceFromVertexXY(pri_vertex);
        dcapv[2] = 10.*lCollection.pion.GetDistanceFromVertex(pri_vertex);

        //still to check the units
        for(int iPar=0; iPar<3; iPar++){
          fHistDCA2dvXY[ctside][iPar]->Fill(dcaXY[iPar]);
          fHistDCA2dvZ[ctside][iPar]->Fill(TMath::Sqrt(dca[iPar]*dca[iPar]-dcaXY[iPar]*dcaXY[iPar]));
          fHistDCA2dv[ctside][iPar]->Fill(dca[iPar]);
          fHistTrackDistance[ctside][iPar]->Fill(dcaPP[iPar]);


          fHistDCA2pvXY[ctside][iPar]->Fill(dcaXYpv[iPar]);
          fHistDCA2pvZ[ctside][iPar]->Fill(TMath::Sqrt(dcapv[iPar]*dcapv[iPar]-dcaXYpv[iPar]*dcaXYpv[iPar]));
          fHistDCA2pv[ctside][iPar]->Fill(dcapv[iPar]);
        }
        fHistVertexChi2[ctside]->Fill(lCollection.hypertriton.GetChi2()/lCollection.hypertriton.GetNDF());
        }
      }
      
      //this part is implemented only for the std vertexer

      if(fAlg==2){
        float pointAngle = lLVhyp.Angle(vRecDecayLenght);
        cospa = std::cos(pointAngle);
        /// compute the DCA of the 3 tracks from the primary and decay vertex
        AliESDVertex lPV(lTruePrimaryVtx, lPrimaryVtxCov, 1., 1000);
        for (int iTrack = 0; iTrack < 3; iTrack++) {
          double dca2dv[2]    = {0.};
          double dca2pv[2]    = {0.};
          double dca2dvcov[3] = {0.};
          double dca2pvcov[3] = {0.};
          lTrack[iTrack]->PropagateToDCA(lDecayVertex, lMagField, 1000., dca2dv, dca2dvcov);
          lTrack[iTrack]->PropagateToDCA(&lPV, lMagField, 1000., dca2pv, dca2pvcov);

          float dcaXYdv = std::abs(dca2dv[0]) * 10.;    // in mm
          float dcaZdv  = std::abs(dca2dv[1]) * 10.;    // in mm
          float dcadv   = Norm(dcaXYdv, dcaZdv) * 10.;  // in mm

          fHistDCA2dvXY[ctside][iTrack]->Fill(dcaXYdv);
          fHistDCA2dvZ[ctside][iTrack]->Fill(dcaZdv);
          fHistDCA2dv[ctside][iTrack]->Fill(dcadv);

          float dcaXYpv = std::abs(dca2pv[0]) * 10.;    // in mm
          float dcaZpv  = std::abs(dca2pv[1]) * 10.;    // in mm
          float dcapv   = Norm(dcaXYpv, dcaZpv) * 10.;  // in mm

          fHistDCA2pvXY[ctside][iTrack]->Fill(dcaXYpv);
          fHistDCA2pvZ[ctside][iTrack]->Fill(dcaZpv);
          fHistDCA2pv[ctside][iTrack]->Fill(dcapv);
        }
        /// compute the track2track distance used in the vertexer
        float pPM[3][3];
        for (int iPerm = 0; iPerm < 3; iPerm++) {
          fHypertritonVertexer.Find2ProngClosestPoint(lTrack[iPerm], lTrack[(iPerm + 1) % 3], lMagField, pPM[iPerm]);
          fHistTrackDistance[ctside][iPerm]->Fill(Point2PointDistance(pPM[iPerm], pPM[(iPerm + 1) % 3]) * 10.);
        }
        double vertexChi2NDF = lDecayVertex->GetChi2perNDF();
        fHistVertexChi2[ctside]->Fill(vertexChi2NDF);
      }
      else if(fAlg==1){
        cospa = recHyp.cosPA;
        /// compute the DCA of the 3 tracks from the primary and decay vertex
        //still not implemented
          fHistDCA2dvXY[ctside][0]->Fill(0);
          fHistDCA2dvXY[ctside][1]->Fill(0);
          fHistDCA2dvXY[ctside][2]->Fill(0);
          fHistDCA2dvZ[ctside][0]->Fill(0);
          fHistDCA2dvZ[ctside][1]->Fill(0);
          fHistDCA2dvZ[ctside][2]->Fill(0);
          fHistDCA2dv[ctside][0]->Fill(10*recHyp.dca_de_sv);
          fHistDCA2dv[ctside][1]->Fill(10*recHyp.dca_pr_sv);
          fHistDCA2dv[ctside][2]->Fill(10*recHyp.dca_pi_sv);
          fHistDCA2pv[ctside][0]->Fill(10*recHyp.dca_de);
          fHistDCA2pv[ctside][1]->Fill(10*recHyp.dca_pr);
          fHistDCA2pv[ctside][2]->Fill(10*recHyp.dca_pi);
          
          for (int iTrack = 0; iTrack < 3; iTrack++) {
            //float dcaXYpv = std::abs(0) * 10.;    // in mm
            //float dcaZpv  = std::abs(0) * 10.;    // in mm
            //float dcapv   = Norm(0, 0) * 10.;  // in mm

            fHistDCA2pvXY[ctside][iTrack]->Fill(0);
            fHistDCA2pvZ[ctside][iTrack]->Fill(0);
          }
          

          fHistTrackDistance[ctside][0]->Fill(10*recHyp.dca_de_pr);
          fHistTrackDistance[ctside][1]->Fill(10*recHyp.dca_de_pi);
          fHistTrackDistance[ctside][2]->Fill(10*recHyp.dca_pr_pi);

          fHistVertexChi2[ctside]->Fill(recHyp.chi2);
        
      }

      //if(lCollection.hypertriton.GetChi2()/lCollection.hypertriton.GetNDF())
      fHistInvMassPtSel[lCharge][2]->Fill(lHypMassRec,lHypPtRec);
      fHistRecPt->Fill(lHypPtGen);
      if(fNclones==0){
        fHistMassResPt[lCharge]->Fill(lHypMassRec-kHypMass,lHypPtGen);
        fHistMassResCt[lCharge]->Fill(lHypMassRec-kHypMass,lHypCtGen);
        fHistPtResPt[lCharge]->Fill(lHypPtRec-lHypPtGen,lHypPtGen);
        fHistCtResCt[lCharge]->Fill(lHypCtRec-lHypCtGen,lHypCtGen);
        fHistCtResCtTrueP[lCharge]->Fill(lHypCtRecTrueP-lHypCtGen,lHypCtGen);
        fHistPResP[lCharge]->Fill(lHypPRec-lHypPGen,lHypPGen);
        fHistXResX[lCharge]->Fill(lHypXRec-lHypXGen,lHypXGen);
        fHistYResY[lCharge]->Fill(lHypYRec-lHypYGen,lHypYGen);
        fHistZResZ[lCharge]->Fill(lHypZRec-lHypZGen,lHypZGen);
      }

      for(int iTrack = 0; iTrack < 3; iTrack++){
        fHistDaughterPt[iTrack][2]->Fill(lDaughterPt[iTrack]);
        fHistNclsITS[iTrack][2]->Fill(lDaughterNclsITS[iTrack]);
        fHistNclsTPC[iTrack][2]->Fill(lDaughterNclsTPC[iTrack]);
        fHistNSigmaTPC[iTrack][2]->Fill(lDaughterNsigmaTPC[iTrack]);
        fHistNSigmaTOF[iTrack][2]->Fill(lDaughterNsigmaTOF[iTrack]);
        fHistDaughterTPCchi2[iTrack][2]->Fill(lDaughterChi2TPC[iTrack]);
        fHistDaughterITSchi2[iTrack][2]->Fill(lDaughterChi2ITS[iTrack]);
      }
      //
      fHistInvMassPt[lCharge][2]->Fill(lHypMassRec,lHypPtRec);

    }
  
    fHistCosPAngle[ctside]->Fill(cospa);
    if(cospa<0.99)
      continue;
    
    /// Efficiency histograms
    int fOldClones = fNclones;
    for(int iCut=0; iCut<kNcuts; iCut++){
      for(int iVar=0; iVar<kNvariations; iVar++)
        if(AcceptCandidate(RCand,iVar,iCut)){
          fNclones=fOldClones;
          if(fNclones==0){
            fHistSingleRecVsCuts[iCut][lCharge]->Fill(lHypPtGen,lHypCtGen,iVar);
            fHistResolutionVsCuts[iCut][lCharge]->Fill(lHypPtRec-lHypPtGen,lHypCtRec-lHypCtGen,iVar);         
            fHistInvMassPt[lCharge][0]->Fill(lHypMassRec,lHypPtRec);
            fNclones++;
          }
          else{
            fHistClonesVsCuts[iCut][lCharge]->Fill(lHypPtGen,lHypCtGen,iVar);
            fHistInvMassPt[lCharge][2]->Fill(lHypMassRec,lHypPtRec);
          }
      }
    }  
  }
  return true;
}

//______________________________________________________________________________
void AliSelectorTaskHyperTriton3VtxPerf::SlaveTerminate(){
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

//______________________________________________________________________________
void AliSelectorTaskHyperTriton3VtxPerf::Terminate(){
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  fTimer->Stop();
  printf("time: %f",fTimer->RealTime());
  fTotTime->SetBinContent(1,fTimer->RealTime());
  fOperations->SetBinContent(1,fNcycles);
  fOperations->SetBinContent(2,fNrec);

  TFile output(Form("%s/%s", fOutputFilePath.Data(), fOutputFileName.Data()), "RECREATE");
  GetOutputList()->Write();
  output.Close();
}

//______________________________________________________________________________
bool AliSelectorTaskHyperTriton3VtxPerf::AcceptCandidate(RCandidate candidate,int iCut = 0, int iVar = 0){
  int CutsSet[3] = {0,0,0};
  CutsSet[iVar] = iCut;
  float dca[2];
  //has tof_de?
  const bool hasTOFout  = candidate.track[0]->GetStatus() & AliVTrack::kTOFout;
  const bool hasTOFtime = candidate.track[0]->GetStatus() & AliVTrack::kTIME;

  if(!(hasTOFout && hasTOFtime)) return false;
  if(TMath::Abs(candidate.NsigmaTPC[0]) > 4.) return false;
  for(int iTrack = 0; iTrack < 3; iTrack++){
    //cut on dca
    candidate.track[iTrack]->GetImpactParameters(dca[0], dca[1]);
    double dcaNorm = std::hypot(dca[0], dca[1]);
    if(dcaNorm<0.05) return false;
    if(!fESDtrackCuts->AcceptTrack(&*candidate.track[iTrack])) return false;
    //cut on NsigmaTPC
    if(candidate.NsigmaTPC[iTrack] > kCuts[0][CutsSet[0]][iTrack]) return false;
    //cut on NclusterTPC
    if(candidate.track[iTrack]->GetTPCNcls() < kCuts[1][CutsSet[1]][iTrack]) return false;
    //cut on NclusterITS
    //if(candidate.track[iTrack]->GetITSNcls() < kCuts[2][CutsSet[2]][iTrack]) return false;
  }
  return true;
}



bool AliSelectorTaskHyperTriton3VtxPerf::KFVertexer(AliESDtrack* kTrack [], RParticles &kResult, double pver[], float kfchi2[]){

  float kMasses[3]{kDeuMass,kPMass,kPiMass};
  double posmom[6],cov[21];
  KFParticle helper[3];

  for (int iT=0; iT < 3; iT++) {
    kTrack[iT]->GetXYZ(posmom);
    kTrack[iT]->GetPxPyPz(posmom+3);
    kTrack[iT]->GetCovarianceXYZPxPyPz(cov);
    helper[iT].Create(posmom,cov,kTrack[iT]->Charge(),kMasses[iT]);
    helper[iT].Chi2() = kTrack[iT]->GetTPCchi2();
    helper[iT].NDF() = kTrack[iT]->GetNumberOfTPCClusters() * 2;
  }
  
  KFParticle oneCandidate;
  oneCandidate.AddDaughter(helper[2]);
  if (kTrack[2] == kTrack[1] || kTrack[1]->Charge() * kTrack[0]->Charge() < 0)
    return false;
  KFParticle twoCandidate{oneCandidate};
  
  twoCandidate.AddDaughter(helper[1]);

  
  float chi2_deuprot = twoCandidate.GetChi2() / twoCandidate.GetNDF();
  kfchi2[0] = chi2_deuprot;
  if (chi2_deuprot > kMaxKFchi2[0])
    return false;

  if (kTrack[0] == kTrack[2] || kTrack[2]->Charge() * kTrack[1]->Charge() > 0) 
    return false;

  KFParticle hyperTriton{twoCandidate};
  hyperTriton.AddDaughter(helper[0]);
  float chi2_3prongs = hyperTriton.GetChi2() / hyperTriton.GetNDF();
  kfchi2[1] = chi2_3prongs;
  if (chi2_3prongs > kMaxKFchi2[1])
    return false;
  //a copy to compute the chi2_topology
  
  kResult.deuteron = helper[0];
  kResult.proton = helper[1];
  kResult.pion = helper[2];
  kResult.hypertriton = hyperTriton;

  KFParticle vertPart{hyperTriton};

  /*
  double pvPos[3], pvCov[6];
  fPrimaryVertex->GetXYZ(pvPos);
  fPrimaryVertex->GetCovarianceMatrix(pvCov);

  KFPVertex kfPVertex;
  kfPVertex.SetXYZ(pvPos[0],pvPos[1],pvPos[2]);
  kfPVertex.SetCovarianceMatrix(pvCov[0],pvCov[1],pvCov[2],pvCov[3],pvCov[4],pvCov[5]);
  kfPVertex.SetChi2(fPrimaryVertex->GetChi2());
  kfPVertex.SetNDF(fPrimaryVertex->GetNDF());
  kfPVertex.SetNContributors(fPrimaryVertex->GetNContributors());
  
  KFParticle prodVertex{kfPVertex};
  vertPart.SetProductionVertex(prodVertex);

  float chi2_topology = vertPart.GetChi2() / vertPart.GetNDF();
  kfchi2[2] = chi2_topology;
  if (chi2_topology > kMaxKFchi2[2])
    return false;
  */

  float dca_de_pr = helper[0].GetDistanceFromParticle(helper[1]);
  float dca_de_pi = helper[0].GetDistanceFromParticle(helper[2]);
  float dca_pr_pi = helper[1].GetDistanceFromParticle(helper[2]);

  if(dca_de_pr > 2. || dca_de_pi > 2. || dca_pr_pi > 2.)
    return false;
  
  return true;
}


bool AliSelectorTaskHyperTriton3VtxPerf::O2Vertexer(AliESDtrack* kTrack [],RHyperTritonO2 &recHyp, double pvPos[], float bz , float decVert[]){
  
  o2::vertexing::DCAFitter3 fVertexer;
  fVertexer.setBz(bz);
  //fVertexer.setMaxR(100.);
  fVertexer.setMaxChi2(10.);
  o2::track::TrackParCov* helper[3] = {nullptr};
  for (int iT=0; iT < 3; iT++)
    helper[iT] = (o2::track::TrackParCov*)((AliExternalTrackParam*)kTrack[iT]);
  
  //the next line is useless
  if(helper[0] == helper[1] || helper[0] == helper[2] || helper[1] == helper[2])
    return false;

  int nVert{0};
  try {
    nVert = fVertexer.process(*helper[0], *helper[1], *helper[2]);
  }
  catch (std::runtime_error& e) {
    std::cout << "MyException caught" << std::endl;
    std::cout << e.what() << std::endl;
    return false;
  }
  if (nVert) {
    auto vert = fVertexer.getPCACandidate();
    auto& deuTrack = fVertexer.getTrack(0);
    auto& prTrack = fVertexer.getTrack(1);
    auto& piTrack = fVertexer.getTrack(2);

    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> ldeu{deuTrack.Pt(), deuTrack.Eta(), deuTrack.Phi(), kDeuMass};
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> lpro{prTrack.Pt(), prTrack.Eta(), prTrack.Phi(), kPMass};
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> lpi{piTrack.Pt(), piTrack.Eta(), piTrack.Phi(), kPiMass};
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> hypertriton{ldeu + lpro + lpi};

    const float mass = hypertriton.mass();
    ROOT::Math::XYZVectorF decayl{(float)(vert[0] - pvPos[0]), (float)(vert[1] - pvPos[1]), (float)(vert[2] - pvPos[2])};
    const float totalMom = hypertriton.P();
    const float len = std::sqrt(decayl.Mag2());

    recHyp.ct = len * kHypMass / totalMom; 
    recHyp.candidates = nVert;
    recHyp.pt = hypertriton.pt();
    recHyp.pz = hypertriton.P();
    recHyp.m = mass;
    recHyp.cosPA = hypertriton.Vect().Dot(decayl) / (totalMom * len);

    auto& deuPos = fVertexer.getTrackPos(0); 
    auto& proPos = fVertexer.getTrackPos(1); 
    auto& piPos = fVertexer.getTrackPos(2); 

    recHyp.dca_de_pr = TMath::Hypot(TMath::Hypot(deuPos[0] - proPos[0], deuPos[1] - proPos[1]), deuPos[2] - proPos[2]);
    recHyp.dca_de_pi = TMath::Hypot(TMath::Hypot(deuPos[0] - piPos[0], deuPos[1] - piPos[1]), deuPos[2] - piPos[2]);
    recHyp.dca_pr_pi = TMath::Hypot(TMath::Hypot(proPos[0] - piPos[0], proPos[1] - piPos[1]), proPos[2] - piPos[2]);

    recHyp.dca_de_sv = TMath::Hypot(TMath::Hypot(deuPos[0] - vert[0], deuPos[1] - vert[1]), deuPos[2] - vert[2]);
    recHyp.dca_pr_sv = TMath::Hypot(TMath::Hypot(proPos[0] - vert[0], proPos[1] - vert[1]), proPos[2] - vert[2]);
    recHyp.dca_pi_sv = TMath::Hypot(TMath::Hypot(piPos[0] - vert[0], piPos[1] - vert[1]), piPos[2] - vert[2]);

    recHyp.chi2 = fVertexer.getChi2AtPCACandidate();
    decVert[0] = vert[0];
    decVert[1] = vert[1];
    decVert[2] = vert[2];
    return true;

  }
  else 
    return false;
}
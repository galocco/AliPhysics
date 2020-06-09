#include <TChain.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include <TCanvas.h>
#include "../AliAnalysisTaskHyperTriton3VtxPerf.h"
void KFAnalysis(const char* output_name="test.root",bool test=false)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  const float mgen=2.99131;

  TChain mcChain("Hyp3");
  mcChain.AddFile("../HyperTritonTree.root");

  TTreeReader fReader(&mcChain);
  TTreeReaderArray<RCandidate> RHyp = {fReader, "RCandidate"};
  TTreeReaderArray<SHyperTriton3> SHyperVec = {fReader, "SHyperTriton"};
  TTreeReaderArray<int> GenMapVec = {fReader, "SGenRecMap"};
  TTreeReaderValue<REvent> RColl = {fReader, "RCollision"};


  const char lAM[3]{"AM"};
  TH1D* fHistCosPA = new TH1D("fHistCosPA","",1000,0,1);
  TH1D* fHistGenPt= new TH1D("fHistGenPt",";p_{T} [GeV/c];",10,0,10);

  float cut_dca[]={0.05,0.4,0.6,0.8,1.,1.2,10.};
  float cut_chi2[]={1.,1.5,2.,2.5,3.,3.5,50.};
  //inizializzo gli istogrammi
  //leggo i file e riempio gli istogrammi
  int ir=0;
  float pt,ptgen,ctgen;
  int clone=0,charge;
  int it=0;
  while (fReader.Next())
  {
    for (auto &SHyper : SHyperVec)
    {
      fHistGenPt->Fill(SHyper.pt);
    }
    
    vector<int> OldMap;
    for (auto &RCand : RHyp)
    {
      //clone=0
      ir++;
      for(auto rec_index : OldMap){
        if(GenMapVec[ir] == rec_index)
          clone++;
      }
      OldMap.push_back(GenMapVec[ir]);
      charge = (SHyperVec[GenMapVec[ir]].positive) ? 1 : 0;
      
    }

  }
  std::cout<<clone<<" "<<ir<<"\n";
  TFile rfile(output_name,"RECREATE");
  fHistGenPt->Write();   
}
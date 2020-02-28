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
//TODO: the number of cuts is everywhere hardcoded

struct SHyperTriton3KF {
  float px = -999.f;
  float py = -999.f;
  float pz = -999.f;
  float l = -1.f;
  float t = -1.f;
  bool positive = false;
};


struct RHyperTriton3KF {
  float px = -999.f;
  float py = -999.f;
  float pz = -999.f;
  float l = -1.f;
  float r = -1.f;
  float chi2_deuprot = -1.f;
  float chi2_3prongs = -1.f;
  float chi2_topology = -1.f;
  float cosPA = -1.f;
  float m = -1;
  Double32_t dca_de = 2.0; //[0.0,2.0,8]
  Double32_t dca_pr = 2.0; //[0.0,2.0,8]
  Double32_t dca_pi = 2.0; //[0.0,2.0,8]
  Double32_t tpcNsig_de = -4.0; //[-4.0,4.0,8]
  Double32_t tpcNsig_pr = -4.0; //[-4.0,4.0,8]
  Double32_t tpcNsig_pi = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_de = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_pr = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_pi = -4.0; //[-4.0,4.0,8]
  Double32_t dca_de_pr = -4.0; //[0.0,2.0,8]
  Double32_t dca_de_pi = -4.0; //[0.0,2.0,8]
  Double32_t dca_pr_pi = -4.0; //[0.0,2.0,8]
  bool hasTOF_de;
  bool hasTOF_pr;
  bool hasTOF_pi;
  UChar_t tpcClus_de = 0u;
  UChar_t tpcClus_pr = 0u;
  UChar_t tpcClus_pi = 0u;
};

void SetEfficiencyErrors(TH1D* , TH1D* );
void SetRateErrors(TH1D*, TH1D*);
void Rainbow_Plot_Eff(TH2D* [][3][7][2],TH2D* [2],const char* [],float [],TDirectory*);
void MassResolutionPlot(TH3D* [][3][7][2],const char* [],float [],TDirectory*);

void KFAnalysis(const char* output_name="test.root",bool test=false)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  const float mgen=2.99131;
  TChain mcChain("Hyp3KF");
  mcChain.AddFile("HyperTritonTree3KF.root");

  TTreeReader fReader(&mcChain);
  TTreeReaderArray<RHyperTriton3KF> RHyperVec = {fReader, "RHyperTriton"};
  TTreeReaderArray<SHyperTriton3KF> SHyperVec = {fReader, "SHyperTriton"};
  TTreeReaderArray<int> GenMapVec = {fReader, "SGenRecMap"};
  //TTreeReaderValue<RCollision> RColl = {fReader, "RCollision"};


  const int n_cut_dca = 7;
  const int n_cut_chi2 = 7;
  const char *lCutDCA[3]{"DCA_pd", "DCA_ppi","DCA_dpi"};
  const char *lCutChi2[3]{"Chi2_deuprot","Chi2_3prongs","Chi2_topology"};
  const char lAM[3]{"AM"};


  TH1D* fHistChi[3];
  TH1D* fHistDCA[3];
  TH1D* fHistCosPA = new TH1D("fHistCosPA","",1000,0,1);
  for(int iCut=0; iCut<3; iCut++)
  {
    fHistChi[iCut] = new TH1D(Form("fHist%s",lCutChi2[iCut]),";#Chi^2",200,0,30);
    fHistDCA[iCut] = new TH1D(Form("fHist%s",lCutDCA[iCut]),";DCA",200,0,2.2);
  }

  TH2D* fHistGen[2];

  TH1D* fHistGenPt= new TH1D("fHistGenPt",";p_{T} [GeV/c];",10,0,10);
  TH1D* fHistGenPX= new TH1D("fHistGenPX",";p_{X} [GeV/c];",400,-10,10);
  TH1D* fHistGenPY= new TH1D("fHistGenPY",";p_{Y} [GeV/c];",400,-10,10);
  TH1D* fHistGenPZ= new TH1D("fHistGenPZ",";p_{Z} [GeV/c];",400,-10,10);
  //[rec/clone][#kinds of cuts][#number of cuts][#matter]
  TH2D* fHistEffDCA[2][3][n_cut_dca][2];
  TH2D* fHistEffChi2[2][3][n_cut_chi2][2];
  TH3D* fHistMassDCA[2][3][n_cut_dca][2];
  TH3D* fHistMassChi2[2][3][n_cut_chi2][2];
  TH2D* fHistResChi2[2][3][n_cut_chi2];
  TH2D* fHistResDCA[2][3][n_cut_dca];

  float cut_dca[]={0.05,0.4,0.6,0.8,1.,1.2,10.};
  float cut_chi2[]={1.,1.5,2.,2.5,3.,3.5,50.};
  //inizializzo gli istogrammi
  for(int iMat=0; iMat<2; iMat++){
    fHistGen[iMat] = new TH2D(Form("fHistGen_%c",lAM[iMat]),";p_{T} [GeV/c];ct [cm]",10,0,10,10,0,50);
    for(int iCut=0; iCut<n_cut_dca; iCut++)
    {
      for(int iPar=0; iPar<3;iPar++){
        fHistEffDCA[0][iPar][iCut][iMat] = new TH2D(Form("fHistEff_%.1f_%s_%c",cut_dca[iCut],lCutDCA[iPar],lAM[iMat]),";p_{T} [GeV/c];ct [cm]",10,0,10,10,0,50);
        fHistEffDCA[1][iPar][iCut][iMat] = new TH2D(Form("fHistEff_%.1f_%s_clone_%c",cut_dca[iCut],lCutDCA[iPar],lAM[iMat]),";p_{T} [GeV/c];ct [cm]",10,0,10,10,0,50);
      
        fHistMassDCA[0][iPar][iCut][iMat] = new TH3D(Form("fHistMass_%.1f_%s_%c",cut_dca[iCut],lCutDCA[iPar],lAM[iMat]),";p_{T} [GeV/c];ct [cm];#Deltam [GeV/c^2]",10,0,10,10,0,50,100,-0.05,0.05);
        fHistMassDCA[1][iPar][iCut][iMat] = new TH3D(Form("fHistMass_%.1f_%s_clone_%c",cut_dca[iCut],lCutDCA[iPar],lAM[iMat]),";p_{T} [GeV/c];ct [cm];#Deltam [GeV/c^2]",10,0,10,10,0,50,100,-0.05,0.05);
      
        //fHistResDCA[iPar][iCut][iMat] = new TH2D(Form("fHistResDCA_%.1f_%s_%c",cut_dca[iCut],lCutDCA[iPar],lAM[iMat]),";#Deltap_{T} [GeV/c];#Deltact [cm]",100,-2,2,100,-10,10);          
      }
    }
    for(int iCut=0; iCut<n_cut_chi2; iCut++)
    {
      for(int iPar=0; iPar<3;iPar++){
        fHistEffChi2[0][iPar][iCut][iMat] = new TH2D(Form("fHistEff_%.1f_%s_%c",cut_chi2[iCut],lCutChi2[iPar],lAM[iMat]),";p_{T} [GeV/c];ct [cm]",10,0,10,10,0,50);
        fHistEffChi2[1][iPar][iCut][iMat] = new TH2D(Form("fHistEff_%.1f_%s_clone_%c",cut_chi2[iCut],lCutChi2[iPar],lAM[iMat]),";p_{T} [GeV/c];ct [cm]",10,0,10,10,0,50);
        fHistMassChi2[0][iPar][iCut][iMat] = new TH3D(Form("fHistMass_%.1f_%s_%c",cut_chi2[iCut],lCutChi2[iPar],lAM[iMat]),";p_{T} [GeV/c];ct [cm];#Deltam [GeV/c^2]",10,0,10,10,0,50,100,-0.05,0.05);
        fHistMassChi2[1][iPar][iCut][iMat] = new TH3D(Form("fHistMass_%.1f_%s_clone_%c",cut_chi2[iCut],lCutChi2[iPar],lAM[iMat]),";p_{T} [GeV/c];ct [cm];#Deltam [GeV/c^2]",10,0,10,10,0,50,100,-0.05,0.05);
        //fHistResChi2[iPar][iCut][iMat] = new TH2D(Form("fHistResChi2_%.1f_%s_%c",cut_dca[iCut],lCutDCA[iPar],lAM[iMat]),";#Deltap_{T} [GeV/c];#Deltact [cm]",100,-2,2,100,-10,10);          
      }
    }
  }

  //leggo i file e riempio gli istogrammi
  int ir=0;
  float pt,ptgen,ctgen;
  int clone,charge;
  int it=0;
  while (fReader.Next())
  {
    if(test){
      it++;
      if(it==10000)
        break;
    }
    for (auto &SHyper : SHyperVec)
    {
      fHistGen[SHyper.positive]->Fill(TMath::Sqrt(SHyper.px*SHyper.px+SHyper.py*SHyper.py),SHyper.l);
      fHistGenPt->Fill(TMath::Sqrt(SHyper.px*SHyper.px+SHyper.py*SHyper.py));
      fHistGenPX->Fill(SHyper.px);
      fHistGenPY->Fill(SHyper.py);
      fHistGenPZ->Fill(SHyper.pz);
    }
    
      
    ir=0;
    vector<int> OldMap;
    for (auto &RHyper : RHyperVec)
    {
      clone=0;
      for(auto t:OldMap){
        if(GenMapVec[ir]==t)
          clone=1;
      }
      OldMap.push_back(GenMapVec[ir]);
      charge = (SHyperVec[GenMapVec[ir]].positive) ? 1 : 0;
      ptgen = TMath::Sqrt(SHyperVec[GenMapVec[ir]].px*SHyperVec[GenMapVec[ir]].px+SHyperVec[GenMapVec[ir]].py*SHyperVec[GenMapVec[ir]].py);
      ctgen = SHyperVec[GenMapVec[ir]].l;

      ir++;
      double dca[3]={RHyper.dca_de_pr,RHyper.dca_de_pi,RHyper.dca_pr_pi};
      double chi2[3]={RHyper.chi2_deuprot,RHyper.chi2_3prongs,RHyper.chi2_topology}; 
      bool no=false;
      for(int iCut=0; iCut<3; iCut++)
      {
        fHistChi[iCut]->Fill(chi2[iCut]);
        fHistDCA[iCut]->Fill(dca[iCut]);
      }
      fHistCosPA->Fill(RHyper.cosPA);
      //messo come test per il selector
      bool t[5];
      t[0]=true;//RHyper.dca_de>0.05 && RHyper.dca_pr>0.005 && RHyper.dca_pi>0.005;
      t[1]=TMath::Abs(RHyper.tpcNsig_de)<3. && TMath::Abs(RHyper.tpcNsig_pr)<3. && TMath::Abs(RHyper.tpcNsig_pi)<3.;
      t[2]=RHyper.cosPA>0.99;
      t[3]=RHyper.hasTOF_de && TMath::Abs(RHyper.tofNsig_de)<4;
      //t[4]=TMath::Abs(RHyper.tpcClus_de)<3.99 && TMath::Abs(RHyper.tpcClus_pi)<3.99 && TMath::Abs(RHyper.tpcClus_pr)<3.99;
      //t[4]=TMath::Abs(RHyper.tpcClus_de)>70 && TMath::Abs(RHyper.tpcClus_pi)>70 && TMath::Abs(RHyper.tpcClus_pr)>70;

      bool tot = true;
      for(int i=0;i<4;i++){
        tot= tot && t[i];
        //cout<<i<<" "<<t[i]<<endl;
      }
      if(tot){
      pt = TMath::Sqrt(RHyper.px*RHyper.px+RHyper.py*RHyper.py);
      for(int cut=0; cut<n_cut_dca; cut++)
      {
        for(int iPar=0; iPar<3;iPar++)
        {
          if(dca[iPar]<cut_dca[cut]){
            fHistEffDCA[clone][iPar][cut][charge]->Fill(ptgen,ctgen);
            fHistMassDCA[clone][iPar][cut][charge]->Fill(ptgen,ctgen,RHyper.m-mgen);
            //if(clone==0)
              //fHistResDCA[iPar][cut][charge]->Fill(pt-ptgen,RHyper.l-ctgen);
          }
        }
      }
      for(int cut=0; cut<n_cut_chi2; cut++)
      {
        for(int iPar=0; iPar<3;iPar++)
        {
          if(chi2[iPar]<cut_chi2[cut]){
            fHistEffChi2[clone][iPar][cut][charge]->Fill(ptgen,ctgen);
            fHistMassChi2[clone][iPar][cut][charge]->Fill(ptgen,ctgen,RHyper.m-mgen);
            //if(clone==0)
            //  fHistResChi2[iPar][cut][charge]->Fill(pt-ptgen,RHyper.l-ctgen);
          }
        }
      }
      }
    }

  }

  TFile rfile(output_name,"RECREATE");
  const int ndir=5;
  const char *lDir[ndir]={"utils","efficiency_dca","efficiency_chi2","resolution_dca","resolution_chi2"};
  TDirectory* subdir[ndir];
  for(int i=0; i<ndir; i++)
    subdir[i] = rfile.mkdir(lDir[i]);
  subdir[0]->cd();
  //subdir[0]->cd();
  for(int iCut=0; iCut<3; iCut++)
  {
    fHistChi[iCut]->Write();
    fHistDCA[iCut]->Write();
  }
  fHistCosPA->Write();
  //salvare istogrammi di massa
  for(int iMat=0; iMat<2; iMat++){
    for(int iCut=0; iCut<n_cut_dca; iCut++)
    {
      for(int iPar=0; iPar<3;iPar++){
        //fHistResDCA[iPar][iCut][iMat]->Write();
      }
    }
    for(int iCut=0; iCut<n_cut_chi2; iCut++)
    {
      for(int iPar=0; iPar<3;iPar++){
        //fHistResChi2[iPar][iCut][iMat]->Write();
      }
    }
  }

  fHistGenPt->Write();
  fHistGenPX->Write();
  fHistGenPY->Write();
  fHistGenPZ->Write();
  Rainbow_Plot_Eff(fHistEffChi2,fHistGen,lCutChi2,cut_chi2,subdir[2]);
  Rainbow_Plot_Eff(fHistEffDCA,fHistGen,lCutDCA,cut_dca,subdir[1]);
  
  MassResolutionPlot(fHistMassChi2,lCutChi2,cut_chi2,subdir[4]);
  MassResolutionPlot(fHistMassDCA,lCutDCA,cut_dca,subdir[3]);

    
}

void Rainbow_Plot_Eff(TH2D* fHistRec[][3][7][2],TH2D* fHistGen[2],const char* lVar[],float cut[],TDirectory* rfile)
{

  const char *lClone[2]{"Eff", "Clone"};
  const char *lYLabel[2]{"efficiency", "#clone/#gen"};
  const char lAM[3]{"AM"};
  const char *lProj[2]{"Pt", "Ct"};
  const float lYRange[2][2]={{0,1.},{0.,0.015}};
  //rfile->ReOpen("UPDATE");
  
  //proiezioni e divisioni
  TH1D* fHistEff;
  TH1D* fHistGenProj;
  TCanvas Rainbow("","");

  TCanvas MultiRainbow("","",1200,400);
  MultiRainbow.Divide(2,1);

  TLegend *leg=new TLegend(0.2,0.2,0.5,0.5);
  //faccio le proiezioni e calcolo le efficienze
  
  for(int iMat=0; iMat<2;iMat++)
  {
    for(int iPar=0; iPar<3;iPar++)
    {
      //pt or ct
      for(int iProj=0; iProj<2; iProj++)
      {
        //eff or rate
        for(int iHist=0; iHist<2;iHist++)
        {
          //cuts
          for(int iCut=0; iCut<7; iCut++)
          {
            if(iProj==0){
              fHistEff = (TH1D*) fHistRec[iHist][iPar][iCut][iMat]->ProjectionX(Form("fHist%s_%s_%s_%.1f_%c",lClone[iHist],lProj[iProj],lVar[iPar],cut[iCut],lAM[iMat]),1,fHistRec[iHist][iPar][iCut][iMat]->GetNbinsX());  
              fHistGenProj = (TH1D*) fHistGen[iMat]->ProjectionX("fHistGenProj",1,fHistGen[iMat]->GetNbinsX());
            }
            else{
              fHistEff = (TH1D*) fHistRec[iHist][iPar][iCut][iMat]->ProjectionY(Form("fHist%s_%s_%s_%.1f_%c",lClone[iHist],lProj[iProj],lVar[iPar],cut[iCut],lAM[iMat]),1,fHistRec[iHist][iPar][iCut][iMat]->GetNbinsY());  
              fHistGenProj = (TH1D*) fHistGen[iMat]->ProjectionY("fHistGenProj",1,fHistGen[iMat]->GetNbinsY());
            }

            fHistEff->SetTitle(Form("%s < %.1f",lVar[iPar],cut[iCut]));
            fHistEff->SetMarkerStyle(8);
            fHistEff->Divide(fHistGenProj);

            if(iHist==0)
              SetEfficiencyErrors(fHistEff,fHistGenProj);
            else
              SetRateErrors(fHistEff,fHistGenProj);

            fHistEff->GetYaxis()->SetRangeUser(lYRange[iHist][0],lYRange[iHist][1]);
            fHistEff->GetYaxis()->SetTitle(lYLabel[iHist]);

            rfile->cd();
            fHistEff->Write();
            Rainbow.cd();
            if(iCut==0)
              fHistEff->Draw("PMC PLC");
            else
              fHistEff->Draw("SAME PMC PLC");
            //draw in the divided canvas
            
            MultiRainbow.cd(iHist+1);
            if(iCut==0)
              fHistEff->Draw("PMC PLC");
            else
              fHistEff->Draw("SAME PMC PLC");
              
            if(iHist==0)
              leg->AddEntry(fHistEff,fHistEff->GetTitle(),"lp");

          }
          //dca/chi2_taglio
          Rainbow.SetName(Form("RainbowPlot%s_%s_%s_%c",lClone[iHist],lProj[iProj],lVar[iPar],lAM[iMat])); 
          Rainbow.BuildLegend();
          rfile->cd();
          Rainbow.Write();
        
        }
        leg->SetTextFont(132);
        MultiRainbow.SetName(Form("MultiRainbowPlot%s_%s_%c",lProj[iProj],lVar[iPar],lAM[iMat])); 
        MultiRainbow.cd(1);
        leg->Draw();
        MultiRainbow.Write();
        leg->Clear();
      }
    }
  }
}

void MassResolutionPlot(TH3D* fHistRec[][3][7][2],const char* lVar[],float cut[],TDirectory* rfile)
{

  const char *lClone[2]{"Rec", "Clone"};
  const char lAM[3]{"AM"};
  const char *lProj[2]{"Pt", "Ct"};
  const char *lProjSet[2]{"zx", "zy"};
  //rfile->ReOpen("UPDATE");
  
  //proiezioni e divisioni
  TH1D* fHistRes;
  TH2D* fHistProto;
  TH1D* fAverage;
  TH1D* fStdDev;
  TCanvas Rainbow("","");
  TCanvas RainbowAve("","");
  TCanvas RainbowStd("","");

  TCanvas MultiRainbow("","",1200,400);
  MultiRainbow.Divide(2,1);

  TLegend *leg=new TLegend(0.2,0.2,0.5,0.5);
  //faccio le proiezioni e calcolo le efficienze
  
  for(int iMat=0; iMat<2;iMat++)
  {
    for(int iPar=0; iPar<3;iPar++)
    {
      for(int iHist=0; iHist<2; iHist++)
      {
        //pt or ct
        for(int iProj=0; iProj<2; iProj++)
        {
          //cuts
          for(int iCut=0; iCut<7; iCut++)
          {
            fHistProto = (TH2D*) fHistRec[iHist][iPar][iCut][iMat]->Project3D(lProjSet[iProj]);

            if(iProj==0){
              fAverage = new TH1D(Form("fAverage%s_%s_%s_%.1f_%c",lClone[iHist],lProj[iProj],lVar[iPar],cut[iCut],lAM[iMat]),";pT [GeV/c];#mu [GeV/c^2]",fHistProto->GetNbinsX(),0,fHistProto->GetXaxis()->GetBinLowEdge(fHistProto->GetNbinsX()+1));
              fStdDev = new TH1D(Form("fStdDev%s_%s_%s_%.1f_%c",lClone[iHist],lProj[iProj],lVar[iPar],cut[iCut],lAM[iMat]),";pT [GeV/c];#sigma [GeV/c^2]",fHistProto->GetNbinsX(),0,fHistProto->GetXaxis()->GetBinLowEdge(fHistProto->GetNbinsX()+1));
            }
            else{
              fAverage = new TH1D(Form("fAverage%s_%s_%s_%.1f_%c",lClone[iHist],lProj[iProj],lVar[iPar],cut[iCut],lAM[iMat]),";ct [cm];#mu [GeV/c^2]",fHistProto->GetNbinsX(),0,fHistProto->GetXaxis()->GetBinLowEdge(fHistProto->GetNbinsX()+1));
              fStdDev = new TH1D(Form("fStdDev%s_%s_%s_%.1f_%c",lClone[iHist],lProj[iProj],lVar[iPar],cut[iCut],lAM[iMat]),";ct [cm];#sigma [GeV/c^2]",fHistProto->GetNbinsX(),0,fHistProto->GetXaxis()->GetBinLowEdge(fHistProto->GetNbinsX()+1));
            }
            
            for(int iBin=1;iBin<=fHistProto->GetNbinsX();iBin++){
              fHistRes  = (TH1D*) fHistProto->ProjectionY("fHistGenProj",iBin,iBin);
              fAverage->SetBinContent(iBin,fHistRes->GetMean());
              fStdDev->SetBinContent(iBin,fHistRes->GetStdDev());
              fAverage->SetBinError(iBin,fHistRes->GetMeanError());
              fStdDev->SetBinError(iBin,fHistRes->GetStdDevError());
            }

            fAverage->SetTitle(Form("%s < %.1f",lVar[iPar],cut[iCut]));
            fStdDev->SetTitle(Form("%s < %.1f",lVar[iPar],cut[iCut]));

            if(iCut==6){
              float deltaAv = fAverage->GetMaximum()-fAverage->GetMinimum();
              float deltaStd = fStdDev->GetMaximum()-fStdDev->GetMinimum();

              fAverage->GetYaxis()->SetRangeUser(fAverage->GetMinimum()-deltaAv,fAverage->GetMaximum()+deltaAv);
              fStdDev->GetYaxis()->SetRangeUser(fStdDev->GetMinimum()-deltaStd,fStdDev->GetMaximum()+deltaStd);
            }

            fAverage->SetMarkerStyle(8);
            fStdDev->SetMarkerStyle(8);
            rfile->cd();
            fAverage->Write();
            fStdDev->Write();

            RainbowAve.cd();
            if(iCut==0)
              fAverage->Draw("PMC PLC");
            else
              fAverage->Draw("SAME PMC PLC");
            RainbowStd.cd();
            if(iCut==0)
              fStdDev->Draw("PMC PLC");
            else
              fStdDev->Draw("SAME PMC PLC");
            //draw in the divided canvas
            MultiRainbow.cd(1);
            if(iCut==0)
              fAverage->Draw("PMC PLC");
            else
              fAverage->Draw("SAME PMC PLC");
            MultiRainbow.cd(2);
            if(iCut==0)
              fStdDev->Draw("PMC PLC");
            else
              fStdDev->Draw("SAME PMC PLC");
              
            leg->AddEntry(fStdDev,fStdDev->GetTitle(),"lp");
          }

          RainbowAve.SetName(Form("RainbowPlotAverage%s_%s_%s_%c",lClone[iHist],lProj[iProj],lVar[iPar],lAM[iMat])); 
          RainbowAve.BuildLegend();
          RainbowStd.SetName(Form("RainbowPlotStd%s_%s_%s_%c",lClone[iHist],lProj[iProj],lVar[iPar],lAM[iMat])); 
          RainbowStd.BuildLegend();
          rfile->cd();
          RainbowAve.Write();
          RainbowStd.Write();
          leg->SetTextFont(132);
          MultiRainbow.SetName(Form("MultiRainbowPlotMass%s_%s_%s_%c",lClone[iHist],lProj[iProj],lVar[iPar],lAM[iMat])); 
          MultiRainbow.cd(1);
          leg->Draw();
          MultiRainbow.Write();
          leg->Clear();
        }
      }
    }
  }
}

void SetEfficiencyErrors(TH1D* HistEff, TH1D* HistGen){
  for(int iBin=1; iBin<=HistEff->GetNbinsX(); iBin++){
    double gen = HistGen->GetBinContent(iBin);
    double eff = HistEff->GetBinContent(iBin);
    if(gen!=0)
      HistEff->SetBinError(iBin,TMath::Sqrt(eff*(1-eff)/gen));
    else{
      HistEff->SetBinError(iBin,0);
      HistEff->SetBinContent(iBin,1);
    }
  }
}

void SetRateErrors(TH1D* HistRate, TH1D* HistGen){
  for(int iBin=1; iBin<=HistRate->GetNbinsX(); iBin++){
    double gen = HistGen->GetBinContent(iBin);
    double rate = HistRate->GetBinContent(iBin);
    if(gen!=0)
      HistRate->SetBinError(iBin,TMath::Sqrt(rate/gen));
    else{
      HistRate->SetBinError(iBin,0);
      HistRate->SetBinContent(iBin,-1);
    }
  }
}
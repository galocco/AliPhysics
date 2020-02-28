#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <string.h>
#include "Hyp3FindConfig.h"
#include "utils.h"
void ordered_indexes(int ,int []);

void selector_efficiencies(char * input_name = "selector_results.root",char * output_name = "efficiency.root",char* pdf_file = "efficiency_plots.pdf",char * folder_name = "plot_folder")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lProj[2] = {"pT","ct"};
  const char* lProjSet[2] = {"zx","zy"};
  const char* lVarName[3] = {"NsigmaTPC","NclusterTPC","NclusterITS"};
  const char* lMeasure[3] = {"Efficiency","FakeRate","CloneRate"};
  const char* lProjMeas[3] = {"Rec","Fake","Clones"};

  TFile input_file(input_name);
  TFile gen_file("Output.root");
  TH3D* fHistGen3D[2];
  fHistGen3D[0] = (TH3D*) gen_file.Get("Hyp3FindTask_summary/fHistGeneratedPtVsYVsCentralityHypTrit");
  fHistGen3D[1] = (TH3D*) gen_file.Get("Hyp3FindTask_summary/fHistGeneratedPtVsYVsCentralityAntiHypTrit");
  int NHyp3[2];
  NHyp3[0] = 4339520;//fHistGen3D[0]->GetEntries();
  NHyp3[1] = 4339520;//fHistGen3D[1]->GetEntries();
  const float cthyp=7.25;
  TFile output_file(Form("%s/%s",folder_name,output_name),"RECREATE");
  const int ndir=5;
  TDirectory* subdir;
  for(int iVar=0; iVar<3; iVar++){
    for(int iMeas=0; iMeas<3; iMeas++){
      for(int iProj=0; iProj<2; iProj++){
        subdir = output_file.mkdir(Form("%s/%s/%s",lVarName[iVar],lMeasure[iMeas],lProj[iProj]));
      }
    }
  }
  
  TH3D* fHistRecVsCuts = nullptr;
  TH2D* fHistRecProj = nullptr;
  TH1D* fHistRec[3][kNvariations] = {{nullptr}};
  TH2D* fHistGenTot = nullptr;
  TH1D* fHistGen[2] = {nullptr}; //for pt and ct
  fHistGen[0] = new TH1D("fHistGetPt","",10,0,10);
  fHistGen[1] = new TH1D("fHistGetCt","",10,0,50); 
  TCanvas cv("","",800,450);
  cv.Print(Form("%s/%s[",folder_name,pdf_file));
  //array for the indexes of the cuts, with the cuts ordered from the tighter to the looser
  int ordered_cuts[kNvariations];
  //matter or antimatter
  for(int iMat=0; iMat<2; iMat++){
    // histograms of the genereted hypetritons
    //fHistGenTot = (TH2D*) input_file.Get(Form("fHistGen_%c",lAM[iMat]));
    for(int iBin=1; iBin<=10; iBin++){
      fHistGen[0]->SetBinContent(iBin,NHyp3[iMat]/10.);
      fHistGen[1]->SetBinContent(iBin,NHyp3[iMat]*(TMath::Exp(-5.*(iBin-1)/cthyp)-TMath::Exp(-5.*iBin/cthyp)*cthyp/5)); 
    }

    //projection on pt or ct
    for(int iProj=0; iProj<2; iProj++){
      //variable of the cut
      for(int iVar=0; iVar<kNcuts; iVar++){
        //meaurements (efficiency,fake-rate,clone-rate)
        for(int iMeas=0; iMeas<3; iMeas++){
          //get the TH3(pt,ct,cut)
          fHistRecVsCuts = (TH3D*) input_file.Get(Form("fHist%s_%s_%c",lProjMeas[iMeas],lVarName[iVar],lAM[iMat]));
          fHistRecProj = (TH2D*) fHistRecVsCuts->Project3D(lProjSet[iProj]);
          //ordering of the cuts
          ordered_indexes(iVar,ordered_cuts);
          //different cuts
          for(auto iCut : ordered_cuts){
            fHistRec[iMeas][iCut] = (TH1D*)fHistRecProj->ProjectionX(Form("fHist%s_%s_%s_cut_%.1f_%c",lProjMeas[iMeas],lVarName[iVar],lProj[iProj],kCuts[iVar][iCut][0],lAM[iMat]),iCut+1,iCut+1);
            if(lProjMeas[iMeas]=="Rec"){
              fHistRec[iMeas][iCut]->Divide(fHistGen[iProj]);
              SetEfficiencyErrors(fHistRec[iMeas][iCut],fHistGen[iProj]);
              fHistRec[iMeas][iCut]->GetYaxis()->SetTitle("efficiency");
            } else {
              fHistRec[iMeas][iCut]->Divide(fHistGen[iProj]);
              SetRateErrors(fHistRec[iMeas][iCut],fHistGen[iProj]); 
              fHistRec[iMeas][iCut]->GetYaxis()->SetTitle(Form("#%s/#Gen",lProjMeas[iMeas]));
            }
            //function to set some graphic features
            set_style(fHistRec[iMeas][iCut],iCut);
            fHistRec[iMeas][iCut]->SetTitle(Form("cut on %s : %.1f",lVarName[iVar],kCuts[iVar][iCut][0]));
            output_file.cd(Form("%s/%s/%s",lVarName[iVar],lMeasure[iMeas],lProj[iProj]));
            fHistRec[iMeas][iCut]->Write();
          }
          //set the range of the histograms
          set_range(fHistRec[iMeas]);
          //plot of the histograms in a single canvas
          rainbow_plot(fHistRec[iMeas],cv,Form("RainbowPlot_%s_%s_%s_%c",lVarName[iVar],lMeasure[iMeas],lProj[iProj],lAM[iMat]),true);          
          output_file.cd();
          cv.Write();          
          cv.Print(Form("%s/%s",folder_name,pdf_file)); 
          cv.Clear();
        }
        //plot of the three measurements in a single canvas
        multirainbow_plot(fHistRec,cv,Form("MultiRainbowPlot_%s_%s_%c",lVarName[iVar],lProj[iProj],lAM[iMat]),true);
        cv.Print(Form("%s/%s",folder_name,pdf_file));
        cv.Clear();
      }
    }
  }
  cv.Print(Form("%s/%s]",folder_name,pdf_file));
}


///
void compare_efficiencies(char* Std_name="efficiencyStd.root", char* KF_name="efficiencyKF.root", char* O2_name="efficiencyO2.root", char* output_name="efficiency_comparison.root", char* pdf_file="efficiency_comparison.pdf", char* folder_name="nome")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lProj[2] = {"pT","ct"};
  const char* lMeasure[3] = {"Efficiency","FakeRate","CloneRate"};
  const char* lProjMeas[3] = {"Rec","Fake","Clones"};
  const char* lTitle[3] = {"Standard vertexer","O^{2} vertexer","Kalman Filter"};

  TFile* input_file[3];
  input_file[0] = new TFile(Std_name);
  input_file[1] = new TFile(O2_name);
  input_file[2] = new TFile(KF_name);
  TFile output_file(output_name,"RECREATE");
  TH1D* fHist[3][3];//[measurement(eff,fake,clones)][vertexer]

  TCanvas cv("","",800,450);
  cv.Print(Form("%s[",pdf_file));
  //antimatter or matter
  for(int iMat=0; iMat<2; iMat++){
    //pt or ct
    for(int iProj=0; iProj<2; iProj++){
      //compare efficiency, fake rate and clones rate
      for(int iMeas=0; iMeas<3; iMeas++){
        //loop on the three vertexer results
        for(int iFile=0; iFile<3; iFile++){
          fHist[iMeas][iFile] = (TH1D*) input_file[iFile]->Get(Form("NsigmaTPC/%s/%s/fHist%s_NsigmaTPC_%s_cut_%.1f_%c",lMeasure[iMeas],lProj[iProj],lProjMeas[iMeas],lProj[iProj],kCuts[0][0][0],lAM[iMat]));
          set_style(fHist[iMeas][iFile],iFile);
          fHist[iMeas][iFile]->SetTitle(lTitle[iFile]);
        }
        //set y axis range
        set_range(fHist[iMeas]);
        //plots the all the iMeas-measurements in the same canvas
        rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%c",lMeasure[iMeas],lProj[iProj],lAM[iMat]),true);
        
        output_file.cd();
        cv.Write();          
        cv.Print(Form("%s",pdf_file)); 
        cv.Clear();
      }
      //plots the all the measurements in the same canvas
      multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%c",lProj[iProj],lAM[iMat]),true);
      cv.Write();          
      cv.Print(Form("%s",pdf_file)); 
      cv.Clear();    
    }
  }
  cv.Print(Form("%s]",pdf_file));
}

void compare_resolutions(char* Std_name="MassResolutionStd.root", char* KF_name="MassResolutionKF.root", char* O2_name="MassResolutionO2.root", char* output_name="resolution_comparison.root", char* pdf_file="resolution_comparison.pdf", char* folder_name="resolutionAll")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[3] = {"Standard vertexer","O^{2} vertexer","Kalman Filter #chi^{2}/NDF<50"};
  //const char* lTitle[3] = {"No cuts","#chi^{2} <0.5","#chi^{2}<5"};
  
  const char* lMeasure[2] = {"Average","StdDev"};
  const char* lRes[8] = {"Mass_pT","Mass_ct","Pt_pT","P_p","Ct_ct","X_X","Y_Y","Z_Z"};
  const char* lDir[2] = {"mean_value","stddev"};

  TFile* input_file[3];
  input_file[0] = new TFile(Std_name);
  input_file[1] = new TFile(O2_name);
  input_file[2] = new TFile(KF_name);
  TFile output_file(Form("%s/%s",folder_name,output_name),"RECREATE");
  TH1D* fHist[2][3];//[mean/stddev][vertexer]

  TCanvas cv("","",800,450);
  cv.Print(Form("%s/%s[",folder_name,pdf_file));
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    //all the resolutions 
    for(int iRes=0; iRes<8; iRes++){
      //loop on the measurements: mean and stddev
      for(int iMeas=0; iMeas<2; iMeas++){          
        //Kalman filter, O2 ,standard vertexer
        for(int iFile=0; iFile<3; iFile++){
          std::cout<<Form("%s/f%s_%s_%c",lDir[0],lMeasure[0],lRes[1],lAM[0])<<"\n";
          fHist[iMeas][iFile] = (TH1D*) input_file[iFile]->Get(Form("%s/f%s_%s_%c",lDir[iMeas],lMeasure[iMeas],lRes[iRes],lAM[iMat]));
          set_style(fHist[iMeas][iFile],iFile);
          fHist[iMeas][iFile]->SetTitle(lTitle[iFile]);
        }
        set_range(fHist[iMeas],5);
        rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%c",lMeasure[iMeas],lRes[iRes],lAM[iMat]),true);
        output_file.cd();
        cv.Write();          
        cv.Print(Form("%s/%s",folder_name,pdf_file)); 
        cv.Clear();
      }
      multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%c",lRes[iRes],lAM[iMat]),true);
      cv.Write();          
      cv.Print(Form("%s/%s",folder_name,pdf_file)); 
      cv.Clear();    
    }
  }
  cv.Print(Form("%s/%s]",folder_name,pdf_file));
}
void compare_chi2(char* KF_name1="MassResolutionKFChi2_1.root",char* KF_name2="MassResolutionKFChi2_2.root",char* KF_name3="MassResolutionKFChi2_5.root",char* KF_name4="MassResolutionKFChi2_50.root",char* output_name="chi2_comparison.root", char* pdf_file="chi2_comparison.pdf", char* folder_name="chi2")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  
  //gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[4] = {"#chi^{2}/NDF<1","#chi^{2}/NDF<2","#chi^{2}/NDF<5","#chi^{2}/NDF<50"};
  //const char* lTitle[3] = {"No cuts","#chi^{2} <0.5","#chi^{2}<5"};
  
  const char* lMeasure[2] = {"Average","StdDev"};
  const char* lRes[8] = {"Mass_pT","Mass_ct","Pt_pT","P_p","Ct_ct","X_X","Y_Y","Z_Z"};
  const char* lDir[2] = {"mean_value","stddev"};

  TFile* input_file[4];
  input_file[0] = new TFile(KF_name1);
  input_file[1] = new TFile(KF_name2);
  input_file[2] = new TFile(KF_name3);
  input_file[3] = new TFile(KF_name4);
  TFile output_file(output_name,"RECREATE");
  TH1D* fHist[2][4];//[mean/stddev][cuts]

  TLegend leg(0.2,0.2,0.4,0.4);
  leg.SetHeader("Kalman filter","C"); // option "C" allows to center the header
  for(int iFile=0; iFile<4; iFile++){
    TH1D* help = new TH1D("help","",10,0,10);
    set_style(help,iFile);
    leg.AddEntry(help,lTitle[iFile],"lp");  
  } 

  TCanvas cv("","",800,450);
  cv.Print(Form("%s[",pdf_file));
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    //all the resolutions 
    for(int iRes=0; iRes<8; iRes++){
      //loop on the measurements: mean and stddev
      for(int iMeas=0; iMeas<2; iMeas++){          
        //Kalman filter, O2 ,standard vertexer
        for(int iFile=0; iFile<4; iFile++){
          std::cout<<Form("%s/f%s_%s_%c",lDir[0],lMeasure[0],lRes[1],lAM[0])<<"\n";
          fHist[iMeas][iFile] = (TH1D*) input_file[iFile]->Get(Form("%s/f%s_%s_%c",lDir[iMeas],lMeasure[iMeas],lRes[iRes],lAM[iMat]));
          set_style(fHist[iMeas][iFile],iFile);
          fHist[iMeas][iFile]->SetTitle(lTitle[iFile]);
        }
        set_range(fHist[iMeas],5);
        rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%c",lMeasure[iMeas],lRes[iRes],lAM[iMat]),false);
        cv.cd();
        leg.Draw("SAME");
        output_file.cd();
        cv.Write();          
        cv.Print(Form("%s",pdf_file)); 
        cv.Clear();
      }
      multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%c",lRes[iRes],lAM[iMat]),false);
      cv.cd(1);
      leg.Draw("SAME");
      cv.Write();          
      cv.Print(Form("%s",pdf_file)); 
      cv.Clear();    
    }
  }
  cv.Print(Form("%s]",pdf_file));
}

void compare_KF_O2(char* KF_name="MassResolutionKFChi2_1.root", char* O2_name="MassResolutionO2.root", char* output_name="vertexer_comparison1.root", char* pdf_file="vertexer_comparison1.pdf", char* folder_name="KF_O2_1")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[2] = {"O^{2} vertexer","Kalman Filter #chi^{2}/NDF<1"};
  //const char* lTitle[3] = {"No cuts","#chi^{2} <0.5","#chi^{2}<5"};
  
  const char* lMeasure[2] = {"Average","StdDev"};
  const char* lRes[8] = {"Mass_pT","Mass_ct","Pt_pT","P_p","Ct_ct","X_X","Y_Y","Z_Z"};
  const char* lDir[2] = {"mean_value","stddev"};

  TFile* input_file[2];
  input_file[0] = new TFile(O2_name);
  input_file[1] = new TFile(KF_name);
  TFile output_file(Form("%s/%s",folder_name,output_name),"RECREATE");
  TH1D* fHist[2][2];//[mean/stddev][vertexer]

  TCanvas cv("","",800,450);
  cv.Print(Form("%s/%s[",folder_name,pdf_file));
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    //all the resolutions 
    for(int iRes=0; iRes<8; iRes++){
      //loop on the measurements: mean and stddev
      for(int iMeas=0; iMeas<2; iMeas++){          
        //Kalman filter, O2 
        for(int iFile=0; iFile<2; iFile++){
          std::cout<<Form("%s/f%s_%s_%c",lDir[0],lMeasure[0],lRes[1],lAM[0])<<"\n";
          fHist[iMeas][iFile] = (TH1D*) input_file[iFile]->Get(Form("%s/f%s_%s_%c",lDir[iMeas],lMeasure[iMeas],lRes[iRes],lAM[iMat]));
          set_style(fHist[iMeas][iFile],iFile);
          fHist[iMeas][iFile]->SetTitle(lTitle[iFile]);
        }
        set_range(fHist[iMeas],5);
        rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%c",lMeasure[iMeas],lRes[iRes],lAM[iMat]),true);
        output_file.cd();
        cv.Write();          
        cv.Print(Form("%s/%s",folder_name,pdf_file)); 
        cv.Clear();
      }
      multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%c",lRes[iRes],lAM[iMat]),true);
      cv.Write();          
      cv.Print(Form("%s/%s",folder_name,pdf_file)); 
      cv.Clear();    
    }
  }
  cv.Print(Form("%s/%s]",folder_name,pdf_file));
}


void compare_distribution(char* Std_name="MassResolutionStd.root", char* KF_name="MassResolutionKF.root", char* O2_name="MassResolutionO2.root", char* output_name="distr_comparison.root", char* pdf_file="distr_comparison.pdf", char* folder_name="nome")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  
  //gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[3] = {"Standard vertexer","O^{2} vertexer","Kalman Filter"};
  const char* lRes[8] = {"fMass_pT_Int","fPt_pT_Int","fP_p_Int","fCt_ct_Int","fX_X_Int","fY_Y_Int","fZ_Z_Int"};

  TFile* input_file[3];
  input_file[0] = new TFile(Std_name);
  input_file[1] = new TFile(O2_name);
  input_file[2] = new TFile(KF_name);
  TFile output_file(output_name,"RECREATE");
  TH1D* fHist[3][3];//[mean/stddev][vertexer]

  TCanvas cv("","",800,450);
  cv.Print(Form("%s[",pdf_file));
  TPaveLabel *titlea = new TPaveLabel(.11,.95,.35,.99,"Antimatter","brndc");
  TPaveLabel *titlem = new TPaveLabel(.11,.95,.35,.99,"Antimatter","brndc");
  
  TLine* line;
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    //all the resolutions (mass,pt,ct)vs(pt,ct)
    for(int iDist=0; iDist<3; iDist++){
      //Kalman filter, O2 ,standard vertexer
      for(int iFile=0; iFile<3; iFile++){
        fHist[iDist][iFile] = (TH1D*) input_file[iFile]->Get(Form("int_distr/%s_%c",lRes[iDist],lAM[iMat]));
        fHist[iDist][iFile]->Scale(1./fHist[iDist][iFile]->GetMaximum());
        set_style(fHist[iDist][iFile],iFile);
        fHist[iDist][iFile]->GetYaxis()->SetTitle("#frac{counts}{max counts}");
        fHist[iDist][iFile]->SetTitle(lTitle[iFile]);
      }
      set_range(fHist[iDist]);
      
      rainbow_plot(fHist[iDist],cv,Form("MultiRainbowPlot_%c",lAM[iMat]),true);
      cv.cd();
      line = new TLine(0.,0.,0.,1.0);
      line->SetLineColor(2);
      line->SetLineStyle(10);
      line->Draw("SAME");

      if(iMat==0)
        cv.SetTitle("Antimatter");
      else
        cv.SetTitle("Matter");
      
      cv.Write();          
      cv.Print(Form("%s",pdf_file)); 
      cv.Clear(); 
    }   
  }

  //range ct 10 - 20 cm

  for(int iFile=0; iFile<3; iFile++){
    fHist[0][iFile] = (TH1D*) input_file[iFile]->Get(Form("fSlice_10_20_%c",lAM[0]));
    fHist[0][iFile]->Scale(1./fHist[0][iFile]->GetMaximum());
    set_style(fHist[0][iFile],iFile);
    fHist[0][iFile]->GetYaxis()->SetTitle("#frac{counts}{max counts}");
    fHist[0][iFile]->SetTitle(lTitle[iFile]);
  }
  set_range(fHist[0]);
  cv.Print(Form("%s]",pdf_file));

  rainbow_plot(fHist[0],cv,Form("MultiRainbowPlotSlice_%c",lAM[0]),true);
  cv.cd();
  line = new TLine(0.,0.,0.,1.0);
  line->SetLineColor(2);
  line->SetLineStyle(10);
  line->Draw("SAME");
  cv.Write();          
  cv.Print(Form("slice10_20.pdf"));
  cv.Clear(); 

}

void ordered_indexes(int Var,int kCutsOrdered[]){
  float tmp;
  float OriginalOrder[kNvariations];
  float Ordered[kNvariations];

  for(int i=0;i<kNvariations;i++){
    Ordered[i]=kCuts[Var][i][0];
    OriginalOrder[i]=kCuts[Var][i][0];
  }

  for (int j=0;j<kNvariations;j++){
    for (int i=kNvariations-2;i>=j;i--){
      if (Ordered[i]>Ordered[i+1])
      {
        tmp = Ordered[i];
        Ordered[i] = Ordered[i+1];
        Ordered[i+1] = tmp;
      }
    }
  }

  for(int i=0;i<kNvariations;i++){
    for(int j=0;j<kNvariations;j++){
      if(Ordered[i]==OriginalOrder[j])
        kCutsOrdered[i]=j;
    }
  }
}

//function to compare efficiencies
void compare_trees(char * findable_name = "efficiencyKF.root",char * task_name = "KFResults.root",char * pdf_file = "compare_trees.pdf")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile findable_file(findable_name);
  TFile task_file(task_name);
  const char lAM[3]{"AM"};
  const char* lName[2]={"Antimatter","Matter"};
  TCanvas cv("","",800*2,450);
  cv.Print(Form("%s[",pdf_file));
  TH1D* histos[2][2];//[eff (ct-pt)][tree]
  //matter and antimatter
  for(int iMat=0;iMat<2;iMat++){
    histos[0][0] = (TH1D*) task_file.Get(Form("efficiency_chi2/fHistEff_Pt_Chi2_deuprot_50.0_%c",lAM[iMat]));
    histos[0][1] = (TH1D*) findable_file.Get(Form("NsigmaTPC/Efficiency/pT/fHistRec_NsigmaTPC_pT_cut_3.0_%c",lAM[iMat]));
    histos[1][0] = (TH1D*) task_file.Get(Form("efficiency_chi2/fHistEff_Ct_Chi2_deuprot_50.0_%c",lAM[iMat]));
    histos[1][1] = (TH1D*) findable_file.Get(Form("NsigmaTPC/Efficiency/ct/fHistRec_NsigmaTPC_ct_cut_3.0_%c",lAM[iMat]));
    histos[0][0]->SetTitle("Task 3KF");
    histos[0][1]->SetTitle("Findable");
    set_range(histos[0]);
    set_range(histos[1]);
    multirainbow_plot(histos,cv,Form("efficiency_%c",lAM[iMat]),true,"PLC PMC");
    cv.Print(Form("%s",pdf_file));
    cv.Clear();
  }
  cv.Print(Form("%s]",pdf_file));
}

void compare_times(char* Std_name="selector_resultsStd.root", char* KF_name1="selector_resultsKFChi2_50.root", char* KF_name2="selector_resultsKFChi2_2.root", char* O2_name="selector_resultsO2.root", char* output_name="time_comparison.root", char* pdf_file="time_comparison.pdf", char* folder_name="nome")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char* lTitle[4] = {"Standard vertexer","Kalman Filter #chi^{2}/NDF<50","Kalman Filter #chi^{2}/NDF<2","O^{2} vertexer"};

  TFile* input_file[4];
  input_file[0] = new TFile(Std_name);
  input_file[1] = new TFile(KF_name1);
  input_file[2] = new TFile(KF_name2);
  input_file[3] = new TFile(O2_name);

  TLine* line = new TLine(0.,1.,4.,1.);
  line->SetLineColor(2);
  line->SetLineStyle(10);

  TH1D* histTime = nullptr;
  TH1D* histComparison = new TH1D("histComparison",";;#frac{process time}{std vtx time}",4,0,4);

  TCanvas cv("","",800,450);
  for(int iVert=0; iVert<4; iVert++){
    histTime = (TH1D*) input_file[iVert]->Get("fTotTime");
    histComparison->SetBinContent(iVert+1,histTime->GetBinContent(1));
    histComparison->SetBinError(iVert+1,0);
    histComparison->GetXaxis()->SetBinLabel(iVert+1,lTitle[iVert]);
  }

  TGaxis*axis = new TGaxis(4,0.8,4,1.1,0.8*histComparison->GetBinContent(1),1.1*histComparison->GetBinContent(1),510,"+L");
  axis->SetTitle("process time (s)");
  axis->SetLabelFont(0);
  axis->SetTitleFont(0);
  
  histComparison->SetMarkerStyle(8);
  histComparison->SetMarkerColor(4);
  //nomalized to the std vertexer
  histComparison->Scale(1./histComparison->GetBinContent(1));
  histComparison->GetYaxis()->SetRangeUser(0.8,1.1);
  histComparison->GetXaxis()->SetLabelSize(0.05);
  histComparison->Draw("P");

   //draw an axis on the right side
  axis->Draw("same");
  line->Draw("same");
  cv.Print(Form("%s",pdf_file));
}
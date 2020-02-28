#include "utils.h"
//function to check the difference of adding the particles in a different order
void permutation_difference(char* KF_name1="MassResolution.root",char* KF_name2="MassResolution.root",char* output_name="difference.root", char* pdf_file="difference.pdf", char* folder_name="difference")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[2] = {"de-pr-pi","pi-pr-de"};
  const char* lHistName[8] = {"MassRes_pT","MassRes_ct","PtRes_pT","PRes_p","CtRes_ct","XRes_X","YRes_Y","ZRes_Z"};
  const char* lQuant[8] = {"p_{T} (GeV/#it{c})","#it{c}t (cm)","p_{T} (GeV/#it{c})","p (GeV/#it{c})","#it{c}t (cm)","x (cm)","y (cm)","z (cm)"};
  const char* lVar[8] = {"p_{T}","#it{c}t","p_{T}","p","#it{c}t","x","y","z"};
  const char* lUnit[8] = {"GeV/#it{c}","cm","GeV/#it{c}","GeV/#it{c}","cm","cm","cm","cm"};
  const int begin[8]={2,1,2,2,1,1,1,1};
  const int canvas_index[12]={1,2,3,4,1,2,3,4,1,2,3,4};
  TFile* input_file[2];
  input_file[0] = new TFile(KF_name1);
  input_file[1] = new TFile(KF_name2);
  TFile output_file(output_name,"RECREATE");

  TCanvas* cv_slice;
  TH2D* plot2d[2][1];
  TH1D* proj[10][2];

  TCanvas cv("","",800,450);
  
  TPaveText pinfo= TPaveText(0.5,0.5,0.91,0.9,"NDC");
  pinfo.SetBorderSize(0);
  pinfo.SetFillStyle(0);
  pinfo.SetTextAlign(30+3);
  pinfo.SetTextFont(42);
  pinfo.AddText("Kalman filter counts difference:");
  pinfo.AddText("particle order de-pr-pi and pi-pr-de");
      
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    for(int iDist=0;iDist<8;iDist++){
      for(int iFile=0; iFile<2; iFile++){
        plot2d[iFile][0] = (TH2D*) input_file[iFile]->Get(Form("TH2Plot/fHist%s_%c",lHistName[iDist],lAM[iMat]));
        TH1D* test = (TH1D*) plot2d[iFile][0]->ProjectionY("fTest",1,10);
        for(int iBin=begin[iDist];iBin<=plot2d[iFile][0]->GetNbinsY();iBin++){
          proj[iBin][iFile]  = (TH1D*) plot2d[iFile][0]->ProjectionX(Form("fSlice_%s_%i_%i_%c",lHistName[iDist],iBin,iFile,lAM[iMat]),iBin,iBin);
          int low = (int)test->GetBinLowEdge(iBin);
          proj[iBin][iFile]->SetTitle(Form("%i #leq %s < %i (%s)",low,lVar[iDist],(int)test->GetBinLowEdge(2)-(int)test->GetBinLowEdge(1)+low,lUnit[iDist]));
          proj[iBin][iFile]->GetYaxis()->SetTitle("counts");
          set_style(proj[iBin][iFile],iFile);
          output_file.cd();
          proj[iBin][iFile]->Write();          
        }
      }

      cv_slice = new TCanvas("","",800*2,450*2);
      cv_slice->Divide(2,2);

      cv_slice->Print(Form("%s/Difference_%s_%c.pdf[",folder_name,lHistName[iDist],lAM[iMat])); 
      int counter=0;
      for(int iBin=begin[iDist];iBin<=10;iBin++){
        if(canvas_index[iBin-begin[iDist]]==1 && iBin!=begin[iDist]){
          counter++;
          cv_slice->cd(1);
          pinfo.Draw();
          cv_slice->Print(Form("%s/Difference_%s_%c.pdf",folder_name,lHistName[iDist],lAM[iMat])); 
          cv_slice->Delete();
          cv_slice = new TCanvas("","",800*2,450*2);
          cv_slice->Divide(2,2);
        }
        cv_slice->cd(canvas_index[iBin-begin[iDist]]);
        proj[iBin][0]->Add(proj[iBin][1],-1);
        proj[iBin][0]->Draw();
      }
      cv_slice->cd(1);
      pinfo.Draw();
      cv_slice->Print(Form("%s/Difference_%s_%c.pdf",folder_name,lHistName[iDist],lAM[iMat])); 
      cv_slice->Delete();
      cv_slice = new TCanvas("","",800*2,450*2);
      cv_slice->Divide(2,2);
      cv_slice->Write();
      cv_slice->Delete();
      cv_slice->Print(Form("%s/Difference_%s_%c.pdf]",folder_name,lHistName[iDist],lAM[iMat])); 
    }
  }   

}

void permutation_comparison(char* KF_name1="MassResolution.root",char* KF_name2="MassResolution.root",char* output_name="comparison.root", char* pdf_file="comparison.pdf", char* folder_name="comparison")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[2] = {"de-pr-pi","pi-pr-de"};
  const char* lHistName[8] = {"MassRes_pT","MassRes_ct","PtRes_pT","PRes_p","CtRes_ct","XRes_X","YRes_Y","ZRes_Z"};
  const char* lQuant[8] = {"p_{T} (GeV/#it{c})","#it{c}t (cm)","p_{T} (GeV/#it{c})","p (GeV/#it{c})","#it{c}t (cm)","x (cm)","y (cm)","z (cm)"};
  const char* lVar[8] = {"p_{T}","#it{c}t","p_{T}","p","#it{c}t","x","y","z"};
  const char* lUnit[8] = {"GeV/#it{c}","cm","GeV/#it{c}","GeV/#it{c}","cm","cm","cm","cm"};
  const int begin[8]={2,1,2,2,1,1,1,1};
  const int canvas_index[12]={1,2,3,4,1,2,3,4,1,2,3,4};
  TFile* input_file[2];
  input_file[0] = new TFile(KF_name1);
  input_file[1] = new TFile(KF_name2);
  TFile output_file(output_name,"RECREATE");

  TCanvas* cv_slice;
  TH2D* plot2d[2][1];
  TH1D* proj[10][2];

  TCanvas cv("","",800,450);
  
  TLegend leg(0.2,0.2,0.4,0.4);
  leg.SetHeader("Kalman filter order:");
  for(int iFile=0; iFile<2; iFile++){
    TH1D* help = new TH1D("help","",10,0,10);
    set_style(help,iFile);
    leg.AddEntry(help,lTitle[iFile],"lp");  
  } 
  
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    for(int iDist=0;iDist<8;iDist++){
      for(int iFile=0; iFile<2; iFile++){
        plot2d[iFile][0] = (TH2D*) input_file[iFile]->Get(Form("TH2Plot/fHist%s_%c",lHistName[iDist],lAM[iMat]));
        TH1D* test = (TH1D*) plot2d[iFile][0]->ProjectionY("fTest",1,10);
        for(int iBin=begin[iDist];iBin<=plot2d[iFile][0]->GetNbinsY();iBin++){
          proj[iBin][iFile]  = (TH1D*) plot2d[iFile][0]->ProjectionX(Form("fSlice_%s_%i_%i_%c",lHistName[iDist],iBin,iFile,lAM[iMat]),iBin,iBin);
          int low = (int)test->GetBinLowEdge(iBin);
          proj[iBin][iFile]->SetTitle(Form("%i #leq %s < %i (%s)",low,lVar[iDist],(int)test->GetBinLowEdge(2)-(int)test->GetBinLowEdge(1)+low,lUnit[iDist]));
          proj[iBin][iFile]->GetYaxis()->SetTitle("counts difference");
          set_style(proj[iBin][iFile],iFile);
          output_file.cd();
          proj[iBin][iFile]->Write();          
        }

      }

      cv_slice = new TCanvas("","",800*2,450*2);
      cv_slice->Divide(2,2);

      cv_slice->Print(Form("%s/Difference_%s_%c.pdf[",folder_name,lHistName[iDist],lAM[iMat])); 
      int counter=0;
      for(int iBin=begin[iDist];iBin<=10;iBin++){
        if(canvas_index[iBin-begin[iDist]]==1 && iBin!=begin[iDist]){
          counter++;
          cv_slice->cd(1);
          leg.Draw();
          cv_slice->Print(Form("%s/Difference_%s_%c.pdf",folder_name,lHistName[iDist],lAM[iMat])); 
          cv_slice->Delete();
          cv_slice = new TCanvas("","",800*2,450*2);
          cv_slice->Divide(2,2);
        }
        cv_slice->cd(canvas_index[iBin-begin[iDist]]);
        for(int iFile=0;iFile<2;iFile++)
          proj[iBin][iFile]->Draw("SAME");
        set_range(proj[iBin]);
      }
      cv_slice->cd(1);
      leg.Draw();
      cv_slice->Print(Form("%s/Comparison_%s_%c.pdf",folder_name,lHistName[iDist],lAM[iMat])); 
      cv_slice->Delete();
      cv_slice = new TCanvas("","",800*2,450*2);
      cv_slice->Divide(2,2);
      cv_slice->Write();
      cv_slice->Delete();
      cv_slice->Print(Form("%s/Comparison_%s_%c.pdf]",folder_name,lHistName[iDist],lAM[iMat])); 
    }
  }   

}

//function to check the topological difference of the second peak
//list of functions to compare the topology for ctrec>ctgen and ctrec<ctgen
void compare_sides(char* input_name="selector_resultsKF.root", char* output_name="sides_comparison.root", char* pdf_file="sides_comparison.pdf", char* folder_name="side",bool norm=true)
{
  gStyle->SetOptStat(0);  
  gSystem->Exec(Form("mkdir %s",folder_name));

  TFile input_file(input_name);
  const char lRL[3]{"RL"};
  const char *lSpecies[3]{"d", "p", "pi"};
  const char* lSide[2] = {"ct_{rec}>ct_{gen}","ct_{rec}<ct_{gen}"};
  TH1D* histVert[2][2];
  TH1D* hist[7][3][2];
  TFile output_file(Form("%s/%s",folder_name,output_name),"recreate");

  TLegend legend(0.2,0.2,0.4,0.4);
  for(int iFile=0; iFile<2; iFile++){
    TH1D* help = new TH1D("help","",10,0,10);
    set_style(help,iFile);
    legend.AddEntry(help,lSide[iFile],"lp");  
  } 

  for(int iSide = 0; iSide < 2; iSide++){
      histVert[0][iSide] = (TH1D*) input_file.Get(Form("fHistVertexChi2_%c",lRL[iSide]));
      histVert[1][iSide] = (TH1D*) input_file.Get(Form("fCosPointingAngle_%c",lRL[iSide]));
      /// DCA to primary vertex
      for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
        hist[0][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2PrimaryvtxXY_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[1][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2PrimaryvtxZ_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[2][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2Primaryvtx_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[3][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2DecayvtxXY_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[4][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2DecayvtxZ_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[5][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2Decayvtx_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[6][iSpecies][iSide] = (TH1D*) input_file.Get(Form("lTrackDistance_%s-%s_%c",lSpecies[iSpecies], lSpecies[(iSpecies + 1) % 3],lRL[iSide]));

      }
  }
  //TFile output_file(Form("%s/%s",folder_name,input_name),"recreate");
  TCanvas cv("","",800,450);
  cv.Print(Form("%s/%s[",folder_name,pdf_file));
  //not DCA 
  for(int iHist=0; iHist<2; iHist++){
    for(int iSide=0; iSide<2; iSide++){
        set_style(histVert[iHist][iSide],iSide);
        histVert[iHist][iSide]->SetTitle(histVert[iHist][1]->GetName());
        if(norm){
          histVert[iHist][iSide]->Scale(1./histVert[iHist][iSide]->GetEntries());
        }
    }
    set_range(histVert[iHist]);
    histVert[iHist][0]->Draw("");
    histVert[iHist][1]->Draw("SAME");
    legend.Draw();
    cv.Print(Form("%s/%s",folder_name,pdf_file));
    cv.Write();
  }
  //DCA
  for(int iHist=0; iHist<7; iHist++){
    for(int iSpecies=0; iSpecies<3; iSpecies++){
      for(int iSide=0; iSide<2; iSide++){
        set_style(hist[iHist][iSpecies][iSide],iSide);
        hist[iHist][iSpecies][iSide]->SetTitle(hist[iHist][iSpecies][1]->GetName());
        if(norm){
          hist[iHist][iSpecies][iSide]->Scale(1./hist[iHist][iSpecies][iSide]->GetEntries());
        }
      }
      set_range(hist[iHist][iSpecies]);
      hist[iHist][iSpecies][0]->Draw("");
      hist[iHist][iSpecies][1]->Draw("SAME");
      legend.Draw();
      cv.Print(Form("%s/%s",folder_name,pdf_file));
      cv.Write();
    }
  } 
  cv.Print(Form("%s/%s]",folder_name,pdf_file));
}


//this function compare the two distribution dividing the two normalized sides
void divide_sides(char* input_name="selector_resultsKF.root", char* output_name="rates_comparison.root", char* pdf_file="rates_comparison.pdf", char* folder_name="rate")
{
  gStyle->SetOptStat(0);  
  gSystem->Exec(Form("mkdir %s",folder_name));

  TFile input_file(input_name);
  const char lRL[3]{"RL"};
  const char *lSpecies[3]{"d", "p", "pi"};
  const char* lSide[2] = {"ct_{rec}>ct_{gen}","ct_{rec}<ct_{gen}"};
  TH1D* histVert[2][2];
  TH1D* hist[7][3][2];

  for(int iSide = 0; iSide < 2; iSide++){
      histVert[0][iSide] = (TH1D*) input_file.Get(Form("fHistVertexChi2_%c",lRL[iSide]));
      histVert[1][iSide] = (TH1D*) input_file.Get(Form("fCosPointingAngle_%c",lRL[iSide]));
      /// DCA to primary vertex
      for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
        hist[0][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2PrimaryvtxXY_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[1][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2PrimaryvtxZ_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[2][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2Primaryvtx_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[3][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2DecayvtxXY_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[4][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2DecayvtxZ_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[5][iSpecies][iSide] = (TH1D*) input_file.Get(Form("fDCA2Decayvtx_%s_%c", lSpecies[iSpecies],lRL[iSide]));
        hist[6][iSpecies][iSide] = (TH1D*) input_file.Get(Form("lTrackDistance_%s-%s_%c",lSpecies[iSpecies], lSpecies[(iSpecies + 1) % 3],lRL[iSide]));

      }
  }

  //TFile output_file(Form("%s/%s",folder_name,input_name),"recreate");
  TCanvas cv("","",800,450);
  cv.Print(Form("%s/%s[",folder_name,pdf_file));
  //not DCA 
  for(int iHist=0; iHist<2; iHist++){
    for(int iSide=0; iSide<2; iSide++){
      set_style(histVert[iHist][iSide],iSide);
      histVert[iHist][iSide]->SetTitle(histVert[iHist][1]->GetName());
      histVert[iHist][iSide]->Scale(1./histVert[iHist][iSide]->GetEntries());
    }
    histVert[iHist][0]->Divide(histVert[iHist][1]);
    SetRateErrors(histVert[iHist][0],histVert[iHist][1]);
    histVert[iHist][0]->Draw("");
    cv.Print(Form("%s/%s",folder_name,pdf_file));
  }
  //DCA
  for(int iHist=0; iHist<7; iHist++){
    for(int iSpecies=0; iSpecies<3; iSpecies++){
      for(int iSide=0; iSide<2; iSide++){
        set_style(hist[iHist][iSpecies][iSide],iSide);
        hist[iHist][iSpecies][iSide]->SetTitle(hist[iHist][iSpecies][1]->GetName());
        hist[iHist][iSpecies][iSide]->Scale(1./hist[iHist][iSpecies][iSide]->GetEntries());
      }
      hist[iHist][iSpecies][0]->Divide(hist[iHist][iSpecies][1]);
      hist[iHist][iSpecies][0]->Draw("");
      cv.Print(Form("%s/%s",folder_name,pdf_file));
    }
  } 
  cv.Print(Form("%s/%s]",folder_name,pdf_file));
}
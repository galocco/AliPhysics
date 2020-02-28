#include "utils.h"
//function to plot the TH2D varrec-vargen vs vargen
void plot2D_All(char* Std_name="MassResolutionStd.root", char* KF_name="MassResolutionKF.root", char* O2_name="MassResolutionO2.root", char* output_name="plot2d.root", char* folder_name="plot2dAll")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[3] = {"Standard vertexer","O^{2} vertexer","Kalman Filter"};
  const char* lHistName[8] = {"MassRes_pT","MassRes_ct","PtRes_pT","PRes_p","CtRes_ct","XRes_X","YRes_Y","ZRes_Z"};
  const char* lQuant[8] = {"p_{T} (GeV/#it{c})","#it{c}t (cm)","p_{T} (GeV/#it{c})","p (GeV/#it{c})","#it{c}t (cm)","x (cm)","y (cm)","z (cm)"};

  TFile* input_file[3];
  input_file[0] = new TFile(Std_name);
  input_file[1] = new TFile(O2_name);
  input_file[2] = new TFile(KF_name);
  TFile output_file(Form("%s/%s",folder_name,output_name),"RECREATE");
  TCanvas cv("","",800*3,450);
  TCanvas cv_slice("","",800*3,450);
  cv_slice.Divide(4,4);
  TH2D* plot2d[3][1];
  TH1D* proj[10][3];
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    for(int iDist=0;iDist<8;iDist++){
      //Kalman filter, O2 ,standard vertexer
      for(int iFile=0; iFile<3; iFile++){
        std::cout<<Form("TH2Plot/fHist%s_%c",lHistName[iDist],lAM[iMat])<<"\n";
        plot2d[iFile][0] = (TH2D*) input_file[iFile]->Get(Form("TH2Plot/fHist%s_%c",lHistName[iDist],lAM[iMat]));
        for(int x=0;x<10;x++)
          for(int y=0;y<10;y++)
            if(plot2d[iFile][0]->GetBinContent(x,y)<1){
              plot2d[iFile][0]->SetBinContent(x,y,1000);
            }
        
        plot2d[iFile][0]->SetTitle(lTitle[iFile]);
      }

      multirainbow_plot(plot2d,cv,Form("MultiRainbowPlot_%s_%c",lHistName[iDist],lAM[iMat]),false,"colz");
      cv.cd();

      cv.Print(Form("%s/Plot_%s_%c.pdf[",folder_name,lHistName[iDist],lAM[iMat])); 
      cv.Write();          
      cv.Print(Form("%s/Plot_%s_%c.pdf",folder_name,lHistName[iDist],lAM[iMat])); 
      cv.Print(Form("%s/Plot_%s_%c.pdf]",folder_name,lHistName[iDist],lAM[iMat])); 
      cv.Clear(); 
    }
  }   

}

//function to plot the TH2D varrec-vargen vs vargen
void plot2D_KF_O2(char* KF_name="MassResolutionKFChi2_2.root", char* O2_name="MassResolutionO2.root", char* output_name="plot2D_KF_O2.root", char* folder_name="plot2D_KF_O2")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[2] = {"Kalman Filter #chi^{2}/NDF<2","O^{2} vertexer"};
  const char* lHistName[8] = {"MassRes_pT","MassRes_ct","PtRes_pT","PRes_p","CtRes_ct","XRes_X","YRes_Y","ZRes_Z"};
  const char* lQuant[8] = {"p_{T} (GeV/#it{c})","#it{c}t (cm)","p_{T} (GeV/#it{c})","p (GeV/#it{c})","#it{c}t (cm)","x (cm)","y (cm)","z (cm)"};

  TFile* input_file[2];
  input_file[1] = new TFile(O2_name);
  input_file[0] = new TFile(KF_name);
  TFile output_file(Form("%s/%s",folder_name,output_name),"RECREATE");
  TCanvas cv("","",800*2,450);
  TCanvas cv_slice("","",800*2,450);
  cv_slice.Divide(4,4);
  TH2D* plot2d[2][1];
  TH1D* proj[10][2];
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    for(int iDist=0;iDist<8;iDist++){
      //Kalman filter, O2 ,standard vertexer
      for(int iFile=0; iFile<2; iFile++){
        std::cout<<Form("TH2Plot/fHist%s_%c",lHistName[iDist],lAM[iMat])<<"\n";
        plot2d[iFile][0] = (TH2D*) input_file[iFile]->Get(Form("TH2Plot/fHist%s_%c",lHistName[iDist],lAM[iMat]));
        for(int x=0;x<10;x++)
          for(int y=0;y<10;y++)
            if(plot2d[iFile][0]->GetBinContent(x,y)<1){
              plot2d[iFile][0]->SetBinContent(x,y,1);
            }
        
        plot2d[iFile][0]->SetTitle(lTitle[iFile]);
      }

      multirainbow_plot(plot2d,cv,Form("MultiRainbowPlot_%s_%c",lHistName[iDist],lAM[iMat]),false,"colz");
      cv.cd();

      cv.Print(Form("%s/Plot_%s_%c.pdf[",folder_name,lHistName[iDist],lAM[iMat])); 
      cv.Write();          
      cv.Print(Form("%s/Plot_%s_%c.pdf",folder_name,lHistName[iDist],lAM[iMat])); 
      cv.Print(Form("%s/Plot_%s_%c.pdf]",folder_name,lHistName[iDist],lAM[iMat])); 
      cv.Clear(); 
    }
  }   

}

//function to plot the distribution of the resolution bin per bin
void plot_res_slices(bool single=false, char* Std_name="MassResolutionStd.root", char* KF_name="MassResolutionKFChi2_50.root", char* O2_name="MassResolutionO2.root", char* output_name="distr.root", char* pdf_file="distr.pdf", char* folder_name="distr")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[3] = {"Standard vertexer","O^{2} vertexer","Kalman Filter #chi^{2}/NDF<50"};
  const char* lHistName[8] = {"MassRes_pT","MassRes_ct","PtRes_pT","PRes_p","CtRes_ct","XRes_X","YRes_Y","ZRes_Z"};
  const char* lQuant[8] = {"p_{T} (GeV/#it{c})","#it{c}t (cm)","p_{T} (GeV/#it{c})","p (GeV/#it{c})","#it{c}t (cm)","x (cm)","y (cm)","z (cm)"};
  const char* lVar[8] = {"p_{T}","#it{c}t","p_{T}","p","#it{c}t","x","y","z"};
  const char* lUnit[8] = {"GeV/#it{c}","cm","GeV/#it{c}","GeV/#it{c}","cm","cm","cm","cm"};
  const int begin[8]={2,1,2,2,1,1,1,1};
  const int canvas_index[12]={1,2,3,4,1,2,3,4,1,2,3,4};
  TFile* input_file[3];
  input_file[0] = new TFile(Std_name);
  input_file[1] = new TFile(O2_name);
  input_file[2] = new TFile(KF_name);
  TFile output_file(output_name,"RECREATE");
  TCanvas* cv_slice;
  TH2D* plot2d[3][1];
  TH1D* proj[10][3];
  TCanvas cv("","",800,450);
  TLegend leg(0.7,0.7,0.9,0.9);

  for(int iFile=0; iFile<3; iFile++){
    TH1D* help = new TH1D("help","",10,0,10);
    set_style(help,iFile);
    leg.AddEntry(help,lTitle[iFile],"lp");  
  } 
  //cv_slice.Print(Form("%s/Slice_%s_%c.pdf[",folder_name,lHistName[0],lAM[0])); 
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    for(int iDist=0;iDist<8;iDist++){
      //Kalman filter, O2 ,standard vertexer
      //leg->Clear();
      for(int iFile=0; iFile<3; iFile++){
        plot2d[iFile][0] = (TH2D*) input_file[iFile]->Get(Form("TH2Plot/fHist%s_%c",lHistName[iDist],lAM[iMat]));
        TH1D* test = (TH1D*) plot2d[iFile][0]->ProjectionY("fTest",1,10);
        for(int iBin=begin[iDist];iBin<=plot2d[iFile][0]->GetNbinsY();iBin++){
          proj[iBin][iFile]  = (TH1D*) plot2d[iFile][0]->ProjectionX(Form("fSlice_%s_%i_%i_%c",lHistName[iDist],iBin,iFile,lAM[iMat]),iBin,iBin);
          int low = (int)test->GetBinLowEdge(iBin);
          proj[iBin][iFile]->SetTitle(Form("%i #leq %s < %i (%s)",low,lVar[iDist],(int)test->GetBinLowEdge(2)-(int)test->GetBinLowEdge(1)+low,lUnit[iDist]));
          proj[iBin][iFile]->Scale(1./proj[iBin][iFile]->GetMaximum());          
          proj[iBin][iFile]->GetYaxis()->SetTitle("#frac{counts}{max counts}");
          set_style(proj[iBin][iFile],iFile);
          output_file.cd();
          proj[iBin][iFile]->Write();
          
        }
      }
      cv_slice = new TCanvas("","",800*2,450*2);
      if(!single){
        cv_slice->Divide(2,2);
      }
      int counter=0;
      for(int iBin=begin[iDist];iBin<=10;iBin++){
        if(canvas_index[iBin-begin[iDist]]==1 && iBin!=begin[iDist] && !single){
          counter++;
          cv_slice->cd(1);
          leg.Draw();
          cv_slice->Print(Form("%s/Slice_%s_%c_%i.pdf",folder_name,lHistName[iDist],lAM[iMat],counter)); 
          cv_slice->Delete();
          cv_slice = new TCanvas("","",800*2,450*2);
          cv_slice->Divide(2,2);
        }
        if(!single){
          cv_slice->cd(canvas_index[iBin-begin[iDist]]);
        }
        else{
          cv_slice->cd();
        }
        for(int iFile=0; iFile<3; iFile++){
          proj[iBin][iFile]->Draw("SAME");
        }
        if(single){
          counter++;
          leg.Draw();
          cv_slice->Print(Form("%s/Slice_%s_%c_%i.pdf",folder_name,lHistName[iDist],lAM[iMat],counter)); 
          cv_slice->Write();
          cv_slice = new TCanvas("","",800*2,450*2);
        }
      }
      if(!single){
        cv_slice->cd(1);
        leg.Draw();
        cv_slice->Print(Form("%s/Slice_%s_%c_%i.pdf",folder_name,lHistName[iDist],lAM[iMat],counter+1)); 
        cv_slice->Delete();
        cv_slice = new TCanvas("","",800*2,450*2);
        cv_slice->Divide(2,2);
        cv_slice->Write();
      }
      cv_slice->Delete();
    }
  }   

  //cv_slice.Print(Form("%s/Slice_%s_%c.pdf]",folder_name,lHistName[0],lAM[0]));
}


//function to plot the distribution of the resolution bin per bin
void compare_distr_chi2(char* KF_name1="MassResolutionKFChi2_1.root",char* KF_name2="MassResolutionKFChi2_2.root",char* KF_name3="MassResolutionKFChi2_5.root",char* KF_name4="MassResolutionKFChi2_50.root",char* output_name="distr_chi2.root", char* pdf_file="distr_chi2.pdf", char* folder_name="distr_chi2")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[4] = {"#chi^{2}/NDF<1","#chi^{2}/NDF<2","#chi^{2}/NDF<5","#chi^{2}/NDF<50"};
  const char* lHistName[8] = {"MassRes_pT","MassRes_ct","PtRes_pT","PRes_p","CtRes_ct","XRes_X","YRes_Y","ZRes_Z"};
  const char* lQuant[8] = {"p_{T} (GeV/#it{c})","#it{c}t (cm)","p_{T} (GeV/#it{c})","p (GeV/#it{c})","#it{c}t (cm)","x (cm)","y (cm)","z (cm)"};
  const char* lVar[8] = {"p_{T}","#it{c}t","p_{T}","p","#it{c}t","x","y","z"};
  const char* lUnit[8] = {"GeV/#it{c}","cm","GeV/#it{c}","GeV/#it{c})","cm","cm","cm","cm"};
  const int begin[8]={2,1,2,2,1,1,1,1};
  const int canvas_index[12]={1,2,3,4,1,2,3,4,1,2,3,4};
  TFile* input_file[4];
  input_file[0] = new TFile(KF_name1);
  input_file[1] = new TFile(KF_name2);
  input_file[2] = new TFile(KF_name3);
  input_file[3] = new TFile(KF_name4);
  TFile output_file(output_name,"RECREATE");
  TCanvas* cv_slice;
  TH2D* plot2d[4][1];
  TH1D* proj[10][4];
  TCanvas cv("","",800,450);
  TLegend leg(0.6,0.6,0.9,0.9);
  leg.SetHeader("Kalman filter:");
  for(int iFile=0; iFile<4; iFile++){
    TH1D* help = new TH1D("help","",10,0,10);
    set_style(help,iFile);
    leg.AddEntry(help,lTitle[iFile],"lp");  
  } 
  //cv_slice.Print(Form("%s/Slice_%s_%c.pdf[",folder_name,lHistName[0],lAM[0])); 
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    for(int iDist=0;iDist<8;iDist++){
      //Kalman filter, O2 ,standard vertexer
      //leg->Clear();
      for(int iFile=0; iFile<4; iFile++){
        plot2d[iFile][0] = (TH2D*) input_file[iFile]->Get(Form("TH2Plot/fHist%s_%c",lHistName[iDist],lAM[iMat]));
        TH1D* test = (TH1D*) plot2d[iFile][0]->ProjectionY("fTest",1,10);
        for(int iBin=begin[iDist];iBin<=plot2d[iFile][0]->GetNbinsY();iBin++){
          proj[iBin][iFile]  = (TH1D*) plot2d[iFile][0]->ProjectionX(Form("fSlice_%s_%i_%i_%c",lHistName[iDist],iBin,iFile,lAM[iMat]),iBin,iBin);
          int low = (int)test->GetBinLowEdge(iBin);
          proj[iBin][iFile]->SetTitle(Form("%i #leq %s < %i (%s)",low,lVar[iDist],(int)test->GetBinLowEdge(2)-(int)test->GetBinLowEdge(1)+low,lUnit[iDist]));
          proj[iBin][iFile]->Scale(1./proj[iBin][iFile]->GetMaximum());          
          proj[iBin][iFile]->GetYaxis()->SetTitle("#frac{counts}{max counts}");
          set_style(proj[iBin][iFile],iFile);
          output_file.cd();
          proj[iBin][iFile]->Write();
          
        }
      }
      cv_slice = new TCanvas("","",800*2,450*2);
      cv_slice->Divide(2,2);

      //cv_slice->Print(Form("%s/Slice_%s_%c.pdf[",folder_name,lHistName[iDist],lAM[iMat])); 
      int counter=0;
      for(int iBin=begin[iDist];iBin<=10;iBin++){
        if(canvas_index[iBin-begin[iDist]]==1 && iBin!=begin[iDist]){
          counter++;
          cv_slice->cd(1);
          leg.Draw();
          cv_slice->Print(Form("%s/Slice_%s_%c_%i.pdf",folder_name,lHistName[iDist],lAM[iMat],counter)); 
          cv_slice->Delete();
          cv_slice = new TCanvas("","",800*2,450*2);
          cv_slice->Divide(2,2);
        }
        //cv_slice->cd(((iBin-begin[iDist]+1)%5)+counter);
        cv_slice->cd(canvas_index[iBin-begin[iDist]]);
        //std::cout<<((iBin-begin[iDist]+1)%5)+counter<<" "<<((iBin-begin[iDist]+1+counter)%5)<<" "<<iBin<<"\n";
        for(int iFile=0; iFile<4; iFile++){
          proj[iBin][iFile]->Draw("SAME");
          /*
          if(iBin==2){
            cv.cd();
            proj[iBin][iFile]->Draw("SAME")
          }
          */
        }
        //if(iBin==2)
        //  cv_slice->Print(Form("%s/Slice_%s_%c_%i.pdf",folder_name,lHistName[iDist],lAM[iMat],counter));
        //leg->Draw();
        //std::cout<<iBin<<"\n";
      }
      cv_slice->cd(1);
      leg.Draw();
      cv_slice->Print(Form("%s/Slice_%s_%c_%i.pdf",folder_name,lHistName[iDist],lAM[iMat],counter+1)); 
      cv_slice->Delete();
      cv_slice = new TCanvas("","",800*2,450*2);
      cv_slice->Divide(2,2);
      cv_slice->Write();//Print(Form("%s/Slice_%s_%c.pdf",folder_name,lHistName[0],lAM[0])); 
      cv_slice->Delete();
      //cv_slice->Print(Form("%s/Slice_%s_%c.pdf]",folder_name,lHistName[iDist],lAM[iMat])); 
    }
  }   

  //cv_slice.Print(Form("%s/Slice_%s_%c.pdf]",folder_name,lHistName[0],lAM[0]));
}


//function to plot the distribution of the resolution bin per bin
void distr_KF_O2(char* KF_name="MassResolutionKFChi2_2.root",char* O2_name="MassResolutionO2.root",bool single=false, char* output_name="distr_KF_O2.root", char* pdf_file="distr_KF_O2.pdf", char* folder_name="distr_KF_O2")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[2] = {"O^{2} vertexer","Kalman filter #chi^{2}/NDF<2"};
  const char* lHistName[8] = {"MassRes_pT","MassRes_ct","PtRes_pT","PRes_p","CtRes_ct","XRes_X","YRes_Y","ZRes_Z"};
  const char* lQuant[8] = {"p_{T} (GeV/#it{c})","#it{c}t (cm)","p_{T} (GeV/#it{c})","p (GeV/#it{c})","#it{c}t (cm)","x (cm)","y (cm)","z (cm)"};
  const char* lVar[8] = {"p_{T}","#it{c}t","p_{T}","p","#it{c}t","x","y","z"};
  const char* lUnit[8] = {"GeV/#it{c}","cm","GeV/#it{c}","GeV/#it{c})","cm","cm","cm","cm"};
  const int begin[8]={2,1,2,2,1,1,1,1};
  const int canvas_index[12]={1,2,3,4,1,2,3,4,1,2,3,4};
  TFile* input_file[2];
  input_file[0] = new TFile(O2_name);
  input_file[1] = new TFile(KF_name);
  TFile output_file(Form("%s/%s",folder_name,output_name),"RECREATE");

  TCanvas* cv_slice;
  TH2D* plot2d[4][1];
  TH1D* proj[10][4];
  TCanvas cv("","",800,450);
  TLegend leg(0.6,0.6,0.95,0.95);
  for(int iFile=0; iFile<2; iFile++){
    TH1D* help = new TH1D("help","",10,0,10);
    set_style(help,iFile);
    leg.AddEntry(help,lTitle[iFile],"lp");  
  } 
  //cv_slice.Print(Form("%s/Slice_%s_%c.pdf[",folder_name,lHistName[0],lAM[0])); 
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    for(int iDist=0;iDist<8;iDist++){
      //Kalman filter, O2 ,standard vertexer
      for(int iFile=0; iFile<2; iFile++){
        plot2d[iFile][0] = (TH2D*) input_file[iFile]->Get(Form("TH2Plot/fHist%s_%c",lHistName[iDist],lAM[iMat]));
        TH1D* test = (TH1D*) plot2d[iFile][0]->ProjectionY("fTest",1,10);
        for(int iBin=begin[iDist];iBin<=plot2d[iFile][0]->GetNbinsY();iBin++){
          proj[iBin][iFile]  = (TH1D*) plot2d[iFile][0]->ProjectionX(Form("fSlice_%s_%i_%i_%c",lHistName[iDist],iBin,iFile,lAM[iMat]),iBin,iBin);
          int low = (int)test->GetBinLowEdge(iBin);
          proj[iBin][iFile]->SetTitle(Form("%i #leq %s < %i (%s)",low,lVar[iDist],(int)test->GetBinLowEdge(2)-(int)test->GetBinLowEdge(1)+low,lUnit[iDist]));
          proj[iBin][iFile]->Scale(1./proj[iBin][iFile]->GetMaximum());          
          proj[iBin][iFile]->GetYaxis()->SetTitle("#frac{counts}{max counts}");
          set_style(proj[iBin][iFile],iFile);
          output_file.cd();
          proj[iBin][iFile]->Write();
          
        }
      }
      cv_slice = new TCanvas("","",800*2,450*2);
      if(!single){
        cv_slice->Divide(2,2);
      }
      int counter=0;
      for(int iBin=begin[iDist];iBin<=10;iBin++){
        if(canvas_index[iBin-begin[iDist]]==1 && iBin!=begin[iDist] && !single){
          counter++;
          cv_slice->cd(1);
          leg.Draw();
          cv_slice->Print(Form("%s/Slice_%s_%c_%i.pdf",folder_name,lHistName[iDist],lAM[iMat],counter)); 
          cv_slice->Delete();
          cv_slice = new TCanvas("","",800*2,450*2);
          cv_slice->Divide(2,2);
        }
        if(single){
          cv_slice->cd(1);
        }
        else{
          cv_slice->cd(canvas_index[iBin-begin[iDist]]);
        }
        for(int iFile=0; iFile<2; iFile++){
          proj[iBin][iFile]->Draw("SAME");
        }
        if(single){
          counter++;
          cv_slice->Print(Form("%s/Slice_%s_%c_%i.pdf",folder_name,lHistName[iDist],lAM[iMat],counter)); 
          leg.Draw();
        }
      }
      if(!single){
        cv_slice->cd(1);
        leg.Draw();
        cv_slice->Print(Form("%s/Slice_%s_%c_%i.pdf",folder_name,lHistName[iDist],lAM[iMat],counter+1)); 
      }
      cv_slice->Delete();
      cv_slice = new TCanvas("","",800*2,450*2);
      cv_slice->Divide(2,2);
      cv_slice->Write();//Print(Form("%s/Slice_%s_%c.pdf",folder_name,lHistName[0],lAM[0])); 
      cv_slice->Delete();
      //cv_slice->Print(Form("%s/Slice_%s_%c.pdf]",folder_name,lHistName[iDist],lAM[iMat])); 
    }
  }   

  //cv_slice.Print(Form("%s/Slice_%s_%c.pdf]",folder_name,lHistName[0],lAM[0]));
}


#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include "Hyp3FindConfig.h"


void MassResolution(char * input_name = "selector_results.root",char * output_name = "MassResolutio.root",bool fit_slice=true)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptTitle(0);
  TFile input_file(input_name);

  const char  lAM[2] = {'A','M'};
  const char* lProjSet[2] = {"zx","zy"};
  const char* lMeasure[3] = {"Res","Average","StdDev"};
  const char* lRes[8] = {"Mass_pT","Mass_ct","Pt_pT","P_p","Ct_ct","X_X","Y_Y","Z_Z"};
  const char* lHistName[8] = {"MassRes_pT","MassRes_ct","PtRes_pT","PRes_p","CtRes_ct","XRes_X","YRes_Y","ZRes_Z"};
  //res solution of: mass,pt,p,ct,x,y,z
  const char* lUnit[8] = {"(GeV/#it{c}^{2})","(GeV/#it{c}^2)","(GeV/#it{c})","(GeV/#it{c})","(cm)","(cm)","(cm)","(cm)"};
  const char* lVar[8] = {"m","m","p_{T}","p","ct","x","y","z"};
  const char* lQuant[8] = {"p_{T} (GeV/#it{c})","#it{c}t (cm)","p_{T} (GeV/#it{c})","p (GeV/#it{c})","#it{c}t (cm)","x (cm)","y (cm)","z (cm)"};

  TFile output_file(output_name,"RECREATE");
  const char *lDir[6]={"TH2Plot","mean_value","stddev","int_distr","multi_plots","chi2"};
  TDirectory* subdir[6];
  int ndir = (fit_slice) ? 6 : 5;
  for(int i=0; i<ndir; i++)
    subdir[i] = output_file.mkdir(lDir[i]);

  TCanvas MultiCV("","",1800,400);
  MultiCV.Divide(3,1);

  TH2D* fHistProto;
  TH3D* fHistVsCuts;
  TH1D* fMeasure[4];
  TH1D* fChi2;
  TObjArray fSlices;
  //dm vs pt/ct -- dct vs ct -- dpt vs pt
  for(int iRes=0; iRes<8; iRes++){
    //matter or antimat
    for(int iMat=0; iMat<2; iMat++){
      fHistProto = (TH2D*) input_file.Get(Form("fHist%s_%c",lHistName[iRes], lAM[iMat]));

      if(fit_slice){
        fHistProto->FitSlicesX(0, 1, 10, 0, "QNR", &fSlices);
        fMeasure[1] = (TH1D*) fSlices[1];
        fMeasure[2] = (TH1D*) fSlices[2];
        fMeasure[1]->SetName(Form("fAverage_%s_%c",lRes[iRes],lAM[iMat]));
        fMeasure[1]->SetTitle(Form(";%s;mean[%srec-%sgen] %s",lQuant[iRes],lVar[iRes],lVar[iRes],lUnit[iRes]));
        fMeasure[2]->SetName(Form("fStdDev_%s_%c",lRes[iRes],lAM[iMat]));
        fMeasure[2]->SetTitle(Form(";%s;std dev[%srec-%sgen] %s",lQuant[iRes],lVar[iRes],lVar[iRes],lUnit[iRes]));        
        subdir[5]->cd();
        fChi2 = (TH1D*) fSlices[3];
        fChi2->SetName(Form("fChi2_%s_%c",lRes[iRes],lAM[iMat]));
        fChi2->Write();
      }
      else{
        fMeasure[1] = new TH1D(Form("fAverage_%s_%c",lRes[iRes],lAM[iMat]),Form(";%s;mean[%srec-%sgen] %s",lQuant[iRes],lVar[iRes],lVar[iRes],lUnit[iRes]),fHistProto->GetNbinsY(),0,fHistProto->GetYaxis()->GetBinLowEdge(fHistProto->GetNbinsY()+1));
        fMeasure[2] = new TH1D(Form("fStdDev_%s_%c",lRes[iRes],lAM[iMat]),Form(";%s;std dev[%srec-%sgen] %s",lQuant[iRes],lVar[iRes],lVar[iRes],lUnit[iRes]),fHistProto->GetNbinsY(),0,fHistProto->GetYaxis()->GetBinLowEdge(fHistProto->GetNbinsY()+1));
        for(int iBin=1;iBin<=fHistProto->GetNbinsY();iBin++){
          fMeasure[0]  = (TH1D*) fHistProto->ProjectionX(Form("fRes_%s_%c",lRes[iRes],lAM[iMat]),iBin,iBin);
          fMeasure[1]->SetBinContent(iBin,fMeasure[0]->GetMean());
          fMeasure[2]->SetBinContent(iBin,fMeasure[0]->GetStdDev());
          fMeasure[1]->SetBinError(iBin,fMeasure[0]->GetMeanError());
          fMeasure[2]->SetBinError(iBin,fMeasure[0]->GetStdDevError());
        }
      }
      if(iRes==4){
        output_file.cd();
        fMeasure[3]  = (TH1D*) fHistProto->ProjectionX(Form("fSlice_10_20_%c",lAM[iMat]),3,4);
        fMeasure[3]->Write();
      }

      fMeasure[3]  = (TH1D*) fHistProto->ProjectionX(Form("f%s_Int_%c",lRes[iRes],lAM[iMat]),1,fHistProto->GetNbinsY());
      
      subdir[0]->cd();
      fHistProto->Write();
      
      MultiCV.cd(1);
      fHistProto->Draw("colz");
      for(int iMeas=1;iMeas<4;iMeas++){
        fMeasure[iMeas]->SetMarkerStyle(8);
        subdir[iMeas]->cd();
        fMeasure[iMeas]->Write();
        if(iMeas!=3){
          MultiCV.cd(iMeas+1);
          fMeasure[iMeas]->Draw("PMC PLC");
        }
      }

      MultiCV.SetName(Form("MultiPlot%s_%c",lRes[iRes],lAM[iMat]));
      subdir[4]->cd();
      MultiCV.Write();
      
    } 
  }

}

/*
void compare_resolution(char * kf_name = "MassResolutionKF.root",char * std_name = "MassResolutionStd.root",char * output_name = "res_comparison.root")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptTitle(0);
  TFile std_file(std_name);
  TFile kf_file(kf_name);

  const char lAM[3]{"AM"};
  const char *lProj[2]{"pT", "ct"};
  const char *lProjSet[2]{"zx", "zy"};
  const char *lVarName[3]={ "NsigmaTPC","NclusterTPC","NclusterITS"};
  const char *lMeasure[3]={ "Res","Average","StdDev"};
  const char *lRes[3]={ "Mass","Pt","Ct"};
  const char *lUnit[3]={ "(GeV/#it{c}^{2})","(GeV/#it{c})","(cm)"};

  TFile output_file(output_name,"RECREATE");
  output_file.cd();

  TCanvas MultiCV("","",1200,400);
  TCanvas AverageCV("","",600,400);
  TCanvas StdDevCV("","",600,400);
  MultiCV.Divide(2,1);

  TH1D* fMeasure[2][2];

  TLegend *leg=new TLegend(0.2,0.2,0.5,0.5);

  //dm vs pt/ct -- dct vs ct -- dpt vs pt
  for(int iRes=0; iRes<3; iRes++){
    //matter or antimat
    for(int iMat=0; iMat<2; iMat++){
      //projection on pt or ct
      for(int iProj=0; iProj<2; iProj++){
        //variable of the cut
        if((iRes==1 && iProj==1)||(iRes==2&& iProj==0))//to do only dpt/ct vs pt/ct
          continue;
        
        for(int iMeas=1; iMeas<3; iMeas++){
          fMeasure[0][iMeas] = (TH1D*) kf_file.Get(Form("f%s_%s_%s_%c",lMeasure[iMeas],lRes[iRes],lProj[iProj],lAM[iMat]));
          fMeasure[1][iMeas] = (TH1D*) std_file.Get(Form("f%s_%s_%s_%c",lMeasure[iMeas],lRes[iRes],lProj[iProj],lAM[iMat]));
        
          fMeasure[0][iMeas]->SetName(Form("KF_%s_%s_%s_%c",lMeasure[iMeas],lRes[iRes],lProj[iProj],lAM[iMat]));
          fMeasure[1][iMeas]->SetName(Form("Std_%s_%s_%s_%c",lMeasure[iMeas],lRes[iRes],lProj[iProj],lAM[iMat]));
          fMeasure[0][iMeas]->SetTitle("Kalman filter");
          fMeasure[1][iMeas]->SetTitle("Standard vertexer");
        }
        

        for(int iMeas=1;iMeas<3;iMeas++){
          for(int iComp=0;iComp<2;iComp++){
            output_file.cd();
            fMeasure[iComp][iMeas]->SetMarkerStyle(8);
            fMeasure[iComp][iMeas]->Write();
            if(iMeas==1){
              AverageCV.cd();
              if(iComp==0)
                fMeasure[iComp][iMeas]->Draw("PMC PLC");
              else
                fMeasure[iComp][iMeas]->Draw("SAME PMC PLC");
            }
            else{
              StdDevCV.cd();
              if(iComp==0)
                fMeasure[iComp][iMeas]->Draw("PMC PLC");
              else
                fMeasure[iComp][iMeas]->Draw("SAME PMC PLC");
            }

            MultiCV.cd(iMeas);
            if(iComp==0){
              MultiCV.SetLeftMargin(30.);
              fMeasure[iComp][iMeas]->Draw("PMC PLC");
            }
            else
              fMeasure[iComp][iMeas]->Draw("SAME PMC PLC");
            if(iMeas==1)
              leg->AddEntry(fMeasure[iComp][1],fMeasure[iComp][1]->GetTitle(),"lp");
          }
          
        }


        AverageCV.SetName(Form("AverageComparison_%s_%s_%c",lRes[iRes],lProj[iProj],lAM[iMat])); 
        AverageCV.BuildLegend();
        AverageCV.Write();

        StdDevCV.SetName(Form("StdDevComparison_%s_%s_%c",lRes[iRes],lProj[iProj],lAM[iMat])); 
        StdDevCV.BuildLegend();
        StdDevCV.Write();

        leg->SetTextFont(132);
        MultiCV.cd(1);
        leg->Draw();
        MultiCV.SetName(Form("MultiComparison_%s_%s_%c",lRes[iRes],lProj[iProj],lAM[iMat])); 
        //MultiCV.SaveAs(Form("%s/%s/MultiRainbowPlot_%s_%s_%c.pdf",folder_name); 
        output_file.cd();
        MultiCV.Write();
        leg->Clear();

      }
    } 
  }

}
*/
/*
gSystem->Exec(Form("mkdir %s",folder_name));
for(auto VarName : lVarName){
  gSystem->Exec(Form("mkdir %s/%s/",folder_name,VarName));
  for(auto Measure : lMeasure)
    gSystem->Exec(Form("mkdir %s/%s/%s",folder_name,VarName,Measure));
}
*/
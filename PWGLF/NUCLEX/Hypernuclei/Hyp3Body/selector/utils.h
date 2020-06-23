#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

void set_style(TH1D* hist, int color_index){
  if(color_index>=8) color_index=0;
  color_index++;
  hist->SetMarkerStyle(8);
  hist->SetMarkerColor(color_index);
  hist->SetLineColor(color_index);
}

template <typename T,size_t size>
void set_range(const T (&histos)[size], float err_mult=5){
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  float max,min;
  bool first=true;
  for(auto histo: histos){
    for(int bin=1; bin<=histo->GetNbinsX(); bin++){
      if(first){
        max = histo->GetBinContent(bin)+histo->GetBinError(bin)*err_mult;
        min = histo->GetBinContent(bin)-histo->GetBinError(bin)*err_mult;
        first=false;
      }
      if(histo->GetBinContent(bin)+histo->GetBinError(bin)*err_mult > max)
        max = histo->GetBinContent(bin)+histo->GetBinError(bin)*err_mult;
      if(histo->GetBinContent(bin)-histo->GetBinError(bin)*err_mult < min)
        min = histo->GetBinContent(bin)-histo->GetBinError(bin)*err_mult;
    }
  }
  
  for(auto histo: histos){
    histo->GetYaxis()->SetRangeUser(min,max);
  }
}

template <typename T,size_t size>
void rainbow_plot(const T (&histos)[size], TCanvas &cv, const char* canvas_name, bool legend=true, const char* plot_mode=""){
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  cv.SetWindowSize(800,450);
  cv.SetName(canvas_name);
  cv.cd();
  bool first=true;
  //set_range(histos);
  for(auto histo: histos){
    if(first){
      first=false;
      histo->Draw(plot_mode);
    }
    else
      histo->Draw(Form("SAME %s",plot_mode));
  }
  if(legend)
    cv.BuildLegend();
}

template <typename T,size_t size_x,size_t size_y>
void multirainbow_plot(const T (&histos)[size_x][size_y], TCanvas &multi_cv, const char* canvas_name, bool legend=true, const char* plot_mode=""){
  gStyle->SetPadRightMargin(0.2); 
  gStyle->SetPadLeftMargin(0.15);
  multi_cv.SetWindowSize(800*size_x,450);  
  multi_cv.SetName(canvas_name);
  TLegend *leg=new TLegend(0.1,0.7,0.5,0.1);
  multi_cv.Divide(size_x,1);
  bool first=true;
  for(int x=0; x<size_x; x++){
    for(auto histo: histos[x]){
      multi_cv.cd(x+1);
      if(first){
        first=false;
        histo->Draw(plot_mode);
      }
      else
        histo->Draw(Form("SAME %s",plot_mode));
      if(x==0)
        leg->AddEntry(histo,histo->GetTitle(),"lp");
    }
  }

  if(legend){
    leg->SetTextFont(132);
    multi_cv.cd(1);
    leg->Draw();
  }
}
template <typename T,size_t size_x,size_t size_y>
void multirainbow_plotvert(const T (&histos)[size_x][size_y], TCanvas &multi_cv, const char* canvas_name, bool legend=true, const char* plot_mode=""){
  multi_cv.SetWindowSize(800,450*size_x);  
  multi_cv.SetName(canvas_name);
  TLegend *leg=new TLegend(0.2,0.2,0.5,0.5);
  multi_cv.Divide(1,size_x, 0.1, 0.1);
  bool first=true;
  for(int x=0; x<size_x; x++){
    for(auto histo: histos[x]){
      multi_cv.cd(x+1);
      if(first){
        first=false;
        histo->Draw(plot_mode);
      }
      else
        histo->Draw(Form("SAME %s",plot_mode));
      if(x==0)
        leg->AddEntry(histo,histo->GetTitle(),"lp");
    }
  }

  if(legend){
    leg->SetTextFont(132);
    multi_cv.cd(1);
    leg->Draw();
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

bool rejection(TF1* density_func,float variable, float max_value=1){
  float prob = density_func->Eval(variable);
  if(gRandom->Uniform(max_value)<prob)
    return true;
  else 
    return false;
}
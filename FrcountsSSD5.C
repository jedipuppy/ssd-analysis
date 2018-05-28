double GaussFunction(double *x, double *par){
  double xx= x[0];
  double N= par[0];
  double mu= par[1];
  double sigma= par[2];

  return N*exp(-(xx-mu)*(xx-mu)/(2.*sigma*sigma));
}

double GaussFunctionSum(double *x, double *par){
  return GaussFunction(x, &par[0])
    +GaussFunction(x, &par[3]);
}

double TimeFit(double *x, double *par){

  //SSD5はずっとINなので、輸送中ばしばしFrが届いている
  //よって残留Fr(違う条件で届いたFr)は０とみなす
  //また、測定開始時刻≠キャッチャーにFrが届いた時刻
  //よってTで時刻を補正
  
  double t= x[0];
  double n= par[0];//キャッチャーに届くFrのflux
  double lambda= par[1];//崩壊定数
  double Ia= 0.79;//α崩壊確率(210,211Frの平均)
  double Omega=1.67e-3;//立体角
  double T=par[2];//補正時間

  return n*Omega*(1-exp(-(lambda*(t+T))))*Ia;

    }

int FrcountsSSD5(char *filename, int xlower, int xupper){
  gROOT->SetStyle("Plain");
  TCanvas *energy_spectrum=new TCanvas("energy_spectrum","energy_spectrum",0,0,1000,450);
  TCanvas *time_spectrum=new TCanvas("time_spectrum","time_spectrum",0,550,1000,450);
  FILE *fin;

  TH1D *hist=new TH1D("hist",Form("%s",filename),4096,0.5,4096.5);
  hist->GetXaxis()->SetTitle("energy (channel)");
  hist->GetYaxis()->SetTitle("event number");

  TH1D *hist2=new TH1D("hist2",Form("%s",filename),900,0,900);//1sで1bin
  hist2->GetXaxis()->SetTitle("time (sec)");
  hist2->GetYaxis()->SetTitle("event number");

  char buf[512];
  int evt_id, evt_height;
  double evt_time;
  double t_min=0;
  int cycle=0;
  if((fin = fopen(Form("%s.csv",filename), "r"))==NULL){
    printf("file open error!!\n");
    exit(1);
  }
  
  while(fgets(buf,sizeof(buf),fin)!=NULL){
    if(cycle>=39){
      sscanf(buf, "%d,%lf,%d", &evt_id, &evt_time, &evt_height);
      if(evt_time/1E6>t_min){
      hist->Fill(evt_height);
      }

      if(evt_height>xlower && evt_height<xupper){
	hist2->Fill(evt_time/1E6);
      }
    }
    cycle++;
  }

  energy_spectrum->cd();
  gPad->SetLogy();
  hist->GetXaxis()->SetRangeUser(1000, 4000);
  hist->Draw();
  TLine *llower=new TLine(xlower, 0, xlower, 2e4);
  TLine *lupper=new TLine(xupper, 0, xupper, 2e4);
  llower->SetLineColor(kGray+3);
  lupper->SetLineColor(kGray+3);
  llower->Draw();
  lupper->Draw();


  //double par[9]={210Fr amp, mu, sigma, 209Fr amp, mu, sigma, flux, lambda, T};
  double par[9]={1e3, 2700, 10, 1e2, 2700, 10, 2e4, 0.0036, 1e12};

  //以下energy spectrumのフィッティング
  TF1 *energy_fit=new TF1("energy_fit",GaussFunctionSum, 2600, 2800, 6);
  energy_fit->SetLineColor(kRed);
  energy_fit->SetParameters(&par[0]);

  hist->Fit("energy_fit");
  energy_fit->Draw("same");

  energy_fit->GetParameters(&par[0]);

  TF1 *f_gauss=new TF1("f_gauss",GaussFunction, 2440, 2540, 3);//210,211Fr
  TF1 *g_gauss=new TF1("g_gauss",GaussFunction, 2440, 2540, 3);//208,209Fr
  
  f_gauss->SetParameters(&par[0]);
  g_gauss->SetParameters(&par[3]);

  f_gauss->SetLineColor(kBlue);
  g_gauss->SetLineColor(kMagenta-7);

  f_gauss->SetLineWidth(2);
  g_gauss->SetLineWidth(1.5);

  f_gauss->Draw("same");
  g_gauss->Draw("same");

  double gross=sqrt(2*TMath::Pi())*par[0]*par[2];

  double tt=900-t_min;//測定時間
  double Sigma=0.00167;//SSD5立体角
  double Ia_ave=0.79; //210,211Frの平均
  double eflux=gross/tt/Sigma/Ia_ave;
  
  //TText *txt_gross=new TText(3200,1e2, Form("GROSS:%f",gross));
  TText *txt_eflux=new TText(10,10,Form("Flux:%f [pps]",eflux));
	cout << eflux <<"\n";
  //txt_gross->Draw("same");
  txt_eflux->Draw("same");

  //以下time spectrumのフィッティング
  time_spectrum->cd();
  hist2->Draw();
  
  TF1 *time_fit=new TF1("time_fit",TimeFit,0,900,3);
  time_fit->SetLineColor(kRed);
  time_fit->SetParameters(&par[6]);
  time_fit->SetParLimits(0, 0, 10e9);
  time_fit->SetParLimits(1, 0, 1);
  time_fit->SetParLimits(2, 0, 1e5);
	time_fit.SetParameter(0,1e4);
	time_fit.SetParameter(1,0.04);
	time_fit.SetParameter(2,0);
	
	
  hist2->Fit("time_fit");
  time_fit->Draw("same");
  
  TText *txt_ROI=new TText(10,10,Form("ROI:%d ~ %d [ch]", xlower, xupper));
  txt_ROI->Draw("same");
  
  time_fit->GetParameters(&par[6]);
  double tflux=par[6];
  TText *txt_tflux=new TText(20,20,Form("Flux:%f [pps]",tflux));
  txt_tflux->Draw("same");
    

  return 0;
  }

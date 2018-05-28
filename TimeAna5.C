/*
 * usage: 
 *    $ root -l 'TimeAna.C("run25-20170210-XXXXXX-SlotX-InX", <xlower>, <xupper>, <Ialpha>, <Omega>, <flux>)'
 *
 * xlower and xupper define the region of interest for energy.
 */

double Accumulation(double *x, double *par){
  double Ialpha_210  =par[0];
  double Omega   =par[1];
  double flux    =par[2];
  double lambda  =par[3];
  double residual=par[4];
  double t1      =par[5];
  double Ialpha_211 =par[6];
  //
  double t = x[0] - t1;
  //
  double Fr_new;

  if(t<0){
    Fr_new = 0;
  }
  else{
    Fr_new = Ialpha_210*Omega * ( 4.1/(4.1+1.2)*flux * (1-exp(-lambda*t)))
    +Ialpha_211*Omega * ( 1.2/(4.1+1.2)*flux * (1-exp(-lambda*t)));
  }

  return Fr_new
			  + Ialpha_210*Omega *lambda * 4.1/(4.1+1.2)*residual * exp(-lambda*x[0])
        + Ialpha_211*Omega *lambda * 1.2/(4.1+1.2)*residual * exp(-lambda*x[0]);
}


int TimeAna5(string strname, double xlower, double xupper, double Ialpha_210, double Ialpha_211, double Omega, double flux, double t1, double t2){
  const char *filename=strname.c_str();
  gStyle->SetOptStat("nmrei");
  TCanvas *energy_spectrum=new TCanvas("energy_spectrum","energy_spectrum",0,0,500,700);
  energy_spectrum->Divide(1,2);
  FILE *fin;

  TH1D *hist=new TH1D("hist",Form("%s",filename), 4096, 0.5, 4096.5);
  hist->GetXaxis()->SetTitle("energy (channel)");
  hist->GetYaxis()->SetTitle("event number");
  //  TH1D *htime=new TH1D("htime", "time spectrum", t2-t1, t1, t2);
  TH1D *htime=new TH1D("htime", "time spectrum", 2048, 0, 2048);
  htime->GetXaxis()->SetTitle("time (sec)");
  htime->GetYaxis()->SetTitle("event number");

  char buf[512];
  int evt_id, evt_height;
  double evt_time;
  int cycle=0;
  if((fin = fopen(Form("%s.csv",filename), "r"))==NULL){
    printf("file open error!!\n");
    exit(1);
  }
  while(fgets(buf,sizeof(buf),fin)!=NULL){
    if(cycle>=39){
      sscanf(buf, "%d,%lf,%d", &evt_id, &evt_time, &evt_height);
      hist->Fill(evt_height);
      if(evt_height>xlower && evt_height<xupper){

	htime->Fill(evt_time/1E6);
      }
    }
    cycle++;
  }

  energy_spectrum->cd(1);
  energy_spectrum->cd(1)->SetLogy(1);
  hist->Draw();
  TLine *llower=new TLine(xlower, 0, xlower, 1e6);
  llower->SetLineColor(kGray);
  llower->Draw();
  TLine *lupper=new TLine(xupper, 0, xupper, 1e6);
  lupper->SetLineColor(kGray);
  lupper->Draw();

  energy_spectrum->cd(2);
  htime->Draw();
  double par[7];
  TF1 *faccum=new TF1("faccum", Accumulation, t1, t2, 7);
  faccum->FixParameter(0, Ialpha_210); faccum->SetParName(0, "Ialpha_210");
  faccum->FixParameter(1, Omega);  faccum->SetParName(1, "Omega");
  faccum->SetParameter(2, flux);   faccum->SetParName(2, "flux");
  faccum->SetParameter(3, 3e-3);   faccum->SetParName(3, "lambda");
  faccum->SetParameter(4, 0);      faccum->SetParName(4, "residual");
  faccum->SetParLimits(4, 0, 1e9);
  faccum->SetParameter(5, t1);      faccum->SetParName(5, "t1");
  faccum->FixParameter(6, Ialpha_211); faccum->SetParName(6, "Ialpha_211");
  htime->Fit("faccum", "R");

  faccum->GetParameters(par);
  double flux_err = faccum->GetParError(2);
  double lambda_err = faccum->GetParError(3);
  double t1_err = faccum->GetParError(5);
  TText *tflux=new TText(300, htime->GetBinContent(300)*0.8, Form("flux = %.2le", par[2]));
  tflux->Draw();

  energy_spectrum -> Print( Form("%s-%.0lf-TimeAna5.png", filename, t1) );
  energy_spectrum -> Print( Form("%s-%.0lf-TimeAna5.C", filename, t1) );


  FILE *fp ;
  if ((fp = fopen("flux.txt", "a")) == NULL){
    fprintf(stderr, "file open error\n");
    return 0;
  }

  fprintf(fp, "%s, %lf, %lf, %lf, %lf,%lf,%lf, %lf\n", filename, par[2], flux_err, 4.1/(4.1+1.2)*par[2], par[3], lambda_err, par[5], t1_err);
  //fprintf(fp,"test\n");

  fclose(fp);



  energy_spectrum->cd();


  return 0;
}

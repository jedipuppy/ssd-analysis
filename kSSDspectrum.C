int kSSDspectrum(char *filename, int xlower, int xupper){
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("nmreuio");
  TCanvas *energy_spectrum=new TCanvas("energy_spectrum","energy_spectrum",0,0,500,700);
  energy_spectrum->Divide(1,2);
  //  TCanvas *time_spectrum=new TCanvas("time_spectrum", "time_spectrum");
  FILE *fin;

  TH1D *hist=new TH1D("hist",Form("%s",filename), 4096, 0.5, 4096.5);
  hist->GetXaxis()->SetTitle("energy (channel)");
  hist->GetYaxis()->SetTitle("event number");
  TH1D *hist2=new TH1D("hist2", "time spectrum", 4096, 0, 900);
  hist2->GetXaxis()->SetTitle("time (sec)");
  hist2->GetYaxis()->SetTitle("event number");

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
	hist2->Fill(evt_time/1E6);
      }
    }
    cycle++;
  }

  TText *txt=new TText(0, 0, Form("ROI: %d ~ %d GROSS: %.0lf",xlower, xupper, hist->Integral(xlower, xupper)));
  hist->GetYaxis()->SetTitleOffset(1.5);
  energy_spectrum->cd(1);
  hist->Draw();
  txt->Draw();
  //  time_spectrum->cd();
  energy_spectrum->cd(2);
  hist2->Draw();

  return 0;
}

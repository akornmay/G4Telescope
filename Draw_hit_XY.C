void Draw_hit_XY(){

  //FILE *ff = fopen("p280GeV_all.dat","r");
  //FILE *ff = fopen("e-3GeV_all.dat","r");
  //FILE *ff = fopen("e-3GeV_all_noG10.dat","r");
  //FILE *ff = fopen("e-3GeV_all_noG10_fixorder.dat","r");
  //FILE *ff = fopen("e-3GeV_all_fixorder.dat","r");
  FILE *ff = fopen("pi+350MeV.dat","r");

  int id;
  float e, x, y, z;

  TH1F *h_e = new TH1F("h_e","energy", 100, 0., 300);
  TH1F *h_x[8];
  TH1F *h_y[8];

  char hname[20];
  for(int ii=0; ii<8; ii++){
    sprintf(hname,"h_x_%d",ii);
    h_x[ii] = new TH1F(hname, hname, 200, -4., 4.);
    sprintf(hname,"h_y_%d",ii);
    h_y[ii] = new TH1F(hname, hname, 200, -4., 4.);
  }
  
  for(int line=0; line<9400; line++) {
    //for(int line=0; line<10; line++) {
    
    //fscanf(ff, "%f %f %f %f %f", &id, &e, &x, &y, &z);
    fscanf(ff, "%i %f %f %f %f", &id, &e, &x, &y, &z);    
    int ii = -1;
    if(z<2.) ii=0;
    else if (z<4.) ii=1;
    else if (z<5.) ii=2;
    else if (z<7.) ii=3;
    else if (z<8.5) ii=4;
    else if (z<10.) ii=5;
    else if (z<12.) ii=6;
    else if (z<13.) ii=7;
    //printf("id %d, e %f, x %f, y %f, z %f \n", ii, e, x, y, z);
    
    h_e->Fill(e);
    h_x[ii]->Fill(x*10.);
    h_y[ii]->Fill(y*10.);

  }


  TFile *fout = new TFile("Hit_out.root","recreate");
  h_e->Write();
  for(int ii=0; ii<8; ii++){
    h_x[ii]->Write();
    h_y[ii]->Write();
  }
  fout->Close();

  TCanvas *c1 = new TCanvas("c1","",1000,500);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogy(1);
  h_x[7]->SetXTitle("mm");
  h_x[7]->Fit("gaus","","",-1., 1.);
  h_x[7]->Draw();
  c1->cd(2);
  gPad->SetLogy(1);
  h_y[7]->SetXTitle("mm");
  h_y[7]->Fit("gaus","","",-1., 1.);
  h_y[7]->Draw();


}


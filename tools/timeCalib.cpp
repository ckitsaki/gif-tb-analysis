std::map<int, float> map_strip_ti;
std::map<int, float> map_strip_tf;
TH1F* h_chi2_ndf = new TH1F("h_chi2_ndf", "chi2/ndf", 100, 0, 20);
bool calibSM1=true;

void getHistogram(std::string filein, std::vector<TH2F*>* v_histos_tdo, std::vector<TH2F*>* v_histos_relbcid, std::vector<TH2F*>* v_histos_pdo)
{
  TFile* fin = new TFile(filein.c_str());

  if(!fin->IsOpen()) 
  {
    printf("<E> Cannot open input file %s\n", filein.c_str());
    exit(1);
  }
  
  if(calibSM1)
    fin->cd("SM1/strips");
  else fin->cd("SB/strips");

  TDirectory *current_sourcedir = gDirectory;
  TList* list = current_sourcedir->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
 
  for(int i=0; i<16; i++)
  {
    TObject* f = list->At(i);
    TH2F* h = (TH2F*)gDirectory->Get(f->GetName());
    v_histos_tdo->push_back(h);
    i+=3;
    if(i>16) break;
  }

  for(int i=1; i<16; i++)
  {
    TObject* f = list->At(i);
    TH2F* h = (TH2F*)gDirectory->Get(f->GetName());
    v_histos_relbcid->push_back(h);
    i+=3;
    if(i>16) break;
  }

  for(int i=3; i<16; i++)
  {
    TObject* f = list->At(i);
    TH2F* h = (TH2F*)gDirectory->Get(f->GetName());
    v_histos_pdo->push_back(h);
    i+=3;
    if(i>16) break;
  }
  
}

void fermidiracfit(TH1F* h, std::vector<float>* v_tmin, std::vector<float>* v_tmax, std::vector<float>* v_tmin_err, std::vector<float>* v_tmax_err)
{
  TF1 *func = new TF1("func", "[0]/( 1 + TMath::Exp( ([1]-x)/[2]) + TMath::Exp((x-[3])/[4]) ) ",0,260); 
  func->SetParameters(h->GetMaximum()/4, h->GetBinLowEdge(h->FindFirstBinAbove()), 0.5, h->GetBinLowEdge(h->FindLastBinAbove()), 0.5);
  func->SetParLimits(0, 0, h->GetMaximum());
  func->SetParLimits(1, h->GetBinLowEdge(h->FindFirstBinAbove())-5, h->GetBinLowEdge(h->FindFirstBinAbove())+10);
  func->SetParLimits(2, 0, 10);
  func->SetParLimits(3, h->GetBinLowEdge(h->FindFirstBinAbove()), h->GetBinLowEdge(h->FindLastBinAbove())+5);
  func->SetParLimits(4, 0, 10); 
  h->Fit("func","Q");
  h->Fit("func","Q");

  v_tmin->push_back(func->GetParameter(1));
  v_tmax->push_back(func->GetParameter(3));
  v_tmin_err->push_back(func->GetParError(1));
  v_tmax_err->push_back(func->GetParError(3));

  h_chi2_ndf->Fill(func->GetChisquare()/func->GetNDF());

}

int calib(std::string runNumber="", int layer=0)
{
  gROOT->SetBatch(kTRUE);
  std::vector<TH2F*> v_histos_tdo;
  std::vector<TH2F*> v_histos_relbcid;
  std::vector<TH2F*> v_histos_pdo;

  map_strip_ti = map<int,float>();
  map_strip_tf = map<int,float>();
  std::string histo_name = "./histos/"+runNumber+"_histograms.root";
  
  getHistogram(histo_name, &v_histos_tdo, &v_histos_relbcid, &v_histos_pdo);
  
  const int tot_strips = 8192;
  const int tot_layers = 4;
  
  std::vector<TH1F*> v_strip_tdo[tot_layers];
  std::vector<TH1F*> v_relbcid_strip[tot_layers];
  std::vector<TH1F*> v_strip_pdo[tot_layers];
  
  // do the projections over all strips
 // for(int ilayer=0; ilayer<tot_layers; ilayer++){
    TH2F* h_tdo_lay = v_histos_tdo.at(layer);
    TH2F* h_relbcid_lay = v_histos_relbcid.at(layer);
    TH2F* h_pdo_lay = v_histos_pdo.at(layer);
    
    for(int istrip=0; istrip<tot_strips; istrip++)
    {
      TH1F* h_proj_relbcid = (TH1F*)h_relbcid_lay->ProjectionY(Form("relbcid_%i",istrip), istrip+1, istrip+1);
      TH1F* h_proj_tdo= (TH1F*)h_tdo_lay->ProjectionY(Form("tdo_%i",istrip), istrip+1, istrip+1);
      TH1F* h_proj_pdo= (TH1F*)h_pdo_lay->ProjectionY(Form("pdo_%i",istrip), istrip+1, istrip+1);
      h_proj_tdo->Rebin(4);
      v_strip_tdo[layer].push_back(h_proj_tdo);
      v_relbcid_strip[layer].push_back(h_proj_relbcid);
      v_strip_pdo[layer].push_back(h_proj_pdo);
    }
 // } //end for loop on layers

 // std::cout<<v_strip_tdo[layer].size()<<" "<<v_strip_pdo[layer].size()<<std::endl;
  std::vector<float> v_tmin;
  std::vector<float> v_tmax;
  std::vector<float> v_tmin_err;
  std::vector<float> v_tmax_err;
  std::vector<float> v_pdos;
  std::vector<float> v_relbcid;
  std::vector<float> v_meantdo;
  int counter1=0;

  TH1F* h_time = new TH1F("h_time", "time", 100, 0, 200);
  
  for(int i=0; i<v_strip_tdo[layer].size(); i++)
  {
    if(v_strip_tdo[layer].at(i)->GetEntries()==0)
    {
      counter1++;
      continue;
    }
    else{
    fermidiracfit(v_strip_tdo[layer].at(i), &v_tmin, &v_tmax, &v_tmin_err, &v_tmax_err);
    v_pdos.push_back(v_strip_pdo[layer].at(i)->GetMean());
    v_relbcid.push_back(v_relbcid_strip[layer].at(i)->GetMean());
    v_meantdo.push_back(v_strip_tdo[layer].at(i)->GetMean());
  }
    //if(v_chi2.at(i)/v_ndf.at(i) < 5) v_goodfits.push_back(v_strip_tdo[0].at(i));
   // std::cout<<chi2<<std::endl;
  }
 
  for(int itime=0; itime<v_pdos.size(); itime++)
  {
    float tmp_time = v_relbcid.at(itime)*25-(v_meantdo.at(itime)-v_tmin.at(itime))/((v_tmax.at(itime)-v_tmin.at(itime))/25);
    h_time->Fill(tmp_time);
  }   

  TCanvas *c_time = new TCanvas();
  c_time->cd();
  h_time->Draw();
  c_time->SaveAs("./output/time.pdf"); 

  int iinit = 1200;
  int ifinal = 1210;
  if(calibSM1) {
    iinit = 2200;
    ifinal = 2210;
  }
    TCanvas* kanv = new TCanvas();
    kanv->Divide(5,2);
    gStyle->SetOptFit(111111);
    for(int i=iinit; i<ifinal; i++)
    {
      //if(v_strip_tdo[0].at(i)->GetEntries()==0) continue;
      kanv->cd(i+1-iinit); 
      v_strip_tdo[layer].at(i)->Draw();
    }
    kanv->SaveAs("./output/tdos.pdf");

    TCanvas* kanv2 = new TCanvas();
    kanv2->Divide(5,2);
    for(int i=iinit; i<ifinal; i++)
    {
      //if(v_strip_tdo[0].at(i)->GetEntries()==0) continue;
      kanv2->cd(i+1-iinit); 
      v_relbcid_strip[layer].at(i)->Draw();
    }
    kanv2->SaveAs("./output/relbcid.pdf");

    TCanvas* kanv3 = new TCanvas();
    kanv3->Divide(5,2);
    for(int i=iinit; i<ifinal; i++)
    {
      //if(v_strip_tdo[0].at(i)->GetEntries()==0) continue;
      kanv3->cd(i+1-iinit); 
      v_strip_pdo[layer].at(i)->Draw();
    }
    kanv3->SaveAs("./output/pdos.pdf");
    
    TH1F* h_tmin = new TH1F("h_ti","ti",64,0,256);
    TH1F* h_tmax = new TH1F("h_tf","tf",64,0,256);

    TH1F* h_tmin_err = new TH1F("h_tmin_err","tmin error",50,0,100);
    TH1F* h_tmax_err = new TH1F("h_tmax_err","tmax error",50,0,100);

       
    TH1F* h_Dt = new TH1F("h_slope", "slope", 80, 0, 10);
    TH2F* h_ti_channel = new TH2F("h_ti_channel", "ti vs channel", 5120, 0, 5120, 50, 50, 100);
    TH2F* h_tf_channel = new TH2F("h_tf_channel", "tf vs channel", 5120, 0, 5120, 100, 100, 200);
    TH2F* h_Dt_channel = new TH2F("h_slope_channel", "slope vs channel", 5120, 0, 5120, 50, 2.5, 5);
    TH2F* h_meanrelbcid_channel = new TH2F("h_relbcid_channel","<relbcid> vs channel", 5120, 0, 5120, 8, -0.5, 7.5);

    TH2F* h_ti_pdo = new TH2F("h_ti_pdo", "ti vs <pdo>", 1024, 0, 1024, 50, 50, 100);
    TH2F* h_tf_pdo = new TH2F("h_tf_pdo", "tf vs <pdo>", 1024, 0, 1024, 50, 150, 200);
    TH2F* h_Dt_pdo = new TH2F("h_slope_pdo", "slope vs <pdo>", 1024, 0, 1024, 50, 2.5, 5);

    TH2F* h_meanpdo_channel = new TH2F("h_meanpdo_channel", "<pdo> vs channel", 5120, 0, 5120, 1024, 0, 1024);
    TH2F* h_meanrelbcid_meanpdo = new TH2F("h_relbcid_pdo","<relbcid> vs <pdo>", 1024, 0, 1024, 8, -0.5, 7.5);
    int count_reject=0;

    int index=0;
    for(int i=0; i<v_strip_tdo[layer].size(); i++)
    {
      if(v_strip_tdo[layer].at(i)->GetEntries()==0) continue;

      h_tmin->Fill(v_tmin.at(index));
      h_tmax->Fill(v_tmax.at(index));

      h_Dt->Fill((v_tmax.at(index)-v_tmin.at(index))/25.);

      std::string name = v_strip_tdo[layer].at(i)->GetName();
      std::string channel = "";
      if(i<99) channel = name.substr(4,2);
      if(i>99 && i<1000) channel = name.substr(4,3);
      else channel = name.substr(4,4);

      h_ti_channel->Fill(std::stoi( channel ), v_tmin.at(index));
      h_tf_channel->Fill(std::stoi( channel ), v_tmax.at(index));
      h_Dt_channel->Fill(std::stoi( channel ), (v_tmax.at(index)-v_tmin.at(index))/25.);
      h_meanpdo_channel->Fill(std::stoi( channel ), v_pdos.at(index));
      h_meanrelbcid_channel->Fill(std::stoi( channel ), v_relbcid.at(index));

      h_meanrelbcid_meanpdo->Fill(v_pdos.at(index), v_relbcid.at(index));
     
      h_ti_pdo->Fill(v_pdos.at(index), v_tmin.at(index));
      h_tf_pdo->Fill(v_pdos.at(index), v_tmax.at(index));
      h_Dt_pdo->Fill(v_pdos.at(index), (v_tmax.at(index)-v_tmin.at(index))/25.);

      map_strip_ti[std::stoi( channel )] = v_tmin.at(index);
      map_strip_tf[std::stoi( channel )] = v_tmax.at(index);
      index++;
    }

    TCanvas* c10 = new TCanvas();
    c10->cd();
    h_meanrelbcid_channel->SetStats(0);
    h_meanrelbcid_channel->Draw();
    h_meanrelbcid_channel->GetXaxis()->SetTitle("channel");
    h_meanrelbcid_channel->GetYaxis()->SetTitle("<relbcid>");
    for(int ivmm=0; ivmm<80; ivmm++)
    {
      TLine *line2 = new TLine(64*ivmm, -0.5, 64*ivmm, 7.5);
      line2->SetLineWidth(1);
      line2->SetLineStyle(3);
      line2->SetLineColor(kBlue);
      line2->Draw("same");
    }
    for(int immfe=0;immfe<10; immfe++){
    TLine* line = new TLine(512*immfe, -0.5, 512*immfe, 7.5);
    line->SetLineWidth(2.);
    line->SetLineStyle(1);
    line->SetLineColor(kRed);
    line->Draw("same");
    }
    h_meanrelbcid_channel->Draw("same");

    c10->SaveAs("./output/meanrelbcid_vs_channel.pdf");

    TCanvas* c_tmin = new TCanvas();
    TCanvas* c_tmax = new TCanvas();
    c_tmin->cd();
    h_tmin->Draw();
    h_tmin->GetXaxis()->SetTitle("ti [counts]");
    c_tmax->cd();
    h_tmax->Draw();
    h_tmax->GetXaxis()->SetTitle("tf [counts]");
    c_tmin->SaveAs("./output/ti.pdf");
    c_tmax->SaveAs("./output/tf.pdf");

    TCanvas* c_slope = new TCanvas();
    c_slope->cd();
    h_Dt->Draw();
    h_Dt->GetXaxis()->SetTitle("slope [counts/ns]");
    c_slope->SaveAs("./output/slope.pdf");

    TCanvas* c1 = new TCanvas();
    c1->cd();
    h_ti_channel->SetStats(0);
    h_ti_channel->Draw();
    for(int ivmm=0; ivmm<80; ivmm++)
    {
      TLine *line2 = new TLine(64*ivmm, 50, 64*ivmm, 100);
      line2->SetLineWidth(1);
      line2->SetLineStyle(3);
      line2->SetLineColor(kBlue);
      line2->Draw("same");
    }
    for(int immfe=0;immfe<10; immfe++){
    TLine* line = new TLine(512*immfe, 50, 512*immfe, 100);
    line->SetLineWidth(2.);
    line->SetLineStyle(1);
    line->SetLineColor(kRed);
    line->Draw("same");
    }
    h_ti_channel->Draw("same");

    TCanvas* c2 = new TCanvas();
    c2->cd();
    h_tf_channel->SetStats(0);
    h_tf_channel->Draw();
    for(int ivmm=0; ivmm<80; ivmm++)
    {
      TLine *line2 = new TLine(64*ivmm, 100, 64*ivmm, 200);
      line2->SetLineWidth(1);
      line2->SetLineStyle(3);
      line2->SetLineColor(kBlue);
      line2->Draw("same");
    }
    for(int immfe=0;immfe<10; immfe++){
    TLine* line = new TLine(512*immfe, 100, 512*immfe, 200);
    line->SetLineWidth(2.);
    line->SetLineStyle(1);
    line->SetLineColor(kRed);
    line->Draw("same");
    }
    h_tf_channel->Draw("same");

    TCanvas* c3 = new TCanvas();
    c3->cd();
    h_Dt_channel->SetStats(0);
    h_Dt_channel->Draw();
    for(int ivmm=0; ivmm<80; ivmm++)
    {
      TLine *line2 = new TLine(64*ivmm, 2.5, 64*ivmm, 5);
      line2->SetLineWidth(1);
      line2->SetLineStyle(3);
      line2->SetLineColor(kBlue);
      line2->Draw("same");
    }
    for(int immfe=0;immfe<10; immfe++){
    TLine* line = new TLine(512*immfe, 2.5, 512*immfe, 5);
    line->SetLineWidth(2.);
    line->SetLineStyle(1);
    line->SetLineColor(kRed);
    line->Draw("same");
    }
    h_Dt_channel->Draw("same");

    TCanvas* c4 = new TCanvas();
    c4->cd();
    h_meanpdo_channel->SetStats(0);
    h_meanpdo_channel->Draw();
    for(int ivmm=0; ivmm<80; ivmm++)
    {
      TLine *line2 = new TLine(64*ivmm, 0, 64*ivmm, 1024);
      line2->SetLineWidth(1);
      line2->SetLineStyle(3);
      line2->SetLineColor(kBlue);
      line2->Draw("same");
    }
    for(int immfe=0;immfe<10; immfe++){
    TLine* line = new TLine(512*immfe, 0, 512*immfe, 1024);
    line->SetLineWidth(2.);
    line->SetLineStyle(1);
    line->SetLineColor(kRed);
    line->Draw("same");
    }
    h_meanpdo_channel->Draw("same");

    TCanvas* c5 = new TCanvas();
    c5->cd();
    h_ti_pdo->SetStats(0);
    h_ti_pdo->Draw("colz");

    TCanvas* c6 = new TCanvas();
    c6->cd();
    h_tf_pdo->SetStats(0);
    h_tf_pdo->Draw("colz");

    TCanvas* c7 = new TCanvas();
    c7->cd();
    h_Dt_pdo->SetStats(0);
    h_Dt_pdo->Draw("colz");

    TCanvas* c8 = new TCanvas();
    c8->cd();
    h_meanrelbcid_meanpdo->Draw();
    h_meanrelbcid_meanpdo->GetXaxis()->SetTitle("<PDO> [counts]");
    h_meanrelbcid_meanpdo->GetYaxis()->SetTitle("<relbcid> [counts]");

    h_ti_channel->GetXaxis()->SetTitle("channel");
    h_ti_channel->GetYaxis()->SetTitle("ti [counts]");
    h_tf_channel->GetXaxis()->SetTitle("channel");
    h_tf_channel->GetYaxis()->SetTitle("tf [counts]");
    h_Dt_channel->GetXaxis()->SetTitle("channel");
    h_Dt_channel->GetYaxis()->SetTitle("slope [counts/ns]");
    h_meanpdo_channel->GetXaxis()->SetTitle("channel");
    h_meanpdo_channel->GetYaxis()->SetTitle("PDO [counts]");
    h_ti_pdo->GetXaxis()->SetTitle("PDO [counts]");
    h_ti_pdo->GetYaxis()->SetTitle("ti [counts]");
    h_tf_pdo->GetXaxis()->SetTitle("PDO [counts]");
    h_tf_pdo->GetYaxis()->SetTitle("tf [counts]");
    h_Dt_pdo->GetXaxis()->SetTitle("PDO [counts]");
    h_Dt_pdo->GetYaxis()->SetTitle("slope [counts/ns]");

    c1->SaveAs("./output/ti_vs_channel.pdf");
    c2->SaveAs("./output/tf_vs_channel.pdf");
    c3->SaveAs("./output/slope_vs_channel.pdf");
    c4->SaveAs("./output/meanPDO_vs_channel.pdf");
    c5->SaveAs("./output/ti_vs_meanpdo.pdf");
    c6->SaveAs("./output/tf_vs_meanpdo.pdf");
    c7->SaveAs("./output/slope_vs_meanpdo.pdf");
    c8->SaveAs("./output/meanrelbcid_meanpdo.pdf");

    std::string add = "_SB";
    if(calibSM1) add="_SM1";
    std::string fname = "./calib/" + runNumber + "_calib_" + std::to_string(layer) + add +".root";
    TFile *out_file = new TFile(fname.c_str(),"RECREATE"); 
    out_file->cd();
    h_Dt->Write();
    h_ti_channel->Write();
    h_tf_channel->Write();
    h_Dt_channel->Write();
    h_meanpdo_channel->Write();
    h_ti_pdo->Write();
    h_tf_pdo->Write();
    h_Dt_pdo->Write();
    h_tmax->Write();
    h_tmin->Write();
    h_meanrelbcid_channel->Write();
    h_chi2_ndf->Write();

    c1->Write(); 
    c2->Write(); 
    c3->Write(); 
    c4->Write();     
    c5->Write(); 
    c6->Write(); 
    c7->Write(); 
    c8->Write(); 

    std::cout<<"TIME: "<< h_time->GetMean()<<" ns "<<counter1<<" layer "<<layer<<std::endl;

    out_file->Close();

  return 0;
}


//#define cluster_multiplicity_plot
//#define sPdo_plot
//#define cPdo_plot
//#define n_clusters
//#define SM1
//#define SB
//#define hit_occup
//#define hv_scan
#define source_scan
#define correlations
std::string date = "14_07_v0";
std::string runNumbers[7] = {"1657871992", "1657838505", "1657833031", "1657827799", "1657822252", "1657876123", "1658081138"};//{"1658081138", "1657876123", "1657822252", "1657827799", "1657827799", "1657838505", "1657871992"};//{"1657871992", "1657838505", "1657827799", "1657827799", "1657822252", "1657876123", "1657881309", "1657792515"};//{"1657872953","1657837396","1657834337","1657826437","1657824840","1657876754","1657888357"};//{"1657795207", "1657792515", "1657793941","1657807186", "1657812543", "1657813282","1657808871", "1657811493", "1657813950","1657809945", "1657810753", "1657815073"};//{"1654249364", "1654248067", "1654246091", "1654203741", "1654201895", "1654186091", "1654250774", "1654252014"};

std::string runNumbers_SC9[7] = {"1657872953","1657837396","1657834337","1657826437","1657824840","1657876754","1657888357"};
std::string runNumbers_SC15[7] = {"1657916610", "1657916051", "1657910889", "1657896516", "1657892624", "1657892040", "1657889042"};
std::string runNumbers_SC12[7] = {"1657928611", "1657926198", "1657925468", "1657918619", "1657919132", "1657921590", "1657922079"};

std::string filename = "output/HV_output_" + date + ".root";
TFile* output_file = new TFile(filename.c_str(), "RECREATE");

std::string runNumbers_9rms[3] 	= {"1657795207", "1657792515", "1657793941"};
std::string runNumbers_12rms[3] = {"1657807186", "1657812543", "1657813282"};
std::string runNumbers_15rms[3] = {"1657808871", "1657811493", "1657813950"};
std::string runNumbers_18rms[3] = {"1657809945", "1657810753", "1657815073"};

TH1F* getHistogram(const char* filein="1654168418_histograms.root", int layerIndex=0)
{
  TFile* fin = new TFile(filein);
 
  if (!fin->IsOpen()) {
    printf("<E> Cannot open input file %s\n",filein) ;
    exit(1) ;
  }

#if (defined(SM1) && defined(hit_occup) )
  fin->cd("SM1/raw_hits");
#endif

#if (defined(SM1) && (defined(sPdo_plot) || defined(cPdo_plot) || defined(n_clusters) || defined(cluster_multiplicity_plot) ))
  fin->cd("SM1/clusters");
#endif

#if (defined(SB) && (defined(sPdo_plot) || defined(cPdo_plot) || defined(n_clusters) || defined(cluster_multiplicity_plot) ))
  fin->cd("SB/clusters");
#endif

#if (defined(SB) && defined(hit_occup) )
  fin->cd("SB/raw_hits");
#endif

#if (defined(correlations))
  fin->cd("correlations");
#endif

  TDirectory *current_sourcedir = gDirectory;

  TList* list = current_sourcedir->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TObject* f = list->At(layerIndex);
  TH1F* h = (TH1F*)gDirectory->Get(f->GetName());

  return h;
}

float getMPV(std::string filein="histograms1.root", int layerIndex=0, std::string wePlot="")
{
  TCanvas * c = new TCanvas();
  TH1F* histo = getHistogram(filein.c_str(), layerIndex);
  std::string grName = wePlot + "_lay" + to_string(layerIndex);
  histo->SetName(grName.c_str());

  gStyle->SetOptStat(111111);

//#ifdef cluster_multiplicity_plot
 // histo->SetBinContent(2,0);
//#endif
  c->SetLogy();
  histo->Draw("hist");

  output_file->cd();
  c->Write();
  std::cout<<histo->GetMean()<<std::endl;
  return histo->GetMean(); 
}

float getOverflow(const char* filein="histograms1.root", int layerIndex=0)
{
  TH1F* histo = getHistogram(filein, layerIndex);
  float overflow = histo->GetBinContent(1023);

  float result = overflow/histo->Integral();
 // std::cout<<result<<std::endl;
  return result;
}

EColor setColors(int layerIndex)
{
	if(layerIndex==0) return kRed;
	if(layerIndex==1) return kBlue;
	if(layerIndex==2) return kGreen;
	if(layerIndex==3) return kYellow;

	return kOrange;
}

EMarkerStyle setMarkers(int layerIndex)
{
	if(layerIndex==0) return kFullCircle;
	if(layerIndex==1) return kFullSquare;
	if(layerIndex==2) return kFullTriangleUp;
	if(layerIndex==3) return kFullTriangleDown;

	return kFullStar;
}

EMarkerStyle setMarkers(int layerIndex, bool open)
{
	if(!open) {
	if(layerIndex==0) return kFullCircle;
	if(layerIndex==1) return kFullSquare;
	if(layerIndex==2) return kFullTriangleUp;
	if(layerIndex==3) return kFullTriangleDown;
	}
	else {
	if(layerIndex==0) return kOpenCircle;
	if(layerIndex==1) return kOpenSquare;
	if(layerIndex==2) return kOpenTriangleUp;
	if(layerIndex==3) return kOpenTriangleDown;
	}

	return kFullStar;
}

EColor setColors(int rms, bool same)
{
  if(rms==9) return kRed;
  else if(rms==12) return kAzure;
  else if(rms==15) return kGreen;
  else if(rms==18) return kMagenta;
  //if(layerIndex==1) return kBlue;
  //if(layerIndex==2) return kGreen;
  //if(layerIndex==3) return kYellow;

  return kOrange;
}

TGraph* getGraph(std::string plot = "hit_occup_", int lay_index=0, int rms=9)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(111);
  //const int npoints = 3;
  const int npoints = 7;

 // if(rms==9)
 // 	npoints = sizeof(runNumbers_9rms)/sizeof(runNumbers_9rms[0]);
  std::cout<<"NPOINTS: "<<npoints<<std::endl;

  //float currents[npoints] = { 1.887, 1.737, 1.058, 0.891, 0.617, 0.486, 0.317, 0.167, 0.058};
  std::vector<float> v_xaxis;
#if defined(source_scan)
  v_xaxis.push_back(1./1.);
  v_xaxis.push_back(1./1.5);
  v_xaxis.push_back(1./2.2);
//  v_xaxis.push_back(1./3.3);
  v_xaxis.push_back(1./4.6);
 // v_xaxis.push_back(1./6.9);
  v_xaxis.push_back(1./10.);
  v_xaxis.push_back(1./22.);
  v_xaxis.push_back(1./46.);

#endif

#if defined(hv_scan)
  v_xaxis.push_back(530);
  v_xaxis.push_back(520);
  v_xaxis.push_back(490);
#endif

  std::vector<std::string> v_rootfiles;
  for(int ifile=0; ifile< npoints; ifile++) {
    std::cout<<"hi\n";
  	if(rms==9)
  		v_rootfiles.push_back("/eos/user/i/idrivask/Chara/GIF_JULY2022/histos/" +runNumbers_9rms[ifile] + "_histograms.root"); 
  	else if(rms==12)
  		v_rootfiles.push_back("/eos/user/i/idrivask/Chara/GIF_JULY2022/histos/" +runNumbers_12rms[ifile] + "_histograms.root"); 
  	else if (rms==15)
  		v_rootfiles.push_back("/eos/user/i/idrivask/Chara/GIF_JULY2022/histos/" +runNumbers_15rms[ifile] + "_histograms.root"); 
  	else if(rms==18)
  		v_rootfiles.push_back("/eos/user/i/idrivask/Chara/GIF_JULY2022/histos/" +runNumbers_18rms[ifile] + "_histograms.root"); 
    else if(rms==0) // source scan 9rms 
      v_rootfiles.push_back("/eos/user/i/idrivask/Chara/GIF_JULY2022/histos/" +runNumbers_SC9[ifile] + "_histograms.root");
  	else if(rms==1) //source scan 12rms
      v_rootfiles.push_back("/eos/user/i/idrivask/Chara/GIF_JULY2022/histos/" +runNumbers_SC12[ifile] + "_histograms.root");
  	else if(rms==2) //source scan 15rms
      v_rootfiles.push_back("/eos/user/i/idrivask/Chara/GIF_JULY2022/histos/" +runNumbers_SC15[ifile] + "_histograms.root");
  	
  }
 
 std::string yTitle="";

 std::vector<float> v_layers;
 

  for(int ifile = 0; ifile<npoints; ifile++) {

    std::cout<<v_rootfiles[ifile]<<" "<<v_xaxis.at(ifile)<<std::endl;
   	
   	int iter = 0;

#if defined(cluster_multiplicity_plot)
   	iter=1;
#endif
#if defined(sPdo_plot)
   	iter=2;
#endif
#if defined(cPdo_plot)
   	iter=5;
#endif
//    for(int ilayerIndex=0; ilayerIndex<4; ilayerIndex++) {
//    	std::cout<<iter<<std::endl;
   	iter = iter + lay_index*10; std::cout<<"ITER > "<<iter<<std::endl;
#if defined(hit_occup)
    yTitle = "<hit_occup>";
#endif

#if defined(n_clusters)
    yTitle = "<N_clusters>";
#endif

#if defined(cluster_multiplicity_plot)
    yTitle = "<clus_mult>";
#endif

#if defined(sPdo_plot)
    yTitle = "<sPDO>";
#endif

#if defined(cPdo_plot)
    yTitle = "<cPDO>";
#endif

//#if defined(n_clusters)
//    if(ilayerIndex==0) 
//#endif
   
#if defined(hit_occup)
    v_layers.push_back(getMPV(v_rootfiles[ifile], lay_index, plot)/1024); 
#endif

#if (defined(n_clusters) || defined(cluster_multiplicity_plot) || defined(sPdo_plot) || defined(cPdo_plot) )
    v_layers.push_back(getMPV(v_rootfiles[ifile], iter, plot)); 
	//iter += 10;
#endif

 //   }
 }

 TGraph* graph;

 //for(int i=0; i<4; i++) {
 TVectorT<float> v_xaxis_tv(v_xaxis.size(), &v_xaxis[0]);

 TVectorT<float> v_layer_tv(v_layers.size(), &v_layers[0]);

 graph = new TGraph(npoints, &v_xaxis_tv[0], &v_layer_tv[0] );

 std::string name = "gr_" + std::to_string(lay_index);
 graph->SetName(name.c_str());
 std::cout<<graph->GetName()<<std::endl;

 graph->SetLineColor(setColors(lay_index)+2);
 graph->SetMarkerColor(setColors(lay_index)+2);
 graph->SetMarkerStyle(setMarkers(lay_index));

 lay_index++;
//}
//for(int i=0; i<6; i++)
 TCanvas* canvas = new TCanvas();
 canvas->cd();
 canvas->SetGrid();
#if defined(source_scan)
 graph->GetXaxis()->SetTitle("1/AF");
#endif
#if defined(hv_scan)
 graph->GetXaxis()->SetTitle("HV [V]");
#endif

 graph->GetYaxis()->SetTitle(yTitle.c_str());
#if defined(hit_occup)
 graph->GetYaxis()->SetRangeUser(0,0.06);//(1000,1600);
 #if defined(hv_scan)
 graph->GetYaxis()->SetRangeUser(0.,0.006);
 #endif
#endif

#if defined(n_clusters)
 graph->GetYaxis()->SetRangeUser(0.,7);//(1000,1600);
#endif

#if defined(cluster_multiplicity_plot)
 graph->GetYaxis()->SetRangeUser(3,6);//(1000,1600);
#endif

#if defined(sPdo_plot)
 graph->GetYaxis()->SetRangeUser(120,420);//(1000,1600);
#endif

#if defined(cPdo_plot)
 graph->GetYaxis()->SetRangeUser(1000,1900);//(1000,1600);
#endif

 graph->Draw("APL");
// graph[1]->Draw("PL same");
// graph[2]->Draw("PL same");
// graph[3]->Draw("PL same");

/*
 grap_lay0 = (TGraph*)graph[0]->Clone("my_graph_0");
 grap_lay1 = (TGraph*)graph[1]->Clone("my_graph_1");
 grap_lay2 = (TGraph*)graph[2]->Clone("my_graph_2");
 grap_lay3 = (TGraph*)graph[3]->Clone("my_graph_3");

 std::cout<<grap_lay0->GetN()<<" "<<grap_lay0->GetPointY(0)<<std::endl;
*/// graph[4]->Draw("P same");
// graph[5]->Draw("P same");
// graph[6]->Draw("P same");
// graph[7]->Draw("P same");

/*
 TF1* f = new TF1("f", "[1]+[0]*x", 0., .3);

 graph[1]->Fit("f","R");

  
TF1* fline = new TF1("fline", "[1]+[0]*x", 0., 1.);
  fline->SetParameter(0,f->GetParameter(0));
  fline->SetParameter(1,f->GetParameter(1));
  fline->SetLineStyle(10);
  fline->SetLineColor(kBlack);
  fline->Draw("l same");
  TLatex l;
  l.DrawLatex(0.8,0.03,Form("#scale[.5]{#splitline{const = %3.5f #pm %3.5f}{slope = %3.5f #pm %3.5f}}", f->GetParameter(1), f->GetParError(1), f->GetParameter(0), f->GetParError(0)));
*/


 
// graph->Write();
// graph[1]->Write();
// graph[2]->Write();
// graph[3]->Write();
// graph[4]->Write();
// graph[5]->Write();
// graph[6]->Write();
// graph[7]->Write();

// canvas->Write();

 std::string outname = "./output/" + plot + ".pdf";
 canvas->SaveAs(outname.c_str());
 

// output_file->Close();

 return graph;
}

void makePlotsInOneCanvas_sourcescan(std::string final_plot_name="hit_occup")
{
  TGraph* gr_0_9 = getGraph("nclus_0_9", 0, 0);
  TGraph* gr_1_9 = getGraph("nclus_1_9", 1, 0);
  TGraph* gr_2_9 = getGraph("nclus_2_9", 2, 0);
  TGraph* gr_3_9 = getGraph("nclus_3_9", 3, 0);

  TGraph* gr_0_12 = getGraph("nclus_0_12", 0, 1);
  TGraph* gr_1_12 = getGraph("nclus_1_12", 1, 1);
  TGraph* gr_2_12 = getGraph("nclus_2_12", 2, 1);
  TGraph* gr_3_12 = getGraph("nclus_3_12", 3, 1);

  TGraph* gr_0_15 = getGraph("nclus_0_15", 0, 2);
  TGraph* gr_1_15 = getGraph("nclus_1_15", 1, 2);
  TGraph* gr_2_15 = getGraph("nclus_2_15", 2, 2);
  TGraph* gr_3_15 = getGraph("nclus_3_15", 3, 2);


  gr_0_9->SetLineColor(kRed);
  gr_1_9->SetLineColor(kRed);
  gr_2_9->SetLineColor(kRed);
  gr_3_9->SetLineColor(kRed);

  gr_0_9->SetMarkerColor(kRed);
  gr_1_9->SetMarkerColor(kRed);
  gr_2_9->SetMarkerColor(kRed);
  gr_3_9->SetMarkerColor(kRed);

  gr_0_12->SetLineColor(kGreen);
  gr_1_12->SetLineColor(kGreen);
  gr_2_12->SetLineColor(kGreen);
  gr_3_12->SetLineColor(kGreen);

  gr_0_12->SetMarkerColor(kGreen);
  gr_1_12->SetMarkerColor(kGreen);
  gr_2_12->SetMarkerColor(kGreen);
  gr_3_12->SetMarkerColor(kGreen);

  gr_0_15->SetLineColor(kBlue);
  gr_1_15->SetLineColor(kBlue);
  gr_2_15->SetLineColor(kBlue);
  gr_3_15->SetLineColor(kBlue);

  gr_0_15->SetMarkerColor(kBlue);
  gr_1_15->SetMarkerColor(kBlue);
  gr_2_15->SetMarkerColor(kBlue);
  gr_3_15->SetMarkerColor(kBlue);

  TCanvas* canvas = new TCanvas();
  canvas->SetGrid();
  canvas->cd();
  canvas->SetRightMargin(.2);
  gr_0_9->Draw("APL");
  gr_1_9->Draw("PL SAME");
  gr_2_9->Draw("PL SAME");
  gr_3_9->Draw("PL SAME");
  gr_0_12->Draw("PL SAME");
  gr_1_12->Draw("PL SAME");
  gr_2_12->Draw("PL SAME");
  gr_3_12->Draw("PL SAME");
  gr_0_15->Draw("PL SAME");
  gr_1_15->Draw("PL SAME");
  gr_2_15->Draw("PL SAME");
  gr_3_15->Draw("PL SAME");


  auto legend = new TLegend(.8,0.12,1.,0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.);
  #if defined(SM1)
  legend->AddEntry(gr_0_9, "IP-1 #times 9 rms", "p");
  legend->AddEntry(gr_1_9, "IP-2 #times 9 rms", "p");
  legend->AddEntry(gr_2_9, "IP-3 #times 9 rms", "p");
  legend->AddEntry(gr_3_9, "IP-4 #times 9 rms", "p");
  legend->AddEntry(gr_0_12, "IP-1 #times 12 rms", "p");
  legend->AddEntry(gr_1_12, "IP-2 #times 12 rms", "p");
  legend->AddEntry(gr_2_12, "IP-3 #times 12 rms", "p");
  legend->AddEntry(gr_3_12, "IP-4 #times 12 rms", "p");
  legend->AddEntry(gr_0_15, "IP-1 #times 15 rms", "p");
  legend->AddEntry(gr_1_15, "IP-2 #times 15 rms", "p");
  legend->AddEntry(gr_2_15, "IP-3 #times 15 rms", "p");
  legend->AddEntry(gr_3_15, "IP-4 #times 15 rms", "p");
  #endif

 legend->Draw();

 std::string output_name = "./output/" + final_plot_name + ".pdf";
 canvas->SaveAs(output_name.c_str());

}

void makePlotsInOneCanvas(std::string final_plot_name="hit_occup_")
{
	
	//TGraph gr_0_9, gr_1_9, gr_2_9, gr_3_9, gr_0_12, gr_1_12, gr_2_12, gr_3_12;
	//occupancyPlot("test", &gr_0_9, &gr_1_9, &gr_2_9, &gr_3_9);
	TGraph* gr_0_9 = getGraph("nclus_0_9", 0, 9);
	TGraph* gr_1_9 = getGraph("nclus_1_9", 1, 9);
	TGraph* gr_2_9 = getGraph("nclus_2_9", 2, 9);
	TGraph* gr_3_9 = getGraph("nclus_3_9", 3, 9);

	TGraph* gr_0_12 = getGraph("nclus_0_12", 0, 12);
	TGraph* gr_1_12 = getGraph("nclus_1_12", 1, 12);
	TGraph* gr_2_12 = getGraph("nclus_2_12", 2, 12);
	TGraph* gr_3_12 = getGraph("nclus_3_12", 3, 12);

	TGraph* gr_0_15 = getGraph("nclus_0_15", 0, 15);
	TGraph* gr_1_15 = getGraph("nclus_1_15", 1, 15);
	TGraph* gr_2_15 = getGraph("nclus_2_15", 2, 15);
	TGraph* gr_3_15 = getGraph("nclus_3_15", 3, 15);

	TGraph* gr_0_18 = getGraph("nclus_0_18", 0, 18);
	TGraph* gr_1_18 = getGraph("nclus_1_18", 1, 18);
	TGraph* gr_2_18 = getGraph("nclus_2_18", 2, 18);
	TGraph* gr_3_18 = getGraph("nclus_3_18", 3, 18);


  gr_0_9->SetLineColor(kRed);
  gr_0_9->SetMarkerColor(kRed);
  gr_1_9->SetLineColor(kRed);
  gr_1_9->SetMarkerColor(kRed);
  gr_2_9->SetLineColor(kRed);
  gr_2_9->SetMarkerColor(kRed);
  gr_3_9->SetLineColor(kRed);
  gr_3_9->SetMarkerColor(kRed);

	gr_0_12->SetLineColor(setColors(12,true)+2);
 	gr_0_12->SetMarkerColor(setColors(12,true)+2);
 	gr_0_15->SetMarkerStyle(setMarkers(0,true));
 	gr_0_18->SetMarkerStyle(setMarkers(0,true));
 	gr_0_15->SetLineColor(setColors(15,true)+2);
 	gr_0_15->SetMarkerColor(setColors(15,true)+2);
 	gr_0_18->SetLineColor(setColors(18,true)+2);
 	gr_0_18->SetMarkerColor(setColors(18,true)+2);

 	gr_1_12->SetLineColor(setColors(12,true)+2);
 	gr_1_12->SetMarkerColor(setColors(12,true)+2);
 	gr_1_15->SetMarkerStyle(setMarkers(1,true));
 	gr_1_18->SetMarkerStyle(setMarkers(1,true));
 	gr_1_15->SetLineColor(setColors(15,true)+2);
 	gr_1_15->SetMarkerColor(setColors(15,true)+2);
 	gr_1_18->SetLineColor(setColors(18,true)+2);
 	gr_1_18->SetMarkerColor(setColors(18,true)+2);

 	gr_2_12->SetLineColor(setColors(12,true)+2);
 	gr_2_12->SetMarkerColor(setColors(12,true)+2);
 	gr_2_15->SetMarkerStyle(setMarkers(2,true) );
 	gr_2_18->SetMarkerStyle(setMarkers(2,true) );
 	gr_2_15->SetLineColor(setColors(15,true)+2);
 	gr_2_15->SetMarkerColor(setColors(15,true)+2);
 	gr_2_18->SetLineColor(setColors(18,true)+2);
 	gr_2_18->SetMarkerColor(setColors(18,true)+2);

 	gr_3_12->SetLineColor(setColors(12,true)+2);
 	gr_3_12->SetMarkerColor(setColors(12,true)+2);
 	gr_3_15->SetMarkerStyle(setMarkers(3,true));
 	gr_3_18->SetMarkerStyle(setMarkers(3,true));
 	gr_3_15->SetLineColor(setColors(15,true)+2);
 	gr_3_15->SetMarkerColor(setColors(15,true)+2);
 	gr_3_18->SetLineColor(setColors(18,true)+2);
 	gr_3_18->SetMarkerColor(setColors(18,true)+2);



	TCanvas* canvas = new TCanvas();
	canvas->SetGrid();
	canvas->cd();
	canvas->SetRightMargin(.2);
	gr_0_9->Draw("APL");
	gr_1_9->Draw("PL SAME");
	gr_2_9->Draw("PL SAME");
	gr_3_9->Draw("PL SAME");
	gr_0_12->Draw("PL SAME");
	gr_1_12->Draw("PL SAME");
	gr_2_12->Draw("PL SAME");
	gr_3_12->Draw("PL SAME");
	gr_0_15->Draw("PL SAME");
	gr_1_15->Draw("PL SAME");
	gr_2_15->Draw("PL SAME");
	gr_3_15->Draw("PL SAME");
	gr_0_18->Draw("PL SAME");
	gr_1_18->Draw("PL SAME");
	gr_2_18->Draw("PL SAME");
	gr_3_18->Draw("PL SAME");
	//std::cout<<gr_0_9.GetN()<<" "<<gr_1_9.GetN()<<std::endl;

	//std::cout<<gr->GetName()<<std::endl;


	auto legend = new TLegend(.8,0.12,1.,0.9);
 	legend->SetBorderSize(0);
 	legend->SetFillStyle(0);
 	legend->SetTextSize(0.);
	#if defined(SM1)
 	legend->AddEntry(gr_0_9, "IP-1 #times 9 rms", "p");
 	legend->AddEntry(gr_1_9, "IP-2 #times 9 rms", "p");
 	legend->AddEntry(gr_2_9, "IP-3 #times 9 rms", "p");
 	legend->AddEntry(gr_3_9, "IP-4 #times 9 rms", "p");
 	legend->AddEntry(gr_0_12, "IP-1 #times 12 rms", "p");
 	legend->AddEntry(gr_1_12, "IP-2 #times 12 rms", "p");
 	legend->AddEntry(gr_2_12, "IP-3 #times 12 rms", "p");
 	legend->AddEntry(gr_3_12, "IP-4 #times 12 rms", "p");
 	legend->AddEntry(gr_0_15, "IP-1 #times 15 rms", "p");
 	legend->AddEntry(gr_1_15, "IP-2 #times 15 rms", "p");
 	legend->AddEntry(gr_2_15, "IP-3 #times 15 rms", "p");
 	legend->AddEntry(gr_3_15, "IP-4 #times 15 rms", "p");
 	legend->AddEntry(gr_0_18, "IP-1 #times 18 rms", "p");
 	legend->AddEntry(gr_1_18, "IP-2 #times 18 rms", "p");
 	legend->AddEntry(gr_2_18, "IP-3 #times 18 rms", "p");
 	legend->AddEntry(gr_3_18, "IP-4 #times 18 rms", "p");
	#endif

	//#if defined(SB)
 	//legend->AddEntry(graph[0], "SBX-1", "p");
 	//legend->AddEntry(graph[1], "SBY-2", "p");
 	//legend->AddEntry(graph[2], "SBY-3", "p");
 	//legend->AddEntry(graph[3], "SBY-4", "p");
//	#endif

 legend->Draw();

 std::string output_name = "./output/" + final_plot_name + ".pdf";
 canvas->SaveAs(output_name.c_str());

}



EColor setColorsAtten(int af)
{

  if(af==0) return  kBlack; // AF = 1
  if(af==1) return  kRed; // AF = 1.5
  if(af==2) return  kGreen; // AF = 2.2
  if(af==3) return  kBlue; // AF = 4.6
  if(af==4) return  kMagenta; // AF = 10
  if(af==5) return  kCyan;   // AF = 22
  if(af==7) return  kTeal; // AF = 46
  if(af==6) return  kViolet; // NO SOURCE
  if(af==8) return  kPink;
  if(af==9) return  kTeal;
  if(af==10) return kAzure;
  if(af==11) return kSpring;

  return kOrange;
}

void make_plot(std::string sc_range, int histo_id=12) {
    //gStyle->SetOptStat(0);
    const int npoints = sizeof(runNumbers)/sizeof(runNumbers[0]);
    std::string rootfiles[npoints];// = {"./histos/1654249364_histograms.root", "./histos/1654248067_histograms.root", "./histos/1654246091_histograms.root", "./histos/1654203741_histograms.root", "./histos/1654201895_histograms.root", "./histos/1654186091_histograms.root", "./histos/1654250774_histograms.root", "./histos/1654252014_histograms.root"}; 

    for (int i=0; i<npoints; i++) {
     std::string histo_name = "./histos/" + runNumbers[i] + "_histograms.root";
     rootfiles[i] = histo_name;
   }

    TCanvas* canvas_cluster_mult = new TCanvas();
    std::string canvas_name = "canvas_" + sc_range;
    canvas_cluster_mult->SetName(canvas_name.c_str());
    canvas_cluster_mult->cd();
    canvas_cluster_mult->SetGrid();
    TH1F* histos[npoints];


    for(int ifile = 0; ifile<npoints; ifile++) {

 //   for(int ilayerIndex=0; ilayerIndex<8; ilayerIndex++) {
      histos[ifile] = getHistogram(rootfiles[ifile].c_str(), histo_id);

     // histos[ifile]->SetName(name.c_str());
      
      if(ifile==6) histos[ifile]->SetLineColor(setColorsAtten(ifile)+2);
      else if(ifile==7) histos[ifile]->SetLineColor(setColorsAtten(ifile)-8);
      else histos[ifile]->SetLineColor(setColorsAtten(ifile));
      histos[ifile]->SetLineWidth(2);
     // histos[ifile]->SetMarkerColor(setColors(ifile)+2);
     // histos[ifile]->SetMarkerStyle(setMarkers(ifile));
     
     if(ifile ==0 )histos[ifile]->Draw("histE");
     else histos[ifile]->Draw("histE same");
 //   }
    }

//#ifdef cluster_multiplicity_plot
 //   histos[0]->SetBinContent(2,0);
 //   histos[1]->SetBinContent(2,0);
 //   histos[2]->SetBinContent(2,0);
 //   histos[3]->SetBinContent(2,0);
 //   histos[4]->SetBinContent(2,0);
 //   histos[5]->SetBinContent(2,0);
 //   histos[6]->SetBinContent(2,0);
 //   histos[7]->SetBinContent(2,0);
//    histos[8]->SetBinContent(2,0);
//    histos[9]->SetBinContent(2,0);
//    histos[10]->SetBinContent(2,0);
//    histos[11]->SetBinContent(2,0);
//#endif
/*
#if defined(cPdo_plot) || defined(tot_charge)
    histos[0]->RebinX(80);
    histos[1]->RebinX(80);
    histos[2]->RebinX(80);
    histos[3]->RebinX(80);
    histos[4]->RebinX(80);
    histos[5]->RebinX(80);
    histos[6]->RebinX(80);
    histos[7]->RebinX(80);
    histos[8]->RebinX(80);
    histos[9]->RebinX(80);
    histos[10]->RebinX(80);
    histos[11]->RebinX(80);
#endif
*/    
    std::cout<<histos[0]->GetNbinsX()<<std::endl;
    histos[0]->SetTitle("SM1");
    histos[0]->GetYaxis()->SetRangeUser(10e-01, 10e+06 );

   // histos[0]->GetXaxis()->SetRangeUser(0,1024);
    canvas_cluster_mult->SetLogy();
    //canvas_cluster_mult->SetLogx();

    auto legend = new TLegend(0.8,0.1,1.,0.3);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.);
    legend->AddEntry(histos[0], "AF = 1", "f");
    legend->AddEntry(histos[1], "AF = 1.5", "f");
    legend->AddEntry(histos[2], "AF = 2.2", "f");
    legend->AddEntry(histos[3], "AF = 4.6", "f");
    legend->AddEntry(histos[4], "AF = 10", "f");
    legend->AddEntry(histos[5], "AF = 22", "f");
    legend->AddEntry(histos[6], "source off", "f");

    legend->Draw();

   TLatex text1, text2, text3, text4, text5, text6, text7, text8, text9, text10, text11, text12;

    text1.DrawLatex(5,6.*10e+05,Form("#scale[.4]{#color[618]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[0]->GetEntries(),histos[0]->GetMean(),histos[0]->GetStdDev()));
    text2.DrawLatex(5,2.*10e+05,Form("#scale[.4]{#color[402]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[1]->GetEntries(),histos[1]->GetMean(),histos[1]->GetStdDev()));
    text3.DrawLatex(5,6.*10e+04,Form("#scale[.4]{#color[434]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[2]->GetEntries(),histos[2]->GetMean(),histos[2]->GetStdDev()));
    text4.DrawLatex(5,2.*10e+04,Form("#scale[.4]{#color[602]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[3]->GetEntries(),histos[3]->GetMean(),histos[3]->GetStdDev()));
    text5.DrawLatex(5,6.*10e+03,Form("#scale[.4]{#color[418]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[4]->GetEntries(),histos[4]->GetMean(),histos[4]->GetStdDev()));
    text6.DrawLatex(5,2.*10e+03,Form("#scale[.4]{#color[634]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[5]->GetEntries(),histos[5]->GetMean(),histos[5]->GetStdDev()));
    text7.DrawLatex(5,6.*10e+02,Form("#scale[.4]{#color[802]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[6]->GetEntries(),histos[6]->GetMean(),histos[6]->GetStdDev()));
 //   text8.DrawLatex(52000,2.*10e+02,Form("#scale[.4]{#color[882]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[7]->GetEntries(),histos[7]->GetMean(),histos[7]->GetStdDev()));
 //   text9.DrawLatex(52000,6.*10e+01,Form("#scale[.4]{#color[902]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[8]->GetEntries(),histos[8]->GetMean(),histos[8]->GetStdDev()));
 //   text10.DrawLatex(52000,2.*10e+01,Form("#scale[.4]{#color[842]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[9]->GetEntries(),histos[9]->GetMean(),histos[9]->GetStdDev()));
 //   text11.DrawLatex(52000,2.*10e+00,Form("#scale[.4]{#color[862]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[10]->GetEntries(),histos[10]->GetMean(),histos[10]->GetStdDev()));
 //   text12.DrawLatex(52000,6.*10e+00,Form("#scale[.4]{#color[822]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[11]->GetEntries(),histos[11]->GetMean(),histos[11]->GetStdDev()));
    std::string out_name = "SM1_sclusters_" + sc_range + ".root";
    canvas_cluster_mult->SaveAs(out_name.c_str());

}
/*
    text1.DrawLatex(200,5.*10e+05,Form("#scale[.4]{#color[618]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[0]->GetEntries(),histos[0]->GetMean(),histos[0]->GetStdDev()));
    text2.DrawLatex(200,1.5*10e+05,Form("#scale[.4]{#color[402]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[1]->GetEntries(),histos[1]->GetMean(),histos[1]->GetStdDev()));
    text3.DrawLatex(200,4.5*10e+04,Form("#scale[.4]{#color[434]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[2]->GetEntries(),histos[2]->GetMean(),histos[2]->GetStdDev()));
    text4.DrawLatex(200,1.5*10e+04,Form("#scale[.4]{#color[602]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[3]->GetEntries(),histos[3]->GetMean(),histos[3]->GetStdDev()));
    text5.DrawLatex(200,5.*10e+03,Form("#scale[.4]{#color[418]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[4]->GetEntries(),histos[4]->GetMean(),histos[4]->GetStdDev()));
    text6.DrawLatex(200,5.*10e+05,Form("#scale[.4]{#color[634]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[5]->GetEntries(),histos[5]->GetMean(),histos[5]->GetStdDev()));
    text7.DrawLatex(200,1.5*10e+05,Form("#scale[.4]{#color[802]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[6]->GetEntries(),histos[6]->GetMean(),histos[6]->GetStdDev()));
    text8.DrawLatex(200,4.5*10e+04,Form("#scale[.4]{#color[882]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[7]->GetEntries(),histos[7]->GetMean(),histos[7]->GetStdDev()));
    text9.DrawLatex(200,1.5*10e+04,Form("#scale[.4]{#color[902]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[8]->GetEntries(),histos[8]->GetMean(),histos[8]->GetStdDev()));
    text10.DrawLatex(200,5.*10e+03,Form("#scale[.4]{#color[842]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[9]->GetEntries(),histos[9]->GetMean(),histos[9]->GetStdDev()));
    text11.DrawLatex(200,1.5*10e+03,Form("#scale[.4]{#color[862]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[10]->GetEntries(),histos[10]->GetMean(),histos[10]->GetStdDev()));
    text12.DrawLatex(200,1.5*10e+03,Form("#scale[.4]{#color[822]{#splitline{Events = %3.0f}{#splitline{Mean = %3.2f}{Std = %3.2f}}}}",histos[10]->GetEntries(),histos[11]->GetMean(),histos[11]->GetStdDev()));


    canvas_cluster_mult->SaveAs("nclus_final.pdf");
}
*/
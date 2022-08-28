#include "TreeReader.C"

int run(std::string run_number, std::string sector="C14")
{
	gROOT->SetBatch(kTRUE);
	
	std::string filename = "/eos/atlas/atlascerngroupdisk/det-nsw/bb5/cosmics/data/GIF++JUL2022/data_test."+run_number+"._.daq.RAW._lb0000._BB5-"+sector+"-MM-swROD._0001.simple.root";

	TreeReader* treeReader = new TreeReader(filename.c_str(),"nsw");
	std::cout<<"got the tree\n";
	TTree *tree = treeReader->the_tree;
	std::cout<<"set the tree with entries "<<tree->GetEntries()<<std::endl;
	int Nentries = tree->GetEntries();

	TH1F* h_hot_channels_per_layer[8];
	TH1F* h_hot_channels_per_layer_MMFE8[8][16];

	std::string hot_ch_fname = "./hotch/" +run_number + "_hot.root";
	TFile* f_hotch = new TFile(hot_ch_fname.c_str(), "RECREATE");

	for(int iLayer=0; iLayer<8; iLayer++)
	{

		h_hot_channels_per_layer[iLayer] = new TH1F(Form("h_hot_channels_per_layer_%i",iLayer), Form("Layer - %i",iLayer), 8193, 0, 8192);

		for(int iboard = 0 ; iboard< 16; iboard++)
		{
			int min=0;
			int max=511;
			if(iboard==0){ min=0; max = 511;}
			if(iboard==1){ min=512; max = 1023;}
			if(iboard==2){ min=1024; max = 1535;}
			if(iboard==3){ min=1536; max = 2047;}
			if(iboard==4){ min=2048; max = 2559;}
			if(iboard==5){ min=2560; max = 3071;}
			if(iboard==6){ min=3072; max = 3583;}
			if(iboard==7){ min=3584; max = 4095;}
			if(iboard==8){ min=4096; max = 4607;}
			if(iboard==9){ min=4608; max = 5119;}
			if(iboard==10){ min=5120; max = 5631;}
			if(iboard==11){ min=5632; max = 6143;}
			if(iboard==12){ min=6144; max = 6655;}
			if(iboard==13){ min=6656; max = 7167;}
			if(iboard==14){ min=7168; max = 7679;}
			if(iboard==15){ min=7680; max = 8192;}

			h_hot_channels_per_layer_MMFE8[iLayer][iboard] = new TH1F(Form("h_hot_channels_per_layer_MMFE8_%i_%i",iboard,iLayer), Form("MMFE8 - %i",iboard), 512, min, max);
		} 

	}

	// loop over events
	for(int iEvent = 0; iEvent<Nentries; iEvent++)
	{
		tree->GetEvent(iEvent);
		std::vector<unsigned int> *layers = treeReader->layer;
		std::vector<int> *strips = treeReader->strip;
		std::vector<unsigned int> *nhits = treeReader->nhits;
		std::cout<<iEvent<<std::endl;

			for(int iLayer=0; iLayer<8; iLayer++)
			{
				for(unsigned int iHit=0; iHit<nhits->size(); iHit++)
				{
					if(layers->at(iHit)==iLayer) 
					{
							h_hot_channels_per_layer[iLayer]->Fill(strips->at(iHit));
							if(strips->at(iHit) < 512) h_hot_channels_per_layer_MMFE8[iLayer][0]->Fill(strips->at(iHit));
							if(strips->at(iHit) >= 512 && strips->at(iHit)<=1023) h_hot_channels_per_layer_MMFE8[iLayer][1]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 1023 && strips->at(iHit)<=1535) h_hot_channels_per_layer_MMFE8[iLayer][2]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 1535 && strips->at(iHit)<=2047) h_hot_channels_per_layer_MMFE8[iLayer][3]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 2047 && strips->at(iHit)<=2559) h_hot_channels_per_layer_MMFE8[iLayer][4]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 2559 && strips->at(iHit)<=3071) h_hot_channels_per_layer_MMFE8[iLayer][5]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 3071 && strips->at(iHit)<=3583) h_hot_channels_per_layer_MMFE8[iLayer][6]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 3583 && strips->at(iHit)<=4095) h_hot_channels_per_layer_MMFE8[iLayer][7]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 4095 && strips->at(iHit)<=4607) h_hot_channels_per_layer_MMFE8[iLayer][8]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 4607 && strips->at(iHit)<=5119) h_hot_channels_per_layer_MMFE8[iLayer][9]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 5119 && strips->at(iHit)<=5631) h_hot_channels_per_layer_MMFE8[iLayer][10]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 5631 && strips->at(iHit)<=6143) h_hot_channels_per_layer_MMFE8[iLayer][11]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 6143 && strips->at(iHit)<=6655) h_hot_channels_per_layer_MMFE8[iLayer][12]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 6655 && strips->at(iHit)<=7167) h_hot_channels_per_layer_MMFE8[iLayer][13]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 7167 && strips->at(iHit)<=7679) h_hot_channels_per_layer_MMFE8[iLayer][14]->Fill(strips->at(iHit));
							if(strips->at(iHit) > 7679 && strips->at(iHit)<=8192) h_hot_channels_per_layer_MMFE8[iLayer][15]->Fill(strips->at(iHit));
					}
				}
			}
			
			continue; //next event


			
	} // end loop over layers	


		TF1* low_func[8][16];
		TF1* high_func[8][16];
		int count_noisy_channels=0;
		std::string hot_name = "./hotch/hotch_" + run_number + ".txt";
		std::ofstream outfile (hot_name.c_str());
		for (int ilay = 0; ilay < 8; ilay++){
    			for (int iboard = 0; iboard < 16; iboard++){
      				std::vector<int> bin_contents;
      				int min;
					int max;
					if(iboard==0){ min=0; max = 511;}
					if(iboard==1){ min=512; max = 1023;}
					if(iboard==2){ min=1024; max = 1535;}
					if(iboard==3){ min=1536; max = 2047;}
					if(iboard==4){ min=2048; max = 2559;}
					if(iboard==5){ min=2560; max = 3071;}
					if(iboard==6){ min=3072; max = 3583;}
					if(iboard==7){ min=3584; max = 4095;}
					if(iboard==8){ min=4096; max = 4607;}
					if(iboard==9){ min=4608; max = 5119;}
					if(iboard==10){ min=5120; max = 5631;}
					if(iboard==11){ min=5632; max = 6143;}
					if(iboard==12){ min=6144; max = 6655;}
					if(iboard==13){ min=6656; max = 7167;}
					if(iboard==14){ min=7168; max = 7679;}
					if(iboard==15){ min=7680; max = 8192;}
					low_func[ilay][iboard] = new TF1(Form("pol0_%i_%i", ilay, iboard),"pol0(0)", min, max);
					high_func[ilay][iboard] = new TF1(Form("high_pol0_%i_%i", ilay, iboard),"pol0(0)", min, max);

     				for(int istrip = 0; istrip < h_hot_channels_per_layer_MMFE8[ilay][iboard]->GetNbinsX(); istrip++){
       					 int strip_hits = h_hot_channels_per_layer_MMFE8[ilay][iboard]->GetBinContent(istrip+1);

        				 bin_contents.push_back(strip_hits);
     				 }
     			
     			std::sort(bin_contents.begin(), bin_contents.end());
     			bin_contents.erase(std::remove(bin_contents.begin(), bin_contents.end(), 0), bin_contents.end());
     			
     			if(bin_contents.empty() || bin_contents.size()<1) continue;
      			low_func[ilay][iboard]->SetParameter(0, 0.2*bin_contents.at(bin_contents.size()/2));
      			high_func[ilay][iboard]->SetParameter(0, bin_contents.at(bin_contents.size()/2)*3);

      			//here we select the strips that we need to reject
      			for(int rstrip=0; rstrip<h_hot_channels_per_layer_MMFE8[ilay][iboard]->GetNbinsX(); rstrip++ )
      			{
      			//	std::cout<<h_hot_channels_per_layer_MMFE8[ilay][iboard]->GetBinContent(rstrip)<<" "<< low_func[ilay][iboard]->GetParameter(0)<<std::endl;
      				if(h_hot_channels_per_layer_MMFE8[ilay][iboard]->GetBinContent(rstrip+1) < low_func[ilay][iboard]->GetParameter(0) || h_hot_channels_per_layer_MMFE8[ilay][iboard]->GetBinContent(rstrip+1) > high_func[ilay][iboard]->GetParameter(0) ) 
      				{ count_noisy_channels++; outfile << ilay << " "<< (int)h_hot_channels_per_layer_MMFE8[ilay][iboard]->GetBinCenter(rstrip+1) << std::endl;}
      			}
      			
      			}

      			std::cout<<"#of strips which need to be rejected "<<count_noisy_channels<<std::endl; 
      			count_noisy_channels=0;
      		}
      	outfile.close();
		TCanvas *c_hot_tot=new TCanvas();
		c_hot_tot->SetTitle("Check low & high occupancy channels");
		c_hot_tot->Divide(4,2);

		TCanvas *c_hot_boards_layer[8];
	
		for(int iLayer=0; iLayer<8; iLayer++){
			c_hot_tot->cd(iLayer+1);
			gPad->SetLogy();
			h_hot_channels_per_layer[iLayer]->Draw("ep");
			h_hot_channels_per_layer[iLayer]->GetXaxis()->SetTitle("strip number");
			h_hot_channels_per_layer[iLayer]->GetYaxis()->SetTitle("hits");

			c_hot_boards_layer[iLayer] = new TCanvas();
			c_hot_boards_layer[iLayer]->Divide(4,4);

			h_hot_channels_per_layer[iLayer]->Write();
		
			for(int iboard=0; iboard<16; iboard++){
				c_hot_boards_layer[iLayer]->cd(iboard+1);
		 		gPad->SetLogy();
		 		h_hot_channels_per_layer_MMFE8[iLayer][iboard]->Draw("ep");
		 		h_hot_channels_per_layer_MMFE8[iLayer][iboard]->GetXaxis()->SetTitle("strip number");
		 		h_hot_channels_per_layer_MMFE8[iLayer][iboard]->GetYaxis()->SetTitle("hits");
		 		h_hot_channels_per_layer_MMFE8[iLayer][iboard]->GetYaxis()->SetRangeUser(1, 10e+05);
		 		low_func[iLayer][iboard]->SetLineColor(kRed);
		 		high_func[iLayer][iboard]->SetLineColor(kRed);
		 		low_func[iLayer][iboard]->Draw("same");
		 		high_func[iLayer][iboard]->Draw("same");
			}

			h_hot_channels_per_layer[iLayer]->Write();
		}

		c_hot_tot->SaveAs("./hotch/hot_channels.pdf");
		f_hotch->cd();
		
		for(int lay=0; lay<8; lay++)
		{
			std::string name = "./hotch/hot_channels_layer_"+std::to_string(lay)+".pdf";
			c_hot_boards_layer[lay]->SaveAs(name.c_str());
			c_hot_boards_layer[lay]->Write();
		}

		f_hotch->Close();
	
	return 0;
}

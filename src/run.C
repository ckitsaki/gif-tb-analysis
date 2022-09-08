
//#define TIME_CALIBRATION
#include "TreeReader.C"
#include "Track.h"

int run(std::string run_number, std::string sector="C14")
{

	gROOT->SetBatch(kTRUE); 

	int iev = 0;
	TF1* f_gaus = new TF1("fgaus","gaus",0,200);

#ifdef TIME_CALIBRATION
	std::string f_calibSM1_name, f_calibSB_name;
	std::vector<float> v_ti_SM1, v_ti_SB, v_slope_SM1, v_slope_SB;
	for(int i=0; i<4; i++)
	{
		f_calibSM1_name = "./calib/"+run_number+"_calib_" + std::to_string(i) + "_SM1.root";//"./calib/" + run_number + "_calib_" + std::to_string(i) + "_SM1.root";
		TFile* file_tmp = new TFile(f_calibSM1_name.c_str());
		TH1F* h_ti_tmp = (TH1F*)file_tmp->Get("h_ti");
		TH1F* h_slope_tmp = (TH1F*)file_tmp->Get("h_slope");
		h_ti_tmp->Fit("fgaus","Q");
		float ti = f_gaus->GetParameter(1); //get mean from gaussian fit
		v_ti_SM1.push_back(ti);
		h_slope_tmp->Fit("fgaus","Q");
		float slope = f_gaus->GetParameter(1); //get mean from gaussian fit
		v_slope_SM1.push_back(slope);
		file_tmp->Close();
		f_calibSB_name = "./calib/"+run_number+"_calib_" + std::to_string(i) + "_SB.root";//"./calib/" + run_number + "_calib_" + std::to_string(i) + "_SB.root";
		file_tmp = new TFile(f_calibSB_name.c_str());
		h_ti_tmp = (TH1F*)file_tmp->Get("h_ti");
		h_slope_tmp = (TH1F*)file_tmp->Get("h_slope");
		h_ti_tmp->Fit("fgaus","Q");
		ti = f_gaus->GetParameter(1); //get mean from gaussian fit
		v_ti_SB.push_back(ti);
		h_slope_tmp->Fit("fgaus","Q");
		slope = f_gaus->GetParameter(1); //get mean from gaussian fit
		v_slope_SB.push_back(slope);
		file_tmp->Close();	
	}
#endif

	std::string fname = "./histos/" + run_number + "_histograms.root";
	
	TFile *out_file = new TFile(fname.c_str(),"RECREATE"); 
	out_file->mkdir("track_candidates/mult_tracks");
	out_file->mkdir("track_candidates/clean_mult_tracks");
	out_file->mkdir("track_candidates/clean_mult_tracks_morechoices");
	out_file->mkdir("track_candidates/track_displays");
	out_file->mkdir("track_candidates/residuals");
	out_file->mkdir("track_candidates/SM1");
	
	gROOT->cd();

	std::string filename = "/eos/atlas/atlascerngroupdisk/det-nsw/bb5/cosmics/data/GIF++JUL2022/data_test."+run_number+"._.daq.RAW._lb0000._BB5-"+sector+"-MM-swROD._0001.simple.root";//"/eos/atlas/atlascerngroupdisk/det-nsw/bb5/cosmics/data/GIF++JUN2022/data_test."+run_number+"._.daq.RAW._lb0000._BB5-"+sector+"-MM-swROD._0001.simple.root";

	TreeReader* treeReader = new TreeReader(filename.c_str(),"nsw");
	std::cout<<"got the tree\n";
	TTree *tree = treeReader->the_tree;
	std::cout<<"set the tree with entries "<<tree->GetEntries()<<std::endl;
	int Nentries = tree->GetEntries(); 
	
	Histograms* histos = new Histograms();

	int single_strip_counter[4]={0, 0, 0, 0};
	int single_strip_counter_small[4]={0, 0, 0, 0};
	
// Counters
	int counter0=0; //total number of events
	int counter1=0; //events with exactly one cluster on SBY1
	int counter3=0; //reference tracks
	int counter4=0;


	int count_SM1_zero_hits=0;

	int count_SM1_eta_out=0;
	int count_SM1_eta_in=0;
	int count_SM1_stereo_out=0;
	int count_SM1_stereo_in=0;
	int count_SM1_stereo=0;
	int count_SM1_stereo_both=0;

	int count_events_hit_on_both_stereo=0;
	int count_accepted_singletracks_outside_1deg=0;
	int count_accepted_singletracks_inside_1deg=0;
	int count_accepted_multipletracks_inside_1deg=0;
	int count_mult_track_events=0;
	int counter_mult_single_tracks=0;

///////////////////////////////////////////////////////

	int i=1;
	int divEvents = tree->GetEntries()/8;
	int ievent=0;
	int count_events=0;

	int count_scluster_in_1deg=0;
	int count_scluster_out_1deg=0;
	int count_sclusters_morethan1=0;
	int count_sclusters_morethan1_in_1deg=0;
	int count_sclusters_morethan1_out_1deg=0;
	int count_accepted_tracks=0;
	// loop over events
	for(int iEvent = 0; iEvent< Nentries; iEvent++)
	{
		gROOT->cd();

		tree->GetEvent(iEvent);
		std::vector<unsigned int> *layers = treeReader->layer;
		std::vector<int> *strips = treeReader->strip;
		std::vector<unsigned int> *nhits = treeReader->nhits;
		std::vector<unsigned int> *radius = treeReader->radius;
		std::vector<unsigned int> *pdo = treeReader->pdo;
		std::vector<unsigned int> *relbcid = treeReader->relbcid;
		std::vector<unsigned int> *bcid = treeReader->bcid;
		std::vector<unsigned int> *tdo = treeReader->tdo;

		counter0++; // total events

		int nlay_clus = 0;
		Layer* layer[4];
		Layer* l_small[4];
		Layer* trigger = new Layer(3, true);
		trigger->isTrigger();

		for(int iLayer = 0; iLayer<sizeof(layer)/sizeof(layer[0]); iLayer++){
			layer[iLayer] = new Layer(iLayer, true);
			l_small[iLayer] = new Layer(iLayer, false);
		}

		layer[0]->setAlphaBeta(0, 1);
		layer[1]->setAlphaBeta(0, 1);
		layer[2]->setAlphaBeta(0, 1); 
		layer[3]->setAlphaBeta(0, 1); 
		l_small[0]->setAlphaBeta(0, 1);
		l_small[1]->setAlphaBeta(0, 1);
		l_small[2]->setAlphaBeta(0, 1);
		l_small[3]->setAlphaBeta(0, 1); 

		if(iEvent == i*divEvents )
		{
			std::cout<<12.5*i<<"% of the events processed"<<std::endl;
			i++;
		}

		//std::cout<<iEvent<<std::endl;

		for(int iHit = 0; iHit<nhits->size(); iHit++)
		{
			if(layers->at(iHit) == 3 && radius->at(iHit)==6) trigger->AddHitIndex(iHit);
			for(int ilay=0; ilay<sizeof(layer)/sizeof(layer[0]); ilay++) {//loop over layers
				if(layers->at(iHit) == ilay && ( radius->at(iHit)==4 || radius->at(iHit)==5 )) layer[ilay]->AddHitIndex(iHit);
				else if(layers->at(iHit) == ilay && ( radius->at(iHit)==2 || radius->at(iHit)==3 ) ) l_small[ilay]->AddHitIndex(iHit);
				else continue;
			}

		}

		std::vector<int> v_strips_layer[sizeof(layer)/sizeof(layer[0])];
		std::vector<int> v_strips_layer_small[sizeof(l_small)/sizeof(l_small[0])];
		std::vector<int> v_strips_trigger;

		histos->h_raw_hit_trigger->Fill(trigger->getNHits());
		trigger->bookFiredStrips(strips,pdo);
		v_strips_trigger = trigger->getFiredStrips();

		for(int i=0; i<v_strips_trigger.size(); i++) {
			histos->h_strip_index_vs_tdo_trigger->Fill(v_strips_trigger.at(i),tdo->at(trigger->getHitIndex(v_strips_trigger.at(i))));
			histos->h_strip_index_vs_relbcid_trigger->Fill(v_strips_trigger.at(i),relbcid->at(trigger->getHitIndex(v_strips_trigger.at(i))));
			histos->h_strip_index_vs_bcid_trigger->Fill(v_strips_trigger.at(i),bcid->at(trigger->getHitIndex(v_strips_trigger.at(i))));
			histos->h_strip_index_vs_pdo0_trigger->Fill(v_strips_trigger.at(i),pdo->at(trigger->getHitIndex(v_strips_trigger.at(i))));
		}
		
		Cluster* cl_lay[sizeof(layer)/sizeof(layer[0])]; 
		Cluster* cl_lay_small[sizeof(l_small)/sizeof(l_small[0])];


		int total_hits_SM1 = 0; // how many layers have a hit
		int total_hits_TZ = 0;
		
		for(int iLayer = 0; iLayer<sizeof(layer)/sizeof(layer[0]); iLayer++) {
			histos->h_raw_hits[iLayer]->Fill(layer[iLayer]->getNHits()); // all hits before any rejection
			histos->h_raw_hits_small[iLayer]->Fill(l_small[iLayer]->getNHits());

			layer[iLayer]->bookFiredStrips(strips,pdo); // here we reject strips with pdo<64 and we also mask the noisy strips
			l_small[iLayer]->bookFiredStrips(strips,pdo);

			v_strips_layer[iLayer] = layer[iLayer]->getFiredStrips();
			v_strips_layer_small[iLayer] = l_small[iLayer]->getFiredStrips();
			
			for(int i=0; i<v_strips_layer_small[iLayer].size(); i++)
			{
			  // int stripIndex = (v_strips_layer_small[iLayer].at(i)<1536) ? (std::abs(1535-v_strips_layer_small[iLayer].at(i))+1024) : v_strips_layer_small[iLayer].at(i);
			   histos->h_strip_index_SB[iLayer]->Fill(v_strips_layer_small[iLayer].at(i));
			}

			for(int i=0; i<v_strips_layer[iLayer].size(); i++) {
				histos->h_strip_index_vs_tdo[iLayer]->Fill(v_strips_layer[iLayer].at(i),tdo->at(layer[iLayer]->getHitIndex(v_strips_layer[iLayer].at(i))));
				histos->h_strip_index_vs_relbcid[iLayer]->Fill(v_strips_layer[iLayer].at(i),relbcid->at(layer[iLayer]->getHitIndex(v_strips_layer[iLayer].at(i))));
				histos->h_strip_index_vs_bcid[iLayer]->Fill(v_strips_layer[iLayer].at(i),bcid->at(layer[iLayer]->getHitIndex(v_strips_layer[iLayer].at(i))));
				histos->h_strip_index_vs_pdo0[iLayer]->Fill(v_strips_layer[iLayer].at(i),pdo->at(layer[iLayer]->getHitIndex(v_strips_layer[iLayer].at(i))));
			}

			for(int i=0; i<v_strips_layer_small[iLayer].size(); i++) {
				histos->h_strip_index_vs_tdo_small[iLayer]->Fill(v_strips_layer_small[iLayer].at(i),tdo->at(l_small[iLayer]->getHitIndex(v_strips_layer_small[iLayer].at(i))));
				histos->h_strip_index_vs_relbcid_small[iLayer]->Fill(v_strips_layer_small[iLayer].at(i),relbcid->at(l_small[iLayer]->getHitIndex(v_strips_layer_small[iLayer].at(i))));
				histos->h_strip_index_vs_bcid_small[iLayer]->Fill(v_strips_layer_small[iLayer].at(i),bcid->at(l_small[iLayer]->getHitIndex(v_strips_layer_small[iLayer].at(i))));
				histos->h_strip_index_vs_pdo0_small[iLayer]->Fill(v_strips_layer_small[iLayer].at(i),pdo->at(l_small[iLayer]->getHitIndex(v_strips_layer_small[iLayer].at(i))));
			}

			for(int i=0; i<v_strips_layer[iLayer].size(); i++)
			{
				histos->h_strip_index_SM1[iLayer]->Fill(v_strips_layer[iLayer].at(i));
			}
			
			cl_lay[iLayer] = new Cluster(layer[iLayer],pdo,relbcid,tdo,radius);
			cl_lay_small[iLayer] = new Cluster(l_small[iLayer],pdo,relbcid,tdo,radius);

			histos->h_nclusters[iLayer]->Fill(cl_lay[iLayer]->getNClusters2());
			histos->h_nclusters_small[iLayer]->Fill(cl_lay_small[iLayer]->getNClusters2());
			
			int tot_strips[sizeof(layer)/sizeof(layer[0])] = {0, 0, 0, 0};
			int tot_strips_small[sizeof(l_small)/sizeof(l_small[0])] = {0, 0, 0, 0};
			
			std::vector<float> v_lead_cl_charge;
			
			for(int icl = 0; icl<cl_lay[iLayer]->getNClusters2(); icl++) {
				std::vector<float> v_strip_times;
				std::vector<float> v_maxstrip_pdo_allclus;

				if(cl_lay[iLayer]->getTotPdo(icl)==2044){
					continue;
				}

				histos->h_nstrips[iLayer]->Fill(cl_lay[iLayer]->getNStrips(icl));
				histos->h_cl_charge[iLayer]->Fill(cl_lay[iLayer]->getTotPdo(icl));
				v_lead_cl_charge.push_back(cl_lay[iLayer]->getTotPdo(icl));
				
				histos->h_clus_positions[iLayer]->Fill(cl_lay[iLayer]->getPosition(icl));
				histos->h_cl_size_vs_cl_charge[iLayer]->Fill(cl_lay[iLayer]->getNStrips(icl), cl_lay[iLayer]->getTotPdo(icl));

				tot_strips[iLayer] += cl_lay[iLayer]->getNStrips(icl); // total number of strips for all clusters formed in the event
	
				
				for(int istrip = 0; istrip < cl_lay[iLayer]->getNStrips(icl); istrip++)
				{
					//if(cl_lay[iLayer]->getTotPdo(icl)==2044){
					//	histos->h_strip_2044_pdo_SM1[iLayer]->Fill(v_strips_layer[iLayer].at(istrip));
					//}
				#ifdef TIME_CALIBRATION
					float s_time = cl_lay[iLayer]->getStripRelBcid(icl, istrip)*25 - (cl_lay[iLayer]->getStripTdo(icl, istrip) - v_ti_SM1.at(iLayer))/v_slope_SM1.at(iLayer); //time in ns
					v_strip_times.push_back(s_time);
					histos->h_strip_time_in_clus[iLayer]->Fill(s_time);
				#endif
					v_maxstrip_pdo_allclus.push_back(cl_lay[iLayer]->getStripPdo(icl, istrip));
					histos->h_strip_pdo_in_clus[iLayer]->Fill(cl_lay[iLayer]->getStripPdo(icl, istrip));
					histos->h_strip_tdo_in_clus[iLayer]->Fill(cl_lay[iLayer]->getStripTdo(icl, istrip));
					histos->h_strip_relbcid_in_clus[iLayer]->Fill(cl_lay[iLayer]->getStripRelBcid(icl, istrip));
					histos->h_strip_index_vs_pdo[iLayer]->Fill(cl_lay[iLayer]->getStripIndex(icl, istrip), cl_lay[iLayer]->getStripPdo(icl, istrip));
					histos->h_strip_index_SM1_cluster[iLayer]->Fill(v_strips_layer[iLayer].at(istrip));
				}
				histos->h_maxstrip_pdo_allClus[iLayer]->Fill(*max_element(v_maxstrip_pdo_allclus.begin(),v_maxstrip_pdo_allclus.end()));
			#ifdef TIME_CALIBRATION
				float min_time = *min_element(v_strip_times.begin(), v_strip_times.end());
				histos->h_strip_time_in_clus_earliest[iLayer]->Fill(min_time);

				v_strip_times.clear();
				v_strip_times.shrink_to_fit();
				

			#endif
				v_maxstrip_pdo_allclus.clear();
				v_maxstrip_pdo_allclus.shrink_to_fit();
			} // for-loop over clusters SM1 
		
			if(!v_lead_cl_charge.empty()) {
			histos->h_cl_charge_leadCluster[iLayer]->Fill(*max_element(v_lead_cl_charge.begin(),v_lead_cl_charge.end()));
			auto it = find(v_lead_cl_charge.begin(), v_lead_cl_charge.end(), *max_element(v_lead_cl_charge.begin(),v_lead_cl_charge.end()));
			if (it != v_lead_cl_charge.end()) 
    		{
				int max_clus_index = it - v_lead_cl_charge.begin(); // index of the leading cluster
				std::vector<float> v_leadCluster_maxstrip_pdo;
				for(int ileadclus = 0; ileadclus < cl_lay[iLayer]->getNStrips(max_clus_index); ileadclus++)
				{
					histos->h_strip_pdo_in_leadCluster[iLayer]->Fill(cl_lay[iLayer]->getStripPdo(max_clus_index, ileadclus));
					v_leadCluster_maxstrip_pdo.push_back(cl_lay[iLayer]->getStripPdo(max_clus_index, ileadclus));
				}
				histos->h_cl_charge_leadCluster_maxstrip[iLayer]->Fill(*max_element(v_leadCluster_maxstrip_pdo.begin(),v_leadCluster_maxstrip_pdo.end()));
			
				v_leadCluster_maxstrip_pdo.clear();
				v_leadCluster_maxstrip_pdo.shrink_to_fit();
			}

		}
			v_lead_cl_charge.clear();
			v_lead_cl_charge.shrink_to_fit();

			for(int icl = 0; icl<cl_lay_small[iLayer]->getNClusters2(); icl++) {
				std::vector<float> v_strip_times;
				if(cl_lay_small[iLayer]->getTotPdo(icl)==2044) continue;

				//cl_lay_small[iLayer]->removeStrips_small(icl); // here I reject the events which give a strange peak at low pdos
				if(cl_lay_small[iLayer]->getNStrips(icl)==0) continue;
				histos->h_nstrips_small[iLayer]->Fill(cl_lay_small[iLayer]->getNStrips(icl));
				histos->h_cl_charge_small[iLayer]->Fill(cl_lay_small[iLayer]->getTotPdo(icl));
				histos->h_clus_positions_small[iLayer]->Fill(cl_lay_small[iLayer]->getPosition(icl));
				histos->h_cl_size_vs_cl_charge_small[iLayer]->Fill(cl_lay_small[iLayer]->getNStrips(icl), cl_lay_small[iLayer]->getTotPdo(icl));

				tot_strips_small[iLayer] += cl_lay_small[iLayer]->getNStrips(icl); // total number of strips for all clusters formed in the event

				
			for(int istrip = 0; istrip < cl_lay_small[iLayer]->getNStrips(icl); istrip++)
				{
					//if(cl_lay_small[iLayer]->getTotPdo(icl)==2044) {
					//	histos->h_strip_2044_pdo_SB[iLayer]->Fill(v_strips_layer_small[iLayer].at(istrip));
					//}

				#ifdef TIME_CALIBRATION
					//std::cout<<"strip - "<<cl_lay[iLayer]->getCluster(icl).at(istrip)<<" "<<cl_lay[iLayer]->getStripPdo(icl, istrip)<<" "<<cl_lay[iLayer]->getStripTdo(icl, istrip)<<" "<<cl_lay[iLayer]->getStripIndex(icl, istrip)<<std::endl;
					float s_time = cl_lay_small[iLayer]->getStripRelBcid(icl, istrip)*25 - (cl_lay_small[iLayer]->getStripTdo(icl, istrip) - v_ti_SB.at(iLayer))/v_slope_SB.at(iLayer); //time in ns
					//std::cout<<"strip time: "<<s_time<<" ns, strip "<<cl_lay[iLayer]->getStripIndex(icl, istrip)<<" layer "<< iLayer<<std::endl;;
					v_strip_times.push_back(s_time);
					histos->h_strip_time_in_clus_small[iLayer]->Fill(s_time);
				#endif
					histos->h_strip_pdo_in_clus_small[iLayer]->Fill(cl_lay_small[iLayer]->getStripPdo(icl, istrip));
					histos->h_strip_tdo_in_clus_small[iLayer]->Fill(cl_lay_small[iLayer]->getStripTdo(icl, istrip));
					histos->h_strip_relbcid_in_clus_small[iLayer]->Fill(cl_lay_small[iLayer]->getStripRelBcid(icl, istrip));
					histos->h_strip_index_vs_pdo_small[iLayer]->Fill(cl_lay_small[iLayer]->getStripIndex(icl, istrip), cl_lay_small[iLayer]->getStripPdo(icl, istrip));
					histos->h_strip_index_SB_cluster[iLayer]->Fill(v_strips_layer_small[iLayer].at(istrip));
				}
			#ifdef TIME_CALIBRATION
				float min_time = *min_element(v_strip_times.begin(), v_strip_times.end());
				histos->h_strip_time_in_clus_earliest_small[iLayer]->Fill(min_time);

				v_strip_times.clear();
				v_strip_times.shrink_to_fit();
			#endif
				
			} // for-loop over clusters SB 

			if(cl_lay[iLayer]->getNClusters2()>0){
				histos->h_raw_hits_vs_tot_strips[iLayer]->Fill(tot_strips[iLayer],layer[iLayer]->getNHits());
				single_strip_counter[iLayer] += cl_lay[iLayer]->getNSingleStrips();
				//total_hits_SM1++;		
			}

			if(cl_lay_small[iLayer]->getNClusters2()>0){
				histos->h_raw_hits_vs_tot_strips_small[iLayer]->Fill(tot_strips_small[iLayer],l_small[iLayer]->getNHits());
				single_strip_counter_small[iLayer] += cl_lay_small[iLayer]->getNSingleStrips();
				//total_hits_TZ++;		
			}


		} //end ilayer for-loop
//  End Clustering

// Start tracking
		if( !(cl_lay_small[1]->getNClusters2()==1)) continue; // only events which give exactly one cluster on SBY1 are processed		
		counter1++; // count events with exactly one cluster on SBY1 

// Alignment needs some iteration 
// Step 1: align MMFE8. Move the y position by the required number of strips (Layer::convertStripToGlobalY_radius, Layer::convertStripToGlobalY_radius_corr).
// Step 2: fill 2D distributions BSY1-yposition vs yposition. Perform a linear fit and get the slope and constant value. This should be done for SBY2, SBY3 and eta_out. Use for the rest SM1 layers the same parameters as done for eta_out.
// Step 3: fill the residuals SBY1-SBY2, SBY1-SBY3, SBY1-eta_out, SBY1-eta_in, SBY1-stereo_in, SBY1-stereo_out, SBY1-stereo. Find the mean value and use it as an additional shift
// Remember to apply the additional shift (+0.56) for the stereo first coordinate when is about to be used (ex. inside this function Track::fillStereoCluster() ). 
// some more plots may be needed similar to those described here as a cross-check. After each step the new corrected positions at a time need to be used.
// Another useful plot to study is the distance between the reconstructed cluster and the track found. Check for consistency if it peaks at 0.
// Dont forget to adjust the axes limits in Histograms.h
		layer[0]->setAlphaBeta(-405.641-0.16, (1-0.0318));
		layer[1]->setAlphaBeta(-405.641-0.23, (1-0.0318));
		layer[2]->setAlphaBeta(-405.641-4.56-0.027, (1-0.0318));
		layer[3]->setAlphaBeta(-405.641+4.144+0.022, (1-0.0318)); 
		
		l_small[2]->setAlphaBeta(33.000-0.05, (1-0.0448)); 
		l_small[3]->setAlphaBeta(35.778-0.1225, (1-0.0467)); 

		float pos_sby1 = cl_lay_small[1]->getPosition(0); // SBY1 reference chamber position

		if(cl_lay_small[2]->getNClusters2()==1 && cl_lay_small[3]->getNClusters2()==1)
		{
			float pos_sby2 = cl_lay_small[2]->getCorrPosition(0);
			float pos_sby3 = cl_lay_small[3]->getCorrPosition(0);
			histos->h_res_SBY2_SBY3_vs_SBY2->Fill(pos_sby2-pos_sby3, pos_sby2);
			histos->h_res_SBY2_SBY3_vs_SBY3->Fill(pos_sby2-pos_sby3, pos_sby3);
			histos->h_res_SBY2_SBY1_vs_SBY1->Fill(pos_sby2-pos_sby1, pos_sby1);
			histos->h_res_SBY3_SBY1_vs_SBY1->Fill(pos_sby3-pos_sby1, pos_sby1);
			histos->h_sby1_minus_sby2_vs_sby2->Fill(pos_sby2, pos_sby1-pos_sby2);
			histos->h_sby1_minus_sby3_vs_sby3->Fill(pos_sby3, pos_sby1-pos_sby3);
		}

		if(cl_lay[0]->getNClusters2()==1 && cl_lay[1]->getNClusters2()==1)
		{
			float pos_eta_out = cl_lay[0]->getCorrPosition(0);
			float pos_eta_in = cl_lay[1]->getCorrPosition(0);
			histos->h_diffpos_lay0->Fill(pos_eta_in-pos_eta_out, pos_eta_out);  
			histos->h_diffpos_lay1->Fill(pos_eta_in-pos_eta_out, pos_eta_in); 
			histos->h_pos_etain_vs_pos_etaout->Fill(cl_lay[0]->getCorrPosition(0), cl_lay[1]->getCorrPosition(0));	
		}

		if(cl_lay[0]->getNClusters2()==1)
		{
			float pos_eta_out = cl_lay[0]->getCorrPosition(0);
			histos->h_sby1_minus_eta_out_vs_pos_eta_out->Fill(pos_eta_out, pos_sby1-pos_eta_out);
		}
		if(cl_lay[1]->getNClusters2()==1)
		{
			float pos_eta_in = cl_lay[1]->getCorrPosition(0);
			histos->h_sby1_minus_pos_eta_in_vs_pos_eta_in->Fill(pos_eta_in, pos_sby1-pos_eta_in);
		}

		if(cl_lay[2]->getNClusters2()==1 && cl_lay[3]->getNClusters2()==1 && cl_lay[1]->getNClusters2()==1)
		{
			float pos_l1 = cl_lay[1]->getCorrPosition(0); //eta in
			float pos_l2 = cl_lay[2]->getCorrPosition(0); //stereo in
			float pos_l3 = cl_lay[3]->getCorrPosition(0); // stereo out
			float pos_stereo_y = ((pos_l3+pos_l2) / 2*TMath::Cos(1.5*TMath::Pi()/180.));
			float pos_stereo_x = ((pos_l3-pos_l2) / 2*TMath::Sin(1.5*TMath::Pi()/180.));
			histos->h_sby1_minus_stereo_vs_stereo->Fill(pos_sby1, pos_sby1-pos_stereo_y);
			histos->h_sby1_minus_stereo_in_vs_stereo_in->Fill(pos_sby1, pos_sby1-pos_l2);
			histos->h_sby1_minus_stereo_out_vs_stereo_out->Fill(pos_sby1, pos_sby1-pos_l3);
			histos->h_beamProfile->Fill(pos_stereo_y, pos_stereo_x);
		}

		if(cl_lay[2]->getNClusters2()==1 && cl_lay[3]->getNClusters2()==1 && cl_lay[0]->getNClusters2()==1)
		{
			float pos_l0 = cl_lay[0]->getCorrPosition(0);
			float pos_l2 = cl_lay[2]->getCorrPosition(0);
			float pos_l3 = cl_lay[3]->getCorrPosition(0);
			float pos_stereo_y = ((pos_l3+pos_l2) / 2*TMath::Cos(1.5*TMath::Pi()/180.));
			float pos_stereo_x = ((pos_l3-pos_l2) / 2*TMath::Sin(1.5*TMath::Pi()/180.));
			histos->h_diffpos_lay2->Fill( pos_l2, pos_l0-pos_l2);  
			histos->h_diffpos_lay3->Fill( pos_l3, pos_l0-pos_l3); 
			histos->h_diffpos_stereolay1->Fill(pos_sby1-pos_l3, pos_sby1);
			histos->h_diffpos_stereolay0->Fill(pos_sby1-pos_l2, pos_sby1);
		}

// ==========================================   end aligment		

// Track selection
		// accept all events which give at least 1 cluster on each of SBY2 and SBY3 small chambers 
		if(cl_lay_small[2]->getNClusters2()>0 && cl_lay_small[3]->getNClusters2()>0)
		{	
			//std::cout<<"EVENT-"<<iEvent<<std::endl;
			std::vector<TGraphAsymmErrors*> v_cand_tracks;
			Track* track = new Track(cl_lay_small[1], cl_lay_small[2], cl_lay_small[3], histos);
			
			// some counters 
			count_scluster_in_1deg += track->getSuperClustersWithin1DegWindow(); // count events which find the supercluster inside the 1 deg window
			count_scluster_out_1deg += track->getTheOnlySuperClusterOutside1DegWindow(); // count events which find only one supercluster but is outside the 1deg window
			count_sclusters_morethan1 += track->getEventsWithMoreThan1SuperCluster(); // count events which find more than one supercluster 
			count_sclusters_morethan1_in_1deg += track->getEventsWithMoreThan1SuperClusterIn1DegWindow(); // count events which find exactly 1 supercluster inside the 1deg window
			count_sclusters_morethan1_out_1deg += track->getEventsWithMoreThan1SuperClusterOut1DegWindow(); // count events which find superclusters outside the 1 deg window

			//	std::cout<<"Event-"<<iEvent<<std::endl;
			// study separately single and many track events
			// single tracks
			if(!track->tooMany() && !track->rejectTrack() && track->checkTrackCandidate())// && (track->getAngle()<=0.5 && track->getAngle()>=-0.5))
			{
				counter3++; //single track counter
				int count_layers_with_cluster = 0; // counts the number of layers with cluster on SM1	
				count_events++;
				
				TGraphAsymmErrors* gr_track = track->getTrack();
				gr_track->SetName(Form("gr_track_%i",iEvent));
				float tr_slope = track->trackSlope();
				float tr_const = track->trackIntercept();
				
				histos->h_residuals_single_tracks->Fill(gr_track->GetPointY(1) - gr_track->GetPointY(2)); // fill residuals on track for single tracks: Delta(SBY2 - SBY3)
				histos->h_residuals_accepted_tracks->Fill(gr_track->GetPointY(1) - gr_track->GetPointY(2)); // fill residuals on track : Delta(SBY2 - SBY3)

				std::vector<float> v_dclus_to_track[4]; // filled with all distances between track and each cluster found on eta layers
				for(int iLayer = 0; iLayer<4; iLayer++) 
				{
					if(cl_lay[iLayer]->getNClusters2()>0) count_layers_with_cluster++;
					histos->h_nclusters_per_layer_event[iLayer]->Fill(cl_lay[iLayer]->getNClusters2()); 
					for(int icl=0; icl<cl_lay[iLayer]->getNClusters2(); icl++)
					{
						float distance = std::abs(tr_slope*cl_lay[iLayer]->getZPosition() - cl_lay[iLayer]->getCorrPosition(icl) + tr_const)/(TMath::Sqrt(tr_slope*tr_slope + 1));
						v_dclus_to_track[iLayer].push_back(distance);
					}

				}
				std::vector<int> v_index_lay2;
				std::vector<int> v_index_lay3;
				std::vector<float> v_dclus_to_track_stereo;
				if(cl_lay[2]->getNClusters2()>0 && cl_lay[3]->getNClusters2()>0)
				{
					for(int icl2 =0; icl2<cl_lay[2]->getNClusters2(); icl2++)
					{
						float pos_cl2 = cl_lay[2]->getCorrPosition(icl2);
						for(int icl3=0; icl3<cl_lay[3]->getNClusters2(); icl3++)
						{
							v_index_lay2.push_back(icl2);
							v_index_lay3.push_back(icl3);
							float pos_cl3 = cl_lay[3]->getCorrPosition(icl3);
							float pos_stereo_y = +0.56 +((pos_cl3+pos_cl2) / 2*TMath::Cos(1.5*TMath::Pi()/180.)); 
							float pos_stereo_x = ((pos_cl3-pos_cl2) / 2*TMath::Sin(1.5*TMath::Pi()/180.));
							float distance = std::abs((tr_slope*(cl_lay[2]->getZPosition()+cl_lay[3]->getZPosition())/2 - pos_stereo_y + tr_const)/(TMath::Sqrt(tr_slope*tr_slope + 1)));
							v_dclus_to_track_stereo.push_back(distance);
						}
					}
				}

				for(int iLayer=0; iLayer<4; iLayer++) // fill the track graph with SM1 eta layers 
				{
					if(v_dclus_to_track[iLayer].empty()) continue;
					float min_distance = *min_element(v_dclus_to_track[iLayer].begin(), v_dclus_to_track[iLayer].end());
					int index = find(v_dclus_to_track[iLayer].begin(), v_dclus_to_track[iLayer].end(),min_distance) - v_dclus_to_track[iLayer].begin(); // index of the cluster which belongs to the track
					if(iLayer<2) track->fillSM1events(index, cl_lay[iLayer], iLayer, 0);	// eta layers
					else // for the stereo layers
					{
						if(v_dclus_to_track_stereo.empty()) // events with cluster on stereo in or stereo out
						{
							if(iLayer==2) track->checkDistanceFromTrack(0, iLayer, cl_lay[2], index, stereo_in_cut_eff);
							if(iLayer==3) track->checkDistanceFromTrack(0, iLayer, cl_lay[3], index, stereo_out_cut_eff);
						}
					}	
				}

				if(!v_dclus_to_track_stereo.empty())
				{
					count_events_hit_on_both_stereo++;
					// find the closest stereo cluster
					float min_distance = *min_element(v_dclus_to_track_stereo.begin(), v_dclus_to_track_stereo.end());
					int index = find(v_dclus_to_track_stereo.begin(), v_dclus_to_track_stereo.end(), min_distance) - v_dclus_to_track_stereo.begin();
					int ind_lay2 = v_index_lay2.at(index); // the cluster index on layer 2
					int ind_lay3 = v_index_lay3.at(index); // the cluster index on layer 3
					histos->h_d_track_ystereo->Fill(min_distance);
					float pos_stereo_y = 0.56 + ((cl_lay[3]->getCorrPosition(ind_lay3)+cl_lay[2]->getCorrPosition(ind_lay2)) / 2*TMath::Cos(1.5*TMath::Pi()/180.)); // if i add the 0.56 shift here the distribution is not centered at 0 for the yaxis
					float pos_stereo_x = ((cl_lay[3]->getCorrPosition(ind_lay3)-cl_lay[2]->getCorrPosition(ind_lay2)) / 2*TMath::Sin(1.5*TMath::Pi()/180.)); // if i add the 0.56 shift here the distribution is not centered at 0 for the yaxis
					histos->h_beamProfile_ontrack->Fill(pos_stereo_y, pos_stereo_x);
					float pos_stereo_z = (cl_lay[3]->getZPosition() + cl_lay[2]->getZPosition() ) / 2;
					track->checkDistanceFromTrack(0, 2, cl_lay[2], ind_lay2, stereo_in_cut_eff);
					track->checkDistanceFromTrack(0, 3, cl_lay[3], ind_lay3, stereo_out_cut_eff);
					track->checkStereoDistanceFromTrack(0, pos_stereo_y, pos_stereo_z, stereo_cut_eff);
					track->fillStereoCluster(cl_lay[2], cl_lay[3], ind_lay2, ind_lay3, 0);  // fill the graph with the stereo point if both stereo layers give a hit
				}
				
				if(track->getTrack()->GetN()>3) track->setSM1errors(0);
				if(count_layers_with_cluster==0) count_SM1_zero_hits++;
				
				if(track->acceptEtaOut()) count_SM1_eta_out++;
				if(track->acceptEtaIn()) count_SM1_eta_in++;
				if(track->acceptStereo()) count_SM1_stereo++;
				if(track->acceptStereoIn() && track->acceptStereoOut()) count_SM1_stereo_both++;
				if(track->acceptStereoIn()) count_SM1_stereo_in++;
				if(track->acceptStereoOut()) count_SM1_stereo_out++;
		
				if(count_events < 101 && count_layers_with_cluster>0)
				{
					out_file->cd("track_candidates/track_displays");
					gr_track->SetTitle(Form("Event-%i",iEvent));
					gr_track->Write();
				}
				for(int i=0; i<4; i++)
				{
					v_dclus_to_track[i].clear();
					v_dclus_to_track[i].shrink_to_fit();
				}
				v_dclus_to_track_stereo.clear();
				v_dclus_to_track_stereo.shrink_to_fit();
				v_index_lay2.clear();
				v_index_lay2.shrink_to_fit();
				v_index_lay3.clear();
				v_index_lay3.shrink_to_fit();

				histos->h_nlayers_with_cluster->Fill(count_layers_with_cluster);
				histos->h_singletrack_events->Fill(1);
			}

			
			if(track->tooMany() && !track->getCandidateTracks().empty())
			{
				v_cand_tracks = track->getCandidateTracks();
				std::vector<TGraphAsymmErrors*> v_track_graphs;
			
				counter4 += v_cand_tracks.size(); // multiple tracks  
				if(v_cand_tracks.size()==1) counter_mult_single_tracks++; // the fit works for one track among the many 
				
				std::vector<std::vector<float> > v_track_info_ypos;
				std::vector<std::vector<float> > v_track_info_layer;
				std::vector<std::vector<float> > v_track_info_xpos;
				std::vector<std::vector<int> > v_track_info_cl_indices;
				std::vector<float> ypositions;
				std::vector<float> layers; 
				std::vector<float> xpositions;
				std::vector<int> cluster_indices;
				
				for(int itrack=0; itrack<v_cand_tracks.size(); itrack++) // loop over tracks
				{
					//if( !(track->getAngle(itrack)>=-0.5 && track->getAngle(itrack)<=0.5)){
					// counter4--;
					// continue;
					//	}
					if(v_cand_tracks.size()>1)
					{
						count_mult_track_events++;
						histos->h_mult_track_angle->Fill(track->getAngle(itrack));
						histos->h_prob_mult_tracks->Fill(track->getTrackProb(itrack));
						//std::cout<<"Track-"<<itrack<<" "<<track->getTrackProb(itrack)<<" "<<track->getAngle(itrack)<<std::endl; 
					}
					int count_layers_with_cluster = 0;
					int count_mult_track_nlayer =0;
					TGraphAsymmErrors* gr_mult_track = v_cand_tracks.at(itrack);
			 		gr_mult_track->SetName(Form("gr_track_%i_%i",itrack,iEvent));
			 		float tr_slope = track->trackSlope(itrack);
					float tr_const = track->trackIntercept(itrack);

					count_events++;
					histos->h_residuals_mult_tracks->Fill(gr_mult_track->GetPointY(1) - gr_mult_track->GetPointY(2));
			 		histos->h_residuals_accepted_tracks->Fill(gr_mult_track->GetPointY(1) - gr_mult_track->GetPointY(2));
			 		std::vector<float> v_dclus_to_track[4]; // filled with all distances between track and each cluster 
						
					for(int iLayer = 0; iLayer<4; iLayer++) 
					{
						if(v_cand_tracks.size()==1 )
						{
							track->isSingleTrack();
							if(cl_lay[iLayer]->getNClusters2()>0) count_layers_with_cluster++;
							histos->h_nclusters_per_layer_event[iLayer]->Fill(cl_lay[iLayer]->getNClusters2());
						} 
					//	if(v_cand_tracks.size()>1)
					//	{
					//		if(cl_lay[iLayer]->getNClusters2()>0) count_mult_track_nlayer++;
					//	}
					//	histos->h_nclusters_per_layer_event[iLayer]->Fill(cl_lay[iLayer]->getNClusters2());
						// find the clusters that belong on the track
						// first get track params
						for(int icl=0; icl<cl_lay[iLayer]->getNClusters2(); icl++)
						{
							float distance = std::abs(tr_slope*cl_lay[iLayer]->getZPosition() - cl_lay[iLayer]->getCorrPosition(icl) + tr_const)/(TMath::Sqrt(tr_slope*tr_slope + 1));
							v_dclus_to_track[iLayer].push_back(distance);
						}
						
					}
					
					std::vector<int> v_index_lay2;
					std::vector<int> v_index_lay3;
					std::vector<float> v_dclus_to_track_stereo;
					if(cl_lay[2]->getNClusters2()>0 && cl_lay[3]->getNClusters2()>0)
					{
						for(int icl2 =0; icl2<cl_lay[2]->getNClusters2(); icl2++)
						{
							float pos_cl2 = cl_lay[2]->getCorrPosition(icl2);
							for(int icl3=0; icl3<cl_lay[3]->getNClusters2(); icl3++)
							{
								v_index_lay2.push_back(icl2);
								v_index_lay3.push_back(icl3);
								float pos_cl3 = cl_lay[3]->getCorrPosition(icl3);
								float pos_stereo_y = 0.56 + ((pos_cl3+pos_cl2) / 2*TMath::Cos(1.5*TMath::Pi()/180.));
								float pos_stereo_x = ((pos_cl3-pos_cl2) / 2*TMath::Sin(1.5*TMath::Pi()/180.));	
								float distance = std::abs(tr_slope*(cl_lay[2]->getZPosition()+cl_lay[3]->getZPosition())/2 - pos_stereo_y + tr_const)/(TMath::Sqrt(tr_slope*tr_slope + 1));
								v_dclus_to_track_stereo.push_back(distance);
							}
						}
					}
					for(int iLayer=0; iLayer<4; iLayer++) // eta layers
					{
						if(v_dclus_to_track[iLayer].empty()) continue;
						float min_distance = *min_element(v_dclus_to_track[iLayer].begin(), v_dclus_to_track[iLayer].end());
						int index = find(v_dclus_to_track[iLayer].begin(), v_dclus_to_track[iLayer].end(),min_distance) - v_dclus_to_track[iLayer].begin(); // index of the cluster which belongs to the track
						if(iLayer<2) {
						 track->fillSM1events(index, cl_lay[iLayer], iLayer, itrack);	// eta layers
						 if(!track->checkIfSingleTrack())
						 {
						 	ypositions.push_back(cl_lay[iLayer]->getCorrPosition(index));
						 	xpositions.push_back(99999);
						 	layers.push_back(iLayer);
						 	cluster_indices.push_back(index); 
						 }
						 
						}
						else
						{
							if(v_dclus_to_track_stereo.empty())
							{
								if(iLayer==2) 
								{
									if(!track->checkIfSingleTrack())
						 			{
						 				ypositions.push_back(cl_lay[iLayer]->getCorrPosition(index));
						 				xpositions.push_back(99999);
						 				layers.push_back(iLayer);
						 				cluster_indices.push_back(index);  
						 			}
									track->checkDistanceFromTrack(itrack, iLayer, cl_lay[iLayer], index, stereo_in_cut_eff);
								}
								if(iLayer==3)
								{
									if(!track->checkIfSingleTrack())
						 			{
						 				ypositions.push_back(cl_lay[iLayer]->getCorrPosition(index));
						 				xpositions.push_back(99999);
						 				layers.push_back(iLayer); 
						 				cluster_indices.push_back(index); 
						 			}
									track->checkDistanceFromTrack(itrack, iLayer, cl_lay[iLayer], index, stereo_out_cut_eff);
								}
							}
						}		
					}

					if(!v_dclus_to_track_stereo.empty())
					{
						if(track->checkIfSingleTrack()) count_events_hit_on_both_stereo++;
						// find the closest stereo cluster
						float min_distance = *min_element(v_dclus_to_track_stereo.begin(), v_dclus_to_track_stereo.end());
						int index = find(v_dclus_to_track_stereo.begin(), v_dclus_to_track_stereo.end(), min_distance) - v_dclus_to_track_stereo.begin();
						int ind_lay2 = v_index_lay2.at(index); // the cluster index on layer 2
						int ind_lay3 = v_index_lay3.at(index); // the cluster index on layer 3
						
						float pos_stereo_y = 0.56 + ((cl_lay[3]->getCorrPosition(ind_lay3)+cl_lay[2]->getCorrPosition(ind_lay2)) / 2*TMath::Cos(1.5*TMath::Pi()/180.)); 
						float pos_stereo_x = ((cl_lay[3]->getCorrPosition(ind_lay3)-cl_lay[2]->getCorrPosition(ind_lay2)) / 2*TMath::Sin(1.5*TMath::Pi()/180.)); 
						float pos_stereo_z = (cl_lay[3]->getZPosition() + cl_lay[2]->getZPosition() ) / 2;
						if(track->checkIfSingleTrack()) 
						{ 
							histos->h_d_track_ystereo->Fill(min_distance);
							histos->h_beamProfile_ontrack->Fill(pos_stereo_y, pos_stereo_x);
							histos->h_d_track_lay2->Fill(track->Get_Distance_from_Track(cl_lay[2]->getZPosition(), cl_lay[2]->getCorrPosition(ind_lay2), itrack));
							histos->h_d_track_lay3->Fill(track->Get_Distance_from_Track(cl_lay[3]->getZPosition(), cl_lay[3]->getCorrPosition(ind_lay3), itrack));
							float st_dist_from_track = track->Get_Distance_from_Track( (cl_lay[2]->getZPosition()+cl_lay[3]->getZPosition())/2, pos_stereo_y, itrack);
							histos->h_d_track_stereo->Fill(st_dist_from_track);
							if(st_dist_from_track<=stereo_cut_eff && st_dist_from_track>=-stereo_cut_eff)
								histos->h_d_track_stereo_cut->Fill(st_dist_from_track);
						}
						if(!track->checkIfSingleTrack())
						{
						 	ypositions.push_back(pos_stereo_y);
						 	xpositions.push_back(pos_stereo_x);
						 	layers.push_back(-1); 
						}
						track->checkDistanceFromTrack(itrack, 2, cl_lay[2], ind_lay2, stereo_in_cut_eff);
						track->checkDistanceFromTrack(itrack, 3, cl_lay[3], ind_lay3, stereo_out_cut_eff);
						track->checkStereoDistanceFromTrack(itrack, pos_stereo_y, pos_stereo_z, stereo_cut_eff);
						track->fillStereoCluster(cl_lay[2], cl_lay[3], ind_lay2, ind_lay3, itrack);
					}
					
					if(track->getTrack(itrack)->GetN()>3) track->setSM1errors(itrack);
					
					if(v_cand_tracks.size()>1) 
					{
						v_track_info_ypos.push_back(ypositions);
						v_track_info_layer.push_back(layers);
						v_track_info_xpos.push_back(xpositions);
						v_track_info_cl_indices.push_back(cluster_indices);
					}
					
					if(!track->checkIfSingleTrack())
			 		{
			 			v_track_graphs.push_back(gr_mult_track);
			 		}

					if(count_layers_with_cluster==0 && track->checkIfSingleTrack()) count_SM1_zero_hits++;
				
					if(track->checkIfSingleTrack()) counter3++;
					if(track->acceptEtaOut()&& track->checkIfSingleTrack()) count_SM1_eta_out++;
					if(track->acceptEtaIn()&& track->checkIfSingleTrack()) count_SM1_eta_in++;
					if(track->acceptStereo()&& track->checkIfSingleTrack()) count_SM1_stereo++;
					if(track->acceptStereoIn() && track->acceptStereoOut()&& track->checkIfSingleTrack()) count_SM1_stereo_both++;
					if(track->acceptStereoIn() && track->checkIfSingleTrack()) count_SM1_stereo_in++;
					if(track->acceptStereoOut()&& track->checkIfSingleTrack()) count_SM1_stereo_out++;
				
					if(count_events < 101  && count_layers_with_cluster>0)
			 		{
			 			out_file->cd("track_candidates/track_displays");
						gr_mult_track->SetTitle(Form("Event-%i",iEvent));
						gr_mult_track->Write();
			 		}

			 		if(!track->checkIfSingleTrack())
			 		{
			 			if(count_mult_track_events < 101)
			 			{
			 				out_file->cd("track_candidates/mult_tracks");
			 				gr_mult_track->SetTitle(Form("Event-%i",iEvent));
			 				gr_mult_track->Write();
			 			}
			 			track->setCriterium();
			 			track->refreshTrack(); // set layer flags back to false for the new track
			 		}

			 		for(int i=0; i<4; i++)
					{
						v_dclus_to_track[i].clear();
						v_dclus_to_track[i].shrink_to_fit();
					}
					
					v_dclus_to_track_stereo.clear();
					v_dclus_to_track_stereo.shrink_to_fit();
					v_index_lay2.clear();
					v_index_lay2.shrink_to_fit();
					v_index_lay3.clear();
					v_index_lay3.shrink_to_fit();
					ypositions.clear();
					ypositions.shrink_to_fit();
					xpositions.clear();
					xpositions.shrink_to_fit();
					layers.clear();
					layers.shrink_to_fit();
					cluster_indices.clear();
					cluster_indices.shrink_to_fit();

					if(track->checkIfSingleTrack()) histos->h_nlayers_with_cluster->Fill(count_layers_with_cluster);
				}

				if(v_cand_tracks.size()>1) // select the optimal track among the many candidates
				{
					counter3++;
					std::vector<float> v_good_track = track->getSelectionCriteria();
					
					if(!v_good_track.empty()) {
					float min_elem = *min_element( v_good_track.begin(), v_good_track.end());
					int index = find(v_good_track.begin(), v_good_track.end(), min_elem) - v_good_track.begin();

					for(int iLayer=0; iLayer<4; iLayer++)
					{
						histos->h_nclusters_per_layer_event[iLayer]->Fill(cl_lay[iLayer]->getNClusters2());
					}

					int count_layers_with_cluster_mult_track = v_track_info_layer[index].size();

					if( std::find(v_track_info_layer[index].begin(), v_track_info_layer[index].end(), -1) != v_track_info_layer[index].end() )
						count_layers_with_cluster_mult_track++;

					histos->h_nlayers_with_cluster->Fill(count_layers_with_cluster_mult_track);
					if(count_layers_with_cluster_mult_track==0) count_SM1_zero_hits++;
					bool stereoin = false;
					bool stereout = false;
					for(int i=0; i<v_track_info_layer[index].size(); i++)
					{
						float distance_from_track = 9999;
						if(v_track_info_layer[index].at(i)>-1) distance_from_track = track->Get_Distance_from_Track(cl_lay[(int)v_track_info_layer[index].at(i)]->getZPosition(), v_track_info_ypos[index].at(i), index);
						else distance_from_track = track->Get_Distance_from_Track( (cl_lay[2]->getZPosition() + cl_lay[3]->getZPosition())/2, v_track_info_ypos[index].at(i), index);
						
						if(v_track_info_layer[index].at(i)==0) // eta_out
						{
							histos->h_d_track_etaout->Fill(distance_from_track);
							if(distance_from_track<=eta_out_cut_eff && distance_from_track>=-eta_out_cut_eff)
							{
								histos->h_d_track_etaout_cut->Fill(distance_from_track);
								histos->h_cl_charge_on_track[(int)v_track_info_layer[index].at(i)]->Fill(cl_lay[(int)v_track_info_layer[index].at(i)]->getTotPdo(v_track_info_cl_indices[index].at(i)));
								histos->h_nstrips_on_track[(int)v_track_info_layer[index].at(i)]->Fill(cl_lay[(int)v_track_info_layer[index].at(i)]->getNStrips(v_track_info_cl_indices[index].at(i)));
								count_SM1_eta_out++;
							}
						}
						if(v_track_info_layer[index].at(i)==1) // eta_in
						{
							histos->h_d_track_etain->Fill(distance_from_track);
							if(distance_from_track<=eta_in_cut_eff && distance_from_track>=-eta_in_cut_eff)
							{
								histos->h_d_track_etain_cut->Fill(distance_from_track);
								histos->h_cl_charge_on_track[(int)v_track_info_layer[index].at(i)]->Fill(cl_lay[(int)v_track_info_layer[index].at(i)]->getTotPdo(v_track_info_cl_indices[index].at(i)));
								histos->h_nstrips_on_track[(int)v_track_info_layer[index].at(i)]->Fill(cl_lay[(int)v_track_info_layer[index].at(i)]->getNStrips(v_track_info_cl_indices[index].at(i)));
								count_SM1_eta_in++;
							}
						}
						if(v_track_info_layer[index].at(i)==2) // stereo_in
						{
							histos->h_d_track_lay2->Fill(distance_from_track);
							if(distance_from_track<=stereo_in_cut_eff && distance_from_track>=-stereo_in_cut_eff)
							{
								histos->h_d_track_lay2_cut->Fill(distance_from_track);
								histos->h_cl_charge_on_track[(int)v_track_info_layer[index].at(i)]->Fill(cl_lay[(int)v_track_info_layer[index].at(i)]->getTotPdo(v_track_info_cl_indices[index].at(i)));
								histos->h_nstrips_on_track[(int)v_track_info_layer[index].at(i)]->Fill(cl_lay[(int)v_track_info_layer[index].at(i)]->getNStrips(v_track_info_cl_indices[index].at(i)));
								count_SM1_stereo_in++;
								stereoin = true;
							}
						}
						if(v_track_info_layer[index].at(i)==3) // stereo_out
						{
							histos->h_d_track_lay3->Fill(distance_from_track);
							histos->h_d_track_ystereo->Fill(std::abs(distance_from_track));
							if(distance_from_track<=stereo_out_cut_eff && distance_from_track>=-stereo_out_cut_eff)
							{
								histos->h_d_track_lay3_cut->Fill(distance_from_track);
								histos->h_cl_charge_on_track[(int)v_track_info_layer[index].at(i)]->Fill(cl_lay[(int)v_track_info_layer[index].at(i)]->getTotPdo(v_track_info_cl_indices[index].at(i)));
								histos->h_nstrips_on_track[(int)v_track_info_layer[index].at(i)]->Fill(cl_lay[(int)v_track_info_layer[index].at(i)]->getNStrips(v_track_info_cl_indices[index].at(i)));
								count_SM1_stereo_out++;
								stereout = true;
							}
						}
						if(v_track_info_layer[index].at(i)==-1) // stereo
						{
							//std::cout<<distance_from_track<<std::endl;
							count_events_hit_on_both_stereo++;
							histos->h_d_track_stereo->Fill(distance_from_track);
							histos->h_d_track_ystereo->Fill(std::abs(distance_from_track));
							if(distance_from_track<=stereo_cut_eff && distance_from_track>=-stereo_cut_eff)
							{
								//std::cout<<distance_from_track<<std::endl;
								histos->h_d_track_stereo_cut->Fill(distance_from_track);
								histos->h_beamProfile_ontrack->Fill(v_track_info_ypos[index].at(i), v_track_info_xpos[index].at(i));
								count_SM1_stereo++;
							}
						}
					}
					if(stereoin && stereout) count_SM1_stereo_both++;

					if(count_mult_track_events < 101 ) 
					{
						out_file->cd("track_candidates/clean_mult_tracks_morechoices");
			 			v_track_graphs.at(index)->Write();
					}
				}
				else if(v_good_track.empty()) { // no clusters on any SM1 layers
					histos->h_nlayers_with_cluster->Fill(0); 
					for(int iLayer=0; iLayer<4; iLayer++)
					{
						histos->h_nclusters_per_layer_event[iLayer]->Fill(cl_lay[iLayer]->getNClusters2());
					}
				} 		
				
				v_good_track.clear();
				v_good_track.shrink_to_fit();

			} 
				
				for(int i=0; i<v_track_info_ypos.size(); i++)
				{
					v_track_info_ypos[i].clear();
					v_track_info_ypos[i].shrink_to_fit();
					v_track_info_xpos[i].clear();
					v_track_info_xpos[i].shrink_to_fit();
					v_track_info_layer[i].clear();
					v_track_info_layer[i].shrink_to_fit();	
					v_track_info_cl_indices[i].clear();
					v_track_info_cl_indices[i].shrink_to_fit();				
				}

				v_track_graphs.clear();
				v_track_graphs.shrink_to_fit(); 
			}
			count_accepted_tracks += track->countAcceptedTracks();
			count_accepted_singletracks_outside_1deg += track->getAcceptedTracksTheOnlySuperClusterOutside1DegWindow();
			count_accepted_singletracks_inside_1deg += track->getAcceptedSingleTracksInside1DegWindow();
			count_accepted_multipletracks_inside_1deg += track->getAcceptedTracksManySuperClusterInside1DegWindow();
			
			histos->h_Nscluster_in_1deg_sby2_sby3->Fill(track->getSuperClustersWithin1DegWindow());
			histos->h_Nscluster_theonly_out_1deg_sby2_sby3->Fill(track->getTheOnlySuperClusterOutside1DegWindow());
			histos->h_Nscluster_morethan1_in_1deg_sby2_sby3->Fill(track->getEventsWithMoreThan1SuperClusterIn1DegWindow());
			histos->h_Nscluster_morethan1_out_1deg_sby2_sby3->Fill(track->getEventsWithMoreThan1SuperClusterOut1DegWindow());

			v_cand_tracks.clear();
			v_cand_tracks.shrink_to_fit();

			delete track;
		}

// ========================= End tracking ==============================================

		for(int iLayer=0; iLayer<4; iLayer++) 
		{
			delete layer[iLayer];
			delete cl_lay[iLayer];
			delete l_small[iLayer];
			delete cl_lay_small[iLayer];
		}

		delete trigger;

		for(int i=0; i<sizeof(v_strips_layer)/sizeof(v_strips_layer[0]); i++)
		{
			v_strips_layer[i].clear();
			v_strips_layer[i].shrink_to_fit();
			v_strips_layer_small[i].clear();
			v_strips_layer_small[i].shrink_to_fit();
		}
		v_strips_trigger.clear();
		v_strips_trigger.shrink_to_fit();

	} // end event for-loop


	out_file->cd();
	out_file->mkdir("trigger");
	out_file->mkdir("SM1/raw_hits");
	out_file->mkdir("SM1/clusters");
	out_file->mkdir("SM1/strips");
	out_file->mkdir("SB/raw_hits");
	out_file->mkdir("SB/clusters");
	out_file->mkdir("SB/strips");
	out_file->mkdir("time_resolution");
	out_file->mkdir("alignment");
	out_file->cd();

	histos->h_cutflow->SetBinContent(1, counter0);
	histos->h_cutflow->SetBinContent(2, counter1); // 1 cluster on SBY1
	histos->h_cutflow->SetBinContent(3, count_accepted_tracks); // tracks found
	histos->h_cutflow->SetBinContent(4, (counter3) ); // tracks used

	histos->h_cutflow->GetXaxis()->SetBinLabel(1, "total events");
	histos->h_cutflow->GetXaxis()->SetBinLabel(2, "1 clus SBY1");
	histos->h_cutflow->GetXaxis()->SetBinLabel(3, "N candidate tracks");
	histos->h_cutflow->GetXaxis()->SetBinLabel(4, "N reference tracks");


	for(int iLayer=0; iLayer<4; iLayer++) 
	{
		out_file->cd("SM1/raw_hits");
		histos->h_raw_hits[iLayer]->Write();
		histos->h_strip_index_SM1[iLayer]->Write();
		histos->h_strip_index_SM1_cluster[iLayer]->Write();

		histos->h_nstrips[iLayer]->SetBinContent(2,single_strip_counter[iLayer]);
		out_file->cd("SM1/clusters");
		histos->h_nclusters[iLayer]->Write();
		histos->h_nstrips[iLayer]->Write();
		histos->h_strip_pdo_in_clus[iLayer]->Write();
		histos->h_strip_relbcid_in_clus[iLayer]->Write();
		histos->h_strip_tdo_in_clus[iLayer]->Write();
		histos->h_cl_charge[iLayer]->Write();
		histos->h_cl_charge_leadCluster[iLayer]->Write();
		histos->h_cl_charge_leadCluster_maxstrip[iLayer]->Write();
		histos->h_strip_pdo_in_leadCluster[iLayer]->Write();
		histos->h_maxstrip_pdo_allClus[iLayer]->Write();
		histos->h_clus_positions[iLayer]->Write();
		histos->h_strip_index_vs_pdo[iLayer]->Write();
		histos->h_raw_hits_vs_tot_strips[iLayer]->Write();
		histos->h_cl_size_vs_cl_charge[iLayer]->Write();
	#ifdef TIME_CALIBRATION
		histos->h_strip_time_in_clus[iLayer]->Write();
		histos->h_strip_time_in_clus_earliest[iLayer]->Write();
	#endif
		out_file->cd("SM1/strips");
		histos->h_strip_index_vs_tdo[iLayer]->Write();
		histos->h_strip_index_vs_relbcid[iLayer]->Write();
		histos->h_strip_index_vs_bcid[iLayer]->Write();
		histos->h_strip_index_vs_pdo0[iLayer]->Write();
		out_file->cd("track_candidates/SM1");
		histos->h_nclusters_per_layer_event[iLayer]->Write();
		histos->h_nstrips_on_track[iLayer]->Write();
		histos->h_cl_charge_on_track[iLayer]->Write();
	}

	
	out_file->cd("track_candidates/residuals");
	histos->h_residuals_single_tracks->Write();
	histos->h_residuals_mult_tracks->Write();
	histos->h_residuals_accepted_tracks->Write();
	
	out_file->cd("track_candidates");
	histos->h_angle_only_scluster_sby2->Write();
	histos->h_angle_only_scluster_sby3->Write();
	histos->h_dangle_only_scluster_sby2sby3->Write();
	histos->h_angle->Write();
	histos->h_chi2->Write();
	histos->h_chi2ndf->Write();
	histos->h_prob->Write();
	histos->h_nlayers_with_cluster->Write();
	histos->h_d_track_etaout->Write();
	histos->h_d_track_etaout_cut->Write();
	histos->h_d_track_etain->Write();
	histos->h_d_track_etain_cut->Write();
	histos->h_d_track_stereo->Write();
	histos->h_d_track_stereo_cut->Write();
	histos->h_d_track_ystereo->Write();
	histos->h_mult_track_angle->Write();
	histos->h_nlayers_mult_track_events->Write();
	histos->h_prob_mult_tracks->Write();

	histos->h_d_track_lay2->Write();
	histos->h_d_track_lay3->Write();
	histos->h_d_track_lay2_cut->Write();
	histos->h_d_track_lay3_cut->Write();
	histos->h_beamProfile_ontrack->Write();


	out_file->cd("alignment");
	histos->h_res_SBY2_SBY3_vs_SBY2->Write(); 
	histos->h_res_SBY2_SBY3_vs_SBY3->Write();
	histos->h_res_SBY2_SBY1_vs_SBY1->Write();
	histos->h_res_SBY3_SBY1_vs_SBY1->Write();
	histos->h_sby1_minus_sby2_vs_sby2->Write();
	histos->h_sby1_minus_sby3_vs_sby3->Write();
	histos->h_pos_etain_vs_pos_etaout->Write();
	histos->h_diffpos_lay0->Write();
	histos->h_diffpos_lay1->Write();
	histos->h_diffpos_stereolay1->Write();
	histos->h_diffpos_stereolay0->Write();
	histos->h_sby1_minus_eta_out_vs_pos_eta_out->Write();
	histos->h_sby1_minus_pos_eta_in_vs_pos_eta_in->Write();
	histos->h_sby1_minus_stereo_vs_stereo->Write();
	histos->h_sby1_minus_stereo_in_vs_stereo_in->Write();
	histos->h_sby1_minus_stereo_out_vs_stereo_out->Write();	
	histos->h_beamProfile->Write();

	for(int iLayer=0; iLayer<4; iLayer++) 
	{
		out_file->cd("SB/raw_hits");
		histos->h_raw_hits_small[iLayer]->Write();
		histos->h_strip_index_SB[iLayer]->Write();
		histos->h_strip_index_SB_cluster[iLayer]->Write();

		out_file->cd("SB/clusters");
		histos->h_nclusters_small[iLayer]->Write();
		histos->h_nstrips_small[iLayer]->Write();
		histos->h_strip_pdo_in_clus_small[iLayer]->Write();
		histos->h_strip_relbcid_in_clus_small[iLayer]->Write();
		histos->h_strip_tdo_in_clus_small[iLayer]->Write();
		histos->h_cl_charge_small[iLayer]->Write();
		histos->h_clus_positions_small[iLayer]->Write();
		histos->h_strip_index_vs_pdo_small[iLayer]->Write();
		histos->h_raw_hits_vs_tot_strips_small[iLayer]->Write();
		histos->h_cl_size_vs_cl_charge_small[iLayer]->Write();
		histos->h_strip_2044_pdo_SB[iLayer]->Write();
		histos->h_strip_2044_pdo_SM1[iLayer]->Write();
	#ifdef TIME_CALIBRATION
		histos->h_strip_time_in_clus_small[iLayer]->Write();
		histos->h_strip_time_in_clus_earliest_small[iLayer]->Write();
	#endif
		out_file->cd("SB/strips");
		histos->h_strip_index_vs_tdo_small[iLayer]->Write();
		histos->h_strip_index_vs_relbcid_small[iLayer]->Write();
		histos->h_strip_index_vs_bcid_small[iLayer]->Write();
		histos->h_strip_index_vs_pdo0_small[iLayer]->Write();
	}
	out_file->cd("trigger");
	histos->h_strip_index_vs_tdo_trigger->Write();
	histos->h_strip_index_vs_relbcid_trigger->Write();
	histos->h_strip_index_vs_bcid_trigger->Write();
	histos->h_strip_index_vs_pdo0_trigger->Write();
	histos->h_raw_hit_trigger->Write();

	out_file->cd();
	histos->h_cutflow->Write();

#ifdef TIME_CALIBRATION
	out_file->cd("time_resolution");
	histos->h_time_res_SB01->Write();
	histos->h_time_res_SB12->Write();
	histos->h_time_res_SB13->Write();
	histos->h_time_res_SB1SM1lay0->Write();
	histos->h_time_res_SB1SM1lay1->Write();
	histos->h_time_res_SB1SM1lay2->Write();
	histos->h_time_res_SB1SM1lay3->Write();
#endif

	std::cout<<"N single tracks "<<counter3<<std::endl;
	std::cout<<"Eta out (accepted) "<<count_SM1_eta_out<<std::endl;
	std::cout<<"Eta in (accepted) "<<count_SM1_eta_in<<std::endl;
	std::cout<<"Stereo in (accepted) "<<count_SM1_stereo_in<<std::endl;
	std::cout<<"Stereo out (accepted) "<<count_SM1_stereo_out<<std::endl;
	std::cout<<"Stereo (combined) (accepted) "<<count_SM1_stereo<<std::endl;
	std::cout<<"stereo both (accepted) "<<count_SM1_stereo_both<<std::endl;
	std::cout<<"ntracks with 0 hits on SM1 "<<count_SM1_zero_hits<<std::endl;

	out_file->Close();

	delete histos;
	delete treeReader;

	return 0;
}


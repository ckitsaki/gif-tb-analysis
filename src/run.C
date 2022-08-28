
//#define TIME_CALIBRATION
#include "TreeReader.C"
#include "Cluster.h"
#include "Histograms.h"

int run(std::string run_number, std::string sector="C14")
{

	gROOT->SetBatch(kTRUE); 

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
	out_file->mkdir("event_displays/SM1");
	out_file->mkdir("event_displays/SB");
	out_file->mkdir("event_displays/all");
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
	
// Counters for cutflow	
	int counter0=0; //total number of events
	int counter1=0; //select events with hits only on specific radii
///////////////////////////////////////////////////////

	int i=1;
	int divEvents = tree->GetEntries()/8;

	// loop over events
	for(int iEvent = 0; iEvent< Nentries; iEvent++)
	{
		gROOT->cd();	
		bool candidate = false;

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
		
		for(int iLayer = 0; iLayer<sizeof(layer)/sizeof(layer[0]); iLayer++) {
			histos->h_raw_hits[iLayer]->Fill(layer[iLayer]->getNHits()); // all hits before any rejection
			histos->h_raw_hits_small[iLayer]->Fill(l_small[iLayer]->getNHits());

			layer[iLayer]->bookFiredStrips(strips,pdo); // here we reject strips with pdo<64 and we also mask the noisy strips
			l_small[iLayer]->bookFiredStrips(strips,pdo);

			v_strips_layer[iLayer] = layer[iLayer]->getFiredStrips();
			v_strips_layer_small[iLayer] = l_small[iLayer]->getFiredStrips();
			
			for(int i=0; i<v_strips_layer_small[iLayer].size(); i++)
			{
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
				if(cl_lay_small[iLayer]->getNStrips(icl)==0) continue;
				histos->h_nstrips_small[iLayer]->Fill(cl_lay_small[iLayer]->getNStrips(icl));
				histos->h_cl_charge_small[iLayer]->Fill(cl_lay_small[iLayer]->getTotPdo(icl));
				histos->h_clus_positions_small[iLayer]->Fill(cl_lay_small[iLayer]->getPosition(icl));
				histos->h_cl_size_vs_cl_charge_small[iLayer]->Fill(cl_lay_small[iLayer]->getNStrips(icl), cl_lay_small[iLayer]->getTotPdo(icl));

				tot_strips_small[iLayer] += cl_lay_small[iLayer]->getNStrips(icl); // total number of strips for all clusters formed in the event

				
			for(int istrip = 0; istrip < cl_lay_small[iLayer]->getNStrips(icl); istrip++)
				{
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
			}

			if(cl_lay_small[iLayer]->getNClusters2()>0){
				histos->h_raw_hits_vs_tot_strips_small[iLayer]->Fill(tot_strips_small[iLayer],l_small[iLayer]->getNHits());
				single_strip_counter_small[iLayer] += cl_lay_small[iLayer]->getNSingleStrips();
			}
		} //end ilayer for-loop

//  End Clustering

//  checks

		if( !(cl_lay_small[1]->getNClusters2()==1)) continue; // only events which give exactly one cluster on SBY1 are processed		
		counter1++; // count events with exactly one cluster on SBY1 

		
		layer[0]->setAlphaBeta(-409.859, (1-0.0297));
		layer[1]->setAlphaBeta(-409.859, (1-0.0297));
		layer[2]->setAlphaBeta(-409.859, (1-0.0297));
		layer[3]->setAlphaBeta(-409.859, (1-0.0297)); 

		l_small[2]->setAlphaBeta(33.000, (1-0.0448)); 
		l_small[3]->setAlphaBeta(35.778, (1-0.0467));

		float pos_sby1 = cl_lay_small[1]->getPosition(0); // SBY1 reference chamber position

		if(cl_lay_small[2]->getNClusters2()==1 && cl_lay_small[3]->getNClusters2()==1)
		{
			histos->h_res_SBY2_SBY3_vs_SBY2->Fill(cl_lay_small[2]->getCorrPosition(0)-cl_lay_small[3]->getCorrPosition(0), cl_lay_small[2]->getCorrPosition(0));
			histos->h_res_SBY2_SBY3_vs_SBY3->Fill(cl_lay_small[2]->getCorrPosition(0)-cl_lay_small[3]->getCorrPosition(0), cl_lay_small[3]->getCorrPosition(0));
		}

		if(cl_lay[0]->getNClusters2()==1 && cl_lay[1]->getNClusters2()==1)
		{
			histos->h_diffpos_lay0->Fill(cl_lay[1]->getCorrPosition(0)-cl_lay[0]->getCorrPosition(0), cl_lay[0]->getCorrPosition(0));
			histos->h_diffpos_lay1->Fill(cl_lay[1]->getCorrPosition(0)-cl_lay[0]->getCorrPosition(0), cl_lay[1]->getCorrPosition(0));
		}

		if(cl_lay[2]->getNClusters2()==1 && cl_lay[3]->getNClusters2()==1)
		{
			float pos_l2 = cl_lay[2]->getCorrPosition(0);
			float pos_l3 = cl_lay[3]->getCorrPosition(0);

			float pos_l2_uncorr = cl_lay[2]->getPosition(0);
			float pos_l3_uncorr = cl_lay[3]->getPosition(0);

			float pos_stereo_y = ((pos_l3+pos_l2) / 2*TMath::Cos(1.5*TMath::Pi()/180.));
			float pos_stereo_x = ((pos_l3-pos_l2) / 2*TMath::Sin(1.5*TMath::Pi()/180.));
			float pos_stereo_y_uncorr = ((pos_l3_uncorr+pos_l2_uncorr) / 2*TMath::Cos(1.5*TMath::Pi()/180.));
			float pos_stereo_x_uncorr = ((pos_l3_uncorr-pos_l2_uncorr) / 2*TMath::Sin(1.5*TMath::Pi()/180.));

			histos->h_stereo_y_corr_vs_uncorr->Fill(pos_stereo_y, pos_stereo_y_uncorr);
			histos->h_stereo_x_corr_vs_uncorr->Fill(pos_stereo_x, pos_stereo_x_uncorr);

			if(cl_lay[0]->getNClusters2()>0) {
			for(int icl=0; icl<cl_lay[0]->getNClusters2(); icl++)
			{
				float pos_l0 = cl_lay[0]->getCorrPosition(icl);
				histos->h_diffpos_stereolay0->Fill(pos_stereo_y-pos_l0, pos_stereo_y);
				histos->h_bsy1_minus_etaout_vs_stereo_x->Fill(pos_stereo_x,pos_sby1-pos_l0);
			}

			}
			if(cl_lay[1]->getNClusters2()>0) {
				for(int icl=0; icl<cl_lay[1]->getNClusters2(); icl++)
				{
					float pos_l1 = cl_lay[1]->getCorrPosition(icl);

					histos->h_diffpos_stereolay1->Fill(pos_stereo_y-pos_l1, pos_stereo_y);
					histos->h_bsy1_minus_etain_vs_stereo_x->Fill(pos_stereo_x, pos_sby1-pos_l1);

					histos->h_diffpos_stereo_uncorr_lay1->Fill(pos_stereo_y_uncorr, pos_stereo_y_uncorr-pos_l1);
				}
			}
			
			histos->h_bsy1_minus_stereoin_vs_stereo_x->Fill(pos_stereo_x, pos_sby1-pos_l2);
			histos->h_bsy1_minus_stereoout_vs_stereo_x->Fill(pos_stereo_x, pos_sby1-pos_l3);
			
			histos->h_beamProfile->Fill(pos_stereo_y, pos_stereo_x);

			histos->h_beamProfile_uncorr->Fill(pos_stereo_y_uncorr, pos_stereo_x_uncorr); 
			
		}

		if(cl_lay[2]->getNClusters2()==1 && cl_lay[1]->getNClusters2()==1)
		{
			float pos_lay2 = cl_lay[2]->getCorrPosition(0);
			histos->h_diffpos_lay2->Fill(pos_lay2, pos_sby1-pos_lay2);

			histos->h_diffpos_lay2_uncorr->Fill(cl_lay[2]->getStrip_from_clpos(cl_lay[2]->getPosition(0)), cl_lay[1]->getPosition(0)-cl_lay[2]->getPosition(0));
		}

		if(cl_lay[3]->getNClusters2()==1 && cl_lay[1]->getNClusters2()==1)
		{
			float pos_lay3 = cl_lay[3]->getCorrPosition(0);
			histos->h_diffpos_lay3->Fill(pos_lay3, pos_sby1-pos_lay3);

			histos->h_diffpos_lay3_uncorr->Fill(cl_lay[3]->getStrip_from_clpos(cl_lay[3]->getPosition(0)), cl_lay[1]->getPosition(0)-cl_lay[3]->getPosition(0));
		}

		if(cl_lay_small[2]->getNClusters2()>0)
		{
			for(int icl=0; icl<cl_lay_small[2]->getNClusters2(); icl++)
			{
				float pos_sby2 = cl_lay_small[2]->getPosition(icl);
				float pos_sby2_corr = cl_lay_small[2]->getCorrPosition(icl);
				histos->h_sby1_minus_sby2_vs_sby2->Fill(pos_sby2, pos_sby1-pos_sby2);
				histos->h_sby1_minus_sby2_vs_sby2_corr->Fill(pos_sby2_corr, pos_sby1-pos_sby2_corr);
			}
		}
		if(cl_lay_small[3]->getNClusters2()>0)
		{
			for(int icl=0; icl<cl_lay_small[3]->getNClusters2(); icl++)
			{
				float pos_sby3 = cl_lay_small[3]->getPosition(icl);
				float pos_sby3_corr = cl_lay_small[3]->getCorrPosition(icl);
				histos->h_sby1_minus_sby3_vs_sby3->Fill(pos_sby3, pos_sby1-pos_sby3);
				histos->h_sby1_minus_sby3_vs_sby3_corr->Fill(pos_sby3_corr, pos_sby1-pos_sby3_corr);
			}
		}
		if(cl_lay[0]->getNClusters2()==1)
		{
			for(int icl=0; icl<cl_lay[0]->getNClusters2(); icl++)
			{
				float pos_eta_out = cl_lay[0]->getPosition(icl);
				float pos_eta_out_corr = cl_lay[0]->getCorrPosition(icl);
				histos->h_sby1_minus_eta_out_vs_pos_eta_out->Fill(pos_eta_out, pos_sby1-pos_eta_out);
				histos->h_sby1_minus_eta_out_vs_pos_eta_out_corr->Fill(pos_eta_out_corr, pos_sby1-pos_eta_out_corr);
			}
		}
		if(cl_lay[1]->getNClusters2()>0)
		{
			for(int icl=0; icl<cl_lay[1]->getNClusters2(); icl++)
			{
				float pos_eta_in = cl_lay[1]->getPosition(icl);
				float pos_eta_in_corr = cl_lay[1]->getCorrPosition(icl);
				histos->h_sby1_minus_pos_eta_in_vs_pos_eta_in->Fill(pos_eta_in, pos_sby1-pos_eta_in);
				histos->h_sby1_minus_pos_eta_in_vs_pos_eta_in_corr->Fill(pos_eta_in_corr, pos_sby1-pos_eta_in_corr);
			}
		}
		if(cl_lay[2]->getNClusters2()==1 && cl_lay[3]->getNClusters2()==1)
		{
			for(int icl2=0; icl2<cl_lay[2]->getNClusters2(); icl2++)
			{
				float pos_stereo_in = cl_lay[2]->getPosition(icl2);
				float pos_stereo_in_corr = cl_lay[2]->getCorrPosition(icl2);
				for(int icl3=0; icl3<cl_lay[3]->getNClusters2(); icl3++)
				{
					float pos_stereo_out = cl_lay[3]->getCorrPosition(icl3);
					float pos_stereo_out_corr = cl_lay[3]->getCorrPosition(icl3);
					float pos_stereo_corr = ((pos_stereo_in_corr+pos_stereo_out_corr) / 2*TMath::Cos(1.5*TMath::Pi()/180.));
					float pos_stereo = (pos_stereo_in+pos_stereo_out) / 2*TMath::Cos(1.5*TMath::Pi()/180.);
					histos->h_sby1_minus_stereo_vs_stereo->Fill(pos_stereo, pos_sby1-pos_stereo);
					histos->h_sby1_minus_stereo_vs_stereo_corr->Fill(pos_stereo_corr, pos_sby1-pos_stereo_corr);
				}
			}
		}

// end checks

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
	out_file->mkdir("correlations");
	out_file->mkdir("superclusters");
	out_file->mkdir("corrected_positions");
	out_file->mkdir("tracking/SM1");
	out_file->mkdir("mmfe8_misalign");


	out_file->cd();

	histos->h_cutflow->SetBinContent(1, counter0);
	histos->h_cutflow->SetBinContent(2, counter1);

	histos->h_cutflow->GetXaxis()->SetBinLabel(1, "total events");
	histos->h_cutflow->GetXaxis()->SetBinLabel(2, "1 clus BSY1");

	out_file->cd("mmfe8_misalign");
	histos->h_res_SBY2_SBY3_vs_SBY2->Write();
	histos->h_res_SBY2_SBY3_vs_SBY3->Write();
	histos->h_diffpos_lay0->Write();
	histos->h_diffpos_lay1->Write();
	histos->h_beamProfile->Write();
	histos->h_beamProfile_uncorr->Write();
	histos->h_stereo_y_corr_vs_uncorr->Write();
	histos->h_stereo_x_corr_vs_uncorr->Write();
	histos->h_diffpos_stereo_uncorr_lay1->Write();
	histos->h_diffpos_lay2_uncorr->Write();
	histos->h_diffpos_lay3_uncorr->Write();

	histos->h_diffpos_stereolay1->Write();
	histos->h_diffpos_strip_index_lay2index_lay1->Write();
	histos->h_diffpos_strip_index_lay3index_lay1->Write();

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
		histos->h_charge_on_track[iLayer]->Write();

	}

	out_file->cd("correlations");
	histos->h_sby1_minus_eta_out_vs_pos_eta_out->Write();
	histos->h_sby1_minus_pos_eta_in_vs_pos_eta_in->Write();
	histos->h_sby1_minus_stereo_vs_stereo->Write();
	histos->h_sby1_minus_sby2_vs_sby2->Write();
	histos->h_sby1_minus_sby3_vs_sby3->Write();

	histos->h_sby1_minus_eta_out_vs_pos_eta_out_corr->Write();
	histos->h_sby1_minus_pos_eta_in_vs_pos_eta_in_corr->Write();
	histos->h_sby1_minus_stereo_vs_stereo_corr->Write();
	histos->h_sby1_minus_sby2_vs_sby2_corr->Write();
	histos->h_sby1_minus_sby3_vs_sby3_corr->Write();

	histos->h_bsy1_minus_etaout_vs_stereo_x->Write();
	histos->h_bsy1_minus_etain_vs_stereo_x->Write();
	histos->h_bsy1_minus_stereoin_vs_stereo_x->Write();
	histos->h_bsy1_minus_stereoout_vs_stereo_x->Write();

	histos->h_diffpos_lay2->Write();
	histos->h_diffpos_lay3->Write();

	histos->h_sclusters_lay0_lay1_1mm->Write();
	histos->h_sclusters_lay0_lay1_2mm->Write();
	histos->h_sclusters_lay0_lay1_5mm->Write();
	histos->h_sclusters_lay0_lay1_40mm->Write();

	out_file->cd("corrected_positions");
	histos->h_clus_positions_small[1]->Write();
	histos->h_cluster_position_corr_IP1->Write();
	histos->h_cluster_position_corr_IP2->Write();
	histos->h_cluster_position_corr_stereo->Write();
	histos->h_cluster_position_corr_SBY2->Write();
	histos->h_cluster_position_corr_SBY3->Write();

	out_file->cd("superclusters");
	histos->h_Nscluster_in_1deg_sby2_sby3->Write();
	histos->h_Nscluster_theonly_out_1deg_sby2_sby3->Write();
	histos->h_Nscluster_morethan1_in_1deg_sby2_sby3->Write();
	histos->h_Nscluster_morethan1_out_1deg_sby2_sby3->Write();
	histos->h_residual_sby2_sby3->Write();
	histos->h_residual_sby2_sby3_outside1deg_theonly_sclus->Write();
	histos->h_residual_sby2_sby3_out_1deg_many_sclus->Write();
	histos->h_residual_sby2_sby3_in_1deg_many_sclus->Write();
	histos->h_residual_sby2_sby3_rejected->Write();
	histos->h_residual_sby2_sby3_accepted->Write();
	histos->h_residual_sby2_sby3_zero_prob->Write();
	histos->h_nmanyrack_events->Write();
	histos->h_singletrack_events->Write();
	
	out_file->cd("track_candidates/residuals");
	histos->h_residuals_single_tracks->Write();
	histos->h_residuals_mult_tracks->Write();
	histos->h_residuals_accepted_tracks->Write();
	histos->h_residuals_single_tracks_anglecut->Write();
	histos->h_residuals_mult_tracks_anglecut->Write();
	histos->h_residuals_accepted_tracks_anglecut->Write();
	
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
	histos->h_d_track_etain->Write();
	histos->h_d_track_stereoin->Write();
	histos->h_d_track_stereoout->Write();
	histos->h_d_track_ystereo->Write();
	histos->h_mult_track_angle->Write();
	histos->h_nlayers_mult_track_events->Write();
	histos->h_prob_mult_tracks->Write();

	histos->h_d_track_lay2->Write();
	histos->h_d_track_lay3->Write();


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

	out_file->cd("SB/clusters");
	histos->h_cl_pos_sby2_sby3->Write();
	histos->h_cl_pos_sby1_sby2->Write();
	histos->h_cl_pos_sby1_sby3->Write();
	histos->h_nclus_in_5mm_window->Write();
	histos->h_nclus_in_1mm_window->Write();


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

	out_file->Close();

	delete histos;
	delete treeReader;

	
	return 0;
}


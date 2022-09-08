#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

class Histograms{
public:
	TH1F* h_raw_hits[4];
	TH1F* h_raw_hits_small[4];
	TH1F* h_raw_hit_trigger;

	TH1F* h_strip_index_SM1_cluster[4];
	TH1F* h_strip_index_SB_cluster[4];
	TH1F* h_strip_index_SM1[4];
	TH1F* h_strip_index_SB[4];

	TH1F* h_strip_2044_pdo_SB[4];
	TH1F* h_strip_2044_pdo_SM1[4];
	
	TH1F* h_nclusters[4];
	TH1F* h_nclusters_small[4];
	
	TH1F* h_nstrips[4];
	TH1F* h_nstrips_small[4];

	TH1F* h_cl_charge[4];
	TH1F* h_cl_charge_small[4];

	TH1F* h_cl_charge_leadCluster[4];
	TH1F* h_cl_charge_leadCluster_maxstrip[4];
	TH1F* h_maxstrip_pdo_allClus[4];
	TH1F* h_strip_pdo_in_leadCluster[4];

	TH1F* h_strip_pdo_in_clus[4];
	TH1F* h_strip_pdo_in_clus_small[4];

	TH1F* h_strip_tdo_in_clus[4];
	TH1F* h_strip_tdo_in_clus_small[4];

	TH1F* h_strip_relbcid_in_clus[4];
	TH1F* h_strip_relbcid_in_clus_small[4];

	TH1F* h_clus_positions[4];
	TH1F* h_clus_positions_small[4];

	TH1F* h_strip_time_in_clus[4];
	TH1F* h_strip_time_in_clus_small[4];

	TH1F* h_strip_time_in_clus_earliest[4];
	TH1F* h_strip_time_in_clus_earliest_small[4];

	TH1F* h_nlayers_with_cluster;
	TH1F* h_nlayers_mult_track_events;
	TH1F* h_mult_track_angle;
	TH1F* h_nclusters_per_layer_event[4];
	TH1F* h_cl_charge_on_track[4];
	TH1F* h_nstrips_on_track[4];

	TH1F* h_time_res_SB01;
	TH1F* h_time_res_SB1SM1lay0;
	TH1F* h_time_res_SB1SM1lay1;
	TH1F* h_time_res_SB1SM1lay2;
	TH1F* h_time_res_SB1SM1lay3;
	TH1F* h_time_res_SB12;
	TH1F* h_time_res_SB13;
	
	TH1F* h_strip_tdo[4][8192];
	TH1F* h_strip_tdo_small[4][8192];
	
	TH2F* h_raw_hits_vs_tot_strips[4];
	TH2F* h_raw_hits_vs_tot_strips_small[4];

	TH2F* h_cl_size_vs_cl_charge[4];
	TH2F* h_cl_size_vs_cl_charge_small[4];
	
	TH2F* h_strip_index_vs_pdo[4];
	TH2F* h_strip_index_vs_pdo_small[4];

	TH2F* h_strip_index_vs_pdo0[4];
	TH2F* h_strip_index_vs_pdo0_small[4];
	TH2F* h_strip_index_vs_pdo0_trigger;

	TH2F* h_tdo_vs_pdo[4];
	TH2F* h_tdo_vs_pdo_small[4];

	TH2F* h_strip_index_vs_tdo[4];
	TH2F* h_strip_index_vs_tdo_small[4];
	TH2F* h_strip_index_vs_tdo_trigger;

	TH2F* h_strip_index_vs_relbcid[4];
	TH2F* h_strip_index_vs_relbcid_small[4];
	TH2F* h_strip_index_vs_relbcid_trigger;

	TH2F* h_strip_index_vs_bcid[4];
	TH2F* h_strip_index_vs_bcid_small[4];
	TH2F* h_strip_index_vs_bcid_trigger;

	TH2F* h_diffpos_lay2;
	TH2F* h_diffpos_lay3;

	TH2F* h_beamProfile;
	TH2F* h_beamProfile_ontrack;

	TH2F* h_res_SBY2_SBY3_vs_SBY2;
	TH2F* h_res_SBY2_SBY3_vs_SBY3;

	TH2F* h_diffpos_lay0;
	TH2F* h_diffpos_lay1;
	TH2F* h_diffpos_stereolay0;
	TH2F* h_diffpos_stereolay1;

	TH1F* h_angle;
	TH1F* h_chi2ndf;
	TH1F* h_chi2;
	TH1F* h_prob;
	TH1F* h_prob_mult_tracks;
	TH1F* h_angle_only_scluster_sby2;
	TH1F* h_angle_only_scluster_sby3;
	TH1F* h_dangle_only_scluster_sby2sby3;

	TH2F* h_sby1_minus_eta_out_vs_pos_eta_out;
	TH2F* h_sby1_minus_pos_eta_in_vs_pos_eta_in;
	TH2F* h_sby1_minus_stereo_vs_stereo;
	TH2F* h_sby1_minus_stereo_in_vs_stereo_in;
	TH2F* h_sby1_minus_stereo_out_vs_stereo_out;
	TH2F* h_sby1_minus_sby2_vs_sby2;
	TH2F* h_sby1_minus_sby3_vs_sby3;

	TH1F* h_residual_sby2_sby3;
	TH1F* h_residual_sby2_sby3_outside1deg_theonly_sclus;
	TH1F* h_residual_sby2_sby3_out_1deg_many_sclus;
	TH1F* h_residual_sby2_sby3_in_1deg_many_sclus;
	TH1F* h_residual_sby2_sby3_rejected;
	TH1F* h_residual_sby2_sby3_accepted;
	TH1F* h_residual_sby2_sby3_zero_prob;
	TH1F* h_residuals_single_tracks;
	TH1F* h_residuals_mult_tracks;
	TH1F* h_residuals_accepted_tracks;

	TH1F* h_cutflow;

	TH1F* h_Nscluster_in_1deg_sby2_sby3;
	TH1F* h_Nscluster_theonly_out_1deg_sby2_sby3;
	TH1F* h_Nscluster_morethan1_in_and_out_sby2_sby3;
	TH1F* h_Nscluster_morethan1_in_1deg_sby2_sby3;
	TH1F* h_Nscluster_morethan1_out_1deg_sby2_sby3;

	TH1F* h_nmanyrack_events;
	TH1F* h_singletrack_events;

	TH1F* h_d_track_etaout;
	TH1F* h_d_track_etaout_cut;
	TH1F* h_d_track_etain;
	TH1F* h_d_track_etain_cut;
	TH1F* h_d_track_stereo;
	TH1F* h_d_track_stereo_cut;
	TH1F* h_d_track_ystereo;
	TH1F* h_d_track_lay2;
	TH1F* h_d_track_lay2_cut;
	TH1F* h_d_track_lay3;
	TH1F* h_d_track_lay3_cut;

	TH2F* h_res_SBY2_SBY1_vs_SBY1;
	TH2F* h_res_SBY3_SBY1_vs_SBY1;
	TH2F* h_pos_etain_vs_pos_etaout;
	TH2F* h_real_vs_exp_position_etaout;

	Histograms();
	inline void init();

};

inline Histograms::Histograms()
{
	init();
}

inline void Histograms::init()
{
	h_cutflow = new TH1F("h_cutflow", "Cutflow" ,10, 0, 10);

	h_real_vs_exp_position_etaout = new TH2F("h_real_vs_exp_position_etaout", "real vs expected position" ,600, 2000, 5000, 600, 1300, 1900);
	h_real_vs_exp_position_etaout->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_real_vs_exp_position_etaout->GetYaxis()->SetTitle("expected position [mm]");

	h_d_track_etaout = new TH1F("h_d_track_etaout", "Distance from track - eta_out" ,1000, -10, 10);
	h_d_track_etaout->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etaout_cut = new TH1F("h_d_track_etaout_cut", "Distance from track - eta_out" ,1000, -10, 10);
	h_d_track_etaout_cut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etain = new TH1F("h_d_track_etain", "Distance from track - eta_in" ,1000, -10, 10);
	h_d_track_etain->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etain_cut = new TH1F("h_d_track_etain_cut", "Distance from track - eta_in" ,1000, -10, 10);
	h_d_track_etain_cut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_stereo = new TH1F("h_d_track_stereo", "Distance from track - stereo" ,1000, -10, 10);
	h_d_track_stereo->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_stereo_cut = new TH1F("h_d_track_stereo_cut", "Distance from track - stereo" ,1000, -10, 10);
	h_d_track_stereo_cut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay2 = new TH1F("h_d_track_lay2", "Distance from track - stereo_in" ,1000, -10, 10);
	h_d_track_lay2->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay2_cut = new TH1F("h_d_track_lay2_cut", "Distance from track - stereo_in" ,1000, -10, 10);
	h_d_track_lay2_cut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay3 = new TH1F("h_d_track_lay3", "Distance from track - stereo_out" ,1000, -10, 10);
	h_d_track_lay3->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay3_cut = new TH1F("h_d_track_lay3_cut", "Distance from track - stereo_out" ,1000, -10, 10);
	h_d_track_lay3_cut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_ystereo = new TH1F("h_d_track_ystereo", "Distance from track - stereo_y" ,1000, -10, 10);
	h_d_track_ystereo->GetXaxis()->SetTitle("distance [mm]");

	h_nlayers_with_cluster = new TH1F("h_nlayers_with_cluster", "Number of layers (SM1) with reconstructed clusters" ,5, 0, 5);
	h_nlayers_with_cluster->GetXaxis()->SetTitle("Nlayers");

	h_nlayers_mult_track_events = new TH1F("h_nlayers_mult_track_events", "Number of layers (SM1) with reconstructed clusters (>1 tracks)" ,5, 0, 5);
	h_nlayers_mult_track_events->GetXaxis()->SetTitle("Nlayers");

	h_nmanyrack_events = new TH1F("h_nmanyrack_events", "Number of multiple tracks found" ,20, 0, 20);
	h_nmanyrack_events->GetXaxis()->SetTitle("Ntracks");

	h_singletrack_events = new TH1F("h_singletrack_events", "Number of single tracks found" ,3, 0, 3);
	h_singletrack_events->GetXaxis()->SetTitle("Ntracks");
	
	h_Nscluster_in_1deg_sby2_sby3 = new TH1F("h_Nscluster_in_1deg_sby2_sby3", "superclusters in 1deg" ,20, 0, 20);
	h_Nscluster_in_1deg_sby2_sby3->GetXaxis()->SetTitle("Nsclusters");
	h_Nscluster_theonly_out_1deg_sby2_sby3 = new TH1F("h_Nscluster_theonly_out_1deg_sby2_sby3", "the only supercluster out 1deg" ,20, 0, 20);
	h_Nscluster_theonly_out_1deg_sby2_sby3->GetXaxis()->SetTitle("Nsclusters");
	h_Nscluster_morethan1_in_and_out_sby2_sby3 = new TH1F("h_Nscluster_morethan1_in_and_out_sby2_sby3", "events with many superclusters" ,20, 0, 20);
	h_Nscluster_morethan1_in_and_out_sby2_sby3->GetXaxis()->SetTitle("Nsclusters");
	h_Nscluster_morethan1_in_1deg_sby2_sby3 = new TH1F("h_Nscluster_morethan1_in_1deg_sby2_sby3", "events with many superclusters in 1deg" ,20, 0, 20);
	h_Nscluster_morethan1_in_1deg_sby2_sby3->GetXaxis()->SetTitle("Nsclusters");
	h_Nscluster_morethan1_out_1deg_sby2_sby3 = new TH1F("h_Nscluster_morethan1_out_1deg_sby2_sby3", "events with many superclusters out 1deg" ,20, 0, 20);
	h_Nscluster_morethan1_out_1deg_sby2_sby3->GetXaxis()->SetTitle("Nsclusters");

	h_residual_sby2_sby3 = new TH1F("h_residual_sby2_sby3", "residuals SBY2-SBY3 in 1deg (exactly 1 supercluster)" ,1000, -10, 10);
	h_residual_sby2_sby3_outside1deg_theonly_sclus = new TH1F("h_residual_sby2_sby3_outside1deg_theonly_sclus", "residuals SBY2-SBY3 the only supercluster but outside 1deg" ,1000, -10, 10);
	h_residual_sby2_sby3_out_1deg_many_sclus = new TH1F("h_residual_sby2_sby3_out_1deg_many_sclus", "residuals SBY2-SBY3 (events with many superclusters none of them in 1deg)" ,1000, -10, 10); // these tracks are not studied at all
	h_residual_sby2_sby3_in_1deg_many_sclus = new TH1F("h_residual_sby2_sby3_in_1deg_many_sclus", "residuals SBY2-SBY3 (events with many superclusters in 1deg - all possible tracks)" ,1000, -10, 10);
	h_residual_sby2_sby3_rejected = new TH1F("h_residual_sby2_sby3_in_1deg_many_sclus_but_reject", "residuals SBY2-SBY3 (events with many superclusters in 1deg - rejections)" ,1000, -10, 10);
	h_residual_sby2_sby3_accepted = new TH1F("h_residual_sby2_sby3_accepted", "residuals SBY2-SBY3 (accepted events)" ,1000, -10, 10);
	h_residual_sby2_sby3_zero_prob = new TH1F("h_residual_sby2_sby3_zero_prob", "residuals SBY2-SBY3 (tracks Prob=0)" ,1000, -10, 10);
	h_residuals_single_tracks = new TH1F("h_residuals_single_tracks", "residuals SBY2-SBY3 (accepted single tracks)" ,1000, -10, 10);
	h_residuals_mult_tracks = new TH1F("h_residuals_mult_tracks", "residuals SBY2-SBY3 (accepted multiple tracks)" ,1000, -10, 10);
	h_residuals_accepted_tracks = new TH1F("h_residuals_accepted_tracks", "residuals SBY2-SBY3 (all accepted tracks)" ,1000, -10, 10);

	h_raw_hit_trigger = new TH1F("h_nhits_trigger", "trigger", 10, 0, 10);

	h_angle = new TH1F("h_angle", "track angle", 5000, -50, 50);
	h_angle->GetXaxis()->SetTitle("#theta [deg]");

	h_mult_track_angle = new TH1F("h_mult_track_angle", "track angle (>1 tracks)", 5000, -50, 50);
	h_mult_track_angle->GetXaxis()->SetTitle("#theta [deg]");

	h_angle_only_scluster_sby2 = new TH1F("h_angle_only_scluster_sby2", "track angle (the only supercluster outside 1deg window)", 5000, -50, 50);
	h_angle_only_scluster_sby2->GetXaxis()->SetTitle("#theta [deg]");
	h_angle_only_scluster_sby3 = new TH1F("h_angle_only_scluster_sby3", "track angle (the only supercluster outside 1deg window)", 5000, -50, 50);
	h_angle_only_scluster_sby3->GetXaxis()->SetTitle("#theta [deg]");

	h_dangle_only_scluster_sby2sby3 = new TH1F("h_dangle_only_scluster_sby2sby3", "#Delta#theta (the only supercluster outside 1deg window)", 5000, -50, 50);
	h_dangle_only_scluster_sby2sby3->GetXaxis()->SetTitle("#theta^{SBY2} - #theta^{SBY3} [deg]");

	h_chi2ndf = new TH1F("h_chi2ndf", "chi2 / ndf", 150, 0, 15);
	h_chi2 = new TH1F("h_chi2", "chi2", 100, 0, 100);
	h_prob = new TH1F("h_prob", "h_prob", 100, 0, 1);
	h_prob_mult_tracks = new TH1F("h_prob_mult_tracks", "h_prob_mult_tracks", 100, 0, 1);

	h_strip_index_vs_pdo0_trigger = new TH2F("h_strip_index_vs_pdo_trigger", "trigger", 8192, 0, 8192, 1024, 0, 1024);
	h_strip_index_vs_pdo0_trigger->GetXaxis()->SetTitle("strip index");
	h_strip_index_vs_pdo0_trigger->GetYaxis()->SetTitle("strip pdo [adc counts]");

	h_strip_index_vs_tdo_trigger = new TH2F("h_strip_index_vs_tdo_trigger", "trigger", 8192, 0, 8192, 256, 0, 256);
	h_strip_index_vs_relbcid_trigger = new TH2F("h_strip_index_vs_relbcid_trigger", "trigger", 8192, 0, 8192, 8, -0.5, 7.5);
	h_strip_index_vs_bcid_trigger = new TH2F("h_strip_index_vs_bcid_trigger", "trigger", 8192, 0, 8192, 3564, 0, 3564);

	h_time_res_SB01 = new TH1F("h_time_res_SB01","Time resolution SBY1 - SBX0",100, -250, 250);
	h_time_res_SB1SM1lay0 = new TH1F("h_time_res_SB1SM10","Time resolution SBY1 - SM1(0)",100, -250, 250);
	h_time_res_SB1SM1lay1 = new TH1F("h_time_res_SB1SM11","Time resolution SBY1 - SM1(1)",100, -250, 250);
	h_time_res_SB1SM1lay2 = new TH1F("h_time_res_SB1SM12","Time resolution SBY1 - SM1(2)",100, -250, 250);
	h_time_res_SB1SM1lay3 = new TH1F("h_time_res_SB1SM13","Time resolution SBY1 - SM1(3)",100, -250, 250);
	h_time_res_SB12 = new TH1F("h_time_res_SB12","Time resolution SBY1 - SBY2",100, -250, 250);
	h_time_res_SB13 = new TH1F("h_time_res_SB13","Time resolution SBY1 - SBY3",100, -250, 250);

//end here


// histograms used for the alignment

	h_diffpos_lay0 = new TH2F("h_diffpos_lay0", "", 100, -2, 2,100,1800,2400); 
	h_diffpos_lay0->GetXaxis()->SetTitle("#Delta y(IP2-IP1) [mm]");
	h_diffpos_lay0->GetYaxis()->SetTitle("y_{IP1} [mm]");

	h_diffpos_lay1 = new TH2F("h_diffpos_lay1", "", 100, -2, 2,100,1800,2400);
	h_diffpos_lay1->GetXaxis()->SetTitle("#Delta y(IP2-IP1) [mm]");
	h_diffpos_lay1->GetYaxis()->SetTitle("y_{IP2} [mm]");

	h_diffpos_lay2 = new TH2F("h_diffpos_lay2", "IP1-IP3 vs IP3 position ", 600, 1800, 2400, 40000, -20, 20);
	h_diffpos_lay2->GetXaxis()->SetTitle("IP3 cluster position [mm]");
	h_diffpos_lay2->GetYaxis()->SetTitle("IP1-IP3 cluster position [mm]");

	h_diffpos_lay3 = new TH2F("h_diffpos_lay3", "IP1-IP4 vs IP4 position ", 600, 1800, 2400, 40000, -20, 20);
	h_diffpos_lay3->GetXaxis()->SetTitle("IP4 cluster position [mm]");
	h_diffpos_lay3->GetYaxis()->SetTitle("IP1-IP4 cluster position [mm]");

	h_diffpos_stereolay0 = new TH2F("h_diffpos_stereolay0", "", 20000, -20, 20, 600, 1300,1900);
	h_diffpos_stereolay0->GetXaxis()->SetTitle("#Delta y(SBY1 - IP3) [mm]");
	h_diffpos_stereolay0->GetYaxis()->SetTitle("y_{stereo} [mm]");

	h_diffpos_stereolay1 = new TH2F("h_diffpos_stereolay1", "", 20000, -20., 20.,600,1300,1900);
	h_diffpos_stereolay1->GetXaxis()->SetTitle("#Delta y(SBY1 - IP4) [mm]");
	h_diffpos_stereolay1->GetYaxis()->SetTitle("y_{SBY1} [mm]");

	h_res_SBY2_SBY3_vs_SBY2 = new TH2F("h_res_SBY2_SBY3_vs_SBY2", "", 20000, -2, 2,600,1350,1950);
	h_res_SBY2_SBY3_vs_SBY2->GetXaxis()->SetTitle("#Delta y(SBY2-SBY3) [mm]");
	h_res_SBY2_SBY3_vs_SBY2->GetYaxis()->SetTitle("y_{SBY2} [mm]");

	h_res_SBY2_SBY3_vs_SBY3 = new TH2F("h_res_SBY2_SBY3_vs_SBY3", "", 20000, -2, 2,600,1350,1950);
	h_res_SBY2_SBY3_vs_SBY3->GetXaxis()->SetTitle("#Delta y(SBY2-SBY3) [mm]");
	h_res_SBY2_SBY3_vs_SBY3->GetYaxis()->SetTitle("y_{SBY3} [mm]");

	h_res_SBY2_SBY1_vs_SBY1 = new TH2F("h_res_SBY2_SBY1_vs_SBY1", "", 10000, -100, 100,600,1350,1950);
	h_res_SBY2_SBY1_vs_SBY1->GetXaxis()->SetTitle("#Delta y(SBY2-SBY1) [mm]");
	h_res_SBY2_SBY1_vs_SBY1->GetYaxis()->SetTitle("y_{SBY1} [mm]");

	h_res_SBY3_SBY1_vs_SBY1 = new TH2F("h_res_SBY3_SBY1_vs_SBY1", "", 10000, -100, 100,600,1350,1950);
	h_res_SBY3_SBY1_vs_SBY1->GetXaxis()->SetTitle("#Delta y(SBY3-SBY1) [mm]");
	h_res_SBY3_SBY1_vs_SBY1->GetYaxis()->SetTitle("y_{SBY1} [mm]");

	h_beamProfile = new TH2F("h_beamProfile", "", 600, 1300, 2000, 100, -.4, .4);
	h_beamProfile->GetYaxis()->SetTitle("#phi [mm]");
	h_beamProfile->GetXaxis()->SetTitle("#eta [mm]");

	h_beamProfile_ontrack = new TH2F("h_beamProfile_ontrack", "", 600, 1300, 2000, 100, -.4, .4);
	h_beamProfile_ontrack->GetYaxis()->SetTitle("#phi [mm]");
	h_beamProfile_ontrack->GetXaxis()->SetTitle("#eta [mm]");

	h_sby1_minus_eta_out_vs_pos_eta_out = new TH2F("h_sby1_minus_eta_out_vs_pos_eta_out", "SB1Y-IP1 vs IP1 position", 700, 1800, 2500, 700, -700, 0);
	h_sby1_minus_eta_out_vs_pos_eta_out->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_sby1_minus_eta_out_vs_pos_eta_out->GetYaxis()->SetTitle("SB1Y-IP1 cluster position [mm]");

	h_sby1_minus_pos_eta_in_vs_pos_eta_in = new TH2F("h_sby1_minus_pos_eta_in_vs_pos_eta_in", "SB1Y-IP2 vs IP2 position", 700, 1800, 2500, 700, -700, 0);
	h_sby1_minus_pos_eta_in_vs_pos_eta_in->GetXaxis()->SetTitle("IP2 cluster position [mm]");
	h_sby1_minus_pos_eta_in_vs_pos_eta_in->GetYaxis()->SetTitle("SB1Y-IP2 cluster position [mm]");

	h_sby1_minus_stereo_vs_stereo = new TH2F("h_sby1_minus_stereo_vs_stereo", "BS1Y-stereo vs SBY1 position", 700, 1300, 2000, 700, -700, 0);
	h_sby1_minus_stereo_vs_stereo->GetXaxis()->SetTitle("SBY1 cluster position [mm]");
	h_sby1_minus_stereo_vs_stereo->GetYaxis()->SetTitle("SB1Y-stereo cluster position [mm]");

	h_sby1_minus_stereo_in_vs_stereo_in = new TH2F("h_sby1_minus_stereo_in_vs_stereo_in", "SB1Y-stereo_in vs SBY1 position", 700, 1300, 2000, 700, -700, 0);
	h_sby1_minus_stereo_in_vs_stereo_in->GetXaxis()->SetTitle("SBY1 cluster position [mm]");
	h_sby1_minus_stereo_in_vs_stereo_in->GetYaxis()->SetTitle("SBY1-stereo_in cluster position [mm]");

	h_sby1_minus_stereo_out_vs_stereo_out = new TH2F("h_sby1_minus_stereo_out_vs_stereo_out", "SBY1-stereo_out vs SBY1 position", 700, 1300, 2000, 700, -700, 0);
	h_sby1_minus_stereo_out_vs_stereo_out->GetXaxis()->SetTitle("SBY1 cluster position [mm]");
	h_sby1_minus_stereo_out_vs_stereo_out->GetYaxis()->SetTitle("SBY1-stereo_out cluster position [mm]");

	h_sby1_minus_sby2_vs_sby2 = new TH2F("h_sby1_minus_sby2_vs_sby2", "BSY1-BSY2 vs BSY2 position", 700, 1300, 2000, 700, 350, -350);
	h_sby1_minus_sby2_vs_sby2->GetXaxis()->SetTitle("BSY2 cluster position [mm]");
	h_sby1_minus_sby2_vs_sby2->GetYaxis()->SetTitle("BSY1-BSY2 cluster position [mm]");

	h_sby1_minus_sby3_vs_sby3 = new TH2F("h_sby1_minus_sby3_vs_sby3", "BS1Y-BSY3 vs BSY3 position", 700, 1300, 2000, 700, 350, -350);
	h_sby1_minus_sby3_vs_sby3->GetXaxis()->SetTitle("BSY3 cluster position [mm]");
	h_sby1_minus_sby3_vs_sby3->GetYaxis()->SetTitle("BSY1-BSY3 cluster position [mm]");

	h_pos_etain_vs_pos_etaout = new TH2F("h_pos_etain_vs_pos_etaout", "eta_in vs eta_out", 450, 1800, 2250, 450, 1800, 2250);
	h_pos_etain_vs_pos_etaout->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_pos_etain_vs_pos_etaout->GetYaxis()->SetTitle("IP2 cluster position [mm]");

// end alignment

	std::string type = "X";
	for(int ilayer = 0; ilayer<4; ilayer++) {
		if(ilayer>0) type = "Y";

		h_strip_2044_pdo_SM1[ilayer] = new TH1F(Form("h_strip_2044_pdo_SM1_%i",ilayer), Form("SM1 Layer - %i (tot clus PDO=2044)",ilayer), 5120, 0, 5119);
		h_strip_2044_pdo_SM1[ilayer]->GetXaxis()->SetTitle("strip index");

		h_strip_2044_pdo_SB[ilayer] = new TH1F(Form("h_strip_2044_pdo_SB_%i",ilayer), Form("SB Layer - %i (tot clus PDO=2044)",ilayer), 5120, 0, 5119);
		h_strip_2044_pdo_SB[ilayer]->GetXaxis()->SetTitle("strip index");

		h_strip_index_SM1[ilayer] = new TH1F(Form("h_strip_index_SM1_%i",ilayer), Form("SM1 Layer - %i",ilayer), 5120, 0, 5119);
		h_strip_index_SM1[ilayer]->GetXaxis()->SetTitle("strip index");

		h_strip_index_SB[ilayer] = new TH1F(Form("h_strip_index_SB_%i",ilayer), Form("SB Layer - %i",ilayer), 5120, 0, 5119);
		h_strip_index_SB[ilayer]->GetXaxis()->SetTitle("strip index");

		h_strip_index_SM1_cluster[ilayer] = new TH1F(Form("h_strip_index_SM1_cluster_%i",ilayer), Form("SM1 Layer - %i: strips fired after clustering",ilayer), 5120, 0, 5119);
		h_strip_index_SM1_cluster[ilayer]->GetXaxis()->SetTitle("strip index");

		h_strip_index_SB_cluster[ilayer] = new TH1F(Form("h_strip_index_SB_cluster_%i",ilayer), Form("SB Layer - %i: strips fired after clustering",ilayer), 5120, 0, 5119);
		h_strip_index_SB_cluster[ilayer]->GetXaxis()->SetTitle("strip index");

		h_raw_hits[ilayer] = new TH1F(Form("h_nhits_lay_SM1_%i",ilayer), Form("SM1 Layer - %i",ilayer), 100, 0, 300);
		h_raw_hits[ilayer]->GetXaxis()->SetTitle("nhits");

		h_raw_hits_small[ilayer] = new TH1F(Form("h_nhits_small_lay_%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 100, 0, 300);
		h_raw_hits_small[ilayer]->GetXaxis()->SetTitle("nhits");

		h_nclusters[ilayer] = new TH1F(Form("h_nclusters_SM1_lay_%i",ilayer), Form("SM1 Layer - %i",ilayer), 40, 0, 40);
		h_nclusters[ilayer]->GetXaxis()->SetTitle("Nclusters");

		h_nclusters_small[ilayer] = new TH1F(Form("h_nclusters_small_lay_%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 40, 0, 40);
		h_nclusters_small[ilayer]->GetXaxis()->SetTitle("Nclusters");

		h_nstrips[ilayer] = new TH1F(Form("h_nstrips_SM1_lay_%i",ilayer), Form("SM1 Layer - %i",ilayer), 40, 0, 40);
		h_nstrips[ilayer]->GetXaxis()->SetTitle("Nstrips");

		h_nstrips_on_track[ilayer] = new TH1F(Form("h_nstrips_SM1_ontrack_lay_%i",ilayer), Form("SM1 Layer - %i",ilayer), 40, 0, 40);
		h_nstrips_on_track[ilayer]->GetXaxis()->SetTitle("Nstrips");

		h_nstrips_small[ilayer] = new TH1F(Form("h_nstrips_small_lay_%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 40, 0, 40);
		h_nstrips_small[ilayer]->GetXaxis()->SetTitle("Nstrips");

		h_strip_time_in_clus[ilayer] = new TH1F(Form("h_strip_time_in_clus_lay_%i",ilayer), Form("SM1 Layer - %i",ilayer), 100, 0, 200);
		h_strip_time_in_clus[ilayer]->GetXaxis()->SetTitle("time [ns]");

		h_strip_time_in_clus_earliest[ilayer] = new TH1F(Form("h_strip_earliest_time_in_clus_lay_%i",ilayer), Form("SM1 Layer - %i",ilayer), 100, 0, 200);
		h_strip_time_in_clus_earliest[ilayer]->GetXaxis()->SetTitle("time [ns]");

		h_strip_time_in_clus_small[ilayer] = new TH1F(Form("h_strip_time_in_clus_small_lay_%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 100, 0, 200);
		h_strip_time_in_clus_small[ilayer]->GetXaxis()->SetTitle("time [ns]");

		h_strip_time_in_clus_earliest_small[ilayer] = new TH1F(Form("h_strip_earliest_time_in_clus_small_lay_%i",ilayer),  Form("SB%s  - %i",type.c_str(),ilayer), 100, 0, 200);
		h_strip_time_in_clus_earliest_small[ilayer]->GetXaxis()->SetTitle("time [ns]");

		h_strip_pdo_in_clus[ilayer] = new TH1F(Form("h_strip_pdo_SM1_lay_%i",ilayer), Form("SM1 Layer - %i",ilayer), 1024, 0, 1024);
		h_strip_pdo_in_clus[ilayer]->GetXaxis()->SetTitle("strip pdo [adc counts]");

		h_strip_pdo_in_clus_small[ilayer] = new TH1F(Form("h_strip_pdo_small_lay_%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 1024, 0, 1024);
		h_strip_pdo_in_clus_small[ilayer]->GetXaxis()->SetTitle("strip pdo [adc counts]");

		h_strip_index_vs_pdo[ilayer] = new TH2F(Form("h_strip_index_vs_pdo_SM1_lay_%i",ilayer), Form("SM1 Layer - %i",ilayer), 8192, 0, 8192, 1024, 0, 1024);
		h_strip_index_vs_pdo[ilayer]->GetXaxis()->SetTitle("strip index");
		h_strip_index_vs_pdo[ilayer]->GetYaxis()->SetTitle("strip pdo [adc counts]");

		h_strip_index_vs_pdo_small[ilayer] = new TH2F(Form("h_strip_index_vs_pdo_small_lay_%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 8192, 0, 8192, 1024, 0, 1024);
		h_strip_index_vs_pdo_small[ilayer]->GetXaxis()->SetTitle("strip index");
		h_strip_index_vs_pdo_small[ilayer]->GetYaxis()->SetTitle("strip pdo [adc counts]");

		h_strip_index_vs_pdo0[ilayer] =  new TH2F(Form("h_strip_index_vs_pdo0_SM1_lay_%i",ilayer), Form("SM1 Layer - %i",ilayer), 8192, 0, 8192, 1024, 0, 1024);
		h_strip_index_vs_pdo0[ilayer]->GetXaxis()->SetTitle("strip index");
		h_strip_index_vs_pdo0[ilayer]->GetYaxis()->SetTitle("strip pdo [adc counts]");

		h_strip_index_vs_pdo0_small[ilayer] =  new TH2F(Form("h_strip_index_vs_pdo0_small_lay_%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 8192, 0, 8192, 1024, 0, 1024);
		h_strip_index_vs_pdo0_small[ilayer]->GetXaxis()->SetTitle("strip index");
		h_strip_index_vs_pdo0_small[ilayer]->GetYaxis()->SetTitle("strip pdo [adc counts]");

		h_tdo_vs_pdo[ilayer] = new TH2F(Form("h_tdo_vs_pdo_SM1_lay_%i",ilayer), Form("SM1 Layer - %i",ilayer), 256, 0, 256, 1024, 0, 1024);
		h_tdo_vs_pdo[ilayer]->GetXaxis()->SetTitle("strip tdo [adc counts]");
		h_tdo_vs_pdo[ilayer]->GetYaxis()->SetTitle("strip pdo [adc counts]");

		h_tdo_vs_pdo_small[ilayer] = new TH2F(Form("h_tdo_vs_pdo_small_lay_%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 256, 0, 256, 1024, 0, 1024);
		h_tdo_vs_pdo_small[ilayer]->GetXaxis()->SetTitle("strip tdo [adc counts]");
		h_tdo_vs_pdo_small[ilayer]->GetYaxis()->SetTitle("strip pdo [adc counts]");

		h_strip_index_vs_tdo[ilayer] = new TH2F(Form("h_strip_index_vs_tdo_SM1_lay%i",ilayer), Form("SM1 Layer - %i",ilayer), 8192, 0, 8192, 256, 0, 256);
		h_strip_index_vs_relbcid[ilayer] = new TH2F(Form("h_strip_index_vs_relbcid_SM1_lay%i",ilayer), Form("SM1 Layer - %i",ilayer), 8192, 0, 8192, 8, -0.5, 7.5);
		h_strip_index_vs_bcid[ilayer] = new TH2F(Form("h_strip_index_vs_bcid_SM1_lay%i",ilayer), Form("SM1 Layer - %i",ilayer), 8192, 0, 8192, 3564, 0, 3564);

		h_strip_index_vs_tdo_small[ilayer] = new TH2F(Form("h_strip_index_vs_tdo_small_lay%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 8192, 0, 8192, 256, 0, 256);
		h_strip_index_vs_relbcid_small[ilayer] = new TH2F(Form("h_strip_index_vs_relbcid_small_lay%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 8192, 0, 8192, 8, -0.5, 7.5);
		h_strip_index_vs_bcid_small[ilayer] = new TH2F(Form("h_strip_index_vs_bcid_small_lay%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 8192, 0, 8192, 3564, 0, 3564);

		h_strip_relbcid_in_clus[ilayer] = new TH1F(Form("h_strip_relbcid_lay_SM1_%i",ilayer), Form("SM1 Layer - %i",ilayer), 8, -0.5, 7.5);
		h_strip_relbcid_in_clus[ilayer]->GetXaxis()->SetTitle("relative BCID");
		h_strip_tdo_in_clus[ilayer] = new TH1F(Form("h_strip_tdo_lay_SM1_%i",ilayer), Form("SM1 Layer - %i",ilayer),85,0,255);
		h_strip_tdo_in_clus[ilayer]->GetXaxis()->SetTitle("strip tdo");

		h_strip_relbcid_in_clus_small[ilayer] = new TH1F(Form("h_strip_relbcid_lay_small_%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer), 8, -0.5, 7.5);
		h_strip_relbcid_in_clus_small[ilayer]->GetXaxis()->SetTitle("relative BCID");
		h_strip_tdo_in_clus_small[ilayer] = new TH1F(Form("h_strip_tdo_lay_small_%i",ilayer), Form("SB%s  - %i",type.c_str(),ilayer),85,0,255);
		h_strip_tdo_in_clus_small[ilayer]->GetXaxis()->SetTitle("strip tdo");

		h_raw_hits_vs_tot_strips[ilayer] = new TH2F(Form("h_raw_hits_vs_tot_strips_SM1_lay_%i",ilayer), Form("SM1 Layer - %i", ilayer), 40, 0, 40, 15, 0, 50);
		h_raw_hits_vs_tot_strips[ilayer]->GetYaxis()->SetTitle("Nhits");
		h_raw_hits_vs_tot_strips[ilayer]->GetXaxis()->SetTitle("Nstrips (all clusters)");

		h_raw_hits_vs_tot_strips_small[ilayer] = new TH2F(Form("h_raw_hits_vs_tot_strips_small_lay_%i",ilayer), Form("SB%s  - %i",type.c_str(), ilayer), 40, 0, 40, 15, 0, 50);
		h_raw_hits_vs_tot_strips_small[ilayer]->GetYaxis()->SetTitle("Nhits");
		h_raw_hits_vs_tot_strips_small[ilayer]->GetXaxis()->SetTitle("Nstrips (all clusters)");

		h_cl_size_vs_cl_charge[ilayer] = new TH2F(Form("h_cl_size_vs_cl_charge_SM1_lay_%i",ilayer), Form("SM1 Layer - %i", ilayer), 40, 0, 40, 100024, 0, 100024);
		h_cl_size_vs_cl_charge[ilayer]->GetYaxis()->SetTitle("cluster pdo [adc counts]");
		h_cl_size_vs_cl_charge[ilayer]->GetXaxis()->SetTitle("Nstrips (each cluster)");

		h_cl_size_vs_cl_charge_small[ilayer] = new TH2F(Form("h_cl_size_vs_cl_charge_lay_small_%i",ilayer), Form("SB%s  - %i",type.c_str(), ilayer), 40, 0, 40, 100024, 0, 100024);
		h_cl_size_vs_cl_charge_small[ilayer]->GetYaxis()->SetTitle("cluster pdo [adc counts]");
		h_cl_size_vs_cl_charge_small[ilayer]->GetXaxis()->SetTitle("Nstrips (each cluster)");

		h_cl_charge[ilayer] = new TH1F(Form("h_cl_charge_SM1_lay%i",ilayer), Form("SM1 Layer - %i", ilayer), 100024, 0, 100024);
		h_cl_charge[ilayer]->GetXaxis()->SetTitle("cluster pdo [adc counts]");

		h_cl_charge_on_track[ilayer] = new TH1F(Form("h_cl_charge_on_track_lay%i",ilayer), Form("SM1 Layer - %i", ilayer), 100024, 0, 100024);
		h_cl_charge_on_track[ilayer]->GetXaxis()->SetTitle("cluster pdo [adc counts]");

		h_cl_charge_leadCluster[ilayer] = new TH1F(Form("h_cl_charge_leadCluster_SM1_lay%i",ilayer), Form("SM1 Layer - %i", ilayer), 100024, 0, 100024);
		h_cl_charge_leadCluster[ilayer]->GetXaxis()->SetTitle("leading cluster pdo [adc counts]");

		h_cl_charge_leadCluster_maxstrip[ilayer] = new TH1F(Form("h_cl_charge_leadCluster_maxstrip_SM1_lay%i",ilayer), Form("SM1 Layer - %i", ilayer), 1024, 0, 1024);
		h_cl_charge_leadCluster_maxstrip[ilayer]->GetXaxis()->SetTitle("max strip pdo in leading cluster [adc counts]");

		h_maxstrip_pdo_allClus[ilayer] = new TH1F(Form("h_maxstrip_pdo_allClus_SM1_lay%i",ilayer), Form("SM1 Layer - %i", ilayer), 1024, 0, 1024);
		h_maxstrip_pdo_allClus[ilayer]->GetXaxis()->SetTitle("max strip pdo in all clusters [adc counts]");

		h_strip_pdo_in_leadCluster[ilayer] = new TH1F(Form("h_strip_pdo_in_leadCluster_SM1_lay%i",ilayer), Form("SM1 Layer - %i", ilayer), 1024, 0, 1024);
		h_strip_pdo_in_leadCluster[ilayer]->GetXaxis()->SetTitle("strip pdo in leading cluster [adc counts]");

		h_cl_charge_small[ilayer] = new TH1F(Form("h_cl_charge_small_lay%i",ilayer), Form("SB%s  - %i",type.c_str(), ilayer), 100024, 0, 100024);
		h_cl_charge_small[ilayer]->GetXaxis()->SetTitle("cluster pdo [adc counts]");

		h_clus_positions[ilayer] = new TH1F(Form("h_clus_position_SM1_lay%i",ilayer), Form("SM1 Layer - %i", ilayer) ,7000, 0, 7000);
		h_clus_positions_small[ilayer] = new TH1F(Form("h_clus_position_small_lay%i",ilayer), Form("SB%s  - %i",type.c_str(), ilayer) ,7000, 0, 7000);
		
		h_nclusters_per_layer_event[ilayer] = new TH1F(Form("h_nclusters_per_layer_event_%i",ilayer), Form("SM1 Layer - %i", ilayer) ,40, 0, 40);
		
		for(int istrip=0; istrip<8192; istrip++){
			h_strip_tdo[ilayer][istrip] = new TH1F(Form("h_strip_%i_tdo_SM1_lay%i",istrip,ilayer), Form("SM1 strip - %i",istrip),85,0,255);
			h_strip_tdo_small[ilayer][istrip] = new TH1F(Form("h_strip_%i_tdo_small_lay%i",istrip,ilayer), Form("SB strip - %i",istrip),85,0,255);
		}
	}

	
}

#endif
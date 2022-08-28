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
	TH1F* h_charge_on_track[4];

	TH1F* h_time_res_SB01;
	TH1F* h_time_res_SB1SM1lay0;
	TH1F* h_time_res_SB1SM1lay1;
	TH1F* h_time_res_SB1SM1lay2;
	TH1F* h_time_res_SB1SM1lay3;
	TH1F* h_time_res_SB12;
	TH1F* h_time_res_SB13;

	TH1F* h_s_diff_y;
	TH1F* h_s_diff_x;
	TH1F* h_nlay_rec_clus;

	TH2F* h_pos23_vs_pos45_y;
	TH2F* h_pos23_vs_pos45_x;
	
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

	TH2F* h_cl_pos_sby2_sby3;
	TH2F* h_cl_pos_sby1_sby2;
	TH2F* h_cl_pos_sby1_sby3;
	TH2F* h_nclus_in_5mm_window;
	TH2F* h_nclus_in_1mm_window;
	TH1F* h_nclus_sby2_and_sby3_5mm;
	TH1F* h_nclus_sby2_and_sby3_1mm;
	TH1F* h_res_lay01;
	TH1F* h_res_lay01_5mm;
	TH2F* h_pos_lay0_pos_lay1;
	

	TH2F* h_diffpos_strip_index_lay0;
	TH2F* h_diffpos_strip_index_lay1;
	TH2F* h_diffpos_strip_index_lay2index_lay0;
	TH2F* h_diffpos_strip_index_lay3index_lay0;
	TH2F* h_diffpos_strip_index_lay2index_lay1;
	TH2F* h_diffpos_strip_index_lay3index_lay1;

	TH2F* h_diffpos_lay2;
	TH2F* h_diffpos_lay3;

	TH2F* h_beamProfile;
	TH2F* h_beamProfile_small;
	TH2F* h_beamProfile_small_2;
	TH2F* h_beamProfile_small_3;
	
	TH2F* h_diffpos_sby2;
	TH2F* h_diffpos_sby3;
	TH2F* h_res_SBY2_SBY3_vs_SBY2;
	TH2F* h_res_SBY2_SBY3_vs_SBY3;

	TH2F* h_diffpos_lay0;
	TH2F* h_diffpos_lay1;
	TH2F* h_diffpos_stereolay0;
	TH2F* h_diffpos_stereolay1;

	TH2F* h_diffstereolay1_phi;
	TH2F* h_diffstereolay1_phi0;

	TH1F* h_hit_pos_sbx;
	TH1F* h_hit_pos_sby1;
	TH1F* h_hit_pos_sby2;
	TH1F* h_hit_pos_sby3;

	TH1F* h_res_secondCoordinate;

	TH1F* h_angle;
	TH1F* h_chi2ndf;
	TH1F* h_chi2;
	TH1F* h_prob;
	TH1F* h_prob_mult_tracks;
	TH1F* h_angle_only_scluster_sby2;
	TH1F* h_angle_only_scluster_sby3;
	TH1F* h_dangle_only_scluster_sby2sby3;

	TH2F* h_sby3_minus_sby2_vs_sby2;
	TH2F* h_pos_stereo_minus_sby1_vs_sby1;
	TH2F* h_pos_stereo_minus_eta_in_vs_eta_in;
	TH2F* h_pos_eta_out_minus_eta_in_vs_eta_in;

	TH2F* h_sby1_minus_eta_out_vs_eta_out;
	TH2F* h_eta_in_minus_eta_out_vs_eta_out;
	TH2F* h_stereo_minus_eta_out_vs_eta_out;
	TH2F* h_pos_sby2_minus_eta_out_vs_eta_out;
	TH2F* h_pos_sby3_minus_eta_out_vs_eta_out;

	TH2F* h_eta_out_minus_sby1_vs_sby1;
	TH2F* h_stereo_minus_sby1_vs_sby1;
	TH2F* h_eta_in_minus_sby1_vs_sby1;
	TH2F* h_sby2_minus_sby1_vs_sby1;
	TH2F* h_sby3_minus_sby1_vs_sby1;

	TH2F* h_eta_out_minus_sby1_vs_sby1_corr;
	TH2F* h_stereo_minus_sby1_vs_sby1_corr;
	TH2F* h_eta_in_minus_sby1_vs_sby1_corr;
	TH2F* h_sby2_minus_sby1_vs_sby1_corr;
	TH2F* h_sby3_minus_sby1_vs_sby1_corr;

	TH2F* h_sby1_minus_eta_out_vs_pos_eta_out_corr;
	TH2F* h_sby1_minus_pos_eta_in_vs_pos_eta_in_corr;
	TH2F* h_sby1_minus_stereo_vs_stereo_corr;
	TH2F* h_sby1_minus_sby2_vs_sby2_corr;
	TH2F* h_sby1_minus_sby3_vs_sby3_corr;

	TH2F* h_sby1_minus_eta_out_vs_pos_eta_out;
	TH2F* h_sby1_minus_pos_eta_in_vs_pos_eta_in;
	TH2F* h_sby1_minus_stereo_vs_stereo;
	TH2F* h_sby1_minus_sby2_vs_sby2;
	TH2F* h_sby1_minus_sby3_vs_sby3;

	TH2F* h_bsy1_minus_etaout_vs_stereo_x;
	TH2F* h_bsy1_minus_etain_vs_stereo_x;
	TH2F* h_bsy1_minus_stereoin_vs_stereo_x;
	TH2F* h_bsy1_minus_stereoout_vs_stereo_x;


	TH1F* h_cluster_position_corr_IP1;
	TH1F* h_cluster_position_corr_IP2;
	TH1F* h_cluster_position_corr_stereo;
	TH1F* h_cluster_position_corr_SBY2;
	TH1F* h_cluster_position_corr_SBY3;

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

	TH1F* h_residuals_single_tracks_anglecut;
	TH1F* h_residuals_mult_tracks_anglecut;
	TH1F* h_residuals_accepted_tracks_anglecut;

	TH1F* h_cutflow;

	TH1F* h_sclusters_lay0_lay1_1mm;
	TH1F* h_sclusters_lay0_lay1_2mm;
	TH1F* h_sclusters_lay0_lay1_5mm;
	TH1F* h_sclusters_lay0_lay1_40mm;

	TH1F* h_Nscluster_in_1deg_sby2_sby3;
	TH1F* h_Nscluster_theonly_out_1deg_sby2_sby3;
	TH1F* h_Nscluster_morethan1_in_and_out_sby2_sby3;
	TH1F* h_Nscluster_morethan1_in_1deg_sby2_sby3;
	TH1F* h_Nscluster_morethan1_out_1deg_sby2_sby3;

	TH1F* h_nmanyrack_events;
	TH1F* h_singletrack_events;

	TH1F* h_d_track_etaout;
	TH1F* h_d_track_etain;
	TH1F* h_d_track_stereoin;
	TH1F* h_d_track_stereoout;
	TH1F* h_d_track_ystereo;
	TH1F* h_d_track_lay2;
	TH1F* h_d_track_lay3;

	TH2F* h_stereo_y_corr_vs_uncorr;
	TH2F* h_stereo_x_corr_vs_uncorr;
	TH2F* h_diffpos_stereo_uncorr_lay1;
	TH2F* h_beamProfile_uncorr;
	TH2F* h_diffpos_lay2_uncorr;
	TH2F* h_diffpos_lay3_uncorr;

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

	h_d_track_etaout = new TH1F("h_d_track_etaout", "Distance from track - eta_out" ,1000, -10, 10);
	h_d_track_etaout->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etain = new TH1F("h_d_track_etain", "Distance from track - eta_in" ,1000, -10, 10);
	h_d_track_etain->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_stereoin = new TH1F("h_d_track_stereoin", "Distance from track - stereo_in" ,1000, -10, 10);
	h_d_track_stereoin->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay2 = new TH1F("h_d_track_lay2", "Distance from track - stereo_in" ,1000, -10, 10);
	h_d_track_lay2->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay3 = new TH1F("h_d_track_lay3", "Distance from track - stereo_in" ,1000, -10, 10);
	h_d_track_lay3->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_stereoout = new TH1F("h_d_track_stereoout", "Distance from track - stereo_out" ,1000, -10, 10);
	h_d_track_stereoout->GetXaxis()->SetTitle("distance [mm]");

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

	h_sclusters_lay0_lay1_1mm = new TH1F("h_sclusters_lay0_lay1_1mm", "superclusters 1mm" ,20, 0, 20);
	h_sclusters_lay0_lay1_1mm->GetXaxis()->SetTitle("Nsclusters");

	h_sclusters_lay0_lay1_2mm = new TH1F("h_sclusters_lay0_lay1_2mm", "superclusters 2mm" ,20, 0, 20);
	h_sclusters_lay0_lay1_2mm->GetXaxis()->SetTitle("Nsclusters");

	h_sclusters_lay0_lay1_5mm = new TH1F("h_sclusters_lay0_lay1_5mm", "superclusters 5mm" ,20, 0, 20);
	h_sclusters_lay0_lay1_5mm->GetXaxis()->SetTitle("Nsclusters");

	h_sclusters_lay0_lay1_40mm = new TH1F("h_sclusters_lay0_lay1_40mm", "superclusters 40mm" ,20, 0, 20);
	h_sclusters_lay0_lay1_40mm->GetXaxis()->SetTitle("Nsclusters");

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
	h_residuals_single_tracks_anglecut = new TH1F("h_residuals_single_tracks_anglecut", "residuals SBY2-SBY3 |#theta|<0.4deg (accepted single tracks)" ,1000, -10, 10);
	h_residuals_mult_tracks_anglecut = new TH1F("h_residuals_mult_tracks_anglecut", "residuals SBY2-SBY3 |#theta|<0.4deg (accepted multiple tracks)" ,1000, -10, 10);
	h_residuals_accepted_tracks_anglecut = new TH1F("h_residuals_accepted_tracks_anglecut", "residuals SBY2-SBY3 |#theta|<0.4deg (all accepted tracks)" ,1000, -10, 10);

	h_cluster_position_corr_IP1 = new TH1F("h_cluster_position_corr_IP1", "IP1" ,1000, 1000, 2000);
	h_cluster_position_corr_IP2 = new TH1F("h_cluster_position_corr_IP2", "IP2", 1000, 1000, 2000);
	h_cluster_position_corr_stereo = new TH1F("h_cluster_position_corr_stereo", "stereo", 1000, 1000, 2000);
	h_cluster_position_corr_SBY2 = new TH1F("h_cluster_position_corr_SBY2", "SBY2", 1000, 1000, 2000);
	h_cluster_position_corr_SBY3 = new TH1F("h_cluster_position_corr_SBY3", "SBY3", 1000, 1000, 2000);

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

	h_cl_pos_sby2_sby3 = new TH2F("h_cl_pos_sby2_sby3", "Cluster position" ,500, 1300, 1800, 500, 1300, 1800);
	h_cl_pos_sby2_sby3->GetXaxis()->SetTitle("SBY2 [mm]");
	h_cl_pos_sby2_sby3->GetYaxis()->SetTitle("SBY3 [mm]");

	h_cl_pos_sby1_sby2 = new TH2F("h_cl_pos_sby1_sby2", "Cluster position" ,500, 1300, 1800, 500, 1300, 1800);
	h_cl_pos_sby1_sby2->GetXaxis()->SetTitle("SBY1 [mm]");
	h_cl_pos_sby1_sby2->GetYaxis()->SetTitle("SBY2 [mm]");

	h_cl_pos_sby1_sby3 = new TH2F("h_cl_pos_sby1_sby3", "Cluster position" ,500, 1300, 1800, 500, 1300, 1800);
	h_cl_pos_sby1_sby3->GetXaxis()->SetTitle("SBY1 [mm]");
	h_cl_pos_sby1_sby3->GetYaxis()->SetTitle("SBY3 [mm]");

	h_nclus_in_5mm_window = new TH2F("h_nclus_in_5mm_window", "Nclusters in +- 5mm", 40, 0, 40, 40, 0, 40);
	h_nclus_in_5mm_window->GetXaxis()->SetTitle("SBY2");
	h_nclus_in_5mm_window->GetYaxis()->SetTitle("SBY3");

	h_nclus_in_1mm_window = new TH2F("h_nclus_in_1mm_window", "Nclusters in +- 1mm", 40, 0, 40, 40, 0, 40);
	h_nclus_in_1mm_window->GetXaxis()->SetTitle("SBY2");
	h_nclus_in_1mm_window->GetYaxis()->SetTitle("SBY3");


	h_sby2_minus_sby1_vs_sby1 = new TH2F("h_sby2_minus_sby1_vs_sby1", "BS2Y-BS1Y vs BS1Y position", 3000, 1300, 2000, 1000, -350, 350);
	h_sby2_minus_sby1_vs_sby1->GetXaxis()->SetTitle("BS1Y cluster position [mm]");
	h_sby2_minus_sby1_vs_sby1->GetYaxis()->SetTitle("BS2Y-BS1Y cluster position [mm]");

	h_sby2_minus_sby1_vs_sby1_corr = new TH2F("h_sby2_minus_sby1_vs_sby1_corr", "BS2Y-BS1Y vs BS1Y position (corrected)", 3000, 1300, 2000, 1000, -350, 350);
	h_sby2_minus_sby1_vs_sby1_corr->GetXaxis()->SetTitle("BS1Y cluster position [mm]");
	h_sby2_minus_sby1_vs_sby1_corr->GetYaxis()->SetTitle("BS2Y-BS1Y cluster position [mm]");

	h_sby3_minus_sby2_vs_sby2 = new TH2F("h_sby3_minus_sby2_vs_sby2", "BS3Y-BS2Y vs BS2Y position", 3000, 1300, 2000, 1000, -350, 350);
	h_sby3_minus_sby2_vs_sby2->GetXaxis()->SetTitle("BS2Y cluster position [mm]");
	h_sby3_minus_sby2_vs_sby2->GetYaxis()->SetTitle("BS3Y-BS2Y cluster position [mm]");

	h_sby3_minus_sby1_vs_sby1 = new TH2F("h_sby3_minus_sby1_vs_sby1", "BS3Y-BS1Y vs BS1Y position", 3000, 1300, 2000, 1000, -350, 350);
	h_sby3_minus_sby1_vs_sby1->GetXaxis()->SetTitle("BS1Y cluster position [mm]");
	h_sby3_minus_sby1_vs_sby1->GetYaxis()->SetTitle("BS3Y-BS1Y cluster position [mm]");

	h_sby3_minus_sby1_vs_sby1_corr = new TH2F("h_sby3_minus_sby1_vs_sby1_corr", "BS3Y-BS1Y vs BS1Y position (corrected)", 3000, 1300, 2000, 1000, -350, 350);
	h_sby3_minus_sby1_vs_sby1_corr->GetXaxis()->SetTitle("BS1Y cluster position [mm]");
	h_sby3_minus_sby1_vs_sby1_corr->GetYaxis()->SetTitle("BS3Y-BS1Y cluster position [mm]");

	h_eta_out_minus_sby1_vs_sby1 = new TH2F("h_eta_out_minus_sby1_vs_sby1", "IP1-BS1Y vs BS1Y position", 3000, 1300, 2000, 1000, 0, 700);
	h_eta_out_minus_sby1_vs_sby1->GetXaxis()->SetTitle("BS1Y cluster position [mm]");
	h_eta_out_minus_sby1_vs_sby1->GetYaxis()->SetTitle("IP1-BS1Y cluster position [mm]");

	h_eta_out_minus_sby1_vs_sby1_corr = new TH2F("h_eta_out_minus_sby1_vs_sby1_corr", "IP1-BS1Y vs BS1Y position (corrected)", 3000, 1300, 2000, 1000, 0, 700);
	h_eta_out_minus_sby1_vs_sby1_corr->GetXaxis()->SetTitle("BS1Y cluster position [mm]");
	h_eta_out_minus_sby1_vs_sby1_corr->GetYaxis()->SetTitle("IP1-BS1Y cluster position [mm]");

	h_eta_in_minus_sby1_vs_sby1 = new TH2F("h_eta_in_minus_sby1_vs_sby1", "IP2-BS1Y vs BS1Y position", 3000, 1300, 2000, 1000, 0, 700);
	h_eta_in_minus_sby1_vs_sby1->GetXaxis()->SetTitle("BS1Y cluster position [mm]");
	h_eta_in_minus_sby1_vs_sby1->GetYaxis()->SetTitle("IP2-BS1Y cluster position [mm]");

	h_eta_in_minus_sby1_vs_sby1_corr = new TH2F("h_eta_in_minus_sby1_vs_sby1_corr", "IP2-BS1Y vs BS1Y position (corrected)", 3000, 1300, 2000, 1000, 0, 700);
	h_eta_in_minus_sby1_vs_sby1_corr->GetXaxis()->SetTitle("BS1Y cluster position [mm]");
	h_eta_in_minus_sby1_vs_sby1_corr->GetYaxis()->SetTitle("IP2-BS1Y cluster position [mm]");

	h_stereo_minus_sby1_vs_sby1 = new TH2F("h_stereo_minus_sby1_vs_sby1", "stereo-BS1Y vs BS1Y position", 3000, 1300, 2000, 1000, 0, 700);
	h_stereo_minus_sby1_vs_sby1->GetXaxis()->SetTitle("BS1Y cluster position [mm]");
	h_stereo_minus_sby1_vs_sby1->GetYaxis()->SetTitle("stereo-BS1Y cluster position [mm]");
	
	h_stereo_minus_sby1_vs_sby1_corr = new TH2F("h_stereo_minus_sby1_vs_sby1_corr", "stereo-BS1Y vs BS1Y position (corrected)", 3000, 1300, 2000, 1000, 0, 700);
	h_stereo_minus_sby1_vs_sby1_corr->GetXaxis()->SetTitle("BS1Y cluster position [mm]");
	h_stereo_minus_sby1_vs_sby1_corr->GetYaxis()->SetTitle("stereo-BS1Y cluster position [mm]");

	h_pos_stereo_minus_sby1_vs_sby1 = new TH2F("h_pos_stereo_minus_sby1_vs_sby1", "stereo-BS1Y vs BS1Y position", 3000, 1300, 2000, 1000, 0, 700);
	h_pos_stereo_minus_sby1_vs_sby1->GetXaxis()->SetTitle("BS1Y cluster position [mm]");
	h_pos_stereo_minus_sby1_vs_sby1->GetYaxis()->SetTitle("stereo-BS1Y cluster position [mm]");

	h_pos_stereo_minus_eta_in_vs_eta_in = new TH2F("h_pos_stereo_minus_eta_in_vs_eta_in", "stereo-IP2 vs IP2 position", 3000, 1600, 2300, 1000, -350, 350);
	h_pos_stereo_minus_eta_in_vs_eta_in->GetXaxis()->SetTitle("IP2 cluster position [mm]");
	h_pos_stereo_minus_eta_in_vs_eta_in->GetYaxis()->SetTitle("stereo-IP2 cluster position [mm]");

	h_pos_eta_out_minus_eta_in_vs_eta_in = new TH2F("h_pos_eta_out_minus_eta_in_vs_eta_in", "IP1-IP2 vs IP2 position", 3000, 1600, 2300, 1000, -350, 350);
	h_pos_eta_out_minus_eta_in_vs_eta_in->GetXaxis()->SetTitle("IP2 cluster position [mm]");
	h_pos_eta_out_minus_eta_in_vs_eta_in->GetYaxis()->SetTitle("IP1-IP2 cluster position [mm]");

///start here

	h_diffpos_lay2 = new TH2F("h_diffpos_lay2", "BS1Y-IP3 vs IP3 position ", 700, 1300, 2000, 700, -350, 350);
	h_diffpos_lay2->GetXaxis()->SetTitle("IP3 cluster position [mm]");
	h_diffpos_lay2->GetYaxis()->SetTitle("BS1Y-IP3 cluster position [mm]");
	h_diffpos_lay3 = new TH2F("h_diffpos_lay3", "BS1Y-IP4 vs IP4 position ", 700, 1300, 2000, 700, -350, 350);
	h_diffpos_lay3->GetXaxis()->SetTitle("IP4 cluster position [mm]");
	h_diffpos_lay3->GetYaxis()->SetTitle("BS1Y-IP3 cluster position [mm]");

	h_sby1_minus_eta_out_vs_pos_eta_out_corr = new TH2F("h_sby1_minus_eta_out_vs_pos_eta_out_corr", "BS1Y-IP1 vs IP1 position (corrected)", 700, 1300, 2000, 700, -350, 350);
	h_sby1_minus_eta_out_vs_pos_eta_out_corr->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_sby1_minus_eta_out_vs_pos_eta_out_corr->GetYaxis()->SetTitle("BS1Y-IP1 cluster position [mm]");

	h_sby1_minus_pos_eta_in_vs_pos_eta_in_corr = new TH2F("h_sby1_minus_pos_eta_in_vs_pos_eta_in_corr", "BS1Y-IP2 vs IP2 position (corrected)", 700, 1300, 2000, 700, -350, 350);
	h_sby1_minus_pos_eta_in_vs_pos_eta_in_corr->GetXaxis()->SetTitle("IP2 cluster position [mm]");
	h_sby1_minus_pos_eta_in_vs_pos_eta_in_corr->GetYaxis()->SetTitle("BS1Y-IP2 cluster position [mm]");

	h_sby1_minus_stereo_vs_stereo_corr = new TH2F("h_sby1_minus_stereo_vs_stereo_corr", "BS1Y-stereo vs stereo position (corrected)", 700, 1300, 2000, 700, -350, 350);
	h_sby1_minus_stereo_vs_stereo_corr->GetXaxis()->SetTitle("stereo cluster position [mm]");
	h_sby1_minus_stereo_vs_stereo_corr->GetYaxis()->SetTitle("BS1Y-stereo cluster position [mm]");

	h_sby1_minus_sby2_vs_sby2_corr = new TH2F("h_sby1_minus_sby2_vs_sby2_corr", "BS1Y-BSY2 vs BSY2 position (corrected)", 700, 1300, 2000, 700, 350, -350);
	h_sby1_minus_sby2_vs_sby2_corr->GetXaxis()->SetTitle("BSY2 cluster position [mm]");
	h_sby1_minus_sby2_vs_sby2_corr->GetYaxis()->SetTitle("BS1Y-BSY2 cluster position [mm]");

	h_sby1_minus_sby3_vs_sby3_corr = new TH2F("h_sby1_minus_sby3_vs_sby3_corr", "BS1Y-BSY3 vs BSY3 position (corrected)", 700, 1300, 2000, 700, 350, -350);
	h_sby1_minus_sby3_vs_sby3_corr->GetXaxis()->SetTitle("BSY3 cluster position [mm]");
	h_sby1_minus_sby3_vs_sby3_corr->GetYaxis()->SetTitle("BS1Y-BSY3 cluster position [mm]");


	h_sby1_minus_eta_out_vs_pos_eta_out = new TH2F("h_sby1_minus_eta_out_vs_pos_eta_out", "BS1Y-IP1 vs IP1 position", 700, 1800, 2500, 700, -700, 0);
	h_sby1_minus_eta_out_vs_pos_eta_out->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_sby1_minus_eta_out_vs_pos_eta_out->GetYaxis()->SetTitle("BS1Y-IP1 cluster position [mm]");

	h_sby1_minus_pos_eta_in_vs_pos_eta_in = new TH2F("h_sby1_minus_pos_eta_in_vs_pos_eta_in", "BS1Y-IP2 vs IP2 position", 700, 1800, 2500, 700, -700, 0);
	h_sby1_minus_pos_eta_in_vs_pos_eta_in->GetXaxis()->SetTitle("IP2 cluster position [mm]");
	h_sby1_minus_pos_eta_in_vs_pos_eta_in->GetYaxis()->SetTitle("BS1Y-IP2 cluster position [mm]");

	h_sby1_minus_stereo_vs_stereo = new TH2F("h_sby1_minus_stereo_vs_stereo", "BS1Y-stereo vs stereo position", 1000, 1300, 2300, 700, -700, 0);
	h_sby1_minus_stereo_vs_stereo->GetXaxis()->SetTitle("stereo cluster position [mm]");
	h_sby1_minus_stereo_vs_stereo->GetYaxis()->SetTitle("BS1Y-stereo cluster position [mm]");

	h_sby1_minus_sby2_vs_sby2 = new TH2F("h_sby1_minus_sby2_vs_sby2", "BS1Y-BSY2 vs BSY2 position", 700, 1300, 2000, 700, 350, -350);
	h_sby1_minus_sby2_vs_sby2->GetXaxis()->SetTitle("BSY2 cluster position [mm]");
	h_sby1_minus_sby2_vs_sby2->GetYaxis()->SetTitle("BS1Y-BSY2 cluster position [mm]");

	h_sby1_minus_sby3_vs_sby3 = new TH2F("h_sby1_minus_sby3_vs_sby3", "BS1Y-BSY3 vs BSY3 position", 700, 1300, 2000, 700, 350, -350);
	h_sby1_minus_sby3_vs_sby3->GetXaxis()->SetTitle("BSY3 cluster position [mm]");
	h_sby1_minus_sby3_vs_sby3->GetYaxis()->SetTitle("BS1Y-BSY3 cluster position [mm]");

	h_bsy1_minus_etaout_vs_stereo_x = new TH2F("h_bsy1_minus_etaout_vs_stereo_x", "BS1Y-IP1 vs stereo_x position", 700, 1300, 2000, 700, 350, -350);
	h_bsy1_minus_etaout_vs_stereo_x->GetXaxis()->SetTitle("stereo_x cluster position [mm]");
	h_bsy1_minus_etaout_vs_stereo_x->GetYaxis()->SetTitle("BS1Y-IP1 cluster position [mm]");

	h_bsy1_minus_etain_vs_stereo_x = new TH2F("h_bsy1_minus_etain_vs_stereo_x", "BS1Y-IP2 vs stereo_x position", 700, 1300, 2000, 700, 350, -350);
	h_bsy1_minus_etain_vs_stereo_x->GetXaxis()->SetTitle("stereo_x cluster position [mm]");
	h_bsy1_minus_etain_vs_stereo_x->GetYaxis()->SetTitle("BS1Y-IP2 cluster position [mm]");

	h_bsy1_minus_stereoin_vs_stereo_x = new TH2F("h_bsy1_minus_stereoin_vs_stereo_x", "BS1Y-IP3 vs stereo_x position", 700, 1300, 2000, 700, 350, -350);
	h_bsy1_minus_stereoin_vs_stereo_x->GetXaxis()->SetTitle("stereo_x cluster position [mm]");
	h_bsy1_minus_stereoin_vs_stereo_x->GetYaxis()->SetTitle("BS1Y-IP3 cluster position [mm]");

	h_bsy1_minus_stereoout_vs_stereo_x = new TH2F("h_bsy1_minus_stereoout_vs_stereo_x", "BS1Y-IP4 vs stereo_x position", 700, 1300, 2000, 700, 350, -350);
	h_bsy1_minus_stereoout_vs_stereo_x->GetXaxis()->SetTitle("stereo_x cluster position [mm]");
	h_bsy1_minus_stereoout_vs_stereo_x->GetYaxis()->SetTitle("BS1Y-IP4 cluster position [mm]");

//end here
	h_sby1_minus_eta_out_vs_eta_out = new TH2F("h_sby1_minus_eta_out_vs_eta_out", "SBY1-IP1 vs IP1 position", 700, 1600, 2300, 2000, -1000, 1000);
	h_sby1_minus_eta_out_vs_eta_out->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_sby1_minus_eta_out_vs_eta_out->GetYaxis()->SetTitle("SBY1-IP1 cluster position [mm]");

	h_eta_in_minus_eta_out_vs_eta_out = new TH2F("h_eta_in_minus_eta_out_vs_eta_out", "IP2-IP1 vs IP1 position", 700, 1600, 2300, 2000, -1000, 1000);
	h_eta_in_minus_eta_out_vs_eta_out->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_eta_in_minus_eta_out_vs_eta_out->GetYaxis()->SetTitle("IP2-IP1 cluster position [mm]");

	h_stereo_minus_eta_out_vs_eta_out = new TH2F("h_stereo_minus_eta_out_vs_eta_out", "stereo-IP1 vs IP1 position", 700, 1600, 2300, 2000, -1000, 1000);
	h_stereo_minus_eta_out_vs_eta_out->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_stereo_minus_eta_out_vs_eta_out->GetYaxis()->SetTitle("stereo-IP1 cluster position [mm]");

	h_pos_sby2_minus_eta_out_vs_eta_out = new TH2F("h_pos_sby2_minus_eta_out_vs_eta_out", "SBY2-IP1 vs IP1 position", 700, 1600, 2300, 2000, -1000, 1000);
	h_pos_sby2_minus_eta_out_vs_eta_out->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_pos_sby2_minus_eta_out_vs_eta_out->GetYaxis()->SetTitle("SBY2-IP1 cluster position [mm]");

	h_pos_sby3_minus_eta_out_vs_eta_out = new TH2F("h_pos_sby3_minus_eta_out_vs_eta_out", "SBY3-IP1 vs IP1 position", 700, 1600, 2300, 2000, -1000, 1000);
	h_pos_sby3_minus_eta_out_vs_eta_out->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_pos_sby3_minus_eta_out_vs_eta_out->GetYaxis()->SetTitle("SBY3-IP1 cluster position [mm]");





	h_res_lay01 = new TH1F("h_res_lay01", "SM1 #Delta(#eta_{in}-#eta_{out})", 100, -2, 2);
	h_res_lay01->GetXaxis()->SetTitle("#Delta y [mm]");

	h_res_secondCoordinate = new TH1F("h_res_secondCoordinate", "#Delta(#phi_{stereo}-y_{SBX})", 1000, -5000, 5000);
	h_res_secondCoordinate->GetXaxis()->SetTitle("#Delta y [mm]");

	h_res_lay01_5mm = new TH1F("h_res_lay01_5mm", "SM1 #Delta(#eta_{in}-#eta_{out})", 100, -2, 2);
	h_res_lay01_5mm->GetXaxis()->SetTitle("#Delta y [mm]");

	h_pos_lay0_pos_lay1 = new TH2F("h_pos_lay0_pos_lay1", "cluster position IP1 vs IP2", 600, 1700, 2300, 600, 1700, 2300);
	h_pos_lay0_pos_lay1->GetXaxis()->SetTitle("y_{IP1} [mm]");
	h_pos_lay0_pos_lay1->GetYaxis()->SetTitle("y_{IP2} [mm]");

	h_diffpos_lay0 = new TH2F("h_diffpos_lay0", "", 100, -2, 2,600,1350,1950);
	h_diffpos_lay0->GetXaxis()->SetTitle("#Delta y(IP2-IP1) [mm]");
	h_diffpos_lay0->GetYaxis()->SetTitle("y_{IP1} [mm]");
	h_diffpos_lay1 = new TH2F("h_diffpos_lay1", "", 100, -2, 2,600,1350,1950);
	h_diffpos_lay1->GetXaxis()->SetTitle("#Delta y(IP2-IP1) [mm]");
	h_diffpos_lay1->GetYaxis()->SetTitle("y_{IP2} [mm]");

	h_res_SBY2_SBY3_vs_SBY2 = new TH2F("h_res_SBY2_SBY3_vs_SBY2", "", 100, -2, 2,600,1350,1950);
	h_res_SBY2_SBY3_vs_SBY2->GetXaxis()->SetTitle("#Delta y(SB2-SB3) [mm]");
	h_res_SBY2_SBY3_vs_SBY2->GetYaxis()->SetTitle("y_{SB2} [mm]");

	h_res_SBY2_SBY3_vs_SBY3 = new TH2F("h_res_SBY2_SBY3_vs_SBY3", "", 100, -2, 2,600,1350,1950);
	h_res_SBY2_SBY3_vs_SBY3->GetXaxis()->SetTitle("#Delta y(SB2-SB3) [mm]");
	h_res_SBY2_SBY3_vs_SBY3->GetYaxis()->SetTitle("y_{SB3} [mm]");

	h_diffpos_strip_index_lay0 = new TH2F("h_diffpos_strip_index_lay0", "", 100, -2, 2,600,2000,3200);
	h_diffpos_strip_index_lay0->GetXaxis()->SetTitle("#Delta y(IP2-IP1) [mm]");
	h_diffpos_strip_index_lay0->GetYaxis()->SetTitle("strip index (IP1)");
	h_diffpos_strip_index_lay1 = new TH2F("h_diffpos_strip_index_lay1", "", 100, -2, 2,600,2000,3200);
	h_diffpos_strip_index_lay1->GetXaxis()->SetTitle("#Delta y(IP2-IP1) [mm]");
	h_diffpos_strip_index_lay1->GetYaxis()->SetTitle("strip index (IP2)");


	h_stereo_y_corr_vs_uncorr = new TH2F("h_stereo_y_corr_vs_uncorr", "", 600, 1350, 1950, 600,1850,2450);
	h_stereo_y_corr_vs_uncorr->GetXaxis()->SetTitle("y_{stereo}^{corr} [mm]");
	h_stereo_y_corr_vs_uncorr->GetYaxis()->SetTitle("y_{stereo} [mm]");

	h_stereo_x_corr_vs_uncorr = new TH2F("h_stereo_x_corr_vs_uncorr", "", 500, -.5, .5, 500,-.5,.5);
	h_stereo_x_corr_vs_uncorr->GetXaxis()->SetTitle("x_{stereo}^{corr} [mm]");
	h_stereo_x_corr_vs_uncorr->GetYaxis()->SetTitle("x_{stereo} [mm]");


	h_diffpos_stereolay0 = new TH2F("h_diffpos_stereolay0", "", 20000, -2, 2, 600,1350,1950);
	h_diffpos_stereolay0->GetXaxis()->SetTitle("#Delta y(stereo - IP1) [mm]");
	h_diffpos_stereolay0->GetYaxis()->SetTitle("y_{stereo} [mm]");

	h_diffpos_stereolay1 = new TH2F("h_diffpos_stereolay1", "", 20000, -2., 2.,600,1350,1950);
	h_diffpos_stereolay1->GetXaxis()->SetTitle("#Delta y(stereo - IP2) [mm]");
	h_diffpos_stereolay1->GetYaxis()->SetTitle("y_{stereo} [mm]");

	h_diffpos_stereo_uncorr_lay1 = new TH2F("h_diffpos_stereo_uncorr_lay1", "", 600,1800,2500, 100, -10, 10);
	h_diffpos_stereo_uncorr_lay1->GetYaxis()->SetTitle("#Delta y(stereo - IP2) [mm]");
	h_diffpos_stereo_uncorr_lay1->GetXaxis()->SetTitle("y_{stereo} [mm]");

	h_diffstereolay1_phi = new TH2F("h_diffstereolay1_phi", "", 100, -10, 10,600,-300,300);
	h_diffstereolay1_phi->GetXaxis()->SetTitle("#Delta y(IP2 - stereo) [mm]");
	h_diffstereolay1_phi->GetYaxis()->SetTitle("x [mm]");

	h_diffstereolay1_phi0 = new TH2F("h_diffstereolay1_phi0", "", 100, -10, 10,600,0,600);
	h_diffstereolay1_phi0->GetXaxis()->SetTitle("#Delta y(IP2 - IP3) [mm]");
	h_diffstereolay1_phi0->GetYaxis()->SetTitle("y_{IP3} [mm]");

	h_diffpos_lay2_uncorr = new TH2F("h_diffpos_lay2_uncorr", "", 1024,2048,3072,100, -10, 0);
	h_diffpos_lay2_uncorr->GetYaxis()->SetTitle("#Delta y(IP2 - IP3) [mm]");
	h_diffpos_lay2_uncorr->GetXaxis()->SetTitle("IP3 strip index");

	h_diffpos_lay3_uncorr = new TH2F("h_diffpos_lay3_uncorr", "", 1024,2048,3072, 100, 0, 10);
	h_diffpos_lay3_uncorr->GetYaxis()->SetTitle("#Delta y(IP2 - IP4) [mm]");
	h_diffpos_lay3_uncorr->GetXaxis()->SetTitle("IP4 strip index");

	h_diffpos_strip_index_lay2index_lay0 = new TH2F("h_diffpos_strip_index_lay2index_lay0", "", 100, -10, 10,5120,0,5120);
	h_diffpos_strip_index_lay2index_lay0->GetXaxis()->SetTitle("#Delta y(stereo - IP1) [mm]");
	h_diffpos_strip_index_lay2index_lay0->GetYaxis()->SetTitle("strip index (IP3)");

	h_diffpos_strip_index_lay2index_lay1 = new TH2F("h_diffpos_strip_index_lay2index_lay1", "", 100, -10, 10,5120,0,5120);
	h_diffpos_strip_index_lay2index_lay1->GetXaxis()->SetTitle("#Delta y(stereo - IP2) [mm]");
	h_diffpos_strip_index_lay2index_lay1->GetYaxis()->SetTitle("strip index (IP3)");

	h_diffpos_strip_index_lay3index_lay0 = new TH2F("h_diffpos_strip_index_lay3index_lay0", "", 100, -10, 10,5120,0,5120);
	h_diffpos_strip_index_lay3index_lay0->GetXaxis()->SetTitle("#Delta y(stereo - IP1) [mm]");
	h_diffpos_strip_index_lay3index_lay0->GetYaxis()->SetTitle("strip index (IP4)");

	h_diffpos_strip_index_lay3index_lay1 = new TH2F("h_diffpos_strip_index_lay3index_lay1", "", 100, -10, 10,5120,0,5120);
	h_diffpos_strip_index_lay3index_lay1->GetXaxis()->SetTitle("#Delta y(stereo - IP2) [mm]");
	h_diffpos_strip_index_lay3index_lay1->GetYaxis()->SetTitle("strip index (IP4)");

	h_beamProfile = new TH2F("h_beamProfile", "", 500, 1300, 1800, 100, -.4, .4);
	h_beamProfile->GetYaxis()->SetTitle("#phi [mm]");
	h_beamProfile->GetXaxis()->SetTitle("#eta [mm]");

	h_beamProfile_uncorr = new TH2F("h_beamProfile_uncorr", "", 500, 1800, 2500, 100, -.4, .4);
	h_beamProfile_uncorr->GetYaxis()->SetTitle("#phi [mm]");
	h_beamProfile_uncorr->GetXaxis()->SetTitle("#eta [mm]");

	h_beamProfile_small = new TH2F("h_beamProfile_small", "", 800, 1300, 2100, 800, 1300, 2100);
	h_beamProfile_small->GetXaxis()->SetTitle("position_{SBX} [mm]");
	h_beamProfile_small->GetYaxis()->SetTitle("position_{SBY1} [mm]");

	h_beamProfile_small_2 = new TH2F("h_beamProfile_small_2", "", 800, 1300, 2100, 800, 1300, 2100);
	h_beamProfile_small_2->GetXaxis()->SetTitle("position_{SBX} [mm]");
	h_beamProfile_small_2->GetYaxis()->SetTitle("position_{SBY2} [mm]");

	h_beamProfile_small_3 = new TH2F("h_beamProfile_small_3", "", 800, 1300, 2100, 800, 1300, 2100);
	h_beamProfile_small_3->GetXaxis()->SetTitle("position_{SBX} [mm]");
	h_beamProfile_small_3->GetYaxis()->SetTitle("position_{SBY3} [mm]");

	h_diffpos_sby2 = new TH2F("h_diffpos_sby2", "", 100, -5, 5, 450, 1350, 1800);
	h_diffpos_sby2->GetXaxis()->SetTitle("#Delta y(SBY3 - SBY2) [mm]");
	h_diffpos_sby2->GetYaxis()->SetTitle("y_{SBY2} [mm]");

	h_diffpos_sby3 = new TH2F("h_diffpos_sby3", "", 100, -5, 5, 450, 1350, 1800);
	h_diffpos_sby3->GetXaxis()->SetTitle("#Delta y(SBY3 - SBY2) [mm]");
	h_diffpos_sby3->GetYaxis()->SetTitle("y_{SBY3} [mm]");

	h_hit_pos_sbx = new TH1F("h_hit_pos_sbx", "", 5120, 0, 5119);
	h_hit_pos_sbx->GetXaxis()->SetTitle("strip index");

	h_hit_pos_sby1 = new TH1F("h_hit_pos_sby1", "", 5120, 0, 5119);
	h_hit_pos_sby1->GetXaxis()->SetTitle("strip index");

	h_hit_pos_sby2 = new TH1F("h_hit_pos_sby2", "", 5120, 0, 5119);
	h_hit_pos_sby2->GetXaxis()->SetTitle("strip index");

	h_hit_pos_sby3 = new TH1F("h_hit_pos_sby3", "", 5120, 0, 5119);
	h_hit_pos_sby3->GetXaxis()->SetTitle("strip index");
	
/*
	h_nlay_rec_clus = new TH1F("h_nlay_rec_clus", "Number of layers with at least 1 cluster", 8, 0, 8);
	h_nlay_rec_clus->GetXaxis()->SetTitle("N(layers)");
	h_nlay_rec_clus->GetYaxis()->SetTitle("Events");

	h_s_diff_y = new TH1F("h_s_diff_y","stereo layers: y cluster position difference", 50, -500, 500);
	h_pos23_vs_pos45_y = new TH2F("h_pos23_vs_pos45_y","y position stereo-45 vs stereo-23", 100, 900, 3100, 100, -900, 3100);

	h_s_diff_x = new TH1F("h_s_diff_x","stereo layers: x cluster position difference", 100, -1000, 1000);
	h_pos23_vs_pos45_x = new TH2F("h_pos23_vs_pos45_x","x position stereo-45 vs stereo-23", 3000, -3100, 3100, 3000, -3100, 3100);
*/
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

		h_charge_on_track[ilayer] = new TH1F(Form("h_charge_on_track_lay%i",ilayer), Form("SM1 Layer - %i", ilayer), 100024, 0, 100024);
		h_cl_charge[ilayer]->GetXaxis()->SetTitle("cluster pdo [adc counts]");

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

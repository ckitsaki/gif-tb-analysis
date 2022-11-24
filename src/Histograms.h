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

	TH1F* h_clus_positions_corr[4];
	TH1F* h_clus_positions_small_corr[4];

	TH1F* h_clus_positions_corr_ontrack[4];
	TH1F* h_clus_positions_small_corr_ontrack[4];
	TH1F* h_clus_positions_stripIndex_corr_ontrack[4];

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
	TH1F* h_angle_cut;
	TH1F* h_chi2ndf;
	TH1F* h_chi2;
	TH1F* h_prob;
	TH1F* h_prob_mult_tracks;
	TH1F* h_angle_only_scluster_sby2;
	TH1F* h_angle_only_scluster_sby3;
	TH1F* h_dangle_only_scluster_sby2sby3;

	TH1F* h_angle_4points_IP1;
	TH1F* h_chi2ndf_4points_IP1;
	TH1F* h_chi2_4points_IP1;
	TH1F* h_prob_4points_IP1;

	TH1F* h_angle_4points_IP2;
	TH1F* h_chi2ndf_4points_IP2;
	TH1F* h_chi2_4points_IP2;
	TH1F* h_prob_4points_IP2;

	TH1F* h_angle_4points_stereo;
	TH1F* h_chi2ndf_4points_stereo;
	TH1F* h_chi2_4points_stereo;
	TH1F* h_prob_4points_stereo;

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
	

	TH1F* h_d_track_etaout_cut_anglecut;
	TH1F* h_d_track_etain_cut_anglecut;
	TH1F* h_d_track_lay2_cut_anglecut;
	TH1F* h_d_track_lay3_cut_anglecut;
	TH1F* h_d_track_stereo_cut_anglecut;
	TH1F* h_d_track_sby1_anglecut;
	TH1F* h_d_track_sby2_anglecut;
	TH1F* h_d_track_sby3_anglecut;

	TH1F* h_d_track_etaout_4points;
	TH1F* h_d_track_etaout_cut_4points;
	TH1F* h_d_track_etaout_cut_anglecut_4points;

	TH1F* h_d_track_etain_4points;
	TH1F* h_d_track_etain_cut_4points;
	TH1F* h_d_track_etain_cut_anglecut_4points;

	TH1F* h_d_track_stereo_4points;
	TH1F* h_d_track_stereo_cut_4points;
	TH1F* h_d_track_stereo_cut_anglecut_4points;

	TH2F* h_res_SBY2_SBY1_vs_SBY1;
	TH2F* h_res_SBY3_SBY1_vs_SBY1;
	TH2F* h_pos_etain_vs_pos_etaout;
	TH2F* h_real_vs_exp_position_etaout;

	TH1F* h_pos_diff_eta_out_in;
	TH1F* h_pos_diff_eta_out_stereo;
	TH1F* h_pos_diff_eta_in_stereo;

	TH1F* h_pos_diff_eta_out_in_cutangle;
	TH1F* h_pos_diff_eta_out_stereo_cutangle;
	TH1F* h_pos_diff_eta_in_stereo_cutangle;

	TH1F* h_pos_diff_eta_out_in_cutangle_corr;
	TH1F* h_pos_diff_eta_out_stereo_cutangle_corr;
	TH1F* h_pos_diff_eta_in_stereo_cutangle_corr;

	TH2F* h_pos_diff_eta_out_in_vs_angle_cutangle;
	TH2F* h_pos_diff_eta_out_in_vs_angle_cutangle_corr;
	TH2F* h_pos_diff_eta_out_stereo_vs_angle_cutangle;
	TH2F* h_pos_diff_eta_out_stereo_vs_angle_cutangle_corr;
	TH2F* h_pos_diff_eta_in_stereo_vs_angle_cutangle;
	TH2F* h_pos_diff_eta_in_stereo_vs_angle_cutangle_corr;

	TH3F* h_pos_diff_eta_out_in_vs_out_vs_angle_cutangle;
	TH3F* h_pos_diff_eta_out_in_vs_in_vs_angle_cutangle;
	TH3F* h_pos_diff_eta_out_stereo_vs_out_vs_angle_cutangle;
	TH3F* h_pos_diff_eta_out_stereo_vs_stereo_vs_angle_cutangle;
	TH3F* h_pos_diff_eta_in_stereo_vs_in_vs_angle_cutangle;
	TH3F* h_pos_diff_eta_in_stereo_vs_stereo_vs_angle_cutangle;

	TH2F* h_pos_eta_out_vs_angle_cutangle;
	TH2F* h_pos_eta_out_vs_angle_cutangle_corr;
	TH2F* h_pos_stereo_vs_angle_cutangle;
	TH2F* h_pos_stereo_vs_angle_cutangle_corr;
	TH2F* h_pos_eta_in_vs_angle_cutangle;
	TH2F* h_pos_eta_in_vs_angle_cutangle_corr;
	TH2F* h_pos_phi_vs_angle_cutangle;
	TH2F* h_pos_phi_vs_angle_cutangle_corr;

	TH2F* h_pos_eta_out_vs_eta_in_ontrack_cutangle;
	TH2F* h_pos_eta_out_vs_stereo_ontrack_cutangle;
	TH2F* h_pos_eta_in_vs_stereo_ontrack_cutangle;
	//TH2F* h_pos_diff_eta_out_in_vs_evtno_cutangle;
	//TH2F* h_pos_diff_eta_out_stereo_vs_evtno_cutangle;
	//TH2F* h_pos_diff_eta_in_stereo_vs_evtno_cutangle;

	TH2F* h_pos_diff_eta_out_in_vs_eta_in_cutangle;
	TH2F* h_pos_diff_eta_out_in_vs_eta_out_cutangle;

	TH2F* h_d_track_etaout_cut_anglecut_vs_cl_pos;
	TH2F* h_d_track_etain_cut_anglecut_vs_cl_pos;
	TH2F* h_d_track_stereo_cut_anglecut_vs_cl_pos;

	TH2F* h_pos_stereo_vs_angle_cutangle_corrected;

	TH2F* h_align_etaout_ontrack;

	TH2F* h_align_eta_out_ontrack_onxaxis;
	TH2F* h_align_eta_in_ontrack_onxaxis;
	TH2F* h_align_stereo_ontrack_onxaxis;
	TH2F* h_align_SBY2_ontrack_onxaxis;
	TH2F* h_align_SBY3_ontrack_onxaxis;

	TH2F* h_align_eta_out_onxaxis;
	TH2F* h_align_eta_in_onxaxis;
	TH2F* h_align_stereo_onxaxis;
	TH2F* h_align_SBY2_onxaxis;
	TH2F* h_align_SBY3_onxaxis;

	TH2F* h_alignX_phi_vs_SBX0;

	TH2F* h_align_eta_out_SBX0_onxaxis;
	TH2F* h_align_eta_in_SBX0_onxaxis;
	TH2F* h_align_stereo_SBX0_onxaxis;
	TH2F* h_align_SBY2_SBX0_onxaxis;
	TH2F* h_align_SBY3_SBX0_onxaxis;

	TH2F* h_align_eta_out_SBX0_onxaxis_corrected;
	TH2F* h_align_eta_in_SBX0_onxaxis_corrected;
	TH2F* h_align_stereo_SBX0_onxaxis_corrected;
	TH2F* h_align_SBY2_SBX0_onxaxis_corrected;
	TH2F* h_align_SBY3_SBX0_onxaxis_corrected;

	TH2F* h_align_eta_out_SBX0_ontrack_onxaxis;
	TH2F* h_align_eta_in_SBX0_ontrack_onxaxis;
	TH2F* h_align_stereo_SBX0_ontrack_onxaxis;
	TH2F* h_align_SBY2_SBX0_ontrack_onxaxis;
	TH2F* h_align_SBY3_SBX0_ontrack_onxaxis;

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

	h_pos_diff_eta_out_in = new TH1F("h_pos_diff_eta_out_in", "eta_in - eta_out" ,1000, -5, 5);
	h_pos_diff_eta_out_in->GetXaxis()->SetTitle("#Delta(#eta^{in}-#eta^{out}) [mm]");
	
	h_pos_diff_eta_out_stereo = new TH1F("h_pos_diff_eta_out_stereo", "stereo - eta_out" ,1000, -5, 5);
	h_pos_diff_eta_out_stereo->GetXaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{out}) [mm]");
	
	h_pos_diff_eta_in_stereo = new TH1F("h_pos_diff_eta_in_stereo", "stereo - eta_in" ,1000, -5, 5);
	h_pos_diff_eta_in_stereo->GetXaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{out}) [mm]");

	h_pos_diff_eta_out_in_cutangle = new TH1F("h_pos_diff_eta_out_in_cutangle", "eta_in - eta_out" ,1000, -5, 5);
	h_pos_diff_eta_out_in_cutangle->GetXaxis()->SetTitle("#Delta(#eta^{in}-#eta^{out}) [mm]");

	h_pos_diff_eta_out_in_cutangle_corr = new TH1F("h_pos_diff_eta_out_in_cutangle_corr", "eta_in - eta_out (corrected)" ,1000, -5, 5);
	h_pos_diff_eta_out_in_cutangle_corr->GetXaxis()->SetTitle("#Delta(#eta^{in}-#eta^{out}) [mm]");

	h_pos_diff_eta_out_in_vs_angle_cutangle = new TH2F("h_pos_diff_eta_out_in_vs_angle_cutangle", "eta_in - eta_out vs track angle" ,1000, -1, 1, 5000, -5, 5);
	h_pos_diff_eta_out_in_vs_angle_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{in}-#eta^{out}) [mm]");
	h_pos_diff_eta_out_in_vs_angle_cutangle->GetXaxis()->SetTitle("#theta [deg]");

	h_pos_diff_eta_out_in_vs_angle_cutangle_corr = new TH2F("h_pos_diff_eta_out_in_vs_angle_cutangle_corr", "eta_in - eta_out vs track angle" ,1000, -1, 1, 5000, -5, 5);
	h_pos_diff_eta_out_in_vs_angle_cutangle_corr->GetYaxis()->SetTitle("#Delta(#eta^{in}-#eta^{out}) [mm] (corrected)");
	h_pos_diff_eta_out_in_vs_angle_cutangle_corr->GetXaxis()->SetTitle("#theta [deg]");

	h_pos_diff_eta_out_stereo_vs_angle_cutangle = new TH2F("h_pos_diff_eta_out_stereo_vs_angle_cutangle", "stereo - eta_out vs track angle" ,1000, -1, 1, 5000, -5, 5);
	h_pos_diff_eta_out_stereo_vs_angle_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{out}) [mm]");
	h_pos_diff_eta_out_stereo_vs_angle_cutangle->GetXaxis()->SetTitle("#theta [deg]");

	h_pos_diff_eta_out_stereo_vs_angle_cutangle_corr = new TH2F("h_pos_diff_eta_out_stereo_vs_angle_cutangle_corr", "stereo - eta_out vs track angle" ,1000, -1, 1, 5000, -5, 5);
	h_pos_diff_eta_out_stereo_vs_angle_cutangle_corr->GetYaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{out}) [mm] (corrected)");
	h_pos_diff_eta_out_stereo_vs_angle_cutangle_corr->GetXaxis()->SetTitle("#theta [deg]");

	h_pos_diff_eta_in_stereo_vs_angle_cutangle = new TH2F("h_pos_diff_eta_in_stereo_vs_angle_cutangle", "stereo - eta_in vs track angle" ,1000, -1, 1, 5000, -5, 5);
	h_pos_diff_eta_in_stereo_vs_angle_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{in}) [mm]");
	h_pos_diff_eta_in_stereo_vs_angle_cutangle->GetXaxis()->SetTitle("#theta [deg]");

	h_pos_diff_eta_in_stereo_vs_angle_cutangle_corr = new TH2F("h_pos_diff_eta_in_stereo_vs_angle_cutangle_corr", "stereo - eta_in vs track angle " ,1000, -1, 1, 5000, -5, 5);
	h_pos_diff_eta_in_stereo_vs_angle_cutangle_corr->GetYaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{in}) [mm] (corrected)");
	h_pos_diff_eta_in_stereo_vs_angle_cutangle_corr->GetXaxis()->SetTitle("#theta [deg]");


	h_pos_eta_out_vs_angle_cutangle = new TH2F("h_pos_eta_out_vs_angle_cutangle", "eta_out vs track angle" , 600, 1300, 1900, 1000, -1, 1);
	h_pos_eta_out_vs_angle_cutangle->GetXaxis()->SetTitle("#eta^{out} [mm]");
	h_pos_eta_out_vs_angle_cutangle->GetYaxis()->SetTitle("#theta [deg]");

	h_pos_eta_out_vs_angle_cutangle_corr = new TH2F("h_pos_eta_out_vs_angle_cutangle_corr", "eta_out vs track angle " , 600, 1300, 1900, 1000, -1, 1);
	h_pos_eta_out_vs_angle_cutangle_corr->GetXaxis()->SetTitle("#eta^{out} [mm] (corrected)");
	h_pos_eta_out_vs_angle_cutangle_corr->GetYaxis()->SetTitle("#theta [deg]");

	h_pos_eta_in_vs_angle_cutangle = new TH2F("h_pos_eta_in_vs_angle_cutangle", "eta_in vs track angle " , 600, 1300, 1900, 1000, -1, 1);
	h_pos_eta_in_vs_angle_cutangle->GetXaxis()->SetTitle("#eta^{in} [mm]");
	h_pos_eta_in_vs_angle_cutangle->GetYaxis()->SetTitle("#theta [deg]");

	h_pos_eta_in_vs_angle_cutangle_corr = new TH2F("h_pos_eta_in_vs_angle_cutangle_corr", "eta_in vs track angle " , 600, 1300, 1900, 1000, -1, 1);
	h_pos_eta_in_vs_angle_cutangle_corr->GetXaxis()->SetTitle("#eta^{in} [mm] (corrected)");
	h_pos_eta_in_vs_angle_cutangle_corr->GetYaxis()->SetTitle("#theta [deg]");

	h_pos_stereo_vs_angle_cutangle = new TH2F("h_pos_stereo_vs_angle_cutangle", "stereo vs track angle " , 600, 1300, 1900, 1000, -1, 1);
	h_pos_stereo_vs_angle_cutangle->GetXaxis()->SetTitle("#eta^{stereo} [mm]");
	h_pos_stereo_vs_angle_cutangle->GetYaxis()->SetTitle("#theta [deg]");

	h_pos_stereo_vs_angle_cutangle_corrected = new TH2F("h_pos_stereo_vs_angle_cutangle_corrected", "stereo vs track angle " , 600, 1300, 1900, 1000, -1, 1);
	h_pos_stereo_vs_angle_cutangle_corrected->GetXaxis()->SetTitle("#eta^{stereo} [mm] (corrected)");
	h_pos_stereo_vs_angle_cutangle_corrected->GetYaxis()->SetTitle("#theta [deg]");

	h_pos_stereo_vs_angle_cutangle_corr = new TH2F("h_pos_stereo_vs_angle_cutangle_corr", "stereo vs track angle " , 600, 1300, 1900, 1000, -1, 1);
	h_pos_stereo_vs_angle_cutangle_corr->GetXaxis()->SetTitle("#eta^{stereo} [mm] (corrected)");
	h_pos_stereo_vs_angle_cutangle_corr->GetYaxis()->SetTitle("#theta [deg]");

	h_pos_phi_vs_angle_cutangle = new TH2F("h_pos_phi_vs_angle_cutangle", "phi vs track angle " , 4000, -.4, .4, 1000, -1, 1);
	h_pos_phi_vs_angle_cutangle->GetXaxis()->SetTitle("#phi [mm]");
	h_pos_phi_vs_angle_cutangle->GetYaxis()->SetTitle("#theta [deg]");

	h_pos_phi_vs_angle_cutangle_corr = new TH2F("h_pos_phi_vs_angle_cutangle_corr", "phi vs track angle " , 100, -.4, .4, 1000, -1, 1);
	h_pos_phi_vs_angle_cutangle_corr->GetXaxis()->SetTitle("#phi [mm] (corrected)");
	h_pos_phi_vs_angle_cutangle_corr->GetYaxis()->SetTitle("#theta [deg]");

//----->
	h_pos_diff_eta_out_in_vs_out_vs_angle_cutangle = new TH3F("h_pos_diff_eta_out_in_vs_out_vs_angle_cutangle", "eta_in - eta_out vs eta_out vs track angle" , 600, 1300, 1900, 5000, -1., 1., 30, -.4, .4);
	h_pos_diff_eta_out_in_vs_out_vs_angle_cutangle->GetXaxis()->SetTitle("#eta^{out} [mm]");
	h_pos_diff_eta_out_in_vs_out_vs_angle_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{in}-#eta^{out}) [mm]");
	h_pos_diff_eta_out_in_vs_out_vs_angle_cutangle->GetZaxis()->SetTitle("#theta [deg]");

	h_pos_diff_eta_out_in_vs_in_vs_angle_cutangle = new TH3F("h_pos_diff_eta_out_in_vs_in_vs_angle_cutangle", "eta_in - eta_out vs eta_in vs track angle" , 600, 1300, 1900, 5000, -1., 1., 30, -.4, .4);
	h_pos_diff_eta_out_in_vs_in_vs_angle_cutangle->GetXaxis()->SetTitle("#eta^{in} [mm]");
	h_pos_diff_eta_out_in_vs_in_vs_angle_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{in}-#eta^{out}) [mm]");
	h_pos_diff_eta_out_in_vs_in_vs_angle_cutangle->GetZaxis()->SetTitle("#theta [deg]");

	h_pos_diff_eta_out_stereo_vs_out_vs_angle_cutangle = new TH3F("h_pos_diff_eta_out_stereo_vs_out_vs_angle_cutangle", "stereo - eta_out vs eta_out vs track angle" , 600, 1300, 1900, 5000, -1., 1., 30, -.4, .4);
	h_pos_diff_eta_out_stereo_vs_out_vs_angle_cutangle->GetXaxis()->SetTitle("#eta^{out} [mm]");
	h_pos_diff_eta_out_stereo_vs_out_vs_angle_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{out}) [mm]");
	h_pos_diff_eta_out_stereo_vs_out_vs_angle_cutangle->GetZaxis()->SetTitle("#theta [deg]");

	h_pos_diff_eta_out_stereo_vs_stereo_vs_angle_cutangle = new TH3F("h_pos_diff_eta_out_stereo_vs_stereo_vs_angle_cutangle", "stereo - eta_out vs stereo vs track angle" , 600, 1300, 1900, 5000, -1., 1., 30, -.4, .4);
	h_pos_diff_eta_out_stereo_vs_stereo_vs_angle_cutangle->GetXaxis()->SetTitle("#eta^{stereo} [mm]");
	h_pos_diff_eta_out_stereo_vs_stereo_vs_angle_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{out}) [mm]");
	h_pos_diff_eta_out_stereo_vs_stereo_vs_angle_cutangle->GetZaxis()->SetTitle("#theta [deg]");
	
	h_pos_diff_eta_in_stereo_vs_in_vs_angle_cutangle = new TH3F("h_pos_diff_eta_in_stereo_vs_in_vs_angle_cutangle", "stereo - eta_in vs eta_in vs track angle" , 600, 1300, 1900, 5000, -1., 1., 30, -.4, .4);
	h_pos_diff_eta_in_stereo_vs_in_vs_angle_cutangle->GetXaxis()->SetTitle("#eta^{in} [mm]");
	h_pos_diff_eta_in_stereo_vs_in_vs_angle_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{in}) [mm]");
	h_pos_diff_eta_in_stereo_vs_in_vs_angle_cutangle->GetZaxis()->SetTitle("#theta [deg]");

	h_pos_diff_eta_in_stereo_vs_stereo_vs_angle_cutangle = new TH3F("h_pos_diff_eta_in_stereo_vs_stereo_vs_angle_cutangle", "stereo - eta_in vs stereo vs track angle" , 600, 1300, 1900, 5000, -1., 1., 30, -.4, .4);
	h_pos_diff_eta_in_stereo_vs_stereo_vs_angle_cutangle->GetXaxis()->SetTitle("#eta^{stereo} [mm]");
	h_pos_diff_eta_in_stereo_vs_stereo_vs_angle_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{in}) [mm]");
	h_pos_diff_eta_in_stereo_vs_stereo_vs_angle_cutangle->GetZaxis()->SetTitle("#theta [deg]");
	
//----->
	h_pos_eta_out_vs_eta_in_ontrack_cutangle = new TH2F("h_pos_eta_out_vs_eta_in_ontrack_cutangle", "eta_out vs eta_in" ,600, 1300, 1900, 600, 1300, 1900);
	h_pos_eta_out_vs_eta_in_ontrack_cutangle->GetYaxis()->SetTitle("#eta^{out} [mm]");
	h_pos_eta_out_vs_eta_in_ontrack_cutangle->GetXaxis()->SetTitle("#eta^{in} [mm]");

	h_pos_eta_out_vs_stereo_ontrack_cutangle = new TH2F("h_pos_eta_out_vs_stereo_ontrack_cutangle", "eta_out vs stereo" ,600, 1300, 1900, 600, 1300, 1900);
	h_pos_eta_out_vs_stereo_ontrack_cutangle->GetYaxis()->SetTitle("#eta^{out} [mm]");
	h_pos_eta_out_vs_stereo_ontrack_cutangle->GetXaxis()->SetTitle("#eta^{stereo} [mm]");

	h_pos_eta_in_vs_stereo_ontrack_cutangle = new TH2F("h_pos_eta_in_vs_stereo_ontrack_cutangle", "eta_in vs stereo" ,600, 1300, 1900, 600, 1300, 1900);
	h_pos_eta_in_vs_stereo_ontrack_cutangle->GetYaxis()->SetTitle("#eta^{in} [mm]");
	h_pos_eta_in_vs_stereo_ontrack_cutangle->GetXaxis()->SetTitle("#eta^{stereo} [mm]");


	////
/*
	h_pos_diff_eta_out_in_vs_evtno_cutangle = new TH2F("h_pos_diff_eta_out_in_vs_evtno_cutangle", "eta_in - eta_out vs evtno" ,310000, 0, 310000, 1000, -5, 5);
	h_pos_diff_eta_out_in_vs_evtno_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{in}-#eta^{out}) [mm]");
	h_pos_diff_eta_out_in_vs_evtno_cutangle->GetXaxis()->SetTitle("evtno");

	h_pos_diff_eta_out_stereo_vs_evtno_cutangle = new TH2F("h_pos_diff_eta_out_stereo_vs_evtno_cutangle", "stereo - eta_out vs evtno" ,310000, 0, 310000, 1000, -5, 5);
	h_pos_diff_eta_out_stereo_vs_evtno_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{out}) [mm]");
	h_pos_diff_eta_out_stereo_vs_evtno_cutangle->GetXaxis()->SetTitle("evtno");

	h_pos_diff_eta_in_stereo_vs_evtno_cutangle = new TH2F("h_pos_diff_eta_in_stereo_vs_evtno_cutangle", "stereo - eta_in vs evtno" ,310000, 0, 310000, 1000, -5, 5);
	h_pos_diff_eta_in_stereo_vs_evtno_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{in}) [mm]");
	h_pos_diff_eta_in_stereo_vs_evtno_cutangle->GetXaxis()->SetTitle("evtno");
*/
	/////

	h_pos_diff_eta_out_in_vs_eta_in_cutangle = new TH2F("h_pos_diff_eta_out_in_vs_eta_in_cutangle", "eta_in - eta_out" ,600, 1300, 1900 ,1000, -5, 5);
	h_pos_diff_eta_out_in_vs_eta_in_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{in}-#eta^{out}) [mm]");
	h_pos_diff_eta_out_in_vs_eta_in_cutangle->GetXaxis()->SetTitle("IP2 cluster position [mm]");

	h_pos_diff_eta_out_in_vs_eta_out_cutangle = new TH2F("h_pos_diff_eta_out_in_vs_eta_out_cutangle", "eta_in - eta_out" ,600, 1300, 1900 ,1000, -5, 5);
	h_pos_diff_eta_out_in_vs_eta_out_cutangle->GetYaxis()->SetTitle("#Delta(#eta^{in}-#eta^{out}) [mm]");
	h_pos_diff_eta_out_in_vs_eta_out_cutangle->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	
	h_pos_diff_eta_out_stereo_cutangle = new TH1F("h_pos_diff_eta_out_stereo_cutangle", "stereo - eta_out" ,1000, -5, 5);
	h_pos_diff_eta_out_stereo_cutangle->GetXaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{out}) [mm]");

	h_pos_diff_eta_out_stereo_cutangle_corr = new TH1F("h_pos_diff_eta_out_stereo_cutangle_corr", "stereo - eta_out (corrected)" ,1000, -5, 5);
	h_pos_diff_eta_out_stereo_cutangle_corr->GetXaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{out}) [mm]");
	
	h_pos_diff_eta_in_stereo_cutangle = new TH1F("h_pos_diff_eta_in_stereo_cutangle", "stereo - eta_in" ,1000, -5, 5);
	h_pos_diff_eta_in_stereo_cutangle->GetXaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{in}) [mm]");

	h_pos_diff_eta_in_stereo_cutangle_corr = new TH1F("h_pos_diff_eta_in_stereo_cutangle_corr", "stereo - eta_in (corrected)" ,1000, -5, 5);
	h_pos_diff_eta_in_stereo_cutangle_corr->GetXaxis()->SetTitle("#Delta(#eta^{stereo}-#eta^{in}) [mm]");

	h_d_track_etaout = new TH1F("h_d_track_etaout", "Distance from track - eta_out" ,1000, -10, 10);
	h_d_track_etaout->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etaout_cut = new TH1F("h_d_track_etaout_cut", "Distance from track - eta_out" ,1000, -10, 10);
	h_d_track_etaout_cut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_sby1_anglecut = new TH1F("h_d_track_sby1_anglecut", "Distance from track - SBY1" ,1000, -10, 10);
	h_d_track_sby1_anglecut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_sby2_anglecut = new TH1F("h_d_track_sby2_anglecut", "Distance from track - SBY2" ,1000, -10, 10);
	h_d_track_sby2_anglecut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_sby3_anglecut = new TH1F("h_d_track_sby3_anglecut", "Distance from track - SBY3" ,1000, -10, 10);
	h_d_track_sby3_anglecut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etaout_cut_anglecut = new TH1F("h_d_track_etaout_cut_anglecut", "Distance from track - eta_out" ,1000, -10, 10);
	h_d_track_etaout_cut_anglecut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etaout_cut_anglecut_vs_cl_pos = new TH2F("h_d_track_etaout_cut_anglecut_vs_cl_pos", "Distance from track - eta_out vs cluster position" ,600, 1300, 1900, 1000, -10, 10);
	h_d_track_etaout_cut_anglecut_vs_cl_pos->GetYaxis()->SetTitle("distance [mm]");
	h_d_track_etaout_cut_anglecut_vs_cl_pos->GetXaxis()->SetTitle("cluster position [mm]");

	h_d_track_etain_cut_anglecut_vs_cl_pos = new TH2F("h_d_track_etain_cut_anglecut_vs_cl_pos", "Distance from track - eta_in vs cluster position" ,600, 1300, 1900, 1000, -10, 10);
	h_d_track_etain_cut_anglecut_vs_cl_pos->GetYaxis()->SetTitle("distance [mm]");
	h_d_track_etain_cut_anglecut_vs_cl_pos->GetXaxis()->SetTitle("cluster position [mm]");

	h_d_track_stereo_cut_anglecut_vs_cl_pos = new TH2F("h_d_track_stereo_cut_anglecut_vs_cl_pos", "Distance from track - stereo vs cluster position" ,600, 1300, 1900, 1000, -10, 10);
	h_d_track_stereo_cut_anglecut_vs_cl_pos->GetYaxis()->SetTitle("distance [mm]");
	h_d_track_stereo_cut_anglecut_vs_cl_pos->GetXaxis()->SetTitle("cluster position [mm]");

	h_d_track_etaout_4points = new TH1F("h_d_track_etaout_4points", "Distance from track - eta_out" ,10000, -10, 10);
	h_d_track_etaout_4points->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etaout_cut_4points = new TH1F("h_d_track_etaout_cut_4points", "Distance from track - eta_out" ,10000, -10, 10);
	h_d_track_etaout_cut_4points->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etaout_cut_anglecut_4points = new TH1F("h_d_track_etaout_cut_anglecut_4points", "Distance from track - eta_out" ,1000, -10, 10);
	h_d_track_etaout_cut_anglecut_4points->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etain_4points = new TH1F("h_d_track_etain_4points", "Distance from track - eta_in" ,1000, -10, 10);
	h_d_track_etain_4points->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etain_cut_4points = new TH1F("h_d_track_etain_cut_4points", "Distance from track - eta_in" ,1000, -10, 10);
	h_d_track_etain_cut_4points->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etain_cut_anglecut_4points = new TH1F("h_d_track_etain_cut_anglecut_4points", "Distance from track - eta_in" ,1000, -10, 10);
	h_d_track_etain_cut_anglecut_4points->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_stereo_4points = new TH1F("h_d_track_stereo_4points", "Distance from track - stereo" ,1000, -10, 10);
	h_d_track_stereo_4points->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_stereo_cut_4points = new TH1F("h_d_track_stereo_cut_4points", "Distance from track - stereo" ,1000, -10, 10);
	h_d_track_stereo_cut_4points->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_stereo_cut_anglecut_4points = new TH1F("h_d_track_stereo_cut_anglecut_4points", "Distance from track - stereo" ,1000, -10, 10);
	h_d_track_stereo_cut_anglecut_4points->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etain = new TH1F("h_d_track_etain", "Distance from track - eta_in" ,1000, -10, 10);
	h_d_track_etain->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etain_cut = new TH1F("h_d_track_etain_cut", "Distance from track - eta_in" ,1000, -10, 10);
	h_d_track_etain_cut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_etain_cut_anglecut = new TH1F("h_d_track_etain_cut_anglecut", "Distance from track - eta_in" ,1000, -10, 10);
	h_d_track_etain_cut_anglecut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_stereo = new TH1F("h_d_track_stereo", "Distance from track - stereo" ,1000, -10, 10);
	h_d_track_stereo->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_stereo_cut = new TH1F("h_d_track_stereo_cut", "Distance from track - stereo" ,1000, -10, 10);
	h_d_track_stereo_cut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_stereo_cut_anglecut = new TH1F("h_d_track_stereo_cut_anglecut", "Distance from track - stereo" ,1000, -10, 10);
	h_d_track_stereo_cut_anglecut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay2 = new TH1F("h_d_track_lay2", "Distance from track - stereo_in" ,1000, -10, 10);
	h_d_track_lay2->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay2_cut = new TH1F("h_d_track_lay2_cut", "Distance from track - stereo_in" ,1000, -10, 10);
	h_d_track_lay2_cut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay2_cut_anglecut = new TH1F("h_d_track_lay2_cut_anglecut", "Distance from track - stereo_in" ,1000, -10, 10);
	h_d_track_lay2_cut_anglecut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay3 = new TH1F("h_d_track_lay3", "Distance from track - stereo_out" ,1000, -10, 10);
	h_d_track_lay3->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay3_cut = new TH1F("h_d_track_lay3_cut", "Distance from track - stereo_out" ,1000, -10, 10);
	h_d_track_lay3_cut->GetXaxis()->SetTitle("distance [mm]");

	h_d_track_lay3_cut_anglecut = new TH1F("h_d_track_lay3_cut_anglecut", "Distance from track - stereo_out" ,1000, -10, 10);
	h_d_track_lay3_cut_anglecut->GetXaxis()->SetTitle("distance [mm]");

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

	h_angle_4points_IP1 = new TH1F("h_angle_4points_IP1", "track angle", 5000, -50, 50);
	h_angle_4points_IP1->GetXaxis()->SetTitle("#theta [deg]");

	h_angle_4points_IP2 = new TH1F("h_angle_4points_IP2", "track angle", 5000, -50, 50);
	h_angle_4points_IP2->GetXaxis()->SetTitle("#theta [deg]");

	h_angle_4points_stereo = new TH1F("h_angle_4points_stereo", "track angle", 5000, -50, 50);
	h_angle_4points_stereo->GetXaxis()->SetTitle("#theta [deg]");

	h_angle_cut = new TH1F("h_angle_cut", "track angle", 5000, -50, 50);
	h_angle_cut->GetXaxis()->SetTitle("#theta [deg]");

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

	h_chi2ndf_4points_IP1 = new TH1F("h_chi2ndf_4points_IP1", "chi2 / ndf", 150, 0, 15);
	h_chi2_4points_IP1 = new TH1F("h_chi2_4points_IP1", "chi2", 100, 0, 100);
	h_prob_4points_IP1 = new TH1F("h_prob_4points_IP1", "h_prob", 100, 0, 1);

	h_chi2ndf_4points_IP2 = new TH1F("h_chi2ndf_4points_IP2", "chi2 / ndf", 150, 0, 15);
	h_chi2_4points_IP2 = new TH1F("h_chi2_4points_IP2", "chi2", 100, 0, 100);
	h_prob_4points_IP2 = new TH1F("h_prob_4points_IP2", "h_prob", 100, 0, 1);

	h_chi2ndf_4points_stereo = new TH1F("h_chi2ndf_4points_stereo", "chi2 / ndf", 150, 0, 15);
	h_chi2_4points_stereo = new TH1F("h_chi2_4points_stereo", "chi2", 100, 0, 100);
	h_prob_4points_stereo = new TH1F("h_prob_4points_stereo", "h_prob", 100, 0, 1);

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

	h_diffpos_lay0 = new TH2F("h_diffpos_lay0", "", 100, -2, 2,600,1300,1900); 
	h_diffpos_lay0->GetXaxis()->SetTitle("#Delta y(IP2-IP1) [mm]");
	h_diffpos_lay0->GetYaxis()->SetTitle("y_{IP1} [mm]");

	h_diffpos_lay1 = new TH2F("h_diffpos_lay1", "", 100, -2, 2,600,1300,1900);
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

	h_align_eta_out_ontrack_onxaxis = new TH2F("h_align_eta_out_ontrack_onxaxis", "",200, -.4, .4, 10000, -100, 100);
	h_align_eta_out_ontrack_onxaxis->GetXaxis()->SetTitle("#phi [mm]");
	h_align_eta_out_ontrack_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - IP1) [mm]");

	h_align_eta_out_onxaxis = new TH2F("h_align_eta_out_onxaxis", "",200, -.4, .4, 10000, -100, 100);
	h_align_eta_out_onxaxis->GetXaxis()->SetTitle("#phi [mm]");
	h_align_eta_out_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - IP1) [mm]");

	h_align_eta_out_SBX0_onxaxis = new TH2F("h_align_eta_out_SBX0_onxaxis", "", 600, 1350, 1950, 1000, -100, 100);
	h_align_eta_out_SBX0_onxaxis->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_eta_out_SBX0_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - IP1) [mm]");

	h_align_eta_out_SBX0_onxaxis_corrected = new TH2F("h_align_eta_out_SBX0_onxaxis_corrected", "corrected", 600, 1350, 1950, 1000, -100, 100);
	h_align_eta_out_SBX0_onxaxis_corrected->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_eta_out_SBX0_onxaxis_corrected->GetYaxis()->SetTitle("#Delta y (SBY1 - IP1) [mm]");

	h_alignX_phi_vs_SBX0 = new TH2F("h_alignX_phi_vs_SBX0", "", 600, 1350, 1950, 1000, -.5, .5);
	h_alignX_phi_vs_SBX0->GetXaxis()->SetTitle("SBX0 [mm]");
	h_alignX_phi_vs_SBX0->GetYaxis()->SetTitle("#phi [mm]");

	h_align_eta_out_SBX0_ontrack_onxaxis = new TH2F("h_align_eta_out_SBX0_ontrack_onxaxis", "", 600, 1350, 1950, 1000, -100, 100);
	h_align_eta_out_SBX0_ontrack_onxaxis->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_eta_out_SBX0_ontrack_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - IP1) [mm]");

	h_align_eta_in_ontrack_onxaxis = new TH2F("h_align_eta_in_ontrack_onxaxis", "", 200, -.4, .4, 10000, -100, 100);
	h_align_eta_in_ontrack_onxaxis->GetXaxis()->SetTitle("#phi [mm]");
	h_align_eta_in_ontrack_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - IP2) [mm]");

	h_align_eta_in_onxaxis = new TH2F("h_align_eta_in_onxaxis", "", 200, -.4, .4, 10000, -100, 100);
	h_align_eta_in_onxaxis->GetXaxis()->SetTitle("#phi [mm]");
	h_align_eta_in_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - IP2) [mm]");

	h_align_eta_in_SBX0_onxaxis = new TH2F("h_align_eta_in_SBX0_onxaxis", "", 600, 1350, 1950, 1000, -100, 100);
	h_align_eta_in_SBX0_onxaxis->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_eta_in_SBX0_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - IP2) [mm]");

	h_align_eta_in_SBX0_onxaxis_corrected = new TH2F("h_align_eta_in_SBX0_onxaxis_corrected", "corrected", 600, 1350, 1950, 1000, -100, 100);
	h_align_eta_in_SBX0_onxaxis_corrected->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_eta_in_SBX0_onxaxis_corrected->GetYaxis()->SetTitle("#Delta y (SBY1 - IP2) [mm]");

	h_align_eta_in_SBX0_ontrack_onxaxis = new TH2F("h_align_eta_in_SBX0_ontrack_onxaxis", "", 600, 1350, 1950, 1000, -100, 100);
	h_align_eta_in_SBX0_ontrack_onxaxis->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_eta_in_SBX0_ontrack_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - IP2) [mm]");

	h_align_stereo_ontrack_onxaxis = new TH2F("h_align_stereo_ontrack_onxaxis", "", 200, -.4, .4, 10000, -100, 100);
	h_align_stereo_ontrack_onxaxis->GetXaxis()->SetTitle("#phi [mm]");
	h_align_stereo_ontrack_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - #eta^{stereo}) [mm]");

	h_align_stereo_onxaxis = new TH2F("h_align_stereo_onxaxis", "", 200, -.4, .4, 10000, -100, 100);
	h_align_stereo_onxaxis->GetXaxis()->SetTitle("#phi [mm]");
	h_align_stereo_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - #eta^{stereo}) [mm]");

	h_align_stereo_SBX0_onxaxis = new TH2F("h_align_stereo_SBX0_onxaxis", "", 600, 1350, 1950, 1000, -100, 100);
	h_align_stereo_SBX0_onxaxis->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_stereo_SBX0_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - #eta^{stereo}) [mm]");

	h_align_stereo_SBX0_onxaxis_corrected = new TH2F("h_align_stereo_SBX0_onxaxis_corrected", "corrected", 600, 1350, 1950, 1000, -100, 100);
	h_align_stereo_SBX0_onxaxis_corrected->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_stereo_SBX0_onxaxis_corrected->GetYaxis()->SetTitle("#Delta y (SBY1 - #eta^{stereo}) [mm]");

	h_align_stereo_SBX0_ontrack_onxaxis = new TH2F("h_align_stereo_SBX0_ontrack_onxaxis", "", 600, 1350, 1950, 1000, -100, 100);
	h_align_stereo_SBX0_ontrack_onxaxis->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_stereo_SBX0_ontrack_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - #eta^{stereo}) [mm]");

	h_align_SBY2_ontrack_onxaxis = new TH2F("h_align_SBY2_ontrack_onxaxis", "", 200, -.4, .4, 10000, -100, 100);
	h_align_SBY2_ontrack_onxaxis->GetXaxis()->SetTitle("#phi [mm]");
	h_align_SBY2_ontrack_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - SBY2) [mm]");	

	h_align_SBY2_onxaxis = new TH2F("h_align_SBY2_onxaxis", "", 200, -.4, .4, 10000, -100, 100);
	h_align_SBY2_onxaxis->GetXaxis()->SetTitle("#phi [mm]");
	h_align_SBY2_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - SBY2) [mm]");	

	h_align_SBY2_SBX0_onxaxis = new TH2F("h_align_SBY2_SBX0_onxaxis", "", 600, 1350, 1950, 1000, -100, 100);
	h_align_SBY2_SBX0_onxaxis->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_SBY2_SBX0_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - SBY2) [mm]");	

	h_align_SBY2_SBX0_onxaxis_corrected = new TH2F("h_align_SBY2_SBX0_onxaxis_corrected", "corrected", 600, 1350, 1950, 1000, -100, 100);
	h_align_SBY2_SBX0_onxaxis_corrected->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_SBY2_SBX0_onxaxis_corrected->GetYaxis()->SetTitle("#Delta y (SBY1 - SBY2) [mm]");	

	h_align_SBY2_SBX0_ontrack_onxaxis = new TH2F("h_align_SBY2_SBX0_ontrack_onxaxis", "", 600, 1350, 1950, 1000, -100, 100);
	h_align_SBY2_SBX0_ontrack_onxaxis->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_SBY2_SBX0_ontrack_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - SBY2) [mm]");	

	h_align_SBY3_ontrack_onxaxis = new TH2F("h_align_SBY3_ontrack_onxaxis", "", 200, -.4, .4, 10000, -100, 100);
	h_align_SBY3_ontrack_onxaxis->GetXaxis()->SetTitle("#phi [mm]");
	h_align_SBY3_ontrack_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - SBY3) [mm]");	

	h_align_SBY3_onxaxis = new TH2F("h_align_SBY3_onxaxis", "", 200, -.4, .4, 10000, -100, 100);
	h_align_SBY3_onxaxis->GetXaxis()->SetTitle("#phi [mm]");
	h_align_SBY3_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - SBY3) [mm]");	

	h_align_SBY3_SBX0_onxaxis = new TH2F("h_align_SBY3_SBX0_onxaxis", "", 600, 1350, 1950, 1000, -100, 100);
	h_align_SBY3_SBX0_onxaxis->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_SBY3_SBX0_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - SBY3) [mm]");	

	h_align_SBY3_SBX0_onxaxis_corrected = new TH2F("h_align_SBY3_SBX0_onxaxis_corrected", "corrected", 600, 1350, 1950, 1000, -100, 100);
	h_align_SBY3_SBX0_onxaxis_corrected->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_SBY3_SBX0_onxaxis_corrected->GetYaxis()->SetTitle("#Delta y (SBY1 - SBY3) [mm]");	

	h_align_SBY3_SBX0_ontrack_onxaxis = new TH2F("h_align_SBY3_SBX0_ontrack_onxaxis", "", 600, 1350, 1950, 1000, -100, 100);
	h_align_SBY3_SBX0_ontrack_onxaxis->GetXaxis()->SetTitle("SBX0 [mm]");
	h_align_SBY3_SBX0_ontrack_onxaxis->GetYaxis()->SetTitle("#Delta y (SBY1 - SBY3) [mm]");	

	h_beamProfile = new TH2F("h_beamProfile", "", 600, 1300, 2000, 100, -.4, .4);
	h_beamProfile->GetYaxis()->SetTitle("#phi [mm]");
	h_beamProfile->GetXaxis()->SetTitle("#eta [mm]");

	h_beamProfile_ontrack = new TH2F("h_beamProfile_ontrack", "", 600, 1300, 2000, 100, -.4, .4);
	h_beamProfile_ontrack->GetYaxis()->SetTitle("#phi [mm]");
	h_beamProfile_ontrack->GetXaxis()->SetTitle("#eta [mm]");

	h_sby1_minus_eta_out_vs_pos_eta_out = new TH2F("h_sby1_minus_eta_out_vs_pos_eta_out", "SB1Y-IP1 vs IP1 position", 700, 1300, 2000, 700, -350, 350);
	h_sby1_minus_eta_out_vs_pos_eta_out->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_sby1_minus_eta_out_vs_pos_eta_out->GetYaxis()->SetTitle("SB1Y-IP1 cluster position [mm]");

	h_align_etaout_ontrack = new TH2F("h_align_etaout_ontrack", "SB1Y-IP1 vs IP1 position", 700, 1300, 2000, 700, -350, 350);
	h_align_etaout_ontrack->GetXaxis()->SetTitle("IP1 cluster position [mm]");
	h_align_etaout_ontrack->GetYaxis()->SetTitle("SB1Y-IP1 cluster position [mm]");

	h_sby1_minus_pos_eta_in_vs_pos_eta_in = new TH2F("h_sby1_minus_pos_eta_in_vs_pos_eta_in", "SB1Y-IP2 vs IP2 position", 700, 1300, 2000, 700, -350, 350);
	h_sby1_minus_pos_eta_in_vs_pos_eta_in->GetXaxis()->SetTitle("IP2 cluster position [mm]");
	h_sby1_minus_pos_eta_in_vs_pos_eta_in->GetYaxis()->SetTitle("SB1Y-IP2 cluster position [mm]");

	h_sby1_minus_stereo_vs_stereo = new TH2F("h_sby1_minus_stereo_vs_stereo", "BS1Y-stereo vs SBY1 position", 700, 1300, 2000, 700, -350, 350);
	h_sby1_minus_stereo_vs_stereo->GetXaxis()->SetTitle("SBY1 cluster position [mm]");
	h_sby1_minus_stereo_vs_stereo->GetYaxis()->SetTitle("SB1Y-stereo cluster position [mm]");

	h_sby1_minus_stereo_in_vs_stereo_in = new TH2F("h_sby1_minus_stereo_in_vs_stereo_in", "SB1Y-stereo_in vs SBY1 position", 700, 1300, 2000, 700, -350, 350);
	h_sby1_minus_stereo_in_vs_stereo_in->GetXaxis()->SetTitle("SBY1 cluster position [mm]");
	h_sby1_minus_stereo_in_vs_stereo_in->GetYaxis()->SetTitle("SBY1-stereo_in cluster position [mm]");

	h_sby1_minus_stereo_out_vs_stereo_out = new TH2F("h_sby1_minus_stereo_out_vs_stereo_out", "SBY1-stereo_out vs SBY1 position", 700, 1300, 2000, 700, -350, 350);
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
		
		h_clus_positions_corr[ilayer] = new TH1F(Form("h_clus_position_SM1_lay_corr%i",ilayer), Form("SM1 Layer - %i (corrected)", ilayer) ,700, 1300, 2000);
		h_clus_positions_small_corr[ilayer] = new TH1F(Form("h_clus_position_small_lay_corr%i",ilayer), Form("SB%s  - %i (corrected)",type.c_str(), ilayer) ,700, 1300, 2000);

		h_clus_positions_corr_ontrack[ilayer] = new TH1F(Form("h_clus_position_SM1_lay_corr_ontrack%i",ilayer), Form("SM1 Layer - %i (corrected)", ilayer) ,700, 1300, 2000);
		h_clus_positions_small_corr_ontrack[ilayer] = new TH1F(Form("h_clus_position_small_lay_corr_ontrack%i",ilayer), Form("SB%s  - %i (corrected)",type.c_str(), ilayer) ,700, 1300, 2000);
		h_clus_positions_stripIndex_corr_ontrack[ilayer] = new TH1F(Form("h_clus_positions_stripIndex_corr_ontrack%i",ilayer), Form("SM1 Layer - %i (corrected)", ilayer) ,3072, 0, 3072);

		h_nclusters_per_layer_event[ilayer] = new TH1F(Form("h_nclusters_per_layer_event_%i",ilayer), Form("SM1 Layer - %i", ilayer) ,40, 0, 40);
		
		for(int istrip=0; istrip<8192; istrip++){
			h_strip_tdo[ilayer][istrip] = new TH1F(Form("h_strip_%i_tdo_SM1_lay%i",istrip,ilayer), Form("SM1 strip - %i",istrip),85,0,255);
			h_strip_tdo_small[ilayer][istrip] = new TH1F(Form("h_strip_%i_tdo_small_lay%i",istrip,ilayer), Form("SB strip - %i",istrip),85,0,255);
		}
	}

	
}

#endif

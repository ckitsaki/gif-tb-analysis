float eta_out_cut_eff_down = 1.;
float eta_out_cut_eff_up = 1.;

float eta_in_cut_eff_down = 1.;
float eta_in_cut_eff_up = 1.;

float stereo_in_cut_eff_down = 3;
float stereo_in_cut_eff_up = 2.5;

float stereo_out_cut_eff_down = 2.5;
float stereo_out_cut_eff_up = 3;

float stereo_cut_eff_down = 1.;
float stereo_cut_eff_up = 1.;

float track_angle_cut_up = 0.5;
float track_angle_cut_down = 0.5;

// this cut is inserted to make sure that all chambers see the same beam spectrum.
float beam_acceptance_up = 1740;//1750;
float beam_acceptance_down = 1590;//1400;

// need to perform an extra alignment based on the tracking
struct correction{
	float alpha_corr;
	float beta_corr;
};
	
correction struct_eta_out = {1.19,-0.0009};
correction struct_eta_in = {1.19-0.03,-0.0009};
correction struct_stereo = {1.19-0.03, -0.0009};
//correction struct_eta_in = {1.366,-0.001029};
//correction struct_stereo = {1.67, -0.0012};

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
float beam_acceptance_up = 1750;
float beam_acceptance_down = 1500;//1480;//1400;

// need to perform an extra alignment based on the tracking
struct correction{
	double alpha_corr;
	double beta_corr;
};
	
correction struct_eta_out = {1.253,-0.0009485};
correction struct_eta_in = {1.462,-0.001087};
correction struct_stereo = {1.797,-0.001285};

// fix the correlation due to the track angle
correction struct_eta_in_out_angle = {  0.0004275,0.2937};
correction struct_eta_stereo_in_angle = { 0.0005186,0.4003};
correction struct_eta_stereo_out_angle = { 0.0005391,0.7067};


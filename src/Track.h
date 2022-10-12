#ifndef TRACK_H
#define TRACK_H

#include "Cluster.h"
#include "Histograms.h"
#include "Cuts.h"

class Histograms;
class Cluster;
class Track{

public:
	Track(){};
	Track(Cluster* cl_SBY1, Cluster* cl_SBY2, Cluster* cl_SBY3, Histograms* histos); //constructor
	virtual ~Track()=default; // destructor 
	inline void findSuperClustersWithin1DegWindow();
// return counters
	inline int getSuperClustersWithin1DegWindow() {return m_counter_1deg;};
	inline int getTheOnlySuperClusterOutside1DegWindow() {return m_counter_outside1deg_theonly_scluster;};
	inline int getEventsWithMoreThan1SuperCluster() {return m_counter_more_clusters;};
	inline int getEventsWithMoreThan1SuperClusterIn1DegWindow() {return m_counter_more_clusters_in_1degWindow;};
	inline int getEventsWithMoreThan1SuperClusterOut1DegWindow() {return m_counter_more_clusters_out_1degWindow;};
	inline int getAcceptedTracksTheOnlySuperClusterOutside1DegWindow() {return m_count_singleTracks_outside_1deg;};
	inline int getAcceptedSingleTracksInside1DegWindow() {return m_count_singleTracks_inside_1deg;};
	inline int getAcceptedTracksManySuperClusterInside1DegWindow() {return m_count_many_supercluster;};
// =====================================================================================================================
	inline bool checkTrackCandidate() { return m_isTrackCandidate;}; // for single tracks
	inline void findCandidateTrack();
	inline void isTrackCandidate() { m_isTrackCandidate = ( !(m_map_z_y.empty()) ? true : false );  };
	inline void findTheClosestPair();
	inline int countAcceptedTracks() {return m_accepted_tracks;};
	inline TGraphAsymmErrors* getTrack() {return gr_track;}; // for single tracks
	inline void fitTrack();
	inline float getAngle() {return m_angle;}; // for single tracks
	inline float trackSlope() {return m_slope;}; // for single tracks
	inline float trackIntercept() {return m_intercept;}; // for single tracks
	inline bool rejectTrack() {return m_reject_track;};
	inline void tryAnotherPair();
	inline void fillResiduals() { histograms->h_residual_sby2_sby3_in_1deg_many_sclus->Fill( (m_map_z_y[m_cl_sby2->getZPosition()] - m_map_z_y[m_cl_sby3->getZPosition()]) ); histograms->h_residual_sby2_sby3_accepted->Fill( (m_map_z_y[m_cl_sby2->getZPosition()] - m_map_z_y[m_cl_sby3->getZPosition()]) );};
	inline void fillResidualsZeroProb() { histograms->h_residual_sby2_sby3_zero_prob->Fill( (m_map_z_y[m_cl_sby2->getZPosition()] - m_map_z_y[m_cl_sby3->getZPosition()]) ); histograms->h_residual_sby2_sby3_rejected->Fill( (m_map_z_y[m_cl_sby2->getZPosition()] - m_map_z_y[m_cl_sby3->getZPosition()]) );};
	inline bool tooMany() {return m_many_sclusters;}; // flag which indicates if more than one tracks are found.

	inline int getNMultipleCandTracks() { return m_v_gr_candidate_tracks.size();};
	
	inline std::vector<TGraphAsymmErrors*> getCandidateTracks() {return m_v_gr_candidate_tracks;}; // for multiple tracks in the event
	inline TGraphAsymmErrors* getTrack(int itrack) {return m_v_gr_candidate_tracks.at(itrack);}; // for multiple tracks in the event
	inline float trackSlope(int itrack) {return m_v_slopes_tracks.at(itrack);}; // for multiple tracks in the event
	inline float trackIntercept(int itrack) {return m_v_intercept_tracks.at(itrack);}; // for multiple tracks in the event
	inline float getAngle(int itrack) {return m_v_angle_tracks.at(itrack);}; // for multiple tracks in the event
	inline float getTrackProb(int itrack) {return m_v_prob_tracks.at(itrack);}; // for multiple tracks in the event

	inline void fillSM1events(int clus_index, Cluster* clus_lay, int layIndex, int itrack); // fill the SM1 eta-layers
	inline void fillStereoCluster(Cluster* clus_lay2, Cluster* clus_lay3, int ind_cl2, int ind_cl3, int itrack); // fill with the first coordinate of the stereo information (combination of the two clusters found on both stereo layers)
	inline void setSM1errors(int itrack); // set errors to the SM1 points
	
	inline bool acceptEtaOut() {return m_isEtaOut;}; 
	inline bool acceptEtaIn() {return m_isEtaIn;};
	inline bool acceptStereoIn() {return m_isStereoIn;};
	inline bool acceptStereoOut() {return m_isStereoOut;};
	inline bool acceptStereo() {return m_isStereo;};

	inline void refreshTrack() { m_isEtaIn=false; m_isEtaOut=false; m_isStereoIn=false; m_isStereoIn=false; m_isStereo=false; m_v_distances.clear(); m_v_distances.shrink_to_fit(); }; // for multiple tracks
	
	inline void checkDistanceFromTrack(int itrack); // not used
	inline void checkDistanceFromTrack(int itrack, float threshold); // not used
	inline void checkDistanceFromTrack(int itrack, int layer, Cluster* cluster, int index, float threshold_down, float threshold_up); 
	inline void checkStereoDistanceFromTrack(int itrack, float pos_stereo_y, float pos_stereo_z, float stereo_cut_eff_down, float stereo_cut_eff_up);

	inline void etaOutFlag(bool isEtaOut) {m_isEtaOut = isEtaOut;};
	inline void etaInFlag(bool isEtaIn) {m_isEtaIn = isEtaIn;};

	inline void setCriterium(); // this criterium is used to find the best track among the many for the case of multiple track events
	inline std::vector<float> getSelectionCriteria() {return m_v_track_sel_criterium;};

	inline float Get_Distance_from_Track(float z, float y, int itrack); // return the distance of the cluster found on SM1 layer from the track
	inline void isSingleTrack() {m_mult_track_singles = true;};
	inline bool checkIfSingleTrack() {return m_mult_track_singles;};

	// Methods to calculate the efficiency using the inclusive-exclusive method
	inline TGraphAsymmErrors* track4points(int itrack, int ilayer); 
	inline bool etaOut_fired(int itrack) {return m_etaOut_fired.at(itrack);};
	inline bool etaIn_fired(int itrack) {return m_etaIn_fired.at(itrack);};
	inline void etaOut_fired() {m_etaOut_fired.push_back(false);};
	inline void etaIn_fired() {m_etaIn_fired.push_back(false);};

	// Functions to define the efficiency per MMFE8
	inline float extrapolateTrackOnSM1(int itrack);

	TGraphAsymmErrors* gr_track;
	TGraphAsymmErrors* gr_track_4points;
	TGraphAsymmErrors* gr_track_5points;
	Histograms* histograms;

private:
	Cluster* m_cl_sby1;
	Cluster* m_cl_sby2;
	Cluster* m_cl_sby3;
	
	int m_counter_1deg=0;
	int m_counter_outside1deg_theonly_scluster=0;
	int m_counter_more_clusters=0;
	int m_counter_more_clusters_in_1degWindow=0;
	int m_counter_more_clusters_out_1degWindow=0;
	int m_accepted_tracks=0;
	int m_count_singleTracks_outside_1deg=0;
	int m_count_singleTracks_inside_1deg=0;
	int m_count_many_supercluster=0;
	
	bool m_many_sclusters=false;
	bool m_isTrackCandidate=false;
	bool m_reject_track = false;

	bool m_isStereo = false;
	bool m_isEtaIn = false;
	bool m_isEtaOut = false;
	bool m_isStereoIn = false;
	bool m_isStereoOut = false;
	bool m_mult_track_singles = false;

// flags for inclusive-exclusive method
	std::vector<bool> m_etaOut_fired; 
	std::vector<bool> m_etaIn_fired; 
// =====================================

	bool m_accept_eta_out = false;
	bool m_accept_eta_in = false;

	bool m_the_only_supercluster_outside_1deg_window = false;

	float m_angle;
	float m_slope;
	float m_intercept;
	float m_totCharge;
	float m_prob;

	int m_counter;
	int m_ipoint=2;

	std::vector<float> m_v_pos_abs_diff;
	std::vector<float> m_v_pos_diff_ref;
	
	std::vector<TGraphAsymmErrors*> m_v_gr_candidate_tracks;
	std::vector<float> m_v_slopes_tracks;
	std::vector<float> m_v_intercept_tracks;
	std::vector<float> m_v_angle_tracks;
	std::vector<float> m_v_prob_tracks;
	std::vector<float> m_v_track_sel_criterium;
	std::vector<float> m_v_distances;

	std::map<float,float> m_map_z_y;
	std::map<float,float> m_map_z_errory;
	std::multimap<float,float> m_pos_sby2_sby3;
	std::map<float,float> m_map_point_charge;
	
};

inline Track::Track(Cluster* cl_SBY1, Cluster* cl_SBY2, Cluster* cl_SBY3, Histograms* histos)
{
	histograms = histos;
	m_cl_sby1 = cl_SBY1;
	m_cl_sby2 = cl_SBY2;
	m_cl_sby3 = cl_SBY3;
	m_map_z_y = map<float,float>();
	m_map_z_errory = map<float,float>();
	m_pos_sby2_sby3 = multimap<float,float>();
	m_map_point_charge = map<float,float>();
	findSuperClustersWithin1DegWindow();
	isTrackCandidate(); // check if the m_map_z_y is filled or not
	if(m_isTrackCandidate && !m_many_sclusters)
	{
		findCandidateTrack();
		fitTrack();
	}
	if(m_many_sclusters) // if it s a track candidate we start with the closest pair
	{
	//	std::cout<<m_counter_more_clusters_in_1degWindow<<" tracks found \n";
		for(int icand=0; icand<m_pos_sby2_sby3.size(); icand++)
		{
			if(icand!=0) tryAnotherPair();
			
			if(icand==0) {
				findTheClosestPair();
				findCandidateTrack();
				fitTrack();
			}

			m_map_z_y.clear();
			m_map_z_errory.clear();
		}
	}

//	m_v_pos_diff_ref.clear();
//	m_v_pos_abs_diff.clear();
//	m_v_pos_diff_ref.shrink_to_fit();
//	m_v_pos_abs_diff.shrink_to_fit();
//	m_v_gr_rejected_tracks.clear();
//	m_v_gr_rejected_tracks.shrink_to_fit();
//	m_v_gr_accepted_tracks.clear();
//	m_v_gr_accepted_tracks.shrink_to_fit();

	
}

inline void Track::fillSM1events(int clus_index, Cluster* clus_lay, int layIndex, int itrack)
{

	float alpha_corr = 0;
	float beta_corr = 1;
	float beta=0;
	
	if(!m_many_sclusters)
	{
	
		beta = (1/TMath::Sqrt(m_slope*m_slope + 1));
		m_ipoint++;
	
		if(layIndex==0)
		{
			alpha_corr = struct_eta_out.alpha_corr;
			beta_corr = struct_eta_out.beta_corr + beta;
		}
		if(layIndex==1)
		{
			alpha_corr = struct_eta_in.alpha_corr;
			beta_corr = struct_eta_in.beta_corr + beta;
		}
	
		gr_track->SetPointX(m_ipoint, clus_lay->getZPosition());
		gr_track->SetPointY(m_ipoint, alpha_corr + beta_corr*clus_lay->getCorrPosition(clus_index));
	}
	else 
	{
		
		beta = 1/TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1);
		if(getTrack(itrack)->GetN()<4) m_ipoint=2;
		m_ipoint++;
		
		if(itrack>0 && getTrack(itrack)->GetN()<4){
		 m_map_point_charge.clear();
		// std::cout<<m_ipoint<<std::endl;
		}
		if(layIndex==0)
		{
			alpha_corr = struct_eta_out.alpha_corr;
			beta_corr = struct_eta_out.beta_corr + beta;
		}
		if(layIndex==1)
		{
			alpha_corr = struct_eta_in.alpha_corr;
			beta_corr = struct_eta_in.beta_corr + beta;
		}
	
		getTrack(itrack)->SetPointX(m_ipoint, clus_lay->getZPosition());
		getTrack(itrack)->SetPointY(m_ipoint, alpha_corr + beta_corr*clus_lay->getCorrPosition(clus_index));	
	}
	m_map_point_charge[m_ipoint] = clus_lay->getTotPdo(clus_index);
	m_totCharge += clus_lay->getTotPdo(clus_index);
//	std::cout<<"POINTS: "<<gr_track->GetPointX(m_ipoint)<<" "<<gr_track->GetPointY(m_ipoint)<<std::endl;
	if(layIndex==0)
	{
		checkDistanceFromTrack(itrack, layIndex, clus_lay, clus_index, eta_out_cut_eff_down, eta_out_cut_eff_up);
		m_etaOut_fired.push_back(true);
	} 
	if(layIndex==1)
	{
		checkDistanceFromTrack(itrack, layIndex, clus_lay, clus_index, eta_in_cut_eff_down, eta_in_cut_eff_up);
		m_etaIn_fired.push_back(true);
	} 
}

inline float Track::extrapolateTrackOnSM1(int itrack)
{
	float expected_position_on_SM1 = 9999.;

	if(!m_many_sclusters) expected_position_on_SM1 = m_slope*(4*2415+(13.835+30.615+46.985+63.765))/4 + m_intercept;
	else expected_position_on_SM1 = trackSlope(itrack)*(4*2415+(13.835+30.615+46.985+63.765))/4 + trackIntercept(itrack);

	return expected_position_on_SM1;
}

inline TGraphAsymmErrors* Track::track4points(int itrack, int ilayer)
{
	gr_track_4points = new TGraphAsymmErrors();
	gr_track_4points->GetXaxis()->SetTitle("z [mm]");
	gr_track_4points->GetYaxis()->SetTitle("y [mm]");
	bool check = false;
	if(ilayer>0) check=true; 
	
	if(!m_many_sclusters)
	{
		gr_track_4points->SetPointX(0, gr_track->GetPointX(0));
		gr_track_4points->SetPointX(1, gr_track->GetPointX(1));
		gr_track_4points->SetPointX(2, gr_track->GetPointX(2));

		gr_track_4points->SetPointY(0, gr_track->GetPointY(0));
		gr_track_4points->SetPointY(1, gr_track->GetPointY(1));
		gr_track_4points->SetPointY(2, gr_track->GetPointY(2));

		gr_track_4points->SetPointError(0, gr_track->GetErrorXlow(0), gr_track->GetErrorXhigh(0), gr_track->GetErrorYlow(0), gr_track->GetErrorYhigh(0));
		gr_track_4points->SetPointError(1, gr_track->GetErrorXlow(1), gr_track->GetErrorXhigh(1), gr_track->GetErrorYlow(1), gr_track->GetErrorYhigh(1));
		gr_track_4points->SetPointError(2, gr_track->GetErrorXlow(2), gr_track->GetErrorXhigh(2), gr_track->GetErrorYlow(2), gr_track->GetErrorYhigh(2));

		if(gr_track->GetN()>3 )
		{
			if(check && m_etaOut_fired.at(itrack))
			{
				gr_track_4points->SetPointX(3, gr_track->GetPointX(4));
				gr_track_4points->SetPointY(3, gr_track->GetPointY(4));
				gr_track_4points->SetPointError(3, gr_track->GetErrorXlow(4), gr_track->GetErrorXhigh(4), gr_track->GetErrorYlow(4), gr_track->GetErrorYhigh(4));
			}
			else
			{
				gr_track_4points->SetPointX(3, gr_track->GetPointX(3));
				gr_track_4points->SetPointY(3, gr_track->GetPointY(3));
				gr_track_4points->SetPointError(3, gr_track->GetErrorXlow(3), gr_track->GetErrorXhigh(3), gr_track->GetErrorYlow(3), gr_track->GetErrorYhigh(3));
			}
		}
	}
	else
	{
		gr_track_4points->SetPointX(0, getTrack(itrack)->GetPointX(0));
		gr_track_4points->SetPointX(1, getTrack(itrack)->GetPointX(1));
		gr_track_4points->SetPointX(2, getTrack(itrack)->GetPointX(2));
		
		gr_track_4points->SetPointError(0, getTrack(itrack)->GetErrorXlow(0), getTrack(itrack)->GetErrorXhigh(0), getTrack(itrack)->GetErrorYlow(0), getTrack(itrack)->GetErrorYhigh(0));
		gr_track_4points->SetPointError(1, getTrack(itrack)->GetErrorXlow(1), getTrack(itrack)->GetErrorXhigh(1), getTrack(itrack)->GetErrorYlow(1), getTrack(itrack)->GetErrorYhigh(1));
		gr_track_4points->SetPointError(2, getTrack(itrack)->GetErrorXlow(2), getTrack(itrack)->GetErrorXhigh(2), getTrack(itrack)->GetErrorYlow(2), getTrack(itrack)->GetErrorYhigh(2));

		gr_track_4points->SetPointY(0, getTrack(itrack)->GetPointY(0));
		gr_track_4points->SetPointY(1, getTrack(itrack)->GetPointY(1));
		gr_track_4points->SetPointY(2, getTrack(itrack)->GetPointY(2));

		if(getTrack(itrack)->GetN()>3)
		{
			if(check && m_etaOut_fired.at(itrack))
			{
				gr_track_4points->SetPointX(3, getTrack(itrack)->GetPointX(4));
				gr_track_4points->SetPointY(3, getTrack(itrack)->GetPointY(4));
				gr_track_4points->SetPointError(3, getTrack(itrack)->GetErrorXlow(4), getTrack(itrack)->GetErrorXhigh(4), getTrack(itrack)->GetErrorYlow(4), getTrack(itrack)->GetErrorYhigh(4));
			}
			else
			{
				gr_track_4points->SetPointX(3, getTrack(itrack)->GetPointX(3));
				gr_track_4points->SetPointY(3, getTrack(itrack)->GetPointY(3));
				gr_track_4points->SetPointError(3, getTrack(itrack)->GetErrorXlow(3), getTrack(itrack)->GetErrorXhigh(3), getTrack(itrack)->GetErrorYlow(3), getTrack(itrack)->GetErrorYhigh(3));
			}
			
		}
	}

	TF1* f_track_4points = new TF1("f_track_4points", "[0] + [1]*x", 0, 3800);
	gr_track_4points->Fit("f_track_4points","Q");

	float constant = f_track_4points->GetParameter(0);
	float slope = f_track_4points->GetParameter(1);
	float ndf = f_track_4points->GetNDF();
	float chi2 = f_track_4points->GetChisquare();
	float prob = f_track_4points->GetProb();
	float angle = (180./TMath::Pi())*TMath::ATan(slope);

	if(prob<0.05)
		gr_track_4points->Fit("f_track","Q");

	if(ilayer==0)
	{
		histograms->h_chi2ndf_4points_IP1->Fill(chi2/ndf);
		histograms->h_chi2_4points_IP1->Fill(chi2);
		histograms->h_angle_4points_IP1->Fill(m_angle);
		histograms->h_prob_4points_IP1->Fill(prob);
	}
	if(ilayer==1)
	{
		histograms->h_chi2ndf_4points_IP2->Fill(chi2/ndf);
		histograms->h_chi2_4points_IP2->Fill(chi2);
		histograms->h_angle_4points_IP2->Fill(m_angle);
		histograms->h_prob_4points_IP2->Fill(prob);
	}
	

	if(m_etaOut_fired.at(itrack) && !check)
	{
		float distance_from_track = (slope*gr_track_4points->GetPointX(3) - gr_track_4points->GetPointY(3) + constant)/(TMath::Sqrt(slope*slope + 1));
		histograms->h_d_track_etaout_4points->Fill(distance_from_track);
		if(distance_from_track<=eta_out_cut_eff_up && distance_from_track>=-eta_out_cut_eff_down)
		{
			histograms->h_d_track_etaout_cut_4points->Fill(distance_from_track);
			if(angle>=-track_angle_cut_down && angle<=track_angle_cut_up) histograms->h_d_track_etaout_cut_anglecut_4points->Fill(distance_from_track);
		}
	}

	if(m_etaIn_fired.at(itrack) && check)
	{
		float distance_from_track = 99999.;
		if(!m_etaOut_fired.at(itrack)) distance_from_track = (slope*gr_track_4points->GetPointX(3) - gr_track_4points->GetPointY(3) + constant)/(TMath::Sqrt(slope*slope + 1));
		else distance_from_track = (slope*gr_track_4points->GetPointX(4) - gr_track_4points->GetPointY(4) + constant)/(TMath::Sqrt(slope*slope + 1));
		
		histograms->h_d_track_etain_4points->Fill(distance_from_track);
		if(distance_from_track<=eta_out_cut_eff_up && distance_from_track>=-eta_out_cut_eff_down)
		{
			histograms->h_d_track_etain_cut_4points->Fill(distance_from_track);
			if(angle>=-track_angle_cut_down && angle<=track_angle_cut_up) histograms->h_d_track_etain_cut_anglecut_4points->Fill(distance_from_track);
		}
	}

	return gr_track_4points;
} 

inline void Track::setSM1errors(int itrack)
{
	float totalCharge = m_totCharge;
	for(auto& ipoint : m_map_point_charge)
	{
		if(!m_many_sclusters)
			gr_track->SetPointError(ipoint.first, 10, 10, /*ipoint.second/(2*totalCharge)*//*ipoint.second/(2*totalCharge)*/0.4, 0.4/*ipoint.second/(2*totalCharge)*//*ipoint.second/(2*totalCharge)*/);
		else 
			getTrack(itrack)->SetPointError(ipoint.first, 10, 10, /*ipoint.second/(2*totalCharge)*//*ipoint.second/(2*totalCharge)*/0.4, 0.4/*ipoint.second/(2*totalCharge)*//*ipoint.second/(2*totalCharge)*/);
	}
}

inline float Track::Get_Distance_from_Track(float z, float y, int itrack)
{
	if(!m_many_sclusters)
		return (m_slope*z - y + m_intercept)/(TMath::Sqrt(m_slope*m_slope + 1));
	else
		return (trackSlope(itrack)*z - y + trackIntercept(itrack))/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));

	return 9999.;
}

inline void Track::fillStereoCluster(Cluster* clus_lay2, Cluster* clus_lay3, int ind_cl2, int ind_cl3, int itrack)
{
	float alpha_corr = 0;
	float beta_corr = 1;
	float beta=0;

	m_ipoint++;
	if(!m_many_sclusters)
	{
		beta = (1/TMath::Sqrt(m_slope*m_slope + 1));
		alpha_corr = struct_stereo.alpha_corr;
		beta_corr = struct_stereo.beta_corr + beta;
		gr_track->SetPointX(m_ipoint, (clus_lay2->getZPosition() + clus_lay3->getZPosition() )/2 );
		float stereo_pos = (+0.56 + ( (clus_lay2->getCorrPosition(ind_cl2)+clus_lay3->getCorrPosition(ind_cl3)) / 2*TMath::Cos(1.5*TMath::Pi()/180.)));
		gr_track->SetPointY(m_ipoint, alpha_corr + beta_corr*stereo_pos);
	}
	else {
		beta = 1/TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1);
		alpha_corr = struct_stereo.alpha_corr;
		beta_corr = struct_stereo.beta_corr + beta;
		getTrack(itrack)->SetPointX(m_ipoint, (clus_lay2->getZPosition() + clus_lay3->getZPosition() )/2 );
		float stereo_pos = (+0.56 + ( (clus_lay2->getCorrPosition(ind_cl2)+clus_lay3->getCorrPosition(ind_cl3)) / 2*TMath::Cos(1.5*TMath::Pi()/180.)));
		getTrack(itrack)->SetPointY(m_ipoint, alpha_corr + beta_corr*stereo_pos);
	}
	m_map_point_charge[m_ipoint] = (clus_lay2->getTotPdo(ind_cl2)+clus_lay3->getTotPdo(ind_cl3));
	m_totCharge += clus_lay2->getTotPdo(ind_cl2);
	m_totCharge += clus_lay3->getTotPdo(ind_cl3);
}

inline void Track::setCriterium()
{
	float criterium = 0;
	if(m_v_distances.empty()) return;
	for(int i=0; i<m_v_distances.size(); i++)
	{
		criterium += std::abs(m_v_distances.at(i));
	}
	criterium = TMath::Sqrt(criterium)/m_v_distances.size();
	//std::cout<<criterium<<" "<<m_v_distances.size()<<std::endl;
	m_v_track_sel_criterium.push_back(criterium);
}

inline void Track::checkDistanceFromTrack(int itrack, int layer, Cluster* cluster, int index, float threshold_down, float threshold_up)
{
	float distance_from_track = 9999.;
	float pos_sby1 = m_cl_sby1->getPosition(0);

	float alpha_corr = 0;
	float beta_corr = 1;
	float beta=0;
	if (!m_many_sclusters){
	 beta = (1/TMath::Sqrt(m_slope*m_slope + 1));
	}
	else{
		beta = 1/TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1);
	}
	
	if(layer==0)
	{
		alpha_corr = struct_eta_out.alpha_corr;
		beta_corr = struct_eta_out.beta_corr + beta;
	}
	if(layer==1)
	{
		alpha_corr = struct_eta_in.alpha_corr;
		beta_corr = struct_eta_in.beta_corr + beta;
	}
	
	if(!m_many_sclusters) distance_from_track = (m_slope*cluster->getZPosition() - (alpha_corr + beta_corr*cluster->getCorrPosition(index)) + m_intercept)/(TMath::Sqrt(m_slope*m_slope + 1));
	else distance_from_track = (trackSlope(itrack)*cluster->getZPosition() - (alpha_corr + beta_corr*cluster->getCorrPosition(index)) + trackIntercept(itrack))/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));

	if(m_many_sclusters && !m_mult_track_singles && layer<2) m_v_distances.push_back(distance_from_track);
	
	if(layer==0)
	{
		if(!m_many_sclusters || (m_many_sclusters && m_mult_track_singles)) histograms->h_d_track_etaout->Fill(distance_from_track);	
		if(distance_from_track<=threshold_up && distance_from_track>=-threshold_down && (cluster->getCorrPosition(index)>=beam_acceptance_down && cluster->getCorrPosition(index)<=beam_acceptance_up) ) {
		 m_isEtaOut=true;
		 if(!m_many_sclusters && (m_angle>=-track_angle_cut_down && m_angle<=track_angle_cut_up)){
		  histograms->h_d_track_etaout_cut_anglecut->Fill(distance_from_track);
		  histograms->h_d_track_etaout_cut_anglecut_vs_cl_pos->Fill((alpha_corr + beta_corr*cluster->getCorrPosition(index)), distance_from_track);
		  histograms->h_align_etaout_ontrack->Fill((alpha_corr + beta_corr*cluster->getCorrPosition(index)), pos_sby1-(alpha_corr + beta_corr*cluster->getCorrPosition(index)));
		}
		 
		 if((m_many_sclusters && m_mult_track_singles) && (getAngle(0)>=-track_angle_cut_down && getAngle(0)<=track_angle_cut_up)){
		  histograms->h_d_track_etaout_cut_anglecut->Fill(distance_from_track);
		  histograms->h_d_track_etaout_cut_anglecut_vs_cl_pos->Fill((alpha_corr + beta_corr*cluster->getCorrPosition(index)), distance_from_track);
		  histograms->h_align_etaout_ontrack->Fill((alpha_corr + beta_corr*cluster->getCorrPosition(index)), pos_sby1-(alpha_corr + beta_corr*cluster->getCorrPosition(index)));
		}
		 
		 if(!m_many_sclusters || (m_many_sclusters && m_mult_track_singles)) {
		  histograms->h_clus_positions_corr_ontrack[layer]->Fill((alpha_corr + beta_corr*cluster->getCorrPosition(index)));
		  histograms->h_clus_positions_stripIndex_corr_ontrack[layer]->Fill(cluster->getStrip_from_clpos((alpha_corr + beta_corr*cluster->getCorrPosition(index))));
		  histograms->h_d_track_etaout_cut->Fill(distance_from_track);
		  histograms->h_cl_charge_on_track[layer]->Fill(cluster->getTotPdo(index));
		  histograms->h_nstrips_on_track[layer]->Fill(cluster->getNStrips(index));
		  }
		}
		else m_isEtaOut = false;
	}  

	if(layer==1)
	{
		if(!m_many_sclusters || (m_many_sclusters && m_mult_track_singles)) histograms->h_d_track_etain->Fill(distance_from_track);
		if(distance_from_track<=threshold_up && distance_from_track>=-threshold_down && (cluster->getCorrPosition(index)>=beam_acceptance_down && cluster->getCorrPosition(index)<=beam_acceptance_up)) {
		 m_isEtaIn=true;
		 if(!m_many_sclusters && (m_angle>=-track_angle_cut_down && m_angle<=track_angle_cut_up))
		 {
		 	 histograms->h_d_track_etain_cut_anglecut->Fill(distance_from_track);
		 	 histograms->h_d_track_etain_cut_anglecut_vs_cl_pos->Fill(alpha_corr + beta_corr*cluster->getCorrPosition(index), distance_from_track);
		 }
		 if((m_many_sclusters && m_mult_track_singles) && (getAngle(0)>=-track_angle_cut_down && getAngle(0)<=track_angle_cut_up)) 
		 {
		 	histograms->h_d_track_etain_cut_anglecut->Fill(distance_from_track);
		 	histograms->h_d_track_etain_cut_anglecut_vs_cl_pos->Fill(alpha_corr + beta_corr*cluster->getCorrPosition(index), distance_from_track);
		 }
		  
		 if(!m_many_sclusters || (m_many_sclusters && m_mult_track_singles)) {
		 	histograms->h_clus_positions_corr_ontrack[layer]->Fill(alpha_corr + beta_corr*cluster->getCorrPosition(index));
		 	histograms->h_clus_positions_stripIndex_corr_ontrack[layer]->Fill(cluster->getStrip_from_clpos(cluster->getCorrPosition(index)));
		 	histograms->h_d_track_etain_cut->Fill(distance_from_track);
		 	histograms->h_cl_charge_on_track[layer]->Fill(cluster->getTotPdo(index));
		 	histograms->h_nstrips_on_track[layer]->Fill(cluster->getNStrips(index));
		 } 
		}
		else m_isEtaIn = false;
	}

	if(layer==2)
	{
		if(!m_many_sclusters || (m_many_sclusters && m_mult_track_singles)) histograms->h_d_track_lay2->Fill(distance_from_track);
		if(distance_from_track<=threshold_up && distance_from_track>=-threshold_down && (cluster->getCorrPosition(index)>=beam_acceptance_down && cluster->getCorrPosition(index)<=beam_acceptance_up)) {
		 m_isStereoIn=true;
		 if(!m_many_sclusters && (m_angle>=-track_angle_cut_down && m_angle<=track_angle_cut_up)) histograms->h_d_track_lay2_cut_anglecut->Fill(distance_from_track);
		 
		 if((m_many_sclusters && m_mult_track_singles) && (getAngle(0)>=-track_angle_cut_down && getAngle(0)<=track_angle_cut_up)) histograms->h_d_track_lay2_cut_anglecut->Fill(distance_from_track);
		 
		 if(!m_many_sclusters || (m_many_sclusters && m_mult_track_singles))
		 {
		 	histograms->h_clus_positions_corr_ontrack[layer]->Fill(cluster->getCorrPosition(index));
		 	histograms->h_clus_positions_stripIndex_corr_ontrack[layer]->Fill(cluster->getStrip_from_clpos(cluster->getCorrPosition(index)));
		 	histograms->h_d_track_lay2_cut->Fill(distance_from_track);
		 	histograms->h_cl_charge_on_track[layer]->Fill(cluster->getTotPdo(index));
		 	histograms->h_nstrips_on_track[layer]->Fill(cluster->getNStrips(index));
		 } 
		}
		else m_isStereoIn = false;
	}

	if(layer==3)
	{
		if(!m_many_sclusters || (m_many_sclusters && m_mult_track_singles)) histograms->h_d_track_lay3->Fill(distance_from_track);
		if(distance_from_track<=threshold_up && distance_from_track>=-threshold_down && (cluster->getCorrPosition(index)>=beam_acceptance_down && cluster->getCorrPosition(index)<=beam_acceptance_up) ) {
		 m_isStereoOut=true;
		 
		 if(!m_many_sclusters && (m_angle>=-track_angle_cut_down && m_angle<=track_angle_cut_up)) histograms->h_d_track_lay3_cut_anglecut->Fill(distance_from_track);
		 
		 if((m_many_sclusters && m_mult_track_singles) && (getAngle(0)>=-track_angle_cut_down && getAngle(0)<=track_angle_cut_up)) histograms->h_d_track_lay3_cut_anglecut->Fill(distance_from_track);
		 
		 if(!m_many_sclusters || (m_many_sclusters && m_mult_track_singles)) 
		 {
		 	histograms->h_clus_positions_corr_ontrack[layer]->Fill(cluster->getCorrPosition(index));
		 	histograms->h_clus_positions_stripIndex_corr_ontrack[layer]->Fill(cluster->getStrip_from_clpos(cluster->getCorrPosition(index)));
		 	histograms->h_d_track_lay3_cut->Fill(distance_from_track);
		 	histograms->h_cl_charge_on_track[layer]->Fill(cluster->getTotPdo(index));
		 	histograms->h_nstrips_on_track[layer]->Fill(cluster->getNStrips(index));
		 } 
		}
		else m_isStereoOut = false;
	}	
}

inline void Track::checkStereoDistanceFromTrack(int itrack, float pos_stereo_y, float pos_stereo_z, float stereo_cut_eff_down, float stereo_cut_eff_up)
{
	float alpha_corr = 0;
	float beta_corr = 1;
	float beta=0;
	if(!m_many_sclusters) beta = (1/TMath::Sqrt(m_slope*m_slope + 1));
	else beta = 1/TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1);
	
	alpha_corr = struct_stereo.alpha_corr;
	beta_corr = struct_stereo.beta_corr + beta;

	float distance_from_track = 9999.;
	if(!m_many_sclusters ) distance_from_track = (m_slope*pos_stereo_z - (alpha_corr+beta_corr*pos_stereo_y) + m_intercept)/(TMath::Sqrt(m_slope*m_slope + 1));
	else distance_from_track = (trackSlope(itrack)*pos_stereo_z - (alpha_corr+beta_corr*pos_stereo_y) + trackIntercept(itrack))/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));

	if(m_many_sclusters && !m_mult_track_singles) m_v_distances.push_back(distance_from_track);

	if(!m_many_sclusters ) histograms->h_d_track_stereo->Fill(distance_from_track);
	
	if(distance_from_track<=stereo_cut_eff_up && distance_from_track>=-stereo_cut_eff_down && pos_stereo_y>=beam_acceptance_down && pos_stereo_y<=beam_acceptance_up)
	{
		m_isStereo = true;
		if(!m_many_sclusters) {
		 histograms->h_d_track_stereo_cut->Fill(distance_from_track);
		 histograms->h_d_track_stereo_cut_anglecut_vs_cl_pos->Fill(alpha_corr+beta_corr*pos_stereo_y, distance_from_track);
		}
		if(!m_many_sclusters && (m_angle>=-track_angle_cut_down && m_angle<=track_angle_cut_up)) 
		{
			histograms->h_d_track_stereo_cut_anglecut->Fill(distance_from_track);
			histograms->h_d_track_stereo_cut_anglecut_vs_cl_pos->Fill(alpha_corr+beta_corr*pos_stereo_y, distance_from_track);
		}
	}
	else m_isStereo = false;

}


inline void Track::checkDistanceFromTrack(int itrack, float threshold)
{
	float distance_from_track=9999.;
	if(m_isEtaOut)
	{
		distance_from_track = (m_slope*gr_track->GetPointX(3) - gr_track->GetPointY(3) + m_intercept)/(TMath::Sqrt(m_slope*m_slope + 1));
		if(m_many_sclusters)
		{
			distance_from_track = (trackSlope(itrack)*getTrack(itrack)->GetPointX(3) - getTrack(itrack)->GetPointY(3) + trackIntercept(itrack))/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));
			m_v_distances.push_back(distance_from_track);
		}	
		//histograms->h_d_track_etaout->Fill(distance_from_track);
		if( distance_from_track <= threshold && distance_from_track >= -threshold)
		{
			m_isEtaOut = true;
			//if(threshold==2.)
			histograms->h_d_track_etaout->Fill(distance_from_track);
			histograms->h_real_vs_exp_position_etaout->Fill(gr_track->GetPointX(3), m_slope*gr_track->GetPointX(3)+m_intercept);
		}
		else m_isEtaOut = false;
	}
	if(m_isEtaIn)
	{
		if(m_isEtaOut)
		{
			if(!m_many_sclusters)
				distance_from_track = (m_slope*gr_track->GetPointX(4) - gr_track->GetPointY(4) + m_intercept)/(TMath::Sqrt(m_slope*m_slope + 1));
			else {
				distance_from_track = (trackSlope(itrack)*getTrack(itrack)->GetPointX(4) - getTrack(itrack)->GetPointY(4) + trackIntercept(itrack))/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));
				m_v_distances.push_back(distance_from_track);
			}

		}
		else 
		{
			if(!m_many_sclusters)
				distance_from_track = (m_slope*gr_track->GetPointX(3) - gr_track->GetPointY(3) + m_intercept)/(TMath::Sqrt(m_slope*m_slope + 1));
			else {
				distance_from_track = (trackSlope(itrack)*getTrack(itrack)->GetPointX(3) - getTrack(itrack)->GetPointY(3) + trackIntercept(itrack))/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));
				m_v_distances.push_back(distance_from_track);
			}

		}
		
		if(distance_from_track <= threshold && distance_from_track>=-threshold)
		{
			m_isEtaIn = true;
			//if(threshold==2.)
			histograms->h_d_track_etain->Fill(distance_from_track);
		}
		else m_isEtaIn = false;
	}
	
	if(m_isStereo)
	{
		if(m_isEtaOut && m_isEtaIn)
		{
			if(!m_many_sclusters)
				distance_from_track = (m_slope*gr_track->GetPointX(5) - gr_track->GetPointY(5) + m_intercept)/(TMath::Sqrt(m_slope*m_slope + 1));
			else {
				distance_from_track = (trackSlope(itrack)*getTrack(itrack)->GetPointX(5) - getTrack(itrack)->GetPointY(5) + trackIntercept(itrack))/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));
				m_v_distances.push_back(distance_from_track);
			}
		}
		if( (m_isEtaOut && !m_isEtaIn) || (!m_isEtaOut && m_isEtaIn) )
		{
			if(!m_many_sclusters)
				distance_from_track = (m_slope*gr_track->GetPointX(4) - gr_track->GetPointY(4) + m_intercept)/(TMath::Sqrt(m_slope*m_slope + 1));
			else {
				distance_from_track = (trackSlope(itrack)*getTrack(itrack)->GetPointX(4) - getTrack(itrack)->GetPointY(4) + trackIntercept(itrack))/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));
				m_v_distances.push_back(distance_from_track);
			}
		}
		if( !m_isEtaOut && !m_isEtaIn)
		{
			if(!m_many_sclusters)
				distance_from_track = (m_slope*gr_track->GetPointX(3) - gr_track->GetPointY(3) + m_intercept)/(TMath::Sqrt(m_slope*m_slope + 1));
			else
			{
				distance_from_track = (trackSlope(itrack)*getTrack(itrack)->GetPointX(3) - getTrack(itrack)->GetPointY(3) + trackIntercept(itrack))/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));
				m_v_distances.push_back(distance_from_track);
			}
		}

		if(distance_from_track <= threshold && distance_from_track>=-threshold)
		{
			m_isStereo = true;
			//histograms->h_d_track_stereoin->Fill(distance_from_track);
		}
		else m_isStereo = false;
	}

}

inline void Track::checkDistanceFromTrack(int itrack) // deprecated - not used
{
	if(m_isEtaOut)
	{
		float distance_from_track = m_slope*gr_track->GetPointX(3) - gr_track->GetPointY(3) + m_intercept/(TMath::Sqrt(m_slope*m_slope + 1));
		if(m_many_sclusters)
		{
			distance_from_track = trackSlope(itrack)*getTrack(itrack)->GetPointX(3) - getTrack(itrack)->GetPointY(3) + trackIntercept(itrack)/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));

		}
		histograms->h_d_track_etaout->Fill(distance_from_track);
		if( distance_from_track <= 2. && distance_from_track >= -2.)
			m_accept_eta_out = true;
		else m_accept_eta_out = false;
	}
	if(m_isEtaIn)
	{
		float distance_from_track=9999.;
		if(m_isEtaOut)
		{
			if(!m_many_sclusters)
				distance_from_track = m_slope*gr_track->GetPointX(4) - gr_track->GetPointY(4) + m_intercept/(TMath::Sqrt(m_slope*m_slope + 1));
			else
				distance_from_track = trackSlope(itrack)*getTrack(itrack)->GetPointX(4) - getTrack(itrack)->GetPointY(4) + trackIntercept(itrack)/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));
		}
		else {
			if(!m_many_sclusters)
				distance_from_track = m_slope*gr_track->GetPointX(3) - gr_track->GetPointY(3) + m_intercept/(TMath::Sqrt(m_slope*m_slope + 1));
			else 
			{
				distance_from_track = trackSlope(itrack)*getTrack(itrack)->GetPointX(3) - getTrack(itrack)->GetPointY(3) + trackIntercept(itrack)/(TMath::Sqrt(trackSlope(itrack)*trackSlope(itrack) + 1));
			}

		}
		histograms->h_d_track_etain->Fill(distance_from_track);
		if(distance_from_track <= 2. && distance_from_track>=-2.)
			m_accept_eta_in = true;
		else m_accept_eta_in = false;
	}

}

inline void Track::findSuperClustersWithin1DegWindow()
{
	float pos_sby1 = m_cl_sby1->getPosition(0);
	if(m_cl_sby2->getNClusters2()==1 && m_cl_sby3->getNClusters2()==1 ) //check the events with exactly one cluster on SBY2 and SBY3
	{
		m_many_sclusters = false;
		float pos_sby2 = m_cl_sby2->getCorrPosition(0);
		float pos_sby3 = m_cl_sby3->getCorrPosition(0);
		float totalCharge = m_cl_sby1->getTotPdo(0) + m_cl_sby2->getTotPdo(0) + m_cl_sby3->getTotPdo(0);
		if((pos_sby2>=beam_acceptance_down && pos_sby2<=beam_acceptance_up) && (pos_sby3>=beam_acceptance_down && pos_sby3<=beam_acceptance_up))
		{
			m_totCharge = totalCharge;
			m_map_z_y[m_cl_sby1->getZPosition()] = pos_sby1;
			m_map_z_errory[m_cl_sby1->getZPosition()] = (totalCharge - m_cl_sby1->getTotPdo(0) ) / totalCharge; 
			float SBY2_within_1deg = TMath::ATan( (180./TMath::Pi())* (pos_sby1-pos_sby2)/( std::abs(m_cl_sby2->getZPosition()-m_cl_sby1->getZPosition())) );
			float SBY3_within_1deg = TMath::ATan( (180./TMath::Pi())* (pos_sby1-pos_sby3)/( std::abs(m_cl_sby3->getZPosition()-m_cl_sby1->getZPosition())) );
			if(std::abs(SBY2_within_1deg)<=1. && std::abs(SBY3_within_1deg)<=1. )
			{
				m_counter_1deg++;
				histograms->h_residual_sby2_sby3->Fill(pos_sby2-pos_sby3);
				histograms->h_residual_sby2_sby3_accepted->Fill(pos_sby2-pos_sby3);
				m_map_z_y[m_cl_sby2->getZPosition()] = pos_sby2;
				m_map_z_y[m_cl_sby3->getZPosition()] = pos_sby3;
				m_map_z_errory[m_cl_sby2->getZPosition()] = (totalCharge - m_cl_sby2->getTotPdo(0) ) / totalCharge;
				m_map_z_errory[m_cl_sby3->getZPosition()] = (totalCharge - m_cl_sby3->getTotPdo(0) ) / totalCharge;
			}
			else {
				m_counter_outside1deg_theonly_scluster++;
				histograms->h_residual_sby2_sby3_outside1deg_theonly_sclus->Fill(pos_sby2-pos_sby3);
				histograms->h_residual_sby2_sby3_accepted->Fill(pos_sby2-pos_sby3);
				m_map_z_y[m_cl_sby2->getZPosition()] = pos_sby2;
				m_map_z_y[m_cl_sby3->getZPosition()] = pos_sby3;
				m_map_z_errory[m_cl_sby2->getZPosition()] = (totalCharge - m_cl_sby2->getTotPdo(0) ) / totalCharge;
				m_map_z_errory[m_cl_sby3->getZPosition()] = (totalCharge - m_cl_sby3->getTotPdo(0) ) / totalCharge;
				histograms->h_angle_only_scluster_sby2->Fill(SBY2_within_1deg);
				histograms->h_angle_only_scluster_sby3->Fill(SBY3_within_1deg);
				histograms->h_dangle_only_scluster_sby2sby3->Fill(SBY2_within_1deg-SBY3_within_1deg);
				m_the_only_supercluster_outside_1deg_window = true;
			}
		}

	}
	else {
		m_counter_more_clusters++; // more clusters on the SBY2 and SBY3

		for(int icl2=0; icl2<m_cl_sby2->getNClusters2(); icl2++)
		{
			float pos_sby2 = m_cl_sby2->getCorrPosition(icl2);
			if((pos_sby2>=beam_acceptance_down && pos_sby2<=beam_acceptance_up))
			{
				float SBY2_within_1deg = TMath::ATan( (180./TMath::Pi())* (pos_sby1-pos_sby2)/( std::abs(m_cl_sby2->getZPosition()-m_cl_sby1->getZPosition())) );
			
				for(int icl3=0; icl3<m_cl_sby3->getNClusters2(); icl3++)
				{
					float pos_sby3 = m_cl_sby3->getCorrPosition(icl3);
					if((pos_sby3>=beam_acceptance_down && pos_sby3<=beam_acceptance_up))
					{
						float SBY3_within_1deg = TMath::ATan( (180./TMath::Pi())* (pos_sby1-pos_sby3)/( std::abs(m_cl_sby3->getZPosition()-m_cl_sby1->getZPosition())) );

						if(std::abs(SBY2_within_1deg)<=1. && std::abs(SBY3_within_1deg)<=1. )
						{
							m_pos_sby2_sby3.insert(std::pair<float,float>(pos_sby2,pos_sby3));// [pos_sby2] = pos_sby3;
							m_counter_more_clusters_in_1degWindow++;
							m_many_sclusters = true;
						}
						else {
							m_counter_more_clusters_out_1degWindow++;
							histograms->h_residual_sby2_sby3_out_1deg_many_sclus->Fill(pos_sby2-pos_sby3);
							histograms->h_residual_sby2_sby3_rejected->Fill(pos_sby2-pos_sby3);
						}
					} // end geometrical acceptance cut for SBY3
				} // end SBY3 for-loop
			}  // end geometrical acceptance cut for SBY2
		} // end SBY2 for-loop
	} // end else
}

inline void Track::findTheClosestPair() // among the various superclusters I start with the one with the smallest position difference
{
	for(const auto& scluster : m_pos_sby2_sby3)
	{
		float pos_sby2 = scluster.first;
		float pos_sby3 = scluster.second;
		m_v_pos_abs_diff.push_back(std::abs(pos_sby2-pos_sby3));
	}

	m_counter = m_v_pos_abs_diff.size();
	float min = *min_element(m_v_pos_abs_diff.begin(), m_v_pos_abs_diff.end());
//	std::cout<<"min ->"<<min<<std::endl;
	int index = find(m_v_pos_abs_diff.begin(),m_v_pos_abs_diff.end(),min) - m_v_pos_abs_diff.begin(); 
	auto iter = m_pos_sby2_sby3.begin();
	std::advance(iter, index);
	m_v_pos_diff_ref.push_back(min);
	m_v_pos_abs_diff.at(index) = 100;
	
	int index_cl2 = 0;
	int index_cl3 = 0;
	
	for(int icl2 =0; icl2 < m_cl_sby2->getNClusters2(); icl2++)
	{
		if(m_cl_sby2->getCorrPosition(icl2) == iter->first)
		{
			index_cl2 = icl2;
			break;
		}
	}
	for(int icl3 =0; icl3 < m_cl_sby3->getNClusters2(); icl3++)
	{
		if(m_cl_sby3->getCorrPosition(icl3) == iter->second)
		{
			index_cl3 = icl3;
			break;
		}
	}

	float totalCharge = m_cl_sby1->getTotPdo(0) + m_cl_sby2->getTotPdo(index_cl2) + m_cl_sby3->getTotPdo(index_cl3);
	m_totCharge = totalCharge;	
	m_map_z_y[m_cl_sby1->getZPosition()] = m_cl_sby1->getPosition(0);
	m_map_z_errory[m_cl_sby1->getZPosition()] = (totalCharge - m_cl_sby1->getTotPdo(0) ) / totalCharge;
	m_map_z_y[m_cl_sby2->getZPosition()] = m_cl_sby2->getCorrPosition(index_cl2);
	m_map_z_y[m_cl_sby3->getZPosition()] = m_cl_sby3->getCorrPosition(index_cl3);
	m_map_z_errory[m_cl_sby2->getZPosition()] = (totalCharge - m_cl_sby2->getTotPdo(index_cl2) ) / totalCharge;
	m_map_z_errory[m_cl_sby3->getZPosition()] = (totalCharge - m_cl_sby3->getTotPdo(index_cl3) ) / totalCharge;
}

inline void Track::tryAnotherPair()
{
//	std::cout<<"try other pair "<<m_v_pos_abs_diff.size()<<std::endl;
	m_reject_track = false;
	m_map_z_y.clear();
	m_map_z_errory.clear();

	float min = *min_element(m_v_pos_abs_diff.begin(), m_v_pos_abs_diff.end());
	//std::cout<<"min"<<min<<std::endl;
	int index = find(m_v_pos_abs_diff.begin(),m_v_pos_abs_diff.end(),min) - m_v_pos_abs_diff.begin(); 
	auto iter = m_pos_sby2_sby3.begin();
	std::advance(iter, index);
	if(min!=100)
	{
		m_v_pos_diff_ref.push_back(min);
		m_v_pos_abs_diff.at(index) = 100;
		m_counter--;
	}
	
//	std::cout<<"ref size after"<<m_v_pos_diff_ref.size()<<std::endl;
	if(min!=100) {
		int index_cl2 = 0;
		int index_cl3 = 0;
		
		for(int icl2 =0; icl2 < m_cl_sby2->getNClusters2(); icl2++)
		{
			if(m_cl_sby2->getCorrPosition(icl2) == iter->first)
			{
				index_cl2 = icl2;
				break;
			}
		}
		for(int icl3 =0; icl3 < m_cl_sby3->getNClusters2(); icl3++)
		{
			if(m_cl_sby3->getCorrPosition(icl3) == iter->second)
			{
				index_cl3 = icl3;
				break;
			}
		}
		//std::cout<<index_cl2<<" "<<index_cl3<<std::endl;
		float totalCharge = m_cl_sby1->getTotPdo(0) + m_cl_sby2->getTotPdo(index_cl2) + m_cl_sby3->getTotPdo(index_cl3);
		m_totCharge = totalCharge;
		m_map_z_y[m_cl_sby1->getZPosition()] = m_cl_sby1->getPosition(0);
		m_map_z_errory[m_cl_sby1->getZPosition()] = (totalCharge - m_cl_sby1->getTotPdo(0) ) / totalCharge;
		m_map_z_y[m_cl_sby2->getZPosition()] = m_cl_sby2->getCorrPosition(index_cl2);
		m_map_z_y[m_cl_sby3->getZPosition()] = m_cl_sby3->getCorrPosition(index_cl3);
		m_map_z_errory[m_cl_sby2->getZPosition()] = (totalCharge - m_cl_sby2->getTotPdo(index_cl2) ) / totalCharge;
		m_map_z_errory[m_cl_sby3->getZPosition()] = (totalCharge - m_cl_sby3->getTotPdo(index_cl3) ) / totalCharge;

		findCandidateTrack();
		fitTrack();
	}
	else{
		m_reject_track = true;
	}

}

inline void Track::findCandidateTrack() 
{
	int npoints = m_map_z_y.size();
	float v_z[m_map_z_y.size()];
	float v_y[m_map_z_y.size()];
	float v_z_error[m_map_z_y.size()];
	float v_y_error[m_map_z_y.size()];
	//if(m_many_sclusters) std::cout<<"candidate track has "<<npoints<<" points\n";
	int ip = 0;
	
	for(const auto& ipoint: m_map_z_y) 
	{
		v_z[ip] = ipoint.first;
		v_y[ip] = ipoint.second;
		ip++;
    }
    
    ip=0;
    
    for(const auto& ierr: m_map_z_errory) 
	{
		if(ip!=0)
			v_z_error[ip] = 10.;
		else v_z_error[ip] = 0.;
		v_y_error[ip] = 0.4;//ierr.second/2;
		ip++;
	}
	gr_track = new TGraphAsymmErrors(npoints, v_z, v_y, v_z_error, v_z_error, v_y_error, v_y_error);
	gr_track->GetXaxis()->SetTitle("z [mm]");
	gr_track->GetYaxis()->SetTitle("y [mm]");

}

inline void Track::fitTrack()
{
	TF1* f_track = new TF1("f_track", "[0] + [1]*x", 0, 3800);
	gr_track->Fit("f_track","Q");
	
	float constant = f_track->GetParameter(0);
	float slope = f_track->GetParameter(1);
	float ndf = f_track->GetNDF();
	float chi2 = f_track->GetChisquare();
	float prob = f_track->GetProb();

	m_slope = slope;
	m_intercept = constant;
	m_angle = (180./TMath::Pi())*TMath::ATan(slope);
	m_prob = prob;

	if(prob<0.05){
		gr_track->Fit("f_track","Q");
		prob = f_track->GetProb();
		ndf = f_track->GetNDF();
		chi2 = f_track->GetChisquare();
		slope = f_track->GetParameter(1);
		constant = f_track->GetParameter(0);
		m_slope = slope;
		m_intercept = constant;
		m_angle = (180./TMath::Pi())*TMath::ATan(slope);
		m_prob = prob;
		
		if(prob<0.05) 
		{
	 		fillResidualsZeroProb();
	 		m_reject_track = true;
	 	}
	}


	if(!m_reject_track)
	{
		m_accepted_tracks++;
		if(m_many_sclusters){
		 m_v_gr_candidate_tracks.push_back(gr_track);
		 m_v_slopes_tracks.push_back(m_slope);
		 m_v_intercept_tracks.push_back(m_intercept);
		 m_v_prob_tracks.push_back(m_prob);
		 m_v_angle_tracks.push_back(m_angle);
		 m_count_many_supercluster++;
		 fillResiduals();
		}

		if(m_the_only_supercluster_outside_1deg_window) m_count_singleTracks_outside_1deg++;
		if(!m_many_sclusters && !m_the_only_supercluster_outside_1deg_window) m_count_singleTracks_inside_1deg++;
		histograms->h_chi2ndf->Fill(chi2/ndf);
		histograms->h_chi2->Fill(chi2);
		histograms->h_angle->Fill(m_angle);
		histograms->h_prob->Fill(prob);
	}
	else m_reject_track = true;
}
#endif

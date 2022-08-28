#ifndef CLUSTER_H
#define CLUSTER_H

#include "Layer.h"

#define nholes 2 //consecutive holes

class Layer;
class Cluster{
public:
	Cluster(Layer* layer, std::vector<unsigned int> *pdos, std::vector<unsigned int> *relbcid, std::vector<unsigned int> *tdo, std::vector<unsigned int> *radii); //constructor
	virtual ~Cluster()=default; // destructor 
	inline void sortAndSelectStrips(std::vector<int> v_fired_strips);
	inline int getNClusters2(){return m_v_clusters.size();};
	inline void removeDublicates(std::vector<int> &v);
	inline int getNSingleStrips(){return m_count_single_strips;};
	inline void fillClusters(std::vector<int> v_sorted_fired_strips);
	inline void fillClusterSpecs(std::vector<unsigned int> *pdos, std::vector<unsigned int> *relbcid, std::vector<unsigned int> *tdo, std::vector<unsigned int> *radii);
	inline std::vector<std::vector<int> > getClusters(){return m_v_clusters;};
	inline std::vector<int> getCluster(int subcluster) {return m_v_clusters[subcluster];}; 
	inline float getTotPdo(int subcluster);
	inline int getStripIndex(int subcluster, int strip) {return m_v_clusters[subcluster].at(strip);};
	inline float getStripPdo(int subcluster, int strip) {return m_v_pdos[subcluster].at(strip);};
	inline float getStripTdo(int subcluster, int strip) {return m_v_tdo[subcluster].at(strip);};
	inline int getStripRadius(int subcluster, int strip) {return m_v_radii[subcluster].at(strip);};
	inline float getStripRelBcid(int subcluster, int strip) {return m_v_relbcid[subcluster].at(strip);};
	inline float getPosition(int subcluster);
	inline float getCorrPosition(int subcluster);
	inline float getStrip_from_clpos(float cluster_position);
	inline int getNStrips(int subcluster){return m_v_clusters[subcluster].size();};
	inline void removeStrips_small(int subcluster);
	inline float getZPosition() {return this_layer->convertLayerToGlobalZ();};
	Layer* this_layer;
	
private:
	std::vector<int> m_v_fired_strips;
	std::vector<int> m_v_sorted_strips;
	std::vector<std::vector<int> > m_v_clusters;
	std::vector<std::vector<int> > m_v_pdos;
	std::vector<std::vector<int> > m_v_relbcid;
	std::vector<std::vector<int> > m_v_tdo;
	std::vector<std::vector<int> > m_v_radii;

	int m_count_single_strips;
	int m_nclusters;

};

inline Cluster::Cluster(Layer* layer, std::vector<unsigned int> *pdos, std::vector<unsigned int> *relbcid, std::vector<unsigned int> *tdo, std::vector<unsigned int> *radii)
{
	this_layer = layer;
	m_v_fired_strips = this_layer->getFiredStrips();

	if(!m_v_fired_strips.empty()) sortAndSelectStrips(m_v_fired_strips);
	
	if(!m_v_sorted_strips.empty()) {
		fillClusters(m_v_sorted_strips);
		fillClusterSpecs(pdos, relbcid, tdo, radii);
	}

}

inline void Cluster::fillClusterSpecs(std::vector<unsigned int> *pdos, std::vector<unsigned int> *relbcid, std::vector<unsigned int> *tdo, std::vector<unsigned int> *radii)
{
	std::vector<int> v_tmp_pdo;
	std::vector<int> v_tmp_tdo;
	std::vector<int> v_tmp_relbcid;
	std::vector<int> v_tmp_radii;

	for(int iclus = 0; iclus < m_v_clusters.size(); iclus++) {
		for(int istrip = 0; istrip < m_v_clusters[iclus].size(); istrip++)
		{
			int index = this_layer->map_fired_strips_to_indices.find(m_v_clusters[iclus].at(istrip))->second ;			
			v_tmp_pdo.push_back(pdos->at(index));
			v_tmp_relbcid.push_back(relbcid->at(index));
			v_tmp_tdo.push_back(tdo->at(index));
			v_tmp_radii.push_back(radii->at(index));
		}
		m_v_pdos.push_back(v_tmp_pdo);
		m_v_relbcid.push_back(v_tmp_relbcid);
		m_v_tdo.push_back(v_tmp_tdo);
		m_v_radii.push_back(v_tmp_radii);
		v_tmp_pdo.clear();
		v_tmp_tdo.clear();
		v_tmp_relbcid.clear();
		v_tmp_radii.clear();
	}
}

inline float Cluster::getTotPdo(int subcluster)
{
	int totCharge = 0;
	for(int istrip = 0; istrip < m_v_pdos[subcluster].size(); istrip++)
	{
		totCharge += m_v_pdos[subcluster].at(istrip);
	}
	return totCharge;
}

inline void Cluster::removeStrips_small(int subcluster)
{
	int totCharge = getTotPdo(subcluster);
	int single_strip_before = 0;
	if( this_layer->getLayerIndex()==0 && totCharge<=280)
	{
		m_v_clusters[subcluster].clear();
		m_v_pdos[subcluster].clear();
		m_v_relbcid[subcluster].clear();
		m_v_tdo[subcluster].clear();
		m_v_radii[subcluster].clear();
		m_v_clusters[subcluster].shrink_to_fit();
		m_v_pdos[subcluster].shrink_to_fit();
		m_v_relbcid[subcluster].shrink_to_fit();
		m_v_tdo[subcluster].shrink_to_fit();
		m_v_radii[subcluster].shrink_to_fit();
	}
	else if( this_layer->getLayerIndex()==2 && totCharge<=170)
	{
		m_v_clusters[subcluster].clear();
		m_v_pdos[subcluster].clear();
		m_v_relbcid[subcluster].clear();
		m_v_tdo[subcluster].clear();
		m_v_radii[subcluster].clear();
		m_v_clusters[subcluster].shrink_to_fit();
		m_v_pdos[subcluster].shrink_to_fit();
		m_v_relbcid[subcluster].shrink_to_fit();
		m_v_tdo[subcluster].shrink_to_fit();
		m_v_radii[subcluster].shrink_to_fit();
	}
}


inline float Cluster::getPosition(int subcluster)
{
	float sumpdo = 0; 
	int totCharge = getTotPdo(subcluster);

	for(int istrip = 0; istrip < m_v_clusters[subcluster].size(); istrip++){
	if(m_v_radii[subcluster].at(istrip)==4 || m_v_radii[subcluster].at(istrip)==5) sumpdo += m_v_pdos[subcluster].at(istrip)*this_layer->convertStripToGlobalY_radius(m_v_clusters[subcluster].at(istrip), m_v_radii[subcluster].at(istrip));
	else 
	 sumpdo += m_v_pdos[subcluster].at(istrip)*this_layer->convertStripToGlobalY(m_v_clusters[subcluster].at(istrip));//, getStripRadius(subcluster,m_v_clusters[subcluster].at(istrip)));

	}
	return sumpdo/totCharge;
}

inline float Cluster::getCorrPosition(int subcluster)
{
	float sumpdo = 0; 
	int totCharge = getTotPdo(subcluster);
	for(int istrip = 0; istrip < m_v_clusters[subcluster].size(); istrip++){
	if(m_v_radii[subcluster].at(istrip)==4 || m_v_radii[subcluster].at(istrip)==5) sumpdo += m_v_pdos[subcluster].at(istrip)*this_layer->convertStripToGlobalY_radius_corr(m_v_clusters[subcluster].at(istrip), m_v_radii[subcluster].at(istrip));
	else
	 sumpdo += m_v_pdos[subcluster].at(istrip)*this_layer->convertStripToGlobalY_corr(m_v_clusters[subcluster].at(istrip));

	}
	return sumpdo/totCharge;
}

inline float Cluster::getStrip_from_clpos(float cluster_position) // this function works only for events with one cluster on the under study layer
{
	if(this_layer->isSM1()) {
		if(m_v_radii[0].at(0)==4) {
			if(this_layer->getLayerIndex()==0 || this_layer->getLayerIndex()==1) return (cluster_position-895-35.3)/0.425;
			else if(this_layer->getLayerIndex()==2 || this_layer->getLayerIndex()==3) return (cluster_position-895-35.3)/0.425;
		}
		else if(m_v_radii[0].at(0)==5) {
			if(this_layer->getLayerIndex()==0) return (cluster_position-895-35.3)/0.425;
			else if(this_layer->getLayerIndex()==1) return (cluster_position-895-35.3)/0.425;
			else if(this_layer->getLayerIndex()==2 || this_layer->getLayerIndex()==3) return (cluster_position-895-35.3)/0.425;
		}
	}
	return 0;
}

inline void Cluster::removeDublicates(std::vector<int> &v)
{
    auto end = v.end();
    for (auto it = v.begin(); it != end; ++it) {
        end = std::remove(it + 1, end, *it);
    }
 
    v.erase(end, v.end());
}

inline void Cluster::sortAndSelectStrips(std::vector<int> v_fired_strips)
{
	std::vector<int> sorted_strips;

	sort(v_fired_strips.begin(), v_fired_strips.end());

	for(int iStrip = v_fired_strips.size()-1; iStrip--;)
	{
		if(std::abs(v_fired_strips.at(iStrip+1)-v_fired_strips.at(iStrip))>nholes+1) continue;
		sorted_strips.push_back(v_fired_strips.at(iStrip+1));
		sorted_strips.push_back(v_fired_strips.at(iStrip));
	}
	if(!sorted_strips.empty()) removeDublicates(sorted_strips);

	m_count_single_strips = (v_fired_strips.size()-sorted_strips.size());

	sort(sorted_strips.begin(), sorted_strips.end());

	m_v_sorted_strips = sorted_strips;
}

inline void Cluster::fillClusters(std::vector<int> v_sorted_fired_strips)
{
	m_nclusters = 1;
	std::vector<int> v_pos;
	std::vector<int> v_tmp_strips_in_cluster;

	for(int istrip = 0; istrip<v_sorted_fired_strips.size()-1; istrip++)
	{
		if(std::abs(v_sorted_fired_strips.at(istrip+1)-v_sorted_fired_strips.at(istrip))>nholes+1)
		{
			m_nclusters++;
			v_pos.push_back(istrip+1);
		}
	}
	v_pos.push_back(v_sorted_fired_strips.size());

	int it=0;

	for(int iclus=0; iclus<m_nclusters; iclus++)
	{
		for(int istrip = it; istrip<v_pos.at(iclus); istrip++) v_tmp_strips_in_cluster.push_back(v_sorted_fired_strips.at(istrip));
		m_v_clusters.push_back(v_tmp_strips_in_cluster);
		v_tmp_strips_in_cluster.clear();
		it = v_pos.at(iclus);
	}

}

#endif

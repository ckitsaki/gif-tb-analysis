#ifndef LAYER_H
#define LAYER_H

class Layer{
public:
	Layer(int layerIndex, bool isSM1/*, std::vector<int> v_hot*/); //constructor
	virtual ~Layer()=default; // destructor 
	inline void setLayerIndex(int layerIndex) { m_layerIndex = layerIndex;};
	inline int getLayerIndex() { return m_layerIndex;};
	inline void AddHitIndex(int iHit) { m_v_hit_indices.push_back(iHit);};
	inline int getNHits(){ return m_v_hit_indices.size();};
	inline void bookFiredStrips(std::vector<int> *strips, std::vector<unsigned int> *pdo);
	inline std::vector<int> getFiredStrips() {return m_v_fired_strips;};
	inline bool isDisconnectedStrip(int strip);
	inline float convertStripToGlobalY(int strip);
	inline float convertStripToGlobalY_corr(int strip);
	inline float convertStripToGlobalY_radius(int strip, int radius);
	inline float convertStripToGlobalY_radius_corr(int strip, int radius);
	inline float convertLayerToGlobalZ();
	inline int getHitIndex(int strip) { return map_fired_strips_to_indices.find(strip)->second;};
	inline void isTrigger(){m_isTrigger=true;};
	inline bool isSM1(){return m_isSM1;};
	inline void setAlphaBeta(float alpha, float beta) { m_alpha = alpha; m_beta = beta;};
	inline float getAlpha() {return m_alpha;};
	inline float getBeta() {return m_beta;};

	//inline void addStripsToReject(std::vector<int> v_hot) { for(int hots:v_hot) m_v_strips_to_reject.push_back(hots);};
	std::map<int, int> map_fired_strips_to_indices;

private:
	std::vector<int> m_v_hit_indices;
	std::vector<int> m_v_fired_strips;

	//std::vector<int> m_v_strips_to_reject;
	bool m_isSM1 = false;
	bool m_isTrigger = false;
	int m_layerIndex;
	float m_alpha;
	float m_beta;
};

inline Layer::Layer(int layerIndex, bool isSM1 /*, std::vector<int> v_hot*/)
{
	setLayerIndex(layerIndex);
	//addStripsToReject(v_hot);
	m_isSM1 = isSM1;
	map_fired_strips_to_indices = map<int,int>();
}

inline float Layer::convertLayerToGlobalZ() // return center of drift gaps
{
	float lever_arm_0 = 2153.;
	float lever_arm_1 = 790.;
	float SM1_width = 78.;
	if(!m_isSM1) 
	{
		if(m_layerIndex==0) return 0; 
		if(m_layerIndex==1) return 262;
		if(m_layerIndex==2) return 262+lever_arm_0+SM1_width+lever_arm_1;
		if(m_layerIndex==3) return 262+lever_arm_0+SM1_width+lever_arm_1+136;
	}
	else{
	if(m_layerIndex==0) return lever_arm_0+262+13.835;
	if(m_layerIndex==1) return lever_arm_0+262+30.615;
	if(m_layerIndex==2) return lever_arm_0+262+46.985;
	if(m_layerIndex==3) return lever_arm_0+262+63.765;
	}

	return 0;
}

inline float Layer::convertStripToGlobalY_radius(int strip, int radius) 
{

	if(radius==4) {
	if(m_layerIndex==0) //eta layers
	{
		// SM1
		return 895+35.3+(strip)*0.425; 

	}
	else if(m_layerIndex==1)
	{
		// SM1
		return 895+35.3+(strip)*0.425; 
	}
	else if(m_layerIndex==2) //stereo layers
	{
		// SM1
		return 895+35.3+(strip)*0.425; 
	}
	else if(m_layerIndex==3) //stereo layers
	{
		return 895+35.3+(strip)*0.425;
	}
	} //end radius 4
	else if(radius==5) {
		if(m_layerIndex==0 ) //eta layers
	{
		// SM1
		return 895+35.3+(strip)*0.425; 
	}
	else if(m_layerIndex==1)
	{
		// SM1
		return 895+35.3+(strip+1)*0.425; 
	}

	if(m_layerIndex==2 ) //stereo layers
	{
		return 895+35.3+(strip-1)*0.425;
	}
	else if (m_layerIndex==3)
	{
		if(strip >= 2560 && strip <=2663) return 895+35.3+(strip+1)*0.425;
		return 895+35.3+(strip)*0.425;
	}
	
	}//end radius 5
	return 0;

}

inline float Layer::convertStripToGlobalY_radius_corr(int strip, int radius)
{
	if(radius==4) {
	if(m_layerIndex==0) //eta layers
	{
		// SM1
		return (m_alpha + (895+35.3+(strip)*0.425)*m_beta ); 
	}
	else if(m_layerIndex==1)
	{
		// SM1
		return (m_alpha + (895+35.3+(strip)*0.425)*m_beta ); 
	}
	else if(m_layerIndex==2) //stereo layers
	{
		// SM1
		return (m_alpha + (895+35.3+(strip)*0.425)*m_beta ); 
	}
	else if(m_layerIndex==3) //stereo layers
	{
	    return (m_alpha + (895+35.3+(strip)*0.425)*m_beta );
	}
	} //end radius 4
	else if(radius==5) {
		if(m_layerIndex==0 ) //eta layers
	{
		// SM1
		return (m_alpha + (895+35.3+(strip)*0.425)*m_beta ); 
	}
	else if(m_layerIndex==1)
	{
		// SM1
		return (m_alpha + (895+35.3+(strip+1)*0.425)*m_beta );
	}

	if(m_layerIndex==2 ) //stereo layers
	{
		// SM1
		return (m_alpha + (895+35.3+(strip-1)*0.425)*m_beta ); 
	}
	else if (m_layerIndex==3)
	{
		if(strip >= 2560 && strip <=2663) return (m_alpha + (895+35.3+(strip+1)*0.425)*m_beta);
	    else return (m_alpha + (895+35.3+(strip)*0.425)*m_beta );
	}
	
	}//end radius 5
	return 0;

}

inline float Layer::convertStripToGlobalY_corr(int strip) 
{

	if(!m_isSM1) {
		if(strip<1536) return (m_alpha + (895+35.3+(std::abs(1535-strip)+1024)*0.45)*m_beta);
		else return (m_alpha + (895+35.3+(strip)*0.45)*m_beta); 
	}
	return 0;
}


inline float Layer::convertStripToGlobalY(int strip) 
{
	
// SM
/////////////////
if(m_isSM1) 
{
	return (895+35.3+(strip)*0.425);
}

// SB
/////////////////
else {
	if(strip<1536) return (895+35.3+(std::abs(1535-strip)+1024)*0.45);
	else return (895+35.3+(strip)*0.45); 
}
	return 0;
}

inline bool Layer::isDisconnectedStrip(int strip)
{
	if(!m_isSM1)
	{
		//if(m_layerIndex==1 && (strip==1666 || strip==1152 || strip==1153 || strip==1155 || strip==1156 || strip==1538 || strip==1540 || strip==1539 || strip==1833 || strip==1664 || strip==1665 || strip==1667 || strip==1670 || strip==1671 || strip==1673 || strip==1674 || strip==1675 || strip==1676 || strip==1677 || strip==1683 || strip==1684 || strip==1687 || strip==1688 || strip==1689 || strip==1626 || strip==1627 || strip==1624 || strip==1623 || strip==1610 || strip==1611 || strip==1613 || strip==1614 || strip==1615 || strip==1617 || strip==1618) ) return true;
		if(m_layerIndex==1 && (strip==1833 || strip==1999)) return true;
		else if(m_layerIndex==2 && (strip==2043 || strip==1783 || strip==1796 || strip==1797) ) return true;
		//else if(m_layerIndex==2 && ( (strip>=1781 && strip<=1797) || strip==2044 || strip==1915 || strip==1913 || strip==1911 || strip==1910 || strip==1900 || strip==1898 || strip==1897 || strip==1361 || strip==1362 || strip==1359 || strip==1287 || strip==1288 || strip==1270 || strip==1269 || strip==1263 || strip==1261 || strip==1024 || strip==1026 || strip==1027 || strip==1029 || strip==1034 || strip==1035 || strip==1037 || strip==1039 || strip==1040 || strip==1063 || strip==1064 || strip==1067 || strip==1069 || strip==1070 || strip==1072 || strip==1120 || strip==1122)) return true;
		//else if(m_layerIndex==3 && ( (strip>=1601 && strip<=1607) || strip==1045 || strip==1096 || strip==1098 || strip==1254 || strip==1261 || strip==1702 || strip==1704 || strip==1717 || strip==1523 || strip==1718 || strip==1721 || strip==1729 || strip==1730 || strip==1768 || strip==1770 || strip==1779 || strip==1781 || strip==1783 || strip==1790 || strip==1791 || strip==1795 || strip==1797)) return true;
		else if(m_layerIndex==3 && strip==1522) return true;
	}
	else if(m_isSM1)
	{
		//if(strip==1023 || strip==1024 || strip==2047 || strip==2048 || strip==3071 || strip==3072 || strip==4095 || strip==4096) return true; //strips at the interconnection points
		if(m_layerIndex==0 && (strip==2813 || strip==2568 || strip==2569 || strip==2572) ) return true;
		else if(m_layerIndex==1 && (strip==2153 || strip==2154 || strip==2551 || strip==2556 || strip==2924 || strip==2948 || strip==2950)) return true;//(strip==2950 || strip==2948 || strip==2924 || strip==2551 || strip==2556 || strip==2153 || strip==2154)) return true;
		//else if(m_layerIndex==2 && (strip==2559 || strip==2560)) return true;
		else if(!m_isTrigger && m_layerIndex==3 && (strip==2091 || strip==2227 || strip==2659 || strip==2661)) return true;
		//if(m_layerIndex==0 && ((strip>=2568 && strip<=2579) || strip==2682 || strip==2684 || (strip>=2656 && strip<=2659) || strip==2569 || strip==2572 || strip==2813 || strip==2069 || strip==2094 || strip==2084 || strip==2085 || strip==2081 || strip==2150 || strip==2216)) return true;
		//else if(m_layerIndex==1 && (strip==2927 || strip==2928 || strip==2930 || strip==2545 || strip==2549 || strip==2560 || strip==2153 || strip==2154 || strip==2552 || strip==2557 || strip==2559 || strip==2758 || strip==2760 || strip==2814 || strip==2816 || strip==2924 || strip==2926 || strip==2948 || strip==2950)) return true;
		//else if(m_layerIndex==2 && (strip==2739 || strip==2737)) return true;
		//else if(!m_isTrigger && m_layerIndex==3 && ( (strip>=2665 && strip<=2671) || strip==2866 || strip==2659 || strip==2661 || strip==2868 || strip==2869)) return true;
	}

	return false;
}

inline void Layer::bookFiredStrips(std::vector<int> *strips, std::vector<unsigned int> *pdo)
{	
	for(int ind=0; ind<m_v_hit_indices.size(); ind++)
	{
		int index = m_v_hit_indices.at(ind);
	//	if(std::find(m_v_strips_to_reject.begin(), m_v_strips_to_reject.end(), strips->at(index)) != m_v_strips_to_reject.end()) continue;
		if(isDisconnectedStrip(strips->at(index))) continue;
		//if(!m_isSM1 && m_layerIndex==0 && isDisconnectedStrip(strips->at(index))) continue;
		if(pdo->at(index) > 64 )  // cut on PDO
		{
			//std::cout<<isDisconnectedStrip(strips->at(index))<<" "<<pdo->at(index)<<" "<<strips->at(index)<<" "<<m_layerIndex<<std::endl;
			m_v_fired_strips.push_back(strips->at(index));
			map_fired_strips_to_indices[strips->at(index)] = index;			
		} //end if pdo
		else continue;

	}//end for loop
}
#endif

#ifndef TreeReader_h
#define TreeReader_h

class TreeReader {
public :
   TTree          *the_tree;   

   // Declaration of leaf types
   UInt_t          guid;
   UInt_t          run_number;
   UInt_t          lumiblock;
   UInt_t          sourceId;
   Int_t           triggerCalibrationKey;
   std::vector<unsigned int> *linkId;
   std::vector<unsigned int> *fragmentSize;
   std::vector<unsigned int> *level1Id;
   std::vector<unsigned int> *bcid;
   std::vector<unsigned int> *orbit;
   std::vector<unsigned int> *timeset;
   std::vector<unsigned int> *checksum;
   std::vector<unsigned int> *nhits;
   std::vector<unsigned int> *level0Id;
   std::vector<unsigned int> *missing_data_flags;
   std::vector<unsigned int> *to_flag;
   std::vector<unsigned int> *error_flag;
   std::vector<unsigned int> *alive;
   std::vector<unsigned int> *vmmid;
   std::vector<unsigned int> *channel;
   std::vector<unsigned int> *pdo;
   std::vector<unsigned int> *tdo;
   std::vector<unsigned int> *relbcid;
   std::vector<unsigned int> *nflag;
   std::vector<unsigned int> *parity;
   std::vector<unsigned int> *layer;
   std::vector<unsigned int> *radius;
   std::vector<int>     *strip;
 //  std::vector<unsigned int> *type;
 //  std::vector<int>     *sector;
 //  std::vector<int>     *technology;

   // List of branches
   TBranch        *b_guid;   //!
   TBranch        *b_run_number;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_sourceId;   //!
   TBranch        *b_triggerCalibrationKey;   //!
   TBranch        *b_linkId;   //!
   TBranch        *b_fragmentSize;   //!
   TBranch        *b_level1Id;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_timeset;   //!
   TBranch        *b_checksum;   //!
   TBranch        *b_nhits;   //!
   TBranch        *b_level0Id;   //!
   TBranch        *b_missing_data_flags;   //!
   TBranch        *b_to_flag;   //!
   TBranch        *b_error_flag;   //!
   TBranch        *b_alive;   //!
   TBranch        *b_vmmid;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_pdo;   //!
   TBranch        *b_tdo;   //!
   TBranch        *b_relbcid;   //!
   TBranch        *b_nflag;   //!
   TBranch        *b_parity;   //!
   TBranch        *b_layer;   //!
   TBranch        *b_radius;   //!
   TBranch        *b_strip;   //!
  // TBranch        *b_type;   //!
  // TBranch        *b_sector;   //!
  // TBranch        *b_technology;   //!

   TreeReader(std::string file_name, std::string tree_name);
   virtual ~TreeReader()=default;
   inline void Init(TTree *tree);
};

#endif

inline void TreeReader::Init(TTree *tree)
{
   // Set object pointer
   linkId = 0;
   fragmentSize = 0;
   level1Id = 0;
   bcid = 0;
   orbit = 0;
   timeset = 0;
   checksum = 0;
   nhits = 0;
   level0Id = 0;
   missing_data_flags = 0;
   to_flag = 0;
   error_flag = 0;
   alive = 0;
   vmmid = 0;
   channel = 0;
   pdo = 0;
   tdo = 0;
   relbcid = 0;
   nflag = 0;
   parity = 0;
   layer = 0;
   radius = 0;
   strip = 0;
  // type = 0;
  // sector = 0;
  // technology = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   the_tree = tree;
   
   the_tree->SetBranchAddress("guid", &guid, &b_guid);
   the_tree->SetBranchAddress("run_number", &run_number, &b_run_number);
   the_tree->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   the_tree->SetBranchAddress("sourceId", &sourceId, &b_sourceId);
   the_tree->SetBranchAddress("triggerCalibrationKey", &triggerCalibrationKey, &b_triggerCalibrationKey);
   the_tree->SetBranchAddress("linkId", &linkId, &b_linkId);
   the_tree->SetBranchAddress("fragmentSize", &fragmentSize, &b_fragmentSize);
   the_tree->SetBranchAddress("level1Id", &level1Id, &b_level1Id);
   the_tree->SetBranchAddress("bcid", &bcid, &b_bcid);
   the_tree->SetBranchAddress("orbit", &orbit, &b_orbit);
   the_tree->SetBranchAddress("timeset", &timeset, &b_timeset);
   the_tree->SetBranchAddress("checksum", &checksum, &b_checksum);
   the_tree->SetBranchAddress("nhits", &nhits, &b_nhits);
   the_tree->SetBranchAddress("level0Id", &level0Id, &b_level0Id);
   the_tree->SetBranchAddress("missing_data_flags", &missing_data_flags, &b_missing_data_flags);
   the_tree->SetBranchAddress("to_flag", &to_flag, &b_to_flag);
   the_tree->SetBranchAddress("error_flag", &error_flag, &b_error_flag);
   the_tree->SetBranchAddress("alive", &alive, &b_alive);
   the_tree->SetBranchAddress("vmmid", &vmmid, &b_vmmid);
   the_tree->SetBranchAddress("channel", &channel, &b_channel);
   the_tree->SetBranchAddress("pdo", &pdo, &b_pdo);
   the_tree->SetBranchAddress("tdo", &tdo, &b_tdo);
   the_tree->SetBranchAddress("relbcid", &relbcid, &b_relbcid);
   the_tree->SetBranchAddress("nflag", &nflag, &b_nflag);
   the_tree->SetBranchAddress("parity", &parity, &b_parity);
   the_tree->SetBranchAddress("layer", &layer, &b_layer);
   the_tree->SetBranchAddress("radius", &radius, &b_radius);
   the_tree->SetBranchAddress("strip", &strip, &b_strip);
 //  the_tree->SetBranchAddress("type", &type, &b_type);
 //  the_tree->SetBranchAddress("sector", &sector, &b_sector);
 //  the_tree->SetBranchAddress("technology", &technology, &b_technology);

}

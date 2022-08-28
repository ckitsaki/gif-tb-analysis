#include "TreeReader.h"

TreeReader::TreeReader(std::string file_name, std::string tree_name)
{
   TFile* f = new TFile(file_name.c_str());
   std::cout<<"Opening "<<file_name<<std::endl;
   the_tree = (TTree*)f->Get(tree_name.data());
   Init(the_tree);
}

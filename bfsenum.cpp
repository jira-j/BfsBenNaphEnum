#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <boost/program_options.hpp>
#include "bfsenum.hpp"

#include <time.h>

int main(int argc, char **argv)
{
  using namespace std;
  namespace po = boost::program_options;  

  po::options_description opt("Allowed options");

  opt.add_options()
    ("c,c", po::value<int>(), ": input # C atoms")
    ("n,n", po::value<int>(), ": input # N atoms")
    ("o,o", po::value<int>(), ": input # O atoms")
    ("h,h", po::value<int>(), ": input # H atoms")
    ("r,r", po::value<int>(), ": input round number")
    ("f,f", po::value<string>(), ": input left-heavy tree file")
    ("output,t", ": output results")
    ("help,p", ": show this help message");

  po::variables_map argmap;
  po::store(po::parse_command_line(argc, argv, opt), argmap);
  po::notify(argmap);

  if (argmap.count("help") 
      or ( (not argmap.count("c")) and(not argmap.count("n")) and(not argmap.count("o")) and (not argmap.count("h")) and (not argmap.count("f")) ) 
   ) {
    std::cerr << "Usage: " << argv[0] << " [option]" << std::endl << opt << std::endl;

    return EXIT_FAILURE;
  }
  if( argmap.count("f") && ( argmap.count("c") || argmap.count("n") || argmap.count("o") || argmap.count("h") ) ){
    std::cerr << "tree file and chemical formula cannot be input at the same time " << std::endl;
    
    return EXIT_FAILURE;
  }

  int num_C = 0, num_N = 0, num_O = 0, num_H = 0;
  size_t round_Number;

  if( argmap.count("f") ){
    string filename = argmap["f"].as<string>();
    vector< map<string, int> > node_list; //vector of 3-tuple (label, multi, parent) of all nodes in the input file
    read_tree_file( filename, node_list );
    
    ChemTreeCenter tree( node_list );
    
    cout<<"==================="<<endl;
    cout<<"the input tree is "<<endl;
    tree.show( false );
        
    size_t num = generate_from_tree( tree, argmap.count("output") );

    cout<<"======Result======="<<endl;
    cout<<"Total #enumerated structures from a given tree:" << num << endl;    
  }else{
    if (argmap.count("c")) {
      num_C = argmap["c"].as<int>();
    }
    if (argmap.count("n")) {
      num_N = argmap["n"].as<int>();
    }
    if (argmap.count("o")) {
      num_O = argmap["o"].as<int>();
    }
    if (argmap.count("h")) {
      num_H = argmap["h"].as<int>();
    }
    if (argmap.count("r")) {
      round_Number = argmap["r"].as<int>();
    }else{
      round_Number = 0;
    }

    vector<int> atom_numbers(4);
    atom_numbers[0] = num_C;
    atom_numbers[1] = num_N;
    atom_numbers[2] = num_O;
    atom_numbers[3] = num_H;
    
    size_t num = generate(atom_numbers, argmap.count("output"), round_Number);
    cout<<"======Result======="<<endl;
    cout<<"Total #enumerated structures from a given chemical formula:" << num << endl;    
  }

  return EXIT_SUCCESS;
}


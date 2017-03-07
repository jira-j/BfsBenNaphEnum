#pragma once 
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <array>
#include <valarray>
#include <string>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <time.h>

typedef int valence_value_type;
typedef int label_value_type; // `signed' is required because `-' is used.
typedef std::vector< std::vector<int> > carbon_position;

bool compare_function(std::vector<int>,std::vector<int>);

const int num_distinct_atoms = 5;
const std::string atomchar[] = {"b","C", "N", "O", "H"};
const std::string bondchar[] = {"X", "", "=", "#"};
const int valence[] = {6, 4, 3, 2, 1};
const std::string input_atomchar[] = {"C", "N", "O", "H"};
const std::string input_bondchar[] = {"X", "", "=", "#"};
const int input_valence[] = {4, 3, 2, 1};
const int max_valence = 6;
const int first_atom_valence_one = 4;
const int num_special_atom =1;

//fix cp_list for nodes with smallest adj_set >=3
const std::vector<int> fix311 = {0,1,2};
const std::vector<int> fix312 = {0,1,3};
const std::vector<int> fix313 = {0,2,4};
const std::vector< std::vector<int> > fix_cp3_1 = {fix311, fix312, fix313};
const std::vector<int> fix321 = {3,4,5};
const std::vector<int> fix322 = {2,4,5};
const std::vector<int> fix323 = {1,3,5};
const std::vector< std::vector<int> > fix_cp3_2 = {fix321, fix322, fix323};
const std::vector<int> fix41 = {0,1,2,3};
const std::vector<int> fix42 = {0,1,2,4};
const std::vector<int> fix43 = {0,1,3,4};
const std::vector< std::vector<int> > fix_cp4 = {fix41, fix42, fix43};

static bool do_print = false;
static int num_except_H = 0;
static int num_H = 0;
static int num_lack_H = 0;
static int num_naph_bond = 0;

static std::vector<int> t_r, t_v;
static std::vector<valence_value_type> location;
static std::vector<valence_value_type> lack_degree;

void print_cp(const std::vector< std::vector<int> > & cp){
  for(size_t i=0; i<cp.size(); i++){
    std::cout<<"[";
    for(size_t j=0; j<cp[i].size(); j++){
      std::cout<<cp[i][j];
      if(j+1!=cp[i].size()) std::cout<<" ";
    }
    std::cout<<"]";
  }
  std::cout<<std::endl;
}

void print_cp(const std::vector<int> & cp){
  for(size_t j=0; j<cp.size(); j++){
    std::cout<<cp[j];
    if(j+1!=cp.size()) std::cout<<" ";
  }
  std::cout<<std::endl;
}


void read_tree_file( std::string filename , std::vector< std::map< std::string, int> > & node_list ){
  using namespace std;
  ifstream ifs(filename);
  string buffer;
  int label, multi, parent;

  while( ifs >> label >> multi >> parent ){
    map<string,int> temp;
    temp.insert( pair< string,int >("label", label) );
    temp.insert( pair< string,int >("multi", multi) );
    temp.insert( pair< string,int >("parent", parent) );
    //tuple< int, int, int > temp(label, multi, parent);
    node_list.push_back( temp );
    getline(ifs, buffer);
  }
}

typedef std::valarray<char> is_ident_type;
//typedef std::vector<char> is_ident_type;
class ChemTreeCenter {

	struct Node {
		label_value_type label;
	//bool is_identical;
		std::array<int, max_valence> children;
		valence_value_type num_children;
		int parent;
		valence_value_type multi; // subtree id during construction of simple trees
	        valence_value_type nth; // start from 0 not 1
	//int depth;
		std::array<int, max_valence> bond_position;
	//      std::vector<int> bond_position;
	};

	std::vector<int> rest_atoms;
public:std::valarray<Node> nodes;
	int num_nodes;
#ifdef CUT_BY_NUM_H
	int num_H_to_be_added;
#endif
	//int deepest_head;

	void init() { // initialize nodes
	//nodes.resize(num_except_H);
		num_nodes = 0;
#ifdef CUT_BY_NUM_H
		num_H_to_be_added = 0;
#endif
	//deepest_head = 1;
	}

  ChemTreeCenter( std::vector< std::map< std::string, int> > &node_list){
    num_nodes = node_list.size();
    nodes.resize(num_nodes);

    //add root node
    auto &root = nodes[0];
    root.label = node_list[0]["label"] - 1;
    root.num_children = 0;
    root.parent = -1;
    root.multi = 1;

    for(size_t index = 1; index < num_nodes; index++){

      auto& node = nodes[index];
      node.label = node_list[index]["label"] - 1;
      node.parent = node_list[index]["parent"] - 1;
      node.multi = node_list[index]["multi"];

      auto& parnode = nodes[ node.parent ];
      const valence_value_type nc = parnode.num_children;

      node.nth = nc;
      node.num_children = 0;

      parnode.children[nc] = index;
      ++(parnode.num_children);
    }
  }

	inline void printseq() const;
	inline void printseq_single() const;
	inline void printsmi(const int i = 0) const;
	inline void printsmi_single(const int i = 0) const;
	inline void printmol(const std::string& filename) const;
       
public:

	ChemTreeCenter(const std::initializer_list<int>& _rest, const int _num_except_H) : rest_atoms(_rest), nodes(std::valarray<Node>(_num_except_H)) {
		init();
	}
	ChemTreeCenter(const std::vector<int>& _rest, const int _num_except_H) : rest_atoms(_rest), nodes(std::valarray<Node>(_num_except_H)) {
		init();
	}

	inline void update_identical(is_ident_type& is_ident, const int i) const; // when a new node i is added
	inline void update_identical_multi(is_ident_type& is_ident, const int i) const ; 
	inline void update_identical_end(is_ident_type& is_ident, const int i) const ; // when a child of a node i is recognized not to be added more
	inline bool add_root(const label_value_type atom_label);
	inline bool add_root_child(is_ident_type& is_ident, const label_value_type atom_label);
	inline bool add_node(is_ident_type& is_ident, const int parenti, const label_value_type atom_label);
	inline void del_root() {
		num_nodes = 0;
		++(rest_atoms[nodes[0].label]);
	}
	inline void del_last_node() {
		--num_nodes;
		--(nodes[nodes[num_nodes].parent].num_children);
		++(rest_atoms[nodes[num_nodes].label]); 
	}
	inline label_value_type begin_atom_label(const is_ident_type& is_ident, const int i) const;
	inline label_value_type begin_atom_label_root() const
	{
		const int nc = num_nodes - 1;
		if (nc > 0) {
			return nodes[nc].label;
		}
		return 0;
	}
	inline valence_value_type max_multi(const is_ident_type& is_ident, const int i) const;
	inline int is_normal(const int deepest_head) const;
	inline bool is_multi_normal(const int v) const;
	inline bool can_be_added(const int i) const
	{
		return (valence[nodes[i].label] - nodes[i].num_children - 1 > 0);
	}
	inline bool can_be_added_root() const
	{
		return (valence[nodes[0].label] > nodes[0].num_children);
	}

	inline bool remain(const label_value_type atom_label) const
	{
		return (rest_atoms[atom_label] > 0);
	}
	inline int starti() const
	{
		return nodes[num_nodes-1].parent;
	}
	inline int get_num_nodes() const
	{
		return num_nodes;
	}
  
        inline int get_num_child(int index){
                return nodes[index].num_children;
        } 
  inline void reset_bond_position(size_t index){
    for(size_t i = 0; i<max_valence; i++)
      nodes[index].bond_position[i] = 0;
  }
	inline bool share_only_root(int i) const // This function is available only during construction of simple trees
	{
		return (nodes[i].multi != nodes[num_nodes - 1].multi);
#if 0
int j = num_nodes - 1;
while (i != j) {
i = nodes[i].parent;
j = nodes[j].parent;
}
return (i == 0);
#endif
}
inline valence_value_type get_center(const int deepest_head) const // This function is available only during construction of simple trees
{
	const valence_value_type subtree = nodes[deepest_head].multi;
	if (subtree != nodes[num_nodes - 1].multi) {
		return 0;
	} else {
		return subtree;
	}
#if 0
int i = deepest_head;
int j = num_nodes - 1;
while (i != j) {
i = nodes[i].parent;
j = nodes[j].parent;
}
if (i == 0) {
return 0;
}
int previ;
do {
previ = i;
i = nodes[i].parent;
} while (i != 0);

return previ;
#endif
}
inline void set_multi_bond(const int i, const valence_value_type multiple) 
{
	nodes[i].multi = multiple;
}
inline void fill_rest_single_bond(const int i)
{
	for (int j = i; j < num_nodes; ++j) {
		nodes[j].multi = 1;
	}
}
inline void calc_lack_degree() const
{
	lack_degree.clear();
	lack_degree.push_back(valence[nodes[0].label] - nodes[0].num_children);
	for(int i = 1; i < num_nodes; ++i){
		lack_degree.push_back(valence[nodes[i].label] - nodes[i].num_children - 1);
	}
}
  inline int get_parent(const int i) const
  {
    return nodes[i].parent;
  }
  inline void print() const;
  inline void print_single() const;


  bool adjacent(int index1, int index2){
    if(nodes[index1].parent != index2 && nodes[index2].parent != index1)
      return false;
    return true;
  }
		bool can_be_added_multi(int i){
			if(i>0){
				//if current node is benzene and parent node is pyridine -> true
			  if(nodes[i].label==0 && nodes[get_parent(i)].label==0){//<num_special_atom){
			    if(get_parent(i)!=0 && nodes[get_parent(i)].multi==2)
			      //if its parent has merge bond it cannot have merge bond with its parent
			      return false;
			    for(size_t j=0; j<nodes[i].nth; j++){
			      //if one of its sibling has merger bond with parent -> it cannot have merge bond with parent
			      if(nodes[nodes[get_parent(i)].children[j]].multi==2)
				return false;
			    }
			    return true;
			  }
			  //current node is special atom and parent node is not benzene -> false
			  if(nodes[i].label<num_special_atom && nodes[get_parent(i)].label==0){
			    return true;
			  }
			  //if current nodes is not special atom and parent is special atom -> false
			  if(nodes[i].label>=num_special_atom && nodes[get_parent(i)].label>=num_special_atom){
			    return true;
			  }
			}
			return false;
		}
		int get_label(int i){
			return nodes[i].label;
		}
                int get_parent(int i){
		        return nodes[i].parent;
		}
		int get_multi(int i){
			return nodes[i].multi;
		}
		void assign_bond_position(int index,std::vector< std::vector<int> >&valid_pos,const std::vector< std::vector<int> >& c_index){
		  for(size_t j=0;j<valid_pos.size();j++){
				for(size_t k=0;k<valid_pos[j].size();k++){
					if(c_index[j][k]!=-1){
						nodes[index].bond_position[c_index[j][k]] = valid_pos[j][k];
					}
				}
			}

		}
		void clear_bond_position(int i){
			for(size_t j=0;j<nodes[i].num_children;j++)
				nodes[i].bond_position[j] = -1;
		}

		bool is_equal(int,int);
		void group_child(std::vector< std::vector<int> >&,int);
		
		
                size_t label_benzene(const int index,  std::vector<std::vector<int> > &chain_collection, std::vector<carbon_position> &result, std::ofstream & output_file);
                size_t label_child_and_next(std::vector<carbon_position> &result,int j,int k,int carbon_pos,int max_child,int index,const carbon_position &c_index, std::vector< std::vector<int> > & chain_collection,const std::vector< std::pair< std::vector<int>,std::vector<int> > > &trisym_chain_collection,int flip_mode,bool & normal,  std::ofstream & output_file);
                bool is_normal_chain(const std::vector< std::vector<int> > & position,int j,int k, int index,const carbon_position &c_index,const std::vector<  std::vector<int> > & chain_collection,const std::vector< std::pair<std::vector<int>,std::vector<int> > > & trisym_chain_collection, int flip_mode);
		void find_chain(std::vector< std::vector<int> > &chain);
                void find_chain_down(int step,int index_step,std::vector<int> &temp,std::vector< std::vector<int> > &temp2);
                bool chain_middle_node_symmetry(std::vector<int> chain, carbon_position & carbon_pos);
                int chain_middle_pair_symmetry(std::vector<int> chain, carbon_position & carbon_pos, int previous_unequal);
                int chain_symmetry(int index1,int index2);
                bool is_tsub_equal(int index1,int index2);
                bool is_tsub_cp_equal(int index1, int index2, int, int, const std::vector<int> &chain);
                bool is_tsub_has_benzene(int index);
                bool is_most_middle_bnode(const std::vector<int> &chain, const int index);
                int chain_end_check(std::vector<int> chain, const carbon_position &carbon_pos,const carbon_position &c_index, const int index);
                bool is_symmetry(int index);//find if this node is symmetry or not (all its child nodes must be assigned carbon_pos), use before is_normal_benzene
		bool is_trisymmetry(std::vector<int> chain1,std::vector<int> chain2);
                bool is_normal_benzene(const std::vector< std::vector<int> > &position,int j,int k,int carbon_pos,int index,int flip_mode,bool & normal);
                bool normal_naphthalene(int index);
		bool is_trisymmetry_redundant(int index_end1,int index_end2,const carbon_position & temp_pos);
		bool is_updown_symmetry(int index,int dealing_child,int first_unassigned_nth);
		bool is_fused_benzene(std::vector<int> chain){
			for(int i=0;i<chain.size()-1;i++){
				if(this->nodes[chain[i]].parent==chain[i+1]){
					if(this->nodes[chain[i]].multi!=2)
						return false;
				}
				if(this->nodes[chain[i+1]].parent==chain[i]){
					if(this->nodes[chain[i+1]].multi!=2)
						return false;
				}
			}
			return true;
		}

		void get_bond_position(carbon_position &cp_child,int index){
			for(size_t i=0; i<cp_child.size(); i++){
				for(size_t j=0; j<cp_child[i].size(); j++){
					if(cp_child[i][j]==-1){
						cp_child[i][j] = cp_child[i][j-1]+1;
					}else{
						cp_child[i][j] = this->nodes[index].bond_position[cp_child[i][j]];
					}
				}
			}
		}
                void get_bond_position(std::vector<int> &cp_child, int index){
			for(size_t i=0; i<cp_child.size(); i++){
			  cp_child[i] = this->nodes[index].bond_position[cp_child[i]];
			}
		}

                inline void write_smiles(const int index, std::ofstream & file, int & num_cycle) const;

		void show( bool show_position) const{
		  printtree(0, show_position);
		}

                    void printtree(int index, bool show_pos) const{
			Node n = nodes[index];
       
			int temp=index;
			while(temp!=0){
				std::cout << "        ";
				temp = nodes[temp].parent;
			}
			std::cout<<" label =" <<atomchar[n.label] <<" bond = "<<n.multi;//<<(n.is_identical?" (identical) ":" ");
			std::cout<<" parent = "<<n.parent;
			if(show_pos &&  index>0 && nodes[n.parent].label == 0 && nodes[n.parent].bond_position.size()>0){
				std::cout<<" c_position ="<<nodes[n.parent].bond_position[n.nth];
				if(n.multi ==2)
					std::cout<<","<<nodes[n.parent].bond_position[n.nth]+1;
			}
			std::cout<<std::endl;

			if(index == this->num_nodes-1){
				return;
			}
			if(n.num_children!=0){
			  printtree(n.children[0], show_pos);
			}

			if(nodes[index+1].parent == n.parent){
			  printtree(index+1, show_pos);   
			}

		}
	
	};

	inline void ChemTreeCenter::update_identical(is_ident_type& is_ident, const int i) const{
		using namespace std;

		location.clear();
		int j = i;
		int last_ident = -1;
		while (is_ident[j] >= (char)0) {
			if (is_ident[j] > (char)0) {
				int k = j - 1; // because of BFS 
				BOOST_REVERSE_FOREACH(const auto& loc, location) {
					k = nodes[k].children[loc];
				}
				if (nodes[k].label != nodes[i].label) {
					is_ident[j] = (char)0;
				} else {
					last_ident = j;
				}
			}
			location.push_back(nodes[j].nth);
			j = nodes[j].parent;
		}

		if (last_ident < 0) {
			last_ident = i;
			is_ident[last_ident] = (char)(-1);
		}
		last_ident = nodes[last_ident].parent;
		while (is_ident[last_ident] == (char)0 && last_ident != -1) {
			is_ident[last_ident] = (char)(-1);
			last_ident = nodes[last_ident].parent;
		}

	};



	inline void ChemTreeCenter::update_identical_end(is_ident_type& is_ident, const int i) const
	{
		using namespace std;
		location.clear();
		int j = i;
		int last_ident = -1;
		while (is_ident[j] >= (char)0) {
			if (is_ident[j] > (char)0) {
				int k = j - 1; 
				BOOST_REVERSE_FOREACH(const valence_value_type loc, location) {
					k = nodes[k].children[loc];
				}
				if(nodes[k].num_children != nodes[i].num_children) {
					is_ident[j] = (char)0;
				} else {
					last_ident = j;
				}
			}

			location.push_back(nodes[j].nth);
			j = nodes[j].parent;
		}

		if (last_ident < 0) {
		  last_ident = i;
		  is_ident[i]=(char)(-1);
		}

		last_ident = nodes[last_ident].parent; 

		while (is_ident[last_ident] == (char)0 && last_ident != -1) {
			is_ident[last_ident] = (char)(-1);
			last_ident = nodes[last_ident].parent;
		}
		
	}

	inline bool ChemTreeCenter::add_root(const label_value_type atom_label)
	{ 
		auto& root = nodes[0];
		root.label = atom_label;
		root.num_children = 0;
		root.parent = -1;
		root.multi = 1;
		num_nodes = 1;
		--(rest_atoms[atom_label]); 
		return num_nodes == num_except_H;
	}

	inline bool ChemTreeCenter::add_root_child(is_ident_type& is_ident, const label_value_type atom_label)
	{
		auto& node = nodes[num_nodes];
		auto& parnode = nodes[0];
		node.label = atom_label;
		const valence_value_type nc = parnode.num_children;
		node.parent = 0;
		node.multi = num_nodes;
		node.nth = nc;
		node.num_children = 0;

		parnode.children[nc] = num_nodes;
		++(parnode.num_children);

		if ((nc > 0) and (nodes[nc].label == atom_label)) {
			is_ident[num_nodes] = (char)1;
		} else {
			is_ident[num_nodes] = (char)(-1);
		}

		--(rest_atoms[atom_label]);
		++num_nodes;
		return num_nodes == num_except_H;
	}

	inline bool ChemTreeCenter::add_node(is_ident_type& is_ident, const int parenti, const label_value_type atom_label)
	{
		auto& node = nodes[num_nodes];
		auto& parnode = nodes[parenti];
		node.label = atom_label;
		const valence_value_type nc = parnode.num_children;
		node.parent = parenti;
		if (parenti == 0) { // subtree id
			node.multi = num_nodes;
		} else {
			node.multi = parnode.multi;
		}
		node.nth = nc;
		node.num_children = 0;

		parnode.children[nc] = num_nodes;
		++(parnode.num_children);

		is_ident[num_nodes] = (char)(nc > 0);
		update_identical(is_ident, num_nodes);

		--(rest_atoms[atom_label]);
		++num_nodes;

		return num_nodes == num_except_H;
	}

	inline label_value_type ChemTreeCenter::begin_atom_label(const is_ident_type& is_ident, const int i) const
	{
		using namespace std;

		int begin_atom = 0;
		const valence_value_type nc = nodes[i].num_children;
		if (nc > 0) {
			begin_atom = nodes[num_nodes-1].label;
		}
		location.clear();
		location.push_back(nc);
		int j = i;
		while (is_ident[j] >= (char)0) {
			if (is_ident[j] > (char)0) {
				int k = j - 1; 

				bool flag = false;
				BOOST_REVERSE_FOREACH(const valence_value_type loc, location) {
					if(nodes[k].num_children > loc) {
						k = nodes[k].children[loc];
					} else {
						begin_atom = first_atom_valence_one;
						flag = true;
						break;
					}
				}
				if(flag) break;
				else{
					if (nodes[k].label > begin_atom) {
						begin_atom = nodes[k].label;
					}
				}
			}
			location.push_back(nodes[j].nth);
			j = nodes[j].parent;
		}

		return begin_atom;
	}

	inline int ChemTreeCenter::is_normal(const int deepest_head) const
	{
		using namespace std;

		const valence_value_type v = get_center(deepest_head);

		if (v == 0) {
			return 1;
		}

		{ // root
			const label_value_type labelrv = nodes[0].label - nodes[v].label;
			if (labelrv != 0) {
				return (int)(labelrv > 0);
			}
		}

// children of root
		t_r.clear();
		t_v.clear();
		valence_value_type ncv = nodes[v].num_children + 1;
		const valence_value_type minv = (v < ncv) ? v : ncv;
		for (int i = 1; i < minv; ++i) {
			const int vc = nodes[v].children[i-1];
			const label_value_type labelrv = nodes[i].label - nodes[vc].label;
			if (labelrv != 0) {
				return (int)(labelrv > 0);
			}
			t_r.push_back(i);
			t_v.push_back(vc);
		}
		const valence_value_type nc0 = nodes[0].num_children;
		const valence_value_type minnc = (nc0 < ncv) ? nc0 : ncv;
		for (int i = v + 1; i <= minnc; ++i) {
			const int vc = nodes[v].children[i-2];
			const label_value_type labelrv = nodes[i].label - nodes[vc].label;
			if (labelrv != 0) {
				return (int)(labelrv > 0);
			}
			t_r.push_back(i);
			t_v.push_back(vc);
		}
		ncv -= nc0;
		if (ncv != 0) {
			return (int)(ncv > 0);
		}

		size_t j = 0;
		while (j < t_r.size()) {
			const auto& nr = nodes[t_r[j]];
			const auto& nv = nodes[t_v[j]];
			const valence_value_type ch_r = nr.num_children;
			const valence_value_type ch_v = nv.num_children;
			const valence_value_type minnc = (ch_r < ch_v) ? ch_r : ch_v;
			for(int i = 0;i < minnc;i++){
				const int nri = nr.children[i];
				const int nvi = nv.children[i];
				const label_value_type labelrv = nodes[nri].label - nodes[nvi].label;
				if (labelrv != 0) {
					return (int)(labelrv > 0);
				}
				t_r.push_back(nri);
				t_v.push_back(nvi);
			}
			if (ch_r != ch_v) {
				return (int)(ch_r < ch_v);
			}
			j++;
		}
		return 2+v;
	}

	inline bool ChemTreeCenter::is_multi_normal(const int v) const
	{
		using namespace std;

		t_r.clear();
		t_v.clear();
//cerr << v << endl;
		const int r = 0;
		for(int i = 1; i < v; ++i) {
			t_r.push_back(i);
		}
		for (int i = v + 1; i <= nodes[r].num_children; ++i){
			t_r.push_back(i);
		}
		for(int i =0;i<nodes[v].num_children;++i){
			t_v.push_back(nodes[v].children[i]);
		}
		size_t j = 0;
		while (j < t_r.size()) {
			const auto& nr = nodes[t_r[j]];
			const auto& nv = nodes[t_v[j]];
			const valence_value_type mvr = nv.multi - nr.multi;
			if (mvr != 0) {
				return (mvr > 0);
			} 
			for(int i = 0; i < nr.num_children; ++i){
				t_r.push_back(nr.children[i]);
				t_v.push_back(nv.children[i]);
			}
			++j;
		}

		return true;
	}

	inline void ChemTreeCenter::update_identical_multi(is_ident_type& is_ident, const int i) const
	{
		using namespace std;

		location.clear();
		int j = i;
		int last_ident = -1;
		while (is_ident[j] >= (char)0) {
			if (is_ident[j] > (char)0) {
				int k = j - 1; 
				BOOST_REVERSE_FOREACH(const valence_value_type loc, location) {
					k = nodes[k].children[loc];
				}
				if(nodes[k].multi != nodes[i].multi) {
					is_ident[j]= (char)0;
				} else {
					last_ident = j;
				}
			}
			location.push_back(nodes[j].nth);
			j = nodes[j].parent;
		}

		if (last_ident < 0) {
			last_ident = i;
			is_ident[last_ident] = (char)(-1);
		}
		last_ident = nodes[last_ident].parent;
		while (is_ident[last_ident] == (char)0) {
			is_ident[last_ident] = (char)(-1);
			last_ident = nodes[last_ident].parent;
		}
	}

	inline valence_value_type ChemTreeCenter::max_multi(const is_ident_type& is_ident, const int i) const
	{
		valence_value_type maxmulti = max_valence;
		location.clear();
		int j = i;
		while (is_ident[j] >= (char)0) {
			if (is_ident[j] > (char)0) {
				int k = j - 1; 
				BOOST_REVERSE_FOREACH(const valence_value_type loc, location) {
					k = nodes[k].children[loc];
				}
				const valence_value_type kmulti = nodes[k].multi;
				if(kmulti < maxmulti)
				{
					maxmulti = kmulti;
					if (maxmulti == 1) break;
				}
			}

			location.push_back(nodes[j].nth);
			j = nodes[j].parent;
		}
		return maxmulti;
	}

void count_ring_benzene(int index_b,int index_c,std::vector<int> &num_bond,std::vector<int> &num_input,std::vector< std::vector<int> > &input,std::vector< std::vector<int> > &bond){
  
	if(num_input[index_c]>=4 && num_bond[0]>=6){
		std::vector<int> temp_input(num_input);
		std::vector<int> temp_bond(num_bond);
		int index_p =  std::find(atomchar, atomchar + num_distinct_atoms , "p") - atomchar;
		//add benzene ring with double bond to benzene ring (naphthalene)
		if((num_input[index_b]>0||num_input[index_p]>0) && num_bond[1]*2<num_input[index_b]){
		  //check num_bond and num_input to limit that only two benzenes can have merge bond
			temp_input[index_b] = num_input[index_b]+1;
			temp_input[index_c] = num_input[index_c]-4;

			temp_bond[0]-=6;
			temp_bond[1]++;

			input.push_back(temp_input);
			bond.push_back(temp_bond);
			count_ring_benzene(index_b,index_c,temp_bond,temp_input,input,bond);
			}
		if(num_input[index_c]>=6 && num_bond[0]>=8 && num_bond[1]==0){
		  //add another benzene ring 
			temp_input = num_input;
			temp_input[index_b] = num_input[index_b]+1;
			temp_input[index_c] = num_input[index_c]-6;
			
			temp_bond = num_bond;
			temp_bond[0]-=8;

			input.push_back(temp_input);
			bond.push_back(temp_bond);
			count_ring_benzene(index_b,index_c,temp_bond,temp_input,input,bond);
		}
	}

}

void ChemTreeCenter::group_child(std::vector< std::vector<int> > &c_index,int i){
	std::vector<int> temp;
	for(size_t j=0;j<nodes[i].num_children;j++){
	  temp.push_back(static_cast<int>(j));
	}

	int num_skip=0;
	for(size_t j=0;j<temp.size();j++){
	  num_skip=0;
	  std::vector<int> a(1,temp[j+num_skip]);
	  if(nodes[nodes[i].children[temp[j+num_skip]]].multi==2){
	    a.push_back(-1);
	  }
	  c_index.push_back(a);

	  for(size_t k=j+1;k<temp.size();k++){
	    if(this->is_equal(nodes[i].children[temp[j]],nodes[i].children[temp[k]])){
	      c_index[j].push_back(temp[k]);
	      num_skip++;
	      if(nodes[nodes[i].children[k]].multi==2){
		c_index[j].push_back(-1);
	      }
	      temp.erase(temp.begin()+k);
	      k--;
	    }else{
	      break;
	    }
	  }  
	}

}

bool ChemTreeCenter::is_equal(int i,int j){
	if(nodes[i].label!=nodes[j].label || nodes[i].multi!=nodes[j].multi || nodes[i].num_children!=nodes[j].num_children)
		return false;
	if(nodes[i].num_children==0)
		return true;
	if(nodes[i].label< num_special_atom && nodes[i].bond_position!=nodes[j].bond_position)
	      return false;

	//bool result = is_equal(nodes[i].children[0],nodes[j].children[0]); 
	for(size_t index=0;index<nodes[i].num_children;index++){
	  if(!is_equal(nodes[i].children[index],nodes[j].children[index])){
	    return false;
	  }
	  //result = result & is_equal(nodes[i].children[index],nodes[j].children[index]);
	}
	return true;
	//return result;
}

int max(std::vector<int> container){
	int result=container[0];
	for(int i=1;i<container.size();i++){
		if(container[i]>result){
			result = container[i];
		}
	}
	if(result<0)
	   result = 0;
	return result;
}

bool found(std::vector<int> container,int val){
  //return T if val is found in container, otherwise return F
  return std::find(container.begin(),container.end(),val)!=container.end();
}

int flip_carbon(int position,int label,const std::vector<int> &smallest_set,int type){
	switch(type){
		case 0:
		  switch(label){
		  case 0:
		    return (6-position)%6;		    
		  case 1:
		    if(position<5){
		      return 4-position;
		    }else{
		      std::cout<<"error carbon position of N in pyridine"<<std::endl;
		      return 5;
		    }
		  }
		break;
	        case 1://is_normal_benzene2
		  switch(smallest_set.size()){
		  case 2:
		    //switch(smallest_set[0]){
		      //case 0:
		      //if(position<=smallest_set[1]){
		    //return (smallest_set[0]+smallest_set[1]-position)%6;
		    //}else{
		    return (6+smallest_set[0]+smallest_set[1]-position)%6;
			//}
		    //break;
		      //}
		    break;
		    /*case 3:
		    if(smallest_set[0]==0){ //has no parent
		      switch(smallest_set[1]){
		      case 1:
			if(position<2){
			  return 1-position;
			}else{
			  return 7-position;
			}
			break;
		      case 2:
			if(smallest_set[2]!=3 && position!=3){				       
			  if(position<3){
			    return 2-position;
			  }else{
			    return 8-position;
			  }
			}else{
			  if(position<4){
			    return 3-position;
			  }else{
			    return 9-position;
			  }
			}
			break;
		      }
		    }else{
		      //if it has parent, leave the checking for redundancy to is_normal_benzene1 and is_normal_benzene 3
		      return position+1; 
		      }
		  case 4:
		    if(position==smallest_set[0]+5 || smallest_set[1]==smallest_set[0]+2)
		      return -1; //always eliminate the position
		    break;
		  case 5:
		    if(position==smallest_set[0]+5)
		      return -1; //always eliminate the position
		      break;*/
		  }
	        case 2: //is_normal_benzene3
		  return 5-position;
		  break;
			 	
		
	}
	return position+1; //always true
} 

bool ChemTreeCenter::is_trisymmetry(std::vector<int> chain1,std::vector<int> chain2){
  if(!is_fused_benzene(chain1) || !is_fused_benzene(chain2)){
    return false;
  }
  if(chain1.size()%2==0){
    return false;
  }

  int center_node = (chain1.size())/2;
  
  if(chain1[center_node]!=chain2[center_node]){
    return false;
  }
  if(chain1[center_node]>chain1[center_node+1] || chain2[center_node]>chain2[center_node+1]){
    return false;
  }

  for(int i=1;i<=(chain1.size()/2)-1;i++){
    Node parent_benzene[3] = {this->nodes[chain1[center_node-i]],this->nodes[chain1[center_node+i]],this->nodes[chain2[center_node+i]]};
    Node child_benzene[3] = {this->nodes[chain1[center_node-i-1]],this->nodes[chain1[center_node+1+i]],this->nodes[chain2[center_node+i+1]]};
    //check for carbon position of each node in the chain, if they are not the same one -> no tri radical symmetry
    if(parent_benzene[0].bond_position[child_benzene[0].nth]!=parent_benzene[1].bond_position[child_benzene[1].nth])
      return false;
    if(parent_benzene[0].bond_position[child_benzene[0].nth]!=parent_benzene[2].bond_position[child_benzene[2].nth])
      return false;
  }
  return true;
}

void flip_trisymmetry(carbon_position (&flip_cp)[3],int i,int range){
  for(int k=0;k<3;k++){
    for(int j=0;j<range;j++){
      flip_cp[k][i][j] = 5-flip_cp[k][i][j];
    } 
    std::sort(flip_cp[k][i].begin(),flip_cp[k][i].begin()+range);
  }
}

void rotate_trisymmetry(carbon_position (&flip_cp)[3],int i){
  std::vector<int> temp = flip_cp[0][i];
  flip_cp[0][i] = flip_cp[1][i];
  flip_cp[1][i] = flip_cp[2][i];
  flip_cp[2][i] = temp;
}

bool greater_than(const carbon_position (&original_cp)[3],const carbon_position (&flip_cp)[3],int i,int range){
  for(int k=0;k<3;k++){
    for(int j=0;j<range;j++){
      if(original_cp[k][i][j]>flip_cp[k][i][j])
	return true;
      if(original_cp[k][i][j]<flip_cp[k][i][j])
	return false;
    }
  }
  return false;
}

bool ChemTreeCenter::is_trisymmetry_redundant(int index_end1,int index_end2,const carbon_position & temp_pos){
  carbon_position cp_end1,cp_end2;
  group_child(cp_end1,index_end1);
  std::sort(cp_end1.begin(),cp_end1.end(),compare_function);
  group_child(cp_end2,index_end2);
  std::sort(cp_end2.begin(),cp_end2.end(),compare_function);

  get_bond_position(cp_end1,index_end1);
  get_bond_position(cp_end2,index_end2);

  carbon_position original_cp[3] = {cp_end1,cp_end2,temp_pos};
  carbon_position flip_cp[3] = {cp_end1,cp_end2,temp_pos};

  for(int i=0;i<temp_pos.size();i++){
    int range = std::find(temp_pos[i].begin(),temp_pos[i].end(),-2)-temp_pos[i].begin();
    if(range >0){
      rotate_trisymmetry(flip_cp,i);
      if(greater_than(original_cp,flip_cp,i,range)){
	return true;
      }
      rotate_trisymmetry(flip_cp,i);
      if(greater_than(original_cp,flip_cp,i,range)){
	return true;
      }
      
      rotate_trisymmetry(flip_cp,i); //with third rotation, flip_cp is back to original_cp 
      flip_trisymmetry(flip_cp,i,range);
      if(greater_than(original_cp,flip_cp,i,range)){
	return true;
      }
      rotate_trisymmetry(flip_cp,i);
      if(greater_than(original_cp,flip_cp,i,range)){
	return true;
      }
      rotate_trisymmetry(flip_cp,i);
      if(greater_than(original_cp,flip_cp,i,range)){
	return true;
      }
    }else{
      break;
    }
  }  
  return false;
};

bool ChemTreeCenter::is_updown_symmetry(int index,int dealing_child,int first_unassigned_nth){
  //dealing_child is index (not nth) of child that is filling carbon_position, so we will not consider it in finding symmetry of next benzene

  carbon_position c_index;
  group_child(c_index,index);
  if(index!=0 && dealing_child>index){
    //if 'index' has parent node, add parent into c_index in the form of '-2'
    std::vector<int> temp(1,-2);
    if(this->nodes[index].multi==2){
      temp.push_back(-1);
    }
    c_index.push_back(temp);
  }
  std::sort(c_index.begin(),c_index.end(),compare_function);
  
  if(this->get_label(dealing_child)!=0 || this->nodes[dealing_child].multi==1){
    if(dealing_child > index){
      return true;
    }
  }

  int cp_dealing_child1 = 0;
  int cp_dealing_child2 = 5;
  
  if(dealing_child>index){
    cp_dealing_child1 = this->nodes[index].bond_position[this->nodes[dealing_child].nth];
    cp_dealing_child2 = cp_dealing_child1+1;
  }

  carbon_position cp(c_index),flip_cp(c_index);
  for(int i=0;i<c_index.size();i++){
    for(int j=0;j<c_index[i].size();j++){
      if(c_index[i][j]==-2){
	cp[i][j] = 0;
      }else if(c_index[i][j]==-1){
	if(c_index[i][j-1]!=-2){ //c_index[i][j-1] is not parent
	  cp[i][j] = this->nodes[index].bond_position[c_index[i][j-1]]+1;
	}else{
	  cp[i][j] = 5;
	}
      }else{
	cp[i][j] = this->nodes[index].bond_position[c_index[i][j]];
      }
      flip_cp[i][j] = (6+cp_dealing_child1+cp_dealing_child2-cp[i][j])%6;
    }
    std::sort(cp[i].begin(),cp[i].end());
    std::sort(flip_cp[i].begin(),flip_cp[i].end());
  }
  
  if(cp!=flip_cp){
    return false;
  }else{
    //check if benzene node 'b' bonding with current node is symmetry or not
    //if b is not symmetry, then current node is not symmetry
    for(int i=0;i<c_index.size();i++){
      for(int j=0;j<c_index[i].size();j++){
	/*if(c_index[i][j]==-2){//c_index is parent
	  int index_of_parent =this->nodes[index].parent;
	  if(!is_updown_symmetry(index_of_parent,index,this->nodes[index_of_parent].num_children)){
	    return false;
	  }
	  }else */
	if(c_index[i][j]!=-1 && c_index[i][j]<first_unassigned_nth && c_index[i][j]!=this->nodes[dealing_child].nth){//c_index is child that is not reserve carbon for naph_bond
	  int index_of_child = this->nodes[index].children[c_index[i][j]];
	  if(this->get_label(index_of_child)==0 && this->nodes[index_of_child].multi==2){
	    if(!is_updown_symmetry(index_of_child,index,this->nodes[index_of_child].num_children)){
	      return false;
	    }
	  }
	}
      }
    }
    return true;
  }
}

int ChemTreeCenter::chain_symmetry(int index1,int index2){
  //return  1 if cp[index1] < cp[index2] , where index1 < index2   
  //return  0 if cp[index1] = cp[index2] , where index1 < index2   
  //return -1 if cp[index1] > cp[index2] , where index1 < index2   

  if(index1 > index2){
    int temp = index1;
    index1 = index2;
    index2 = temp;
  }

  if(this->get_label(index1)!=this->get_label(index2)){
    std::cout<<"Error chain_symmetry (1)"<<std::endl;
    return -2;
  }

  if(this->get_label(index1)!=0){
    return 0;//true;
  }
  /*
  std::vector<int> children1,children2;
  for(int i=0;i<this->nodes[index1].num_children;i++){
    children1.push_back(this->nodes[index1].children[i]);
  }
  for(int i=0;i<this->nodes[index2].num_children;i++){
    children2.push_back(this->nodes[index2].children[i]);
  }
  */
  //bool is_adjacent = false;
  //int index_parent,index_child;
  
  //remove index1/index2 if it is another node's child (consider other adjacent nodes only)
  if(this->nodes[index1].parent==index2){
    std::cout<<"Error nodes cannot be adjacent"<<std::endl;
    /*is_adjacent = true;
    index_parent = index2;
    index_child = index1;
    children2.erase(children2.begin()+this->nodes[index1].nth);*/
  }else if(this->nodes[index2].parent==index1){
    std::cout<<"Error nodes cannot be adjacent"<<std::endl;
    /*is_adjacent = true;
    index_parent = index1;
    index_child = index2;
    children1.erase(children1.begin()+this->nodes[index2].nth);*/
  }

  //check if the number of adjacent node is equal or not 
  //in the case that one of them is parent node of another one, because the parent node must be the root so it cannot have parent node 
  if(nodes[index1].num_children != nodes[index2].num_children ){
    std::cout<<"Error chain_symmetry (2)"<<std::endl;
    return -2;//false;
  }

  //check if all adjacent nodes have the same structure or not
  //for(int i=0;i<children1.size();i++){
  //  if(!is_equal(children1[i],children2[i])){
  //    std::cout<<"Error chain_symmetry (3)"<<std::endl;
  //    return -2;//false;
  //  }
  //}

  if(this->nodes[index1].multi!=this->nodes[index2].multi){
    std::cout<<"Error chain_symmetry (4)"<<std::endl;
    return -2;//false;
  }

  carbon_position cp_index1,cp_index2;

  group_child(cp_index1, index1);
  group_child(cp_index2, index2);
  std::sort(cp_index1.begin(), cp_index1.end(), compare_function);
  std::sort(cp_index2.begin(), cp_index2.end(), compare_function);
  get_bond_position(cp_index1, index1);
  get_bond_position(cp_index2, index2);  
  

  if(this->get_label(index1)==0){
    if(cp_index1 < cp_index2){
      return 1;
    }else if(cp_index1 == cp_index2){
      return 0;
    }else{
      return -1;
    } 
    }else{
      return -1;
  }
  //}else{
  //if index1 and index2 are adjacent and not benzene -> always symetry (no need to check for multi because they bond with each other) 
  //if(this->get_label(index1)!=0)
  //  return -1;//true;
  //}
  
  //parent must be root unless the tree is not left heavy
  //if(index_parent!=0){
  //  return false;
  //}
  
  //if it is not the first child of this label_atom(ex.carbon), it and its parent cannot have the same structure
  //for(int i=this->nodes[index_child].nth-1;i>=0;i--){
  //  if(this->get_label(this->nodes[index_parent].children[i])==this->get_label(index_child)){
  //    return false;
  //  }
  //}
  
  //std::array<int,max_valence> cp_child=this->nodes[index_child].bond_position;
  //std::array<int,max_valence> cp_parent=this->nodes[index_parent].bond_position;
  /*
  carbon_position cp_child,cp_parent;
  group_child(cp_child,index_child);
  group_child(cp_parent,index_parent);
  //delete information of child benzene from cp_parent
  for(int i=0;i<cp_parent.size();i++){
    for(int j=0;j<cp_parent[i].size();j++){
      if(cp_parent[i][j]==this->nodes[index_child].nth){
	if(j==0 && cp_parent[i].size()==this->nodes[index_child].multi){
	  cp_parent.erase(cp_parent.begin()+i);
	}else{
	  cp_parent[i].erase(cp_parent[i].begin()+j);
	  if(this->nodes[index_child].multi==2){
	    cp_parent[i].erase(cp_parent[i].begin()+j+1);
	  }
	}
	break;
      }
    }
  }
  std::sort(cp_child.begin(),cp_child.end(),compare_function);
  std::sort(cp_parent.begin(),cp_parent.end(),compare_function);
  get_bond_position(cp_child,index_child);
  get_bond_position(cp_parent,index_parent);  
  
  if(this->nodes[index_child].multi==1){
    std::cout<<" index1 = "<<index1<<"   index2 = "<<index2<<std::endl;
    int nth = this->nodes[index_child].nth;
    //compare carbon position of child nodes (excluding index1 in cp_parent)
    //parent must be root so can compare by cp_child==cp_parent
    
    if(cp_child.size()!=cp_parent.size()){
      std::cout<<"error!!"<<std::endl;
    }else{
      if(cp_child < cp_parent){
	if(index_child == index1){
	  std::cout<<" return 1"<<std::endl;
	  return 1;//false;
	}else{
	  return -1;
	}
      }
      if(cp_child > cp_parent){
	if(index_child == index1){
	  std::cout<<" return -1"<<std::endl;
	  return -1;
	}else{
	  return 1;
	}
      }
    }
    return 0;//true;
  }else if(this->nodes[index1].multi==2){
    //used to shift carbon position of child (in parent bond_position) to 0
    //not consider the case that more than 2 benzenes merge together ================================
    //if consider that case, cp of naph bond must be included (now only one cp per one atom from cp_child/cp_parent)
    int cp_shift = this->nodes[index_parent].bond_position[this->nodes[index_child].nth];
    
    //or(int i=0;i<this->nodes[index_parent].num_children;i++){
      //cp_parent[i] = (cp_shift-cp_parent[i]+6)%6;
      //}
    for(int i=0;i<cp_parent.size();i++){
      for(int j=0;j<cp_parent[i].size();j++){
	cp_parent[i][j] = (cp_shift-cp_parent[i][j]+6)%6;
      }
      std::sort(cp_parent[i].begin(),cp_parent[i].end());
    }
    
    bool leftright_sym = true;
    //left-right symmetry between 2 benzenes
    
    //for(int i=0;i<this->nodes[index_child].num_children;i++){
    //  int shift=0;
    //  if(i>=this->nodes[index_child].nth){
	//shift=1;
      //}
      //if(cp_child[i]!=cp_parent[i+shift]){
	//leftright_sym = false;
      //}
      //}
    if(cp_child.size()!=cp_parent.size()){
      std::cout<<"error!!"<<std::endl;
    }else{
      if(cp_child!=cp_parent){
	leftright_sym = false;
      }
    }
    if(leftright_sym){
      return true;
    }

    //left-right and up-down symmetry between 2 benzenes
    
    //for(int i=0;i<this->nodes[index_child].num_children;i++){
    //  int shift=0;
    //  if(i>=this->nodes[index_child].nth){
	//shift=1;
      //}
      //if(cp_child[i]!=5-cp_parent[i+shift]){
	//return false;
      //}
      //}
    for(int i=0;i<cp_parent.size();i++){
      for(int j=0;j<cp_parent[i].size();j++){
	cp_parent[i][j] = 5-cp_parent[i][j];
      }
      std::sort(cp_parent[i].begin(),cp_parent[i].end());
    }

    if(cp_child.size()!=cp_parent.size()){
      std::cout<<"error!!"<<std::endl;
    }else{
      if(cp_child!=cp_parent){
	return false;
      }
    }
    return true;
    }*/
  std::cout<<"error in chain_symmetry case:1,2"<<std::endl;
  return -1;//true; 
}

//bool ChemTreeCenter::chain_end_symmetry(int index1,int index2){
bool ChemTreeCenter::is_tsub_equal(int index1,int index2){ 
  if(get_label(index1) != get_label(index2) )
    return false;

  std::vector<int> children1,children2;
  for(int i=0;i<this->nodes[index1].num_children;i++){
    children1.push_back(this->nodes[index1].children[i]);
  }
  
  for(int i=0;i<this->nodes[index2].num_children;i++){
    children2.push_back(this->nodes[index2].children[i]);
  }
  
  bool is_adjacent = false; 
  //is two end-nodes adjacent with each other or not

  //remove index1/index2 if it is another node's child (consider other adjacent nodes only)
  if(this->nodes[index1].parent==index2){
    is_adjacent = true;
    children2.erase(children2.begin()+this->nodes[index1].nth);
  }else if(this->nodes[index2].parent==index1){
    is_adjacent = true;
    children1.erase(children1.begin()+this->nodes[index2].nth);
  }
  
  if(!is_adjacent){
    if(this->nodes[index1].multi!=this->nodes[index2].multi){
      return false;
    }
    if(children1.size()!=children2.size()){
      return false;
    }
    if(nodes[index1].parent != nodes[index2].parent){

      if(!adjacent(nodes[index1].parent, nodes[index2].parent)){
	if(nodes[index1].nth != nodes[index2].nth){
	  return false;
	}
      }else{
	int nth;
	if(nodes[ nodes[index1].parent ].parent == nodes[index2].parent){
	  nth = nodes[ nodes[index1].parent ].nth;
	  int shift = (nth>nodes[index2].nth ? 0:-1);
	  if(nodes[index1].nth != nodes[index2].nth+shift){
	    return false;
	  }
	}else{
	  nth = nodes[ nodes[index2].parent ].nth;
	  int shift = (nth>nodes[index1].nth ? 0:-1);
	  if(nodes[index1].nth+shift != nodes[index2].nth){
	    return false;
	  }
	}
      }
    }
  }else{
    if(index2!=0){
      return false;
    }
    if(children1.size()!=children2.size()){
      return false;
    }
  }

  //check if all adjacent nodes have the same structure or not
  for(int i=0; i<children1.size(); i++){
    if(!is_equal(children1[i], children2[i])){
      return false;
    }
  }
  return true;
}

bool ChemTreeCenter::is_tsub_cp_equal(int index1, int index2, int chain_index1, int chain_index2, const std::vector<int> &chain){
   
  std::vector<int> child_index1,child_index2, cp_index1, cp_index2;
  for(int i=0;i<this->nodes[index1].num_children;i++){
    cp_index1.push_back(i);
    child_index1.push_back(this->nodes[index1].children[i]);
  }
  
  for(int i=0;i<this->nodes[index2].num_children;i++){
    cp_index2.push_back(i);
    child_index2.push_back(this->nodes[index2].children[i]);
  }

  //remove index1/index2 if it is another node's child (consider other adjacent nodes only)
  if(this->nodes[index1].parent == index2){
    child_index2.erase(child_index2.begin() + this->nodes[index1].nth);
    cp_index2.erase(cp_index2.begin() + this->nodes[index1].nth);
  }else if(this->nodes[index2].parent == index1){
    child_index1.erase(child_index1.begin() + this->nodes[index2].nth);
    cp_index1.erase(cp_index1.begin() + this->nodes[index2].nth);
  }
  
  if(get_label(index1)==0 && get_label(index2)==0){
    if(!found(chain, index1) && !found(chain, index2) ){
      get_bond_position(cp_index1, index1);
      get_bond_position(cp_index2, index2);
      
      if(cp_index1 != cp_index2){
	return false;
      }
    }
  }

  bool result = true;
  for(int i=0; i<child_index1.size(); i++){
    if(chain_index1==0 || (child_index1[i] > chain[chain_index1-1] && child_index2[i] > chain[chain_index2+1]) ){
      result = result && is_tsub_cp_equal(child_index1[i], child_index2[i], chain_index1-1, chain_index2+1, chain);
    }
  }
  return result;
}

bool ChemTreeCenter::is_tsub_has_benzene(int index){
  bool result = false;
  for(int i=0; i< nodes[index].num_children; i++){
    if(nodes[ nodes[index].children[i] ].num_children == 0){
      return false;
    }
    if(get_label(nodes[index].children[i]) == 0){
      return true;
    }
    result = result || is_tsub_has_benzene(nodes[index].children[i]);
  }
  return result;
}

//bool ChemTreeCenter::is_chain_end_redundant(std::vector<int> chain,const carbon_position &carbon_pos,const carbon_position &c_index, const int index){
int ChemTreeCenter::chain_end_check(std::vector<int> chain,const carbon_position &carbon_pos,const carbon_position &c_index, const int index){
  //get index of benzene chain and check if both end is redundant or not
  //return  1 if chain[index1] < chain[index2]
  //return  0 if chain[index1] = chain[index2]
  //return -1 if chain[index1] > chain[index2]

  int index1 = chain[0], index2 = chain[chain.size()-1];
  //int index1 = chain[chain.size()-1], index2 = chain[0];
  
  bool is_adjacent = false; 
  //is two end-nodes adjacent with each other or not

  carbon_position c_index1;
  group_child(c_index1, index1);
  carbon_position cp_index1;
  cp_index1 = c_index1;
  this->get_bond_position(cp_index1, index1);

  carbon_position c_index2,cp_index2;
  if(index2 == index){
    c_index2 = c_index;
    cp_index2 = carbon_pos;
  }else{
    group_child(c_index2, index2);
    this->get_bond_position(cp_index2, index2);
  }

  //remove index1/index2 if it is another node's child (consider other adjacent nodes only)
  if(this->nodes[index1].parent==index2 || this->nodes[index2].parent==index1){
    is_adjacent = true;
  }

  if(this->nodes[index1].multi==1 && this->nodes[index2].multi==1){
    int nth = this->nodes[index1].nth;
    int shift = 0;
    //std::cout<<"-------------------------"<<std::endl;
    //std::cout<<" cp of "<<index1<<" (index1) is ";
    //for(int i=0; i<cp_index1.size(); i++){    
    //  for(int j=0; j<cp_index1[i].size(); j++){
    //std::cout<<cp_index1[i][j]<<" ";
    //}
    //}
    //std::cout<<std::endl;
    //std::cout<<" cp of "<<index2<<" (index2) is ";
    //for(int i=0; i<cp_index2.size(); i++){    
    //  for(int j=0; j<cp_index2[i].size(); j++){
    //std::cout<<cp_index2[i][j]<<" ";
    //}
    //}
    //std::cout<<std::endl;

    for(int i=0; i<cp_index2.size(); i++){    
      for(int j=0; j<cp_index2[i].size(); j++){
	if(is_adjacent && c_index2[i][j] == nth){
	  //remove this element from c_index2
	  if(cp_index2[i].size() == 1){
	    cp_index2.erase(cp_index2.begin()+i);
	    c_index2.erase(c_index2.begin()+i);
	    j=-1;
	  }else{
	    cp_index2[i].erase(cp_index2[i].begin()+j);
	    c_index2.erase(c_index2.begin()+i);
	    j--;
	  }
	  continue;
	}

	if(cp_index2[i][j]!=-2){
	  if(is_adjacent && c_index2[i][j]>=nth){
	    shift = 0;
	  }    

	  int co_carbon_pos = cp_index1[i+shift][j];//this->nodes[index2].bond_position[c_index1[i][j]+shift];	 
	  /*std::cout<<"  index1="<<index1<<"  index2="<<index2<<std::endl;
	  std::cout<<" (cp index2) i="<<i<<"   j="<<j<<std::endl;
	  std::cout<<" (co-cp index1) i="<<i+shift<<"   j="<<j<<std::endl;
	  std::cout<<" cp = "<<carbon_pos[i][j]<<"    co_cp  ="<<co_carbon_pos<<"   shift ="<<shift<<std::endl;
	  */
	  if(cp_index2[i][j]>co_carbon_pos){
	    return 1;//true;
	  }else if(cp_index2[i][j]<co_carbon_pos){
	    return -1;//false;
	  }
	}else{
	  break;
	}
      }
    }
    return 0;//true;
    //return false;
  }else if(this->nodes[index1].multi==2){	  
    std::vector< std::vector<int> > new_carbon_pos;
    new_carbon_pos = carbon_pos;
   
    if(is_fused_benzene(chain)){
      //new_carbon_pos is rotate carbon position in reverse direction (1-2-3-4) to (4-3-2-1)
      //necessary only when two benzene is in the same cyclic structure
      for(int i=0;i<carbon_pos.size();i++){
	for(int j=0;j<carbon_pos[i].size();j++){
	  new_carbon_pos[i][j] = 5-carbon_pos[i][j]; 
	}
	std::sort(new_carbon_pos[i].begin(),new_carbon_pos[i].end());
      }
    }
	 
    bool check_newcarbon = true;
    bool check_carbon = true;
    bool leftright_updown_sym = false;
    
    if(is_adjacent){
      leftright_updown_sym = true;
    }else if(is_fused_benzene(chain)){
      int root_index; //position of the top node (node with lowest depth) in the chain
      for(int i=0;i<chain.size();i++){
	if(chain[i]<chain[i+1]){
	  root_index = i;
	  break;
	}
      }
      carbon_position temp_c_index;
      group_child(temp_c_index,chain[root_index]);
      std::sort(temp_c_index.begin(),temp_c_index.end(),compare_function);
      //check if chain is a straight chain of benzene rings or not and can be flipped in both left-right and up-down direction
      if(abs(this->nodes[chain[root_index]].bond_position[this->nodes[chain[root_index-1]].nth]-this->nodes[chain[root_index]].bond_position[this->nodes[chain[root_index+1]].nth] )==3 && temp_c_index[0].size()%2==0  && chain[root_index]==0){
	//root_index has even child and no parent
	  
	bool no_odd_child = true;
	for(int i=root_index-1;i>0;i--){
	  temp_c_index.clear();
	  group_child(temp_c_index,chain[i]);
	  if(this->nodes[chain[i]].bond_position[this->nodes[chain[i-1]].nth]!=2){
	    no_odd_child = false;	     
	    break;
	  }
	  for(int j=0;j<temp_c_index.size();j++){
	    if(temp_c_index[j].size()%2!=0){
	      no_odd_child=false;
	      break;
	    }
	  }
	}
	if(no_odd_child){
	  leftright_updown_sym = true;
	}
	/* check only first half because the first half and the latter half are the same (from chain_symmetry)
	if(no_odd_child){
	  for(int i=root_index+1;i<chain.size()-1;i++){
	    temp_c_index.clear();
	    group_child(temp_c_index,chain[i]);
	    if(temp_c_index[0].size()%2!=0 || this->nodes[chain[i]].bond_position[this->nodes[chain[i+1]].nth]!=2){
	      no_odd_child = false;
	      break;
	    }
	  }
	  if(no_odd_child){
	    leftright_updown_sym = true;
	  }
	  }*/
      }
    }

    carbon_position index2_c_index; //store nth of benzene's child
    
    if(is_adjacent){
      group_child(index2_c_index,index2);
      std::sort(index2_c_index.begin(),index2_c_index.end(),compare_function);

      int nth = this->nodes[index1].nth;
      int cp_shift = 6+this->nodes[index2].bond_position[nth];
      //delete information of index1 (child benzene) from naph_c_index
      //and convert carbon position of nodes adjacent to index2 to be between 1-4 based on cp_shift
      for(int i=0;i<index2_c_index.size();i++){
	for(int j=0;j<index2_c_index[i].size();j++){
	  if(index2_c_index[i][j]==this->nodes[index1].nth){
	    if(j==0 && index2_c_index[i].size()==2){
	      index2_c_index.erase(index2_c_index.begin()+i);
	      j=-1;//make j at new iteration =0
	    }else{
	      index2_c_index[i].erase(index2_c_index[i].begin()+j);
	      index2_c_index[i].erase(index2_c_index[i].begin()+j+1);
	      j--;
	    }
	  }else{
	    if(carbon_pos[i][j]!=-2){
	      if(index2_c_index[i][j]!=-1){
		if(this->nodes[index2].bond_position[nth]!=0 && this->nodes[index2].bond_position[nth]<3){
		  index2_c_index[i][j] = (cp_shift-this->nodes[index2].bond_position[index2_c_index[i][j]])%6;
		}else{
		  index2_c_index[i][j] = (this->nodes[index2].bond_position[index2_c_index[i][j]]+5-this->nodes[index2].bond_position[nth])%6;
		}
	      }else{
		if(this->nodes[index2].bond_position[nth]!=0 && this->nodes[index2].bond_position[nth]<3){
		  index2_c_index[i][j] = (cp_shift-(this->nodes[index2].bond_position[index2_c_index[i][j-1]]+1))%6;
		}else{
		  index2_c_index[i][j] = (this->nodes[index2].bond_position[index2_c_index[i][j-1]]+1+5-this->nodes[index2].bond_position[nth])%6;
		}
	      }
	    }else{
	      if(j!=0){
		return -1;//false;
	      }
	    }
	  }
	}
	std::sort(index2_c_index[i].begin(),index2_c_index[i].end());
      }
    }/*else{ //not adjacent
      for(int i=0;i<carbon_pos.size();i++){
	for(int j=0;j<carbon_pos[i].size();j++){
	  this->nodes[index2].bond_position[c_index[i][j]];
	}
	}
	}*/

    for(int i=0;i<new_carbon_pos.size();i++){    
      for(int j=0;j<new_carbon_pos[i].size();j++){
	if(carbon_pos[i][j]!=-2){
	  if(leftright_updown_sym){  
	    //their common parent does not have parent, has only even unique childs
	    if(check_newcarbon && check_carbon){
	      if(new_carbon_pos[i][j]>index2_c_index[i][j] && carbon_pos[i][j]>index2_c_index[i][j]){
		return -1;//false;
	      }
	      if(new_carbon_pos[i][j]<index2_c_index[i][j] || carbon_pos[i][j]<index2_c_index[i][j]){
		return 1;//true;
	      }
		
	      //when reach here either new_carbon_pos or carbon_pos (but not both) must equal to index2.bond_position[....]
	      if(new_carbon_pos[i][j]==index2_c_index[i][j]){
		//if new_carbon_pos is equal, we will compare only new carbon_pos in the next round 
		check_carbon = false;
	      }else{
		check_newcarbon = false;
	      }
	    }else if(check_carbon && !check_newcarbon){
	      if(carbon_pos[i][j]>index2_c_index[i][j]){
		return -1;//false;
	      }
	      if(carbon_pos[i][j]<index2_c_index[i][j]){
		return 1;//true;
	      } 
	    }else if(!check_carbon && check_newcarbon){
	      if(new_carbon_pos[i][j]>index2_c_index[i][j]){
		return -1;//false;
	      }
	      if(new_carbon_pos[i][j]<index2_c_index[i][j]){
		return 1;//true;
	      }
	    }
		  
	  }else{
	    if(new_carbon_pos[i][j]>index2_c_index[i][j]){
	      return -1;//false;
	    }
	    if(new_carbon_pos[i][j]<index2_c_index[i][j]){
	      return 1;//true;
	    }
	  }
    
	}//carbon_pos[i][j]!=-2	  
      }//for j of carbon_pos[i][j]
    }//for i of carbon_pos[i][j]
    return 0;//false;

  }//multi==2
    
  return 0;//false;
}
/*
bool ChemTreeCenter::chain_middle_node_symmetry(std::vector<int> chain, carbon_position & carbon_pos){
  //return true if chain[middle-1] and chain[middle+1] can be fliped with each other

	int middle = chain.size()/2;
	
	if(chain[middle]>chain[middle+1] || chain[middle]>chain[middle-1]){
	  //if middle of chain is not parent of its adjacent node in chain then it cannot be symmetry
	  return false;
	}
	if(this->nodes[chain[middle+1]].multi!=this->nodes[chain[middle-1]].multi){
	  return false;
	}
	if(this->get_label(chain[middle])!=0){
	  return true;
	}

	std::vector< std::vector<int> > c_index,flip_carbon_pos; //order of children having the same structure are grouped together 
	//std::vector< std::vector<int> > carbon_pos;
	group_child(c_index,chain[middle]);

	if(chain[middle]!=0){
	  std::vector<int> temp (1,-2);
	  std::vector<int> temp_cp (1, 0);
	  if(this->get_label(chain[middle])==0 && this->nodes[chain[middle]].multi==2){
	    temp.push_back(-1);
	    temp_cp.push_back(5);
	  }
	  c_index.insert(c_index.begin(), temp);
	  carbon_pos.insert(carbon_pos.begin(), temp_cp);
	}
	std::sort(c_index.begin(), c_index.end(), compare_function);
	std::sort(carbon_pos.begin(), carbon_pos.end(), compare_function);

	if(chain[middle]==0 && c_index.size()==1 && c_index[0].size()==2 ){
	  //if its adjacent node has only two benzene with single bond and not have parent, it is always symmetry
	  return true;
	}

	flip_carbon_pos = carbon_pos;
	//initialize carbon_position of corresponding child node
	//for(int i=0;i<c_index.size();i++){
	  //for(int j=0;j<c_index[i].size();j++){
	    //if(c_index[i][j]==-2){
	      //carbon_pos[i][j] = 0;
	    //}else if(c_index[i][j]!=-1){
	      //carbon_pos[i][j] = this->nodes[chain[middle]].bond_position[c_index[i][j]];
	    //}else{
	      //if(c_index[i][j-1]!=-2){
		//carbon_pos[i][j] = carbon_pos[i][j-1]+1;
	      //}else{
		//carbon_pos[i][j] = 5;
	      //}
	    //}
	  //}
	//}

	if(c_index[0].size()==1 && c_index[0][0]!=this->nodes[chain[middle+1]].nth && c_index[0][0]!=this->nodes[chain[middle-1]].nth){
	  for(size_t i=0; i < carbon_pos.size(); i++){
	    for(size_t j=0; j < carbon_pos[i].size(); j++){
	      flip_carbon_pos[i][j] = (6-carbon_pos[i][j])%6;
	    }
	    std::sort(flip_carbon_pos[i].begin(),flip_carbon_pos[i].end());
	  }
	  if(flip_carbon_pos == carbon_pos){
	    return true;
	  }else{ 
	    return false; 
	  }    
	}

	float c_pos1;//=this->nodes[chain[middle]].bond_position[this->nodes[chain[middle-1]].nth];
	float c_pos2;//=this->nodes[chain[middle]].bond_position[this->nodes[chain[middle+1]].nth];
	std::pair< bool, bool> set (false,false);

	for(size_t i=0; i<c_index.size(); i++){
	  for(size_t j=0; j<c_index[i].size(); j++){
	    if(c_index[i][j] == nodes[chain[middle-1]].nth){
	      c_pos1 = carbon_pos[i][j];
	      set.first = true;
	    }else if(c_index[i][j] == nodes[chain[middle+1]].nth){
	      c_pos2 = carbon_pos[i][j];
	      set.second = true;
	    }
	  }
	}
	if(!set.first || !set.second){
	  std::cout<<"Error: carbon position cannot be found properly"<<std::endl;
	}

	//check symmetry of benzene based on position of two adjacent benzene
	switch(abs(c_pos1-c_pos2)){
	case 1:
	case 2:
	case 4:
	case 5:
	  {  
	    int sum_pos = static_cast<int>(c_pos1+c_pos2)%6;
	    if(this->nodes[chain[middle+1]].multi==2 && this->nodes[chain[middle-1]].multi==2){
	      sum_pos+=1;
	    }
	    for(int i=0;i<c_index.size();i++){
	      for(int j=0;j<c_index[i].size();j++){
		flip_carbon_pos[i][j] = (sum_pos+6-carbon_pos[i][j])%6;
	      }
	      std::sort(flip_carbon_pos[i].begin(),flip_carbon_pos[i].end());
	    }

	    if(flip_carbon_pos == carbon_pos)
	      return true;
	    else 
	      return false;     
	    break;
	  }
	case 3:
	  {
	    //flip by left-right reflection
	    float min = std::min(c_pos1,c_pos2);
	    if(this->nodes[chain[middle+1]].multi==2 && this->nodes[chain[middle-1]].multi==2){
	      min+=0.5;
	    }
	    for(int i=0;i<c_index.size();i++){
	      for(int j=0;j<c_index[i].size();j++){
		flip_carbon_pos[i][j] = static_cast<int>(2*min+6-carbon_pos[i][j])%6;
	      }
	      std::sort(flip_carbon_pos[i].begin(),flip_carbon_pos[i].end());
	    }
	    if(flip_carbon_pos == carbon_pos){
	      return true;
	    }
	    //flip by up-down reflection
	    int sum_pos = static_cast<int>(c_pos1+c_pos2)%6;
	    if(this->nodes[chain[middle+1]].multi==2 && this->nodes[chain[middle-1]].multi==2){
	      sum_pos+=1;
	    }
	    for(int i=0;i<c_index.size();i++){
	      for(int j=0;j<c_index[i].size();j++){
		if(carbon_pos[i][j]>sum_pos){
		  flip_carbon_pos[i][j] = sum_pos+6-carbon_pos[i][j];
		}else{
		  flip_carbon_pos[i][j] = sum_pos-carbon_pos[i][j];
		}
	      }
	      std::sort(flip_carbon_pos[i].begin(),flip_carbon_pos[i].end());
	    }
	    if(flip_carbon_pos == carbon_pos){
	      return true;
	    }
	    //flip by both left-right and up-down reflection
	    for(int i=0;i<c_index.size();i++){
	      for(int j=0;j<c_index[i].size();j++){
		flip_carbon_pos[i][j] = (carbon_pos[i][j]+3)%6;
	      }
	      std::sort(flip_carbon_pos[i].begin(),flip_carbon_pos[i].end());
	    }
	    if(flip_carbon_pos == carbon_pos){
	      return true;
	    }else{
	      return false;
	    }
	    break;
	  }  
	}
	std::cout<<"error in middle_chain_symmetry"<<std::endl;
	return false;
	}*/

int ChemTreeCenter::chain_middle_pair_symmetry(std::vector<int> chain, carbon_position &carbon_pos, int previous_unequal){
  int index1 = chain[chain.size()/2-1];
  int index2 = chain[chain.size()/2];

  //return  1 if cp[index1] < cp[index2] , where index1 < index2   
  //return  0 if cp[index1] = cp[index2] , where index1 < index2   
  //return -1 if cp[index1] > cp[index2] , where index1 < index2   

  if(index1 > index2){
    int temp = index1;
    index1 = index2;
    index2 = temp;
  }

  bool is_adjacent = false;
  int index_parent,index_child;
  
  //remove index1/index2 if it is another node's child (consider other adjacent nodes only)
  if(this->nodes[index1].parent==index2){
    index_parent = index2;
    index_child = index1;
    is_adjacent = true;
  }else if(this->nodes[index2].parent==index1){
    index_parent = index1;
    index_child = index2;
    is_adjacent = true;
  }
  
  if(!is_adjacent){
    carbon_position cp_index1, cp_index2;

    cp_index1 = carbon_pos;
    group_child(cp_index2, index2);
    std::sort(cp_index2.begin(), cp_index2.end(), compare_function);
    get_bond_position(cp_index2, index2);  
    
    if(cp_index1 < cp_index2){
      return 1;
    }else if(cp_index1 == cp_index2){
      return 0;
    }else{
      return -1;
    } 

  }else{
    //index1 and index2 are adjacent to each other

    ///parent must be root unless the tree is not left heavy
    if(index_parent!=0){
      return false;
    }
  
    //if it is not the first child of this label_atom(ex.carbon), it and its parent cannot have the same structure
    //for(int i=this->nodes[index_child].nth-1;i>=0;i--){
    //  if(this->get_label(this->nodes[index_parent].children[i])==this->get_label(index_child)){
    //    return false;
    //  }
    //}
    
    carbon_position cp_child, cp_parent, cindex_parent;
    //index parent must be current node -> set cp_parent to carbon_pos
    cp_parent = carbon_pos;
    group_child(cp_child, index_child);
    group_child(cindex_parent, index_parent);
    std::sort(cindex_parent.begin(), cindex_parent.end(), compare_function);
    int cp_of_child_in_cp_parent;
    
    //delete information of child benzene from cp_parent
    for(size_t i=0; i<cindex_parent.size(); i++){
      for(size_t j=0; j<cindex_parent[i].size(); j++){
	if(cindex_parent[i][j]==this->nodes[index_child].nth){
	  cp_of_child_in_cp_parent = cp_parent[i][j];

	  if(j==0 && cindex_parent[i].size() == this->nodes[index_child].multi){
	    cp_parent.erase(cp_parent.begin()+i);
	    cindex_parent.erase(cindex_parent.begin()+i);
	  }else{
	    cp_parent[i].erase(cp_parent[i].begin()+j);
	    cindex_parent[i].erase(cindex_parent[i].begin()+j);
	    if(this->nodes[index_child].multi==2){
	      //already remove the first one 
	      cp_parent[i].erase(cp_parent[i].begin()+j);
	      cindex_parent[i].erase(cindex_parent[i].begin()+j);
	    }
	  }
	  break;
	}
      }
    }
    std::sort(cp_child.begin(),cp_child.end(),compare_function);
    get_bond_position(cp_child,index_child);  
  
    if(this->nodes[index_child].multi==1){
      int nth = this->nodes[index_child].nth;
      //compare carbon position of child nodes (excluding index1 in cp_parent)
      //parent must be root so can compare by cp_child==cp_parent
    
      if(cp_child.size()!=cp_parent.size()){
	std::cout<<"error!!"<<std::endl;
      }else{
	if(cp_child < cp_parent){
	  if(index_child == index1){
	    return 1;//false;
	  }else{
	    return -1;
	  }
	}
	if(cp_child > cp_parent){
	  if(index_child == index1){
	    return -1;
	  }else{
	    return 1;
	  }
	}
      }
      return 0;//true;
  
    }else if(this->nodes[index_child].multi==2){
      
      if(previous_unequal == -1){
	//if chain[middle-1] and chain[middle+1] is not normal (cp[middle-1] < cp[middle+1], where index of node[middle-1] > index of node[middle+1]) -> not normal
	return -1;
      }
      
      if(previous_unequal == 1){
	//if chain[middle-1] and chain[middle+1] is normal (cp[middle-1] > cp[middle+1], where index of node[middle-1] > index of node[middle+1]) -> normal
	return 1;
      }
     
      int cp_shift = cp_of_child_in_cp_parent;//this->nodes[index_parent].bond_position[this->nodes[index_child].nth];
      /*std::cout<<"cp_parent = ";
      print_cp(cp_parent);

      for(int i=0;i<cp_parent.size();i++){
	for(int j=0;j<cp_parent[i].size();j++){
	  cp_parent[i][j] = (cp_shift-cp_parent[i][j]+6)%6;
	}
	std::sort(cp_parent[i].begin(),cp_parent[i].end());
      }     
      std::cout<<"cp_parent(after) = ";
      print_cp(cp_parent);*/
      
      bool check_newcarbon = true;
      bool check_carbon = true;
      bool check_leftright_updown_sym = true;
      
      carbon_position naph_cp_parent = cp_parent; //store carbon_position of child benzene in naph system
    
      //group_child(index_child_c_index, index_child);
      //std::sort(index_child_c_index.begin(), index_child_c_index.end(), compare_function);

      //int nth = this->nodes[index1].nth;
      cp_shift = 6 + cp_of_child_in_cp_parent;//this->nodes[index2].bond_position[nth];
      
      // convert carbon position of nodes adjacent to index2 to be between 1-4 based on cp_shift
      for(size_t i=0; i<naph_cp_parent.size(); i++){
	for(size_t j=0; j<naph_cp_parent[i].size(); j++){
	  /*if(index2_c_index[i][j]==this->nodes[index1].nth){
	    if(j==0 && index2_c_index[i].size()==2){
	      index2_c_index.erase(index2_c_index.begin()+i);
	      j=-1;//make j at new iteration =0
	    }else{
	      index2_c_index[i].erase(index2_c_index[i].begin()+j);
	      index2_c_index[i].erase(index2_c_index[i].begin()+j+1);
	      j--;
	    }
	  }else{*/
	  if(naph_cp_parent[i][j]!=-2){
	    //if(naph_cp_parent[i][j]!=-1){
	    //cp_parent does not contain -1 
	    if(cp_of_child_in_cp_parent != 0 && cp_of_child_in_cp_parent < 3){
	      naph_cp_parent[i][j] = (cp_shift-cp_parent[i][j])%6;
	    }else{
	      naph_cp_parent[i][j] = (cp_parent[i][j] + 5 - cp_of_child_in_cp_parent)%6;
	    }
	    /*}else{
	      if(cp_of_child_in_cp_parent != 0 && cp_of_child_in_cp_parent < 3){
	      naph_cp_parent[i][j] = (cp_shift - (cp_parent[i][j-1]+1) )%6;
	      }else{
	      naph_cp_parent[i][j] = (cp_parent[i][j-1] + 1 + 5 - cp_of_child_in_cp_parent)%6;
	      }
	      }
	    */
	  }else{
	    if(j!=0){
	      return -1;//false;
	    }
	  }
	  //}
	}
	std::sort(naph_cp_parent[i].begin(), naph_cp_parent[i].end());
      }
       
      std::vector< std::vector<int> > new_cp_child;
      new_cp_child = cp_child;
      
      //new_carbon_pos is rotate carbon position in reverse direction (1-2-3-4) to (4-3-2-1)
      //necessary only when two benzene is in the same cyclic structure
      for(size_t i=0; i<new_cp_child.size(); i++){
	for(size_t j=0; j<new_cp_child[i].size(); j++){
	  new_cp_child[i][j] = 5 - cp_child[i][j]; 
	}
	std::sort(new_cp_child[i].begin(), new_cp_child[i].end());
      }      

      for(size_t i=0; i<cp_child.size(); i++){    
	for(size_t j=0; j<cp_child[i].size(); j++){
	  //if(naph_cp_parent[i][j]!=-2){
	  if(check_leftright_updown_sym){  
	    //their common parent does not have parent, has only even unique childs
	    
	    if(check_newcarbon && check_carbon){
	      if(new_cp_child[i][j] > naph_cp_parent[i][j] && cp_child[i][j] > naph_cp_parent[i][j]){
		return  1;//normal
	      }
	      if(new_cp_child[i][j] < naph_cp_parent[i][j] || cp_child[i][j] < naph_cp_parent[i][j]){
		return -1;//not normal
	      }
		
	      //when reach here either new_carbon_pos or carbon_pos (but not both) must equal to index2.bond_position[....]
	      if(new_cp_child[i][j] == naph_cp_parent[i][j]){
		//if new_carbon_pos is equal, we will compare only new carbon_pos in the next round 
		check_carbon = false;
	      }else{
		check_newcarbon = false;
	      }
	    }else if(check_carbon && !check_newcarbon){
	      if(cp_child[i][j] > naph_cp_parent[i][j]){
		return  1;//normal
	      }
	      if(cp_child[i][j] < naph_cp_parent[i][j]){
		return -1;//not normal
	      } 
	    }else if(!check_carbon && check_newcarbon){
	      if(new_cp_child[i][j] > naph_cp_parent[i][j]){
		return  1;//normal
	      }
	      if(new_cp_child[i][j] < naph_cp_parent[i][j]){
		return -1;//not normal
	      }
	    }
		  
	  }else{//not leftright updown sym
	    if(new_cp_child[i][j] > naph_cp_parent[i][j]){
	      return  1;//normal
	    }
	    if(new_cp_child[i][j] < naph_cp_parent[i][j]){
	      return -1;//not notmal
	    }
	  }
    
	}//carbon_pos[i][j]!=-2	  
      }//for j of carbon_pos[i][j]
    }//for i of carbon_pos[i][j]
    return 0;//false;


      /*
      //used to shift carbon position of child (in parent bond_position) to 0
      //not consider the case that more than 2 benzenes merge together ================================
      //if consider that case, cp of naph bond must be included (now only one cp per one atom from cp_child/cp_parent)
      int cp_shift = cp_of_child_in_cp_parent;//this->nodes[index_parent].bond_position[this->nodes[index_child].nth];

      for(int i=0;i<cp_parent.size();i++){
	for(int j=0;j<cp_parent[i].size();j++){
	  cp_parent[i][j] = (cp_shift-cp_parent[i][j]+6)%6;
	}
	std::sort(cp_parent[i].begin(),cp_parent[i].end());
      }
      
      bool leftright_sym = true;
      //left-right symmetry between 2 benzenes
      
      //for(int i=0; i<this->nodes[index_child].num_children; i++){
      	//if(cp_child[i]!=cp_parent[i]){
	  //leftright_sym = false;
	  //}
	//}
      if(cp_child.size()!=cp_parent.size()){
	std::cout<<"error!!"<<std::endl;
      }else{
	if(cp_child!=cp_parent){
	  leftright_sym = false;
	}
      }
      if(leftright_sym){
	std::cout<<"left-right sym"<<std::endl;
	return 1;
      }
      
      //left-right and up-down symmetry between 2 benzenes
      
      //for(int i=0;i<this->nodes[index_child].num_children;i++){
      //  int shift=0;
      //  if(i>=this->nodes[index_child].nth){
	//shift=1;
      //}
      //if(cp_child[i]!=5-cp_parent[i+shift]){
      //return false;
      //}
      //}
      for(int i=0;i<cp_parent.size();i++){
	for(int j=0;j<cp_parent[i].size();j++){
	  cp_parent[i][j] = 5-cp_parent[i][j];
	}
	std::sort(cp_parent[i].begin(),cp_parent[i].end());
      }

      if(cp_child.size()!=cp_parent.size()){
	std::cout<<"error!!"<<std::endl;
      }else{
	if(cp_child==cp_parent){
	  std::cout<<"left-right up-down sym"<<std::endl;
	  return 1;
	}
	}
    return 1;
    }*/
  }//end of if-else is_adjacent
  std::cout<<"error in chain_symmetry case:1,2"<<std::endl;
  return -1;//true
}

void ChemTreeCenter::find_chain_down(int step,int index_step,std::vector<int> &chain,std::vector< std::vector<int> > &chain_collection){
  //index_step is index of step in temp ex.if step is the first element, index_step = 0
  std::array<int,6> child = this->nodes[step].children;
  
  for(int i=0;i<this->nodes[step].num_children;i++){
    int index_child = this->nodes[step].children[i];
    if(this->nodes[index_child].num_children>0 && !found(chain,index_child) && index_child<chain[0]){ 
      //index_child < step to assure that one chain is checked only once
      chain.push_back(index_child);
      
      find_chain_down(index_child,index_step-1,chain,chain_collection);
      chain.pop_back();
    }
  }
  if(this->get_label(step)==0){
    chain_collection.push_back(chain);
  }		
}

void ChemTreeCenter::find_chain(std::vector< std::vector<int> > &chain_collection){
  for(size_t index = 0; index< this->get_num_nodes(); index++){
    if(this->get_label(index)==0){
      //find list of index of benzene chain 
      std::vector<int> upward_chain;
      upward_chain.push_back(index);
      int step = this->nodes[index].parent;
  
      while(step>=0){
	upward_chain.push_back(step);
	this->find_chain_down(step,upward_chain.size()-1,upward_chain,chain_collection);
	step = this->nodes[step].parent;
      }
    }
  }
}


bool ChemTreeCenter::is_most_middle_bnode(const std::vector<int> &chain, const int index){
  //return true if index is the most middle pair (in the right side) of benzene node of the chain
  for(size_t i=chain.size()/2; i<chain.size(); i++){
    if(this->get_label(chain[i])==0){
      if(chain[i] == index){
	return true;
      }else{
	return false;
      }
    }
  }
  return false;
}


bool ChemTreeCenter::is_normal_chain(const std::vector< std::vector<int> > & position,int j,int k, int index, const carbon_position &c_index,const std::vector< std::vector<int> > & chain_collection,const std::vector<std::pair<std::vector<int>,std::vector<int> > > &trisym_chain_collection, int flip_mode){

	if(chain_collection.size()==0 && trisym_chain_collection.size()==0){
	  return true;
	}

	if(k!=position[j].size()-1){
	  //check for normal chain if all carbon in set j is known (k is the last element of j) only.
	  return true;
	}

	carbon_position temp_pos = position;
	for(size_t i=0; i<chain_collection.size(); i++){
	  
	  //remove non-benzene node from chain_collection[i] and store it in benzene_chain
	  std::vector<int> benzene_chain;
	  for(size_t j=0; j<chain_collection[i].size(); j++){
	    if(get_label(chain_collection[i][j]) == 0){
	      benzene_chain.push_back(chain_collection[i][j]);
	    }
	  }
	  //if index is the most middle benzene node or right node of the middle edge of chain[i] -> find chain_middle symmetry
	 
	  //if(is_most_middle_bnode(chain_collection[i], index) ){
	  if(*std::min_element(benzene_chain.begin(), benzene_chain.end()) == index){

	    /*if(benzene_chain.size()%2 != 0){
	    //edited 11 jan 2016
	      bool is_sym = false;
	      
	      //if the middle node is benzene node, check if it is symmetric or not
	      //but treating subtrees of two child nodes in sympath as they have the same structures
	      //although they have nodes with different carbon position

	      std::vector< std::vector<int> > adj_list;
	      group_child(adj_list,index);
	      std::sort( adj_list.begin(), adj_list.end(), compare_function);

	      for(int k=0; k < adj_list.size(); k++){
		for(int l=0; l < adj_list[k].size(); l++){
		  adj_list[k][l] = this->nodes[index].children[ adj_list[k][l] ];
		}
	      }
	      //adj_list is adjcent node list storing index of adj nodes
	      int index_in_sympath = std::min_element(benzene_chain.begin(), benzene_chain.end()) - benzene_chain.begin();
	      std::vector< std::vector<int> > new_position = position;
	      //new position merge two child nodes together although they have different structures
	      
	      for(int adj_i=0; adj_i < adj_list.size(); adj_i++){	      
		std::vector<int> adj_set = adj_list[adj_i];
		int index_another_node = -1;
		if(found(adj_set, chain_collection[i][index_in_sympath-1]) && !found(adj_set, chain_collection[i][index_in_sympath+1])){
		  index_another_node = chain_collection[i][index_in_sympath + 1];
		}else if(!found(adj_set, chain_collection[i][index_in_sympath-1]) && found(adj_set, chain_collection[i][index_in_sympath+1])){
		  index_another_node = chain_collection[i][index_in_sympath - 1];
		}

		if(index_another_node != -1){
		  for(int adj_j = adj_i+1; adj_j < adj_list.size(); adj_j++){	      		
		    int index_in_adj = std::find(adj_list[adj_j].begin(), adj_list[adj_j].end(), index_another_node) - adj_list[adj_j].begin();
		    if(index_in_adj < adj_list[adj_j].size()){
		      new_position[adj_i].push_back( new_position[adj_j][index_in_adj] );
		      if(adj_list[adj_j].size() == 1){		    		    	      		
			adj_list.erase( adj_list.begin() + adj_j );
			new_position.erase( new_position.begin() + adj_j );
		      }else{ 	
			adj_list[adj_j].erase( adj_list[adj_j].begin() + index_in_adj );
			new_position[adj_j].erase( new_position[adj_j].begin() + index_in_adj );
		      }
		    }
		  }
		  break;
		}
	      }

	      int mode;
	      if(index!=0){
		if(this->nodes[index].multi==2){
		  if(is_updown_symmetry(this->nodes[index].parent,index,this->nodes[index].nth)){
		    mode = 2;
		  }else{
		    mode = -1;
		  }
		}else{
		  mode = 0;
		}
	      }else{
		if(new_position[0].size()==1){
		  mode = 0;
		}else{
		  //if(j>0){
		  //if it is not the first group of child, need reflection due to first group's position
		  mode = 1;
		  //}else{
		  //if it is the first group of child, need left-right reflection
		  //flip_mode = 0;
		  //}
		}
	      }

	      int label = this->get_label(index);
	      std::vector< std::vector<int> > flip_position = new_position;
	      for(int a = 0; a < flip_position.size(); a++){
		for(int b = 0; b < flip_position[a].size(); b++){
		  flip_position[a][b] =  flip_carbon(new_position[a][b],label,new_position[0],mode);
		}
		std::sort(flip_position[a].begin(),flip_position[a].end());
	      }
	      if(flip_position == new_position){
		is_sym = true;
	      }

	      if(mode == 1 && new_position[0].size() == 2 && new_position[0][1] == 3){
		mode = 0;
		
		for(int a = 0; a < flip_position.size(); a++){
		  for(int b = 0; b < flip_position[a].size(); b++){
		    flip_position[a][b] =  flip_carbon(new_position[a][b],label,new_position[0], mode);
		  }
		  std::sort(flip_position[a].begin(),flip_position[a].end());
		}
		
		if(flip_position == new_position){
		  is_sym = true;
		}
	      }
	      
	      if(!is_sym){		
		return true;
		}
	       
	    }*/

	    if(chain_collection[i].size()>3){
	      int middle_left, middle_right;
	      if(chain_collection[i].size()%2==0){
		middle_left = chain_collection[i].size()/2-1;
		middle_right = chain_collection[i].size()/2;
	      }else{
		middle_left = chain_collection[i].size()/2-1;
		middle_right = chain_collection[i].size()/2+1;
	      }

	      //check if carbon position of nodes in the tsub of all pairs is equal or not
	      if(!is_tsub_cp_equal(chain_collection[i][middle_left], chain_collection[i][middle_right], middle_left, middle_right, chain_collection[i]) ){
		return true;
	      }
	    }
	
	    int previous_unequal = 0; // = chain_end_check(chain_collection[i], temp_pos, c_index, index);
	    for(size_t a=0; benzene_chain[benzene_chain.size()-1-a] != index ; a++){
	      //check the symmetry of non-middle benzene in a chain(except for the both end) 
	      int temp = chain_symmetry(benzene_chain[a], benzene_chain[benzene_chain.size()-1-a]);
	      switch(temp){
	      case -1: //cp[index1] > cp[index2], where index1 < index2
	      case  1: //cp[index1] < cp[index2], where index1 < index2
		previous_unequal = temp;
		break;
	      case 0:
		break;
	      }
	    }	      

	    if(benzene_chain.size()%2 == 0){
	      int temp = chain_middle_pair_symmetry(benzene_chain, temp_pos, previous_unequal);
	      switch(chain_middle_pair_symmetry(benzene_chain, temp_pos, previous_unequal) ){
	      case 1:
		return false; //chain[i] < chain[n-i] -> not normal
	      case 0:
		if(previous_unequal == 1) // chain[middle-2] < chain[middle+2] -> not normal
		  return false;
		else  // chain[middle-2] <= chain[middle+2] -> go to the next chain
		  break;
	      case -1:
		break;  //chain[i] > chain[n-i] -> go to the next chain
	      }
	    }else{
	      if(previous_unequal == 1){
		return false;
	      }
	    }
	    //reach here if all pairs are equal to each other -> go to next chain
	    
	  }//current index is not the middle node of chain -> go to next chain
	}
	//if all pairs of all chains are equal -> normal
	
	for(int i=0;i<trisym_chain_collection.size();i++){
	  int last_index = chain_collection[i].size()-1;
	  if(is_trisymmetry_redundant(trisym_chain_collection[i].first[last_index],trisym_chain_collection[i].second[last_index],temp_pos)){
	    return false;
	  }
	
	}
	//if no symmetry occur for every chain -> it is normal_benzene
	//if every symmetry chain is not redundant -> it is normal_benzene
	return true;
}

bool ChemTreeCenter::is_normal_benzene(const std::vector< std::vector<int> > &position,int j,int k,int carbon_pos,int index,int flip_mode,bool &normal){
  int label = this->get_label(index);
  bool recheck = false;
  if(flip_mode==-1){
    return true;
  }
  if(k+1<position[j].size() || label==1){
    return true;
  }

  if(flip_mode==1){
    if(j==0){
      //if the smallest c_index set has size 2 and it is the first group of child, need left-right reflection
      flip_mode = 0;
    }else{
      if(position[0][1]==3){
	//if carbon position of the smallest set are 1 and 3, need left-right reflection
	recheck = true;
      }
    }
  }

  std::vector<int> flip_position;
  int start = 0;
  
  /*if(normal){
    start = std::max(0,j-1);
    }*/
  /*for(int a=start;a<j;a++){
    for(int b=0;b<position[a].size();b++){
      flip_position.push_back(flip_carbon(position[a][b],label,position[0],flip_mode));
    }
    std::sort(flip_position.begin(),flip_position.end());
    for(int b=0;b<flip_position.size();b++){
      if(flip_position[b]<position[a][b]){
	return false;
      }
      if(flip_position[b]>position[a][b]){
	if(recheck){
	  return is_normal_benzene(position,j,k,carbon_pos,index,0,normal);
	}else{
	  if(k==position[j].size()-1 && a>0){
	    normal = true;
	  }
	  return true;
	}
      }
    }
    flip_position.clear();
  }*/
  //check normal form of child group j
  for(int b=0;b<k;b++){
    flip_position.push_back(flip_carbon(position[j][b],label,position[0],flip_mode));
  }
  flip_position.push_back( flip_carbon(carbon_pos,label,position[0],flip_mode) ); //add new carbon position to flip vector
  std::sort(flip_position.begin(),flip_position.end());
  
  /*std::cout<<" position =";
  for(size_t i = 0; i<position[j].size(); i++)
    std::cout<<position[j][i]<<" ";
  std::cout<<std::endl;
  std::cout<<" flip position =";
  for(size_t i = 0; i<flip_position.size(); i++)
    std::cout<<flip_position[i]<<" ";
  std::cout<<std::endl;
  */
  for(int b=0; b < k; b++){
    if(flip_position[b] < position[j][b]){  
      return false;
    }
    if(flip_position[b]>position[j][b]){
      if(recheck){
	if(is_normal_benzene(position,j,k,carbon_pos,index,0,normal)){
	  if(k==position[j].size()-1){
	    if(j>0 || index!=0){
	      normal = true;
	    }
	  }  
	  return true;
	}else{ 
	  return false;
	}
      }else{
	if(k==position[j].size()-1){
	  if(j>0 || index!=0){
	    normal = true;
	  }
	}
	return true;
      }
    }
  }

  //compare last carbon position
  if(flip_position[k] < carbon_pos){
    return false;
  }else if(flip_position[k] > carbon_pos){
    if(recheck){
      if(is_normal_benzene(position,j,k,carbon_pos,index,0,normal)){
	if(k==position[j].size()-1){
	  if(j>0 || index!=0){
	    normal = true;
	  }
	}

	return true;
      }else{
	return false;
      }
    }else{
      if(k==position[j].size()-1){
	if(j>0 ||index!=0){
	  normal = true;
	}
      }
      return true;
    }
  }else{
    //position is equal to flip_position
    if(recheck){
      if(is_normal_benzene(position,j,k,carbon_pos,index,0,normal)){
	if(k==position[j].size()-1){
	  if(j>0 || index!=0){
	    normal = true;
	  }
	}
	return true;
      }else{
	return false;
      }
      //return is_normal_benzene(position,j,k,carbon_pos,index,0,normal);
    }
    return true;
  }
  return true;
}

bool ChemTreeCenter::normal_naphthalene(int index){
  if(nodes[index].num_children==0)
    return true;

  std::vector< std::vector<int> > cp_group;
  group_child( cp_group, index);
  std::sort( cp_group.begin(), cp_group.end(), compare_function);
  get_bond_position(cp_group, index);
  
  //int j = cp_group.size()-1;
  int flip_mode = 2;
  bool normal = false;
  for(size_t j = 0; j<cp_group.size(); j++){
    size_t k = cp_group[j].size()-1;
    int cp = cp_group[j][k];
    bool result = is_normal_benzene(cp_group, j, k, cp, index, flip_mode, normal);
    if(!result){
      return false;
    }
    if(normal){
      return true;
    }
  }

  return true;//result;
}

size_t ChemTreeCenter::label_child_and_next(std::vector<carbon_position> &result, int j, int k, int carbon_pos, int max_child, int index, const carbon_position &c_index, std::vector< std::vector<int> > & chain_collection, const std::vector< std::pair<std::vector<int>, std::vector<int> > > &trisym_chain_collection, int flip_mode, bool &normal, std::ofstream & output_file){
	//j and k are indices of position that will be labeled by carbon position
	//index is index of benzene (used to determine if it is root or not)
        int last_index = result.size()-1;
	result[last_index][j][k] = carbon_pos;		
 	carbon_position current_cp = result[last_index];

	if(j+1==current_cp.size() && k+1==current_cp[j].size()){	  
	  assign_bond_position(index, current_cp, c_index);
	  if(is_normal_chain(current_cp, j, k, index, c_index, chain_collection, trisym_chain_collection, flip_mode)){	    
	    if(get_label(index)==0){ 
	      //check for normal benzene with left-right sym in child benzene of a naphthalene ring
	      for(size_t i=0; i<nodes[index].num_children; i++){
		if(nodes[ nodes[index].children[i] ].multi==2 && is_updown_symmetry(index, nodes[index].children[i], i) ){
		  if(!normal_naphthalene(nodes[index].children[i]) ){
		    return 0;
		  }
		  break;
		}
	      }
	    }
	    return this->label_benzene(index-1, chain_collection, result, output_file);
	  }else{
	    return 0;
	  }
	}
	size_t num=0;
	
	//precomputed carbon_position for first group of child benzene (j==0) when benzene does not have parent (index==0)
	if(j==0 && index==0 && k+1<current_cp[j].size() && this->get_label(this->nodes[index].children[c_index[j][k]])!=0){
	  switch(current_cp[j].size()){
	  case 2:
	    {
	      for(int cp=1;cp<=3;cp++){
		num+=label_child_and_next(result, j, k+1, cp, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
		result[last_index][j][k+1]=-2;
	      }
	    }
	    break;
	  case 3:
	    {
	      switch(k){
	      case 0:
		for(int cp=1;cp<=2;cp++){
		  num+=label_child_and_next(result, j, k+1, cp, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
		  result[last_index][j][k+1]=-2;
		}
		break;
	      case 1:
		if(carbon_pos==1){
		  for(int cp=2;cp<=3;cp++){
		    num+=label_child_and_next(result, j, k+1, cp, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
		    result[last_index][j][k+1]=-2;
		  }
		}else if(carbon_pos==2){
		  num+=label_child_and_next(result, j, k+1, 4, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
		  result[last_index][j][k+1]=-2;
		}else{
		  std::cout<<"error===="<<std::endl;
		}
	      }
	    }
	    break;
	  case 4:
	    {
	      switch(k){
	      case 0:
		{
		  num+=label_child_and_next(result, j, k+1, 1, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
		  result[last_index][j][k+1]=-2;
		}
		break;
	      case 1:
		for(int cp=2;cp<=3;cp++){
		  num+=label_child_and_next(result, j, k+1, cp, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
		  result[last_index][j][k+1]=-2;
		}
		break;
	      case 2:
		for(int cp=carbon_pos+1;cp<=4;cp++){
		  num+=label_child_and_next(result, j, k+1, cp, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
		  result[last_index][j][k+1]=-2;
		}
	      }
	    }
	    break;
	  case 5:
	  case 6:
	    {
	      num+=label_child_and_next(result, j, k+1, k, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
	      result[last_index][j][k+1]=-2;
	    }
	    break;
	  }
	  return num;
	}

	/* need to check for normal_chain -> disable
	//precomputed carbon_position of first group of child of benzene (j==0) when benzene has single bond with parent (multi==1)
	if(j==0 && this->nodes[index].multi==1 && k+1<current_cp[j].size() && this->get_label(this->nodes[index].children[c_index[j][k]])!=0){
	  switch(current_cp[j].size()){
	  case 2:
	    for(int cp = carbon_pos+1;cp<=6-carbon_pos;cp++){
	      num+=label_child_and_next(result,j,k+1,cp,max_child,index,c_index,chain_collection,flip_mode);
	      result[last_index][j][k+1]=-2;
	    }
	    break;
	  case 3:
	    switch(k){
	    case 0:
	      for(int cp=carbon_pos+1;cp<=3;cp++){
		num+=label_child_and_next(result,j,k+1,cp,max_child,index,c_index,chain_collection,flip_mode);
		result[last_index][j][k+1]=-2;
	      }
	      break;
	    case 1:
	      for(int cp=carbon_pos+1;cp<=6-current_cp[j][k-1];cp++){
		num+=label_child_and_next(result,j,k+1,cp,max_child,index,c_index,chain_collection,flip_mode);
		result[last_index][j][k+1]=-2;
	      }
	    }
	    break;
	  case 4:
	    switch(k){
	    case 0:
	      num+=label_child_and_next(result,j,k+1,2,max_child,index,c_index,chain_collection,flip_mode);
	      result[last_index][j][k+1]=-2;
	      break;
	    case 1:
	    case 2:
	      for(int cp=carbon_pos+1;cp<=5 && cp<=carbon_pos+2;cp++){
		num+=label_child_and_next(result,j,k+1,cp,max_child,index,c_index,chain_collection,flip_mode);
		result[last_index][j][k+1]=-2;
	      }
	      break;
	    }
	    break;
	  case 5:
	    num+=label_child_and_next(result,j,k+1,k+1,max_child,index,c_index,chain_collection,flip_mode);
	    result[last_index][j][k+1]=-2;
	    break;
	  }
	  return num;
	  }*/

	if(k+1<current_cp[j].size() && current_cp[j][k+1]==-1){
		// set second carbon position of naphthalene
		int cp = current_cp[j][k]+1; 
		if(cp>=max_child){
		  return 0;
		}

		for(int a=0;a<current_cp.size();a++){
		  if(found(current_cp[a],cp)){
		    return 0;
		  }
		}
		if(normal || is_normal_benzene(current_cp,j,k+1,cp,index,flip_mode,normal)){
		  //if(is_normal_chain(current_cp,j,k+1,cp,index,c_index,chain_collection,trisym_chain_collection)){
		  num+=label_child_and_next(result, j, k+1, cp, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
		  result[last_index][j][k+1]=-1;
		  //}
		}else{
		  return 0;
		}
	}else{
	  bool init_normal = normal;
	  for(int cp=1;cp<max_child;cp++){
	     bool valid = true;
	     normal = init_normal;
	     for(int a=0;a<current_cp.size();a++){
	       if(found(current_cp[a],cp)){
		 valid=false;
		 break;
	       }
	     }
	     if(valid){	
	       if(k+1==current_cp[j].size()){		     		 
		 if(normal || is_normal_benzene(current_cp,j+1,0,cp,index,flip_mode,normal)){
		   //if(this->is_normal_chain(current_cp,j+1,0,cp,index,c_index,chain_collection,trisym_chain_collection)){
		   num+=label_child_and_next(result, j+1, 0, cp, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
		     result[last_index][j+1][0]=-2;
		     //}		
		 }
	       }else{		 
		 if( max(current_cp[j]) < cp ){
		   if(normal || is_normal_benzene(current_cp,j,k+1,cp,index,flip_mode,normal)){
		     //if(this->is_normal_chain(current_cp,j,k+1,cp,index,c_index,chain_collection,trisym_chain_collection)){
		     num+=label_child_and_next(result, j, k+1, cp, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
		     result[last_index][j][k+1]=-2; 
			 //}
		   }
		 }
	       } 
	     }
	   }
	}	
	return num;	
}

bool compare_function(std::vector<int> i,std::vector<int> j){
  return i.size()<j.size();
}

bool chain_compare_function(std::vector<int> chain1, std::vector<int> chain2){
  if(chain1[0] != chain2[0]){
    return chain1[0] > chain2[0];
  }else{
    return chain1.size() < chain2.size();
  }
}

size_t ChemTreeCenter::label_benzene(const int index, std::vector< std::vector<int> > &chain_collection, std::vector<carbon_position> &result, std::ofstream & output_file){
	using namespace std;

	if(index == -1){ //this->num_nodes-1){
	  if(do_print){
	    //std::cout<<"--------------------"<<std::endl;
	    //this->show( true );
	    
	    int num_cycles = 1;
	    write_smiles(0, output_file, num_cycles);
	    output_file<<"\n";
	  }
	  return 1;
	}

	//if(this->get_label(index)==0 && nodes[index].num_children>0){
	if(this->get_label(index)<num_special_atom && nodes[index].num_children>0){
	  //std::vector< std::vector<int> > chain_collection; //store index of benzene chain starting from itself to the the top one	
	  std::vector< std::pair< std::vector<int>, std::vector<int> > > trisym_chain_collection;
	  
	  /*
	  //find all possible chain of benzene at 'index' and store in chain
	  //each chain does not include benzene at the leaf node
	  for(int i=0;i<chain_collection.size();i++){
	    bool skip = false;
	    if(chain_collection[i].size()%2==1 && !chain_middle_symmetry(chain_collection[i])){ 
	      //if the middle node is not symmetry -> goto next chain
	      chain_collection.erase(chain_collection.begin()+i);
	      i--;
	      continue;
	    }

	    for(int a=1;a<chain_collection[i].size()/2;a++){
	      //check the symmetry of non-middle benzene in a chain(except for the both end) 
	      //if one of them is not symmetry -> goto the next chain	 
	      if(!chain_symmetry(chain_collection[i][a],chain_collection[i][chain_collection[i].size()-1-a])){
		chain_collection.erase(chain_collection.begin()+i);
		i--;
		skip = true;
		break;
	      }
	    }
	    if(!skip && !chain_end_symmetry(chain_collection[i][0],chain_collection[i][chain_collection[i].size()-1])){
	      chain_collection.erase(chain_collection.begin()+i);
	      i--;
	      continue;
	    }  
	    if(!skip){
	      for(int j=i+1;j<chain_collection.size();j++){
		if(chain_collection[i].size()==chain_collection[j].size() && is_trisymmetry(chain_collection[i],chain_collection[j])){
		  std::pair< std::vector<int>, std::vector<int> > temp (chain_collection[i],chain_collection[j]);
		  trisym_chain_collection.push_back(temp);
		  chain_collection.erase(chain_collection.begin()+i);
		  chain_collection.erase(chain_collection.begin()+j);
		  i--;
		  if(j==i+1){
		    i--;
		  }
		  break;
		}
	      }
	    }
	    }*/

	        vector< vector<int> > c_index; //order of children having the same structure are grouped together
	        int max_child = 6; 
		size_t num = 0;
		//max_child = maximum carbon position for its child
		 
		//although node[i] is root, after assign one of it child to position 0 => maximum carbon position is 5
		if(nodes[index].multi==2)
		  max_child--;

		//group child of node[i] with same structure together
		//ex. c_index =[ [1 2] [3] ] means 1st and 2nd child have the same structure, while 3rd child has different one
		group_child(c_index,index);

		//sort c_index according to size and valence value
		std::sort(c_index.begin(),c_index.end(),compare_function);

		if(c_index[0].size() >= 3 && index == 0){
		  //if the smalles size of sets in adj_nodes list >=3, assign carbon position directly
		  vector<int> temp1;
		  carbon_position temp2;
		  switch(c_index[0].size()){
		  case 3:
		    for(int i =0; i< fix_cp3_1.size(); i++){
		      temp1 = fix_cp3_1[i];
		      temp2.clear();
		      temp2.push_back(temp1);
		      
		      if(c_index.size() == 2){
			temp1 = fix_cp3_2[i];
			temp2.push_back(temp1);
		      }
		      result.push_back(temp2);
		      num+= this->label_benzene(index-1, chain_collection, result, output_file);
		    }
		    break;
		  case 4:
		    for(int i =0; i< fix_cp4.size(); i++){
		      temp1 = fix_cp4[i];
		      temp2.clear();
		      temp2.push_back(temp1);
		      result.push_back(temp2);
		      num+= this->label_benzene(index-1, chain_collection, result, output_file);
		    }
		    break;
		  case 5:
		    temp1 = {0,1,2,3,4};
		    temp2.push_back(temp1);
		    result.push_back(temp2);
		    num+= this->label_benzene(index-1, chain_collection, result, output_file);
		    break;
		  case 6:
		    temp1 = {0,1,2,3,4,5};
		    temp2.push_back(temp1);
		    result.push_back(temp2);
		    num+= this->label_benzene(index-1, chain_collection, result, output_file);
		    break;
		  }
		  return num;
		}

		carbon_position initial_cp = c_index;
		for(int i=0;i<initial_cp.size();i++){
			for(int j=0;j<initial_cp[i].size();j++){
				if(initial_cp[i][j]!=-1)
					initial_cp[i][j]=-2; // -2 is initialization (-1 is naphthalene bond)
			}
		}
		result.push_back(initial_cp);

		int flip_mode;
		if(index!=0){
		  if(this->nodes[index].multi==2){
		    if(is_updown_symmetry(this->nodes[index].parent,index,this->nodes[index].nth)){
		      flip_mode = 2;
		    }else{
		      flip_mode = -1;
		    }
		  }else{
		    flip_mode = 0;
		  }
		}else{
		  if(initial_cp[0].size()==1){
		    flip_mode = 0;
		  }else{
		    //if(j>0){
		    //if it is not the first group of child, need reflection due to first group's position
		    flip_mode = 1;
		    //}else{
		    //if it is the first group of child, need left-right reflection
		    //flip_mode = 0;
		    //}
		  }
		}
		
		bool normal = false;
		if(index==0){
		  num+= label_child_and_next(result, 0, 0, 0, max_child, index, c_index, chain_collection, trisym_chain_collection, flip_mode, normal, output_file);
		}else{

		  int maximum = max_child;
		  if(this->nodes[index].multi==1) maximum=4;
		  for(int cp=1;cp<maximum;cp++){ 
		    normal = false;

		    if(is_normal_benzene(initial_cp,0,0,cp,index,flip_mode,normal)){
		      //if(is_normal_chain(initial_cp,0,0,cp,index,c_index,chain_collection,trisym_chain_collection)){
		      num+= label_child_and_next(result,0,0,cp,max_child,index,c_index,chain_collection,trisym_chain_collection,flip_mode,normal, output_file);

			result[result.size()-1][0][0] = -2;
			//}
		    }
		  }
		}

		result.pop_back(); //result.erase(result.begin()+index);
		return num;
	}else{ // else of if(index==0 && num_children>0)

	  carbon_position temp;
	  result.push_back(temp);
	  size_t num = this->label_benzene(index-1, chain_collection, result, output_file);
	  result.pop_back();   //erase(result.begin()+index);

	  return num;
	}
}

size_t enum_benzene( ChemTreeCenter& tree, std::ofstream & output_file){

        size_t num = 0; 
	std::vector<int> multi;
	if(num_lack_H == 0 && num_naph_bond==0){
	//store previous information of multi in multi
		for(int index=0;index<tree.get_num_nodes();index++){
			multi.push_back(tree.get_multi(index));
			tree.set_multi_bond(index,1);
		}
	}
	std::vector<carbon_position> result;
	
	std::vector< std::vector<int> > chain_collection; //store index of benzene chain starting from itself to the the top one	
	tree.find_chain(chain_collection);
	
	for(size_t i=0; i<chain_collection.size(); i++){	  
	  bool skip = false;
	  for(size_t a=0; a<chain_collection[i].size()/2; a++){
	    //check the symmetry of non-middle benzene in a chain(except for the both end) 
	    //if one of them is not symmetry -> goto the next chain	 
	    if(!tree.is_tsub_equal( chain_collection[i][a], chain_collection[i][chain_collection[i].size()-1-a] ) ){
	      chain_collection.erase(chain_collection.begin()+i);
	      i--;
	      skip = true;
	      break;
	    }
	  }	  
	  if(!skip && tree.is_tsub_has_benzene(chain_collection[i][0])){
	    chain_collection.erase(chain_collection.begin()+i);
	    i--;
	  }
	}

	std::sort(chain_collection.begin(), chain_collection.end(), chain_compare_function);
	
	num = tree.label_benzene(tree.get_num_nodes()-1, chain_collection, result, output_file);			
	
	for(size_t index = 0; index<tree.get_num_nodes(); index++){
	  tree.reset_bond_position(index);
	}

	if(num_lack_H == 0 && num_naph_bond==0){
	//restore previous information of multi
		for(int index=0;index<multi.size();index++){
			tree.set_multi_bond(index,multi[index]);
		}
	}

	return num;
}

inline size_t set_multi_and_next(ChemTreeCenter& tree, is_ident_type is_ident, const int i, valence_value_type multiple, const int double_bond, const int triple_bond,const int naph_bond, const int normal_type, std::ofstream & output_file) 
{
	using namespace std;
	tree.set_multi_bond(i, multiple);
	if(double_bond==0 && triple_bond==0 && naph_bond==0){
		tree.fill_rest_single_bond(i+1);
		if(normal_type >= 2){
			if(not tree.is_multi_normal(normal_type - 2)) {
				return 0;
			}
		}
		return enum_benzene(tree, output_file);
	}

	--multiple;
	lack_degree[i] -= multiple;
	const int pi = tree.get_parent(i);
	lack_degree[pi] -= multiple;

	tree.update_identical_multi(is_ident, i);

	const int num_nodes = tree.get_num_nodes();

	size_t num = 0;
	for (int k = i + 1; k + triple_bond + double_bond +naph_bond <= num_nodes; ++k) {
		if(tree.can_be_added_multi(k)){ 

			const valence_value_type c_degree = lack_degree[k];
			const int pk = tree.get_parent(k);
			const valence_value_type p_degree = lack_degree[pk];
			const valence_value_type maxmulti = tree.max_multi(is_ident, k);

			if(tree.get_label(k)!=0){

				if ((c_degree >= 1) and (p_degree >= 1) and (maxmulti >= 2)) {
					if (double_bond > 0) {
					  num += set_multi_and_next(tree, is_ident, k, 2, double_bond - 1, triple_bond,naph_bond, normal_type, output_file);
					}
					if ((triple_bond > 0) and (c_degree >= 2) and (p_degree >= 2) and (maxmulti >= 3)) {
					  num += set_multi_and_next(tree, is_ident, k, 3, double_bond, triple_bond - 1,naph_bond, normal_type, output_file );
					}
				}
			}else{ //tree.get_label(k)==0
				if(naph_bond>0 && c_degree>=1 && p_degree>=1 && maxmulti>=2){
				  num += set_multi_and_next(tree, is_ident, k, 2, double_bond, triple_bond, naph_bond-1, normal_type, output_file);
				}
			}
		}//end if(tree.can_be_added_multi(k))

		tree.set_multi_bond(k, 1);
		tree.update_identical_multi(is_ident, k);
	}

	lack_degree[i] += multiple;
	lack_degree[pi] += multiple;

	return num;
}

inline size_t generate_multi(ChemTreeCenter tree, // copy call because replaceing multi
			     is_ident_type& is_ident, const int normal_type, std::ofstream & output_file)
{
  
	using namespace std;

	for (int i = tree.get_parent(num_except_H - 1); i != num_except_H - 1; ++i) {
		tree.update_identical_end(is_ident, i);
	}
	tree.calc_lack_degree();

	const int num_nodes = tree.get_num_nodes();
	size_t num = 0;
	
	for(int triple_bond = 0; triple_bond <= (num_lack_H >> 2); ++triple_bond){
		const int double_bond = (num_lack_H >> 1) - (triple_bond << 1);
		for (int k = 1; k + triple_bond + double_bond + num_naph_bond <= num_nodes; ++k) {

			if(tree.can_be_added_multi(k)){
				const valence_value_type c_degree = lack_degree[k];
				const int pk = tree.get_parent(k);
				const valence_value_type p_degree = lack_degree[pk];
				const valence_value_type maxmulti = tree.max_multi(is_ident, k);

				if(tree.get_label(k)>=num_special_atom){

					if ((c_degree >= 1) and (p_degree >= 1) and (maxmulti >= 2)) {
						if (double_bond > 0) {
						  num += set_multi_and_next(tree, is_ident, k, 2, double_bond - 1, triple_bond,num_naph_bond, normal_type, output_file);
						}
						if ((triple_bond > 0) and (c_degree >= 2) and (p_degree >= 2) and (maxmulti >= 3)) {
						  num += set_multi_and_next(tree, is_ident, k, 3, double_bond, triple_bond - 1,num_naph_bond, normal_type, output_file);
						}
					}

				}else{ //tree.get_label(k)==0
					if(num_naph_bond>0 && c_degree>=1 && p_degree>=1 && maxmulti>=2){
					  num += set_multi_and_next(tree, is_ident, k, 2, double_bond,triple_bond,num_naph_bond-1, normal_type, output_file);
					}
				}
			}//end if(tree.can_be_added_multi(k))

			tree.set_multi_bond(k, 1);  
			tree.update_identical_multi(is_ident, k);
		}
	}
	
	return num;
}

inline size_t add_node_and_next(ChemTreeCenter& tree, is_ident_type is_ident, const int deepest_head, const int parenti, const label_value_type atom_label, std::ofstream & output_file )
{using namespace std;
	if (tree.add_node(is_ident, parenti, atom_label)) {
		const int type = tree.is_normal(deepest_head);
		if (type > 0){
		        if (num_lack_H > 0 || num_naph_bond>0) {
			
			  size_t result = generate_multi(tree, is_ident, type, output_file);
			  return result;
			} else {
			  return enum_benzene(tree, output_file);
			}
		} else {
			return 0;
		}
	}
	size_t num = 0;
	for (int i = tree.starti(); i < deepest_head; ++i){
		if (tree.can_be_added(i)) {
			const label_value_type begin_atom = tree.begin_atom_label(is_ident, i);
			for (label_value_type atomi = begin_atom; atomi != first_atom_valence_one; ++atomi) {
				if (tree.remain(atomi)) {				  
				  num += add_node_and_next(tree, is_ident, deepest_head, i, atomi, output_file);
				  tree.del_last_node();
				}
			}
		}
		tree.update_identical_end(is_ident, i);
	}
	if (tree.share_only_root(deepest_head)) {
		const int num_nodes = tree.get_num_nodes();
		for (int i = deepest_head; i < num_nodes; ++i) {
			if (tree.can_be_added(i)) {
				const label_value_type begin_atom = tree.begin_atom_label(is_ident, i);
				for (label_value_type atomi = begin_atom; atomi != first_atom_valence_one; ++atomi) {
					if (tree.remain(atomi)) {					  
					  num += add_node_and_next(tree, is_ident, num_nodes, i, atomi, output_file);
					  tree.del_last_node();
					}
				}
			}
			tree.update_identical_end(is_ident, i);
		}
	}
	return num;
}

inline size_t add_root_child_and_next(ChemTreeCenter& tree, is_ident_type is_ident, const label_value_type atom_label, std::ofstream & output_file)
{	using namespace std;
  
	if (tree.add_root_child(is_ident, atom_label)) {
		const int type = tree.is_normal(1);
		
		if (type > 0) {
		  if (num_lack_H > 0 || num_naph_bond>0) {	  
		    return generate_multi(tree, is_ident, type, output_file);
		  } else {
		    return enum_benzene(tree, output_file);
		  }
		} else {
			return 0;
		}
	}

	size_t num = 0;
	if (tree.can_be_added_root()) {
		const label_value_type begin_atom = tree.begin_atom_label_root();
		for (label_value_type atomi = begin_atom; atomi != first_atom_valence_one; ++atomi) {
			if (tree.remain(atomi)) {
			  num += add_root_child_and_next(tree, is_ident, atomi, output_file);
			  tree.del_last_node();
			}
		}
	}
	if (tree.share_only_root(1)) { 
		const int num_nodes = tree.get_num_nodes();
	       		for (int i = 1; i < num_nodes; ++i) {
			if (tree.can_be_added(i)) {
				const label_value_type begin_atom = tree.begin_atom_label(is_ident, i);
				for (label_value_type atomi = begin_atom; atomi != first_atom_valence_one; ++atomi) {
					if (tree.remain(atomi)) {
					  num += add_node_and_next(tree, is_ident, num_nodes, i, atomi, output_file);
						tree.del_last_node();
					}
				}
			}
			tree.update_identical_end(is_ident, i);

		}
	}
	
	return num;
}

inline size_t add_root_and_next(ChemTreeCenter& tree, const label_value_type atom_label, std::ofstream & output_file)
{
	if (tree.add_root(atom_label)) {
		if (do_print) {
		  tree.print_single();
		}
		return 1;
	} 

	is_ident_type is_ident(num_except_H);
	is_ident[0] = (char)(-1);

	size_t num = 0;
	for (label_value_type atomi = 0; atomi != first_atom_valence_one; ++atomi) {
		if(tree.remain(atomi)){
		  int temp = 0;
		  if(atomi==0 && tree.get_label(0)==0){
		    temp=1;
		  }
		  num += add_root_child_and_next(tree, is_ident, atomi, output_file);
		  tree.del_last_node(); 
		}
	}
	return num;
}

inline static void chk_H(const std::vector<int>& atom_numbers,const std::vector<int> &bond_numbers)
{
	using namespace std;

	int withoutH_deg = 0;
	int withoutH_size = 0;
	int H_size = 0;
	for(int i = 0;i != first_atom_valence_one;i++){
		withoutH_deg += atom_numbers[i]*valence[i];
		withoutH_size += atom_numbers[i];
	}
	const int needH = withoutH_deg - 2*(withoutH_size-1);
	for(int i = first_atom_valence_one;i != num_distinct_atoms; i++){
		H_size += atom_numbers[i];
	}

	num_except_H = withoutH_size;
//num_lack_H = needH-H_size;//-(bond_numbers[1]*2);
	num_lack_H = bond_numbers[0];
	num_naph_bond = bond_numbers[1];
	num_H = H_size;
	if (num_lack_H % 2 != 0) {
		cerr << "error # hydrogen atoms" << endl;
		exit(EXIT_FAILURE);
	}
}


inline size_t generate_from_tree( ChemTreeCenter &tree, const bool _do_print = false){
  using namespace std;

  do_print = _do_print;
  ofstream output_file;
  output_file.open("output.smi");
  
  size_t num = enum_benzene( tree, output_file );
  
  output_file.close();

  return num;
}

inline size_t generate(const std::vector<int>& atom_numbers, const bool _do_print = false, const size_t round_num = 0)
{
	using namespace std;
	do_print = _do_print;
	vector< vector<int> > input; //input including benzene
	vector< vector<int> > bond; //bond is vector size 2, first element is normal multi, second element is bond between 2 benzene rings
	int lack_valence=2;
	ofstream output_file;
	output_file.open("output.smi");

	for(size_t i=0;i<atom_numbers.size();i++){
		lack_valence+=input_valence[i]*atom_numbers[i];
		lack_valence-=2*atom_numbers[i];
	}

	vector<int> input_number(atom_numbers);
	for(int i=0;i<num_special_atom;i++){
	  input_number.insert(input_number.begin(),0);
	}
	vector<int> num_bond(1,lack_valence);
	num_bond.push_back(0);
	//insert into vector
	input.push_back(input_number);
	bond.push_back(num_bond);

	int index_b = std::find(atomchar, atomchar + num_distinct_atoms, "b") - atomchar;
	int index_c = std::find(atomchar, atomchar + num_distinct_atoms, "C") - atomchar;

	count_ring_benzene(index_b,index_c,num_bond,input_number,input,bond);
	
	//restore input number to the original one
	input_number.clear();
	input_number=atom_numbers;
	
	for(int i=0;i<num_special_atom;i++){
	  input_number.insert(input_number.begin(),0);
	}
	//restore temp_bond
	num_bond.clear();
	num_bond.push_back(lack_valence);
	num_bond.push_back(0);

	size_t num = 0;

	for(size_t i=0; i < input.size(); i++){
	  if(round_num > input.size()){
	    std::cout<<"#round is too large"<<std::endl;
	    return 0;
	  }
	  if(round_num > 0){
	    i = round_num - 1;
	  }
	        cout<<"==================="<<endl;
		cout<<"amount of naphthalene rings: "<<bond[i][1]<<endl;
		int num_b = input[i][0] - 2*bond[i][1];
		cout<<"amount of benzene rings: "<<num_b<<endl;
		for(size_t j=1; j<num_distinct_atoms ; j++){
			cout<<"amount of "<<atomchar[j]<<": "<<input[i][j]<<endl;
		}
		//cout<<"bond[0] ="<<bond[i][0]<<"  bond[1] ="<<bond[i][1]<<endl; 

		t_r.clear();
		t_v.clear();
		location.clear();
		lack_degree.clear();

		chk_H(input[i],bond[i]);
		t_r.reserve(num_except_H>>1);
		t_v.reserve(num_except_H>>1);
		location.reserve(num_except_H>>1);
		lack_degree.reserve(num_except_H);

		ChemTreeCenter tree(input[i], num_except_H);


		for (label_value_type atomi = 0; atomi != first_atom_valence_one; atomi++) {
			if (tree.remain(atomi)){	
			  num += add_root_and_next(tree, atomi, output_file);
			  tree.del_root();
			}
		}

		std::cout<<"accumulated #enumerated structure until "<<i+1<<" round is "<<num<<std::endl;
		if(round_num > 0){
		  break;
		}

	}
	output_file.close();
	return num;
}

inline void ChemTreeCenter::write_smiles(const int index, std::ofstream & file, int & num_cycle) const{
  if( index != 0 && nodes[index].label != 0 ){
    file << bondchar[ nodes[index].multi ];
  }

  if(nodes[index].label == 0){ //if it is a benzene node
    if(nodes[index].multi !=2 ){
      file << "c";
      
      int index_cycle = num_cycle;
      num_cycle++;

      file << index_cycle;

      int next_c = 1;
      if(index == 0){
	//index = 0 has no parent -> search for child node with carbon position 0
	int nth_cp0 = std::find(nodes[index].bond_position.begin(), nodes[index].bond_position.end(), 0) - nodes[index].bond_position.begin();

	if(nth_cp0 < nodes[index].num_children){
	  if(nodes[ nodes[index].children[nth_cp0] ].multi == 2){
	    next_c++;
	    write_smiles(nodes[index].children[nth_cp0], file, num_cycle);
	  }else{	  
	    file << "(";
	    write_smiles(nodes[index].children[nth_cp0], file, num_cycle);
	    file << ")";  
	  }
	}else{
	  std::cout<<"Error: no nodes with carbon position 0"<<std::endl;
	}
      }

      for(int i = next_c; i < 6; i++){
	file << "c";
	int nth_cpi = std::find(nodes[index].bond_position.begin(), nodes[index].bond_position.end(), i) - nodes[index].bond_position.begin();
	
	if(nth_cpi < nodes[index].num_children){
	  if(nodes[ nodes[index].children[nth_cpi] ].multi == 2){
	    i++;
	    write_smiles(nodes[index].children[nth_cpi], file, num_cycle);
	  }else{
	    file << "(";
	    write_smiles(nodes[index].children[nth_cpi], file, num_cycle);
	    file << ")";
	  }
	}
      }
      file << index_cycle;
    }else{ //node[index].multi == 2
      int index_cycle = num_cycle;
      num_cycle++;
 
      file << index_cycle;
      for(int i=1; i < 5; i++){
	file << "c";
	int nth_cpi = std::find(nodes[index].bond_position.begin(), nodes[index].bond_position.end(), i) - nodes[index].bond_position.begin();

	if(nth_cpi < nodes[index].num_children){
	  file << "(";
	  write_smiles(nodes[index].children[nth_cpi], file, num_cycle);
	  file << ")";
	}
      }
      file << "c" << index_cycle;
    }
    
  }else{ // if it is a normal node
    file << atomchar[ nodes[index].label ];

    const int num_child = nodes[index].num_children;
    
    if(num_child > 0){ 
      //all child nodes except the last one -> need "(" and ")"
      for(int i=0; i<num_child-1; i++){
	file << "(";
	write_smiles( nodes[index].children[i], file, num_cycle);
	file << ")";
      }
    
      //the last child node -> add label directly
      write_smiles( nodes[index].children[num_child-1], file, num_cycle);
    }
  }
}

inline void ChemTreeCenter::printsmi(const int i) const
{
	if (i != 0) {
		std::cout << bondchar[nodes[i].multi];
	}
	std::cout << atomchar[nodes[i].label];
	const valence_value_type nc = nodes[i].num_children;
	if (nc > 1) {
		for (valence_value_type v = 0; v != nc; ++v) {
			std::cout << "(";
			printsmi(nodes[i].children[v]);
			std::cout << ")";
		}
	} else if (nc == 1) {
		printsmi(nodes[i].children[0]);
	}
}

inline void ChemTreeCenter::printsmi_single(const int i) const
{
	std::cout << atomchar[nodes[i].label];
	const valence_value_type nc = nodes[i].num_children;
	if (nc > 1) {
		for (valence_value_type v = 0; v != nc; ++v) {
			std::cout << "(";
			printsmi_single(nodes[i].children[v]);
			std::cout << ")";
		}
	} else if (nc == 1) {
		printsmi_single(nodes[i].children[0]);
	}
}

inline void ChemTreeCenter::printseq() const
{
	using namespace std;

	cout << nodes[0].label << nodes[0].num_children;
	for (int i = 1; i!= num_nodes; ++i) {
		cout << nodes[i].label << nodes[i].multi << nodes[i].num_children;
	}
}

inline void ChemTreeCenter::printseq_single() const
{
	using namespace std;

	cout << nodes[0].label << nodes[0].num_children;
	for (int i = 1; i!= num_nodes; ++i) {
		cout << nodes[i].label << nodes[i].num_children;
	}
}

inline void ChemTreeCenter::print() const
{
	printsmi();
	std::cout << std::endl;
}

inline void ChemTreeCenter::print_single() const
{
	printsmi_single(0);
	std::cout << std::endl;
}


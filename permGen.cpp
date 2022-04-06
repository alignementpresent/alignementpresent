#include "permGen.hpp"

using namespace std;
typedef unsigned int uint;

string exec(const char* cmd){
//Exec cmd and grab the stdout
    array<char, 128> buffer;
    string result;
    shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}

vector<vector<int8_t>> getAllSubgraphIsomorphism(uint const indexCentralDiGraph,
												 bool const withSSB){
	//Get all possible isomorphisms from the specified CentralDigraph
	//There are 2027 such digraphs, indexed from 0 to 2026, the graph files are generated using the genCentralDiGraphFiles.sage script

	string cmd = "python preParseOutput.py graphNoSSB" + to_string(indexCentralDiGraph);
	if(withSSB)
		cmd = "python preParseOutput.py graph4444SSB0";
	istringstream s(exec(cmd.c_str()));
	
	vector<vector<int8_t>> v;
	string line;
	while(getline(s, line)){
		istringstream ss(line);
		vector<int8_t> vv(16);
		uint i = 0;
		string val;
		while(getline(ss, val, ',')){
			vv[i] = stoi(val);
			i++;
		}

		//Due to the structure of the graph, one value in the isomorphism isn't determined
		//Just need to find which value is missing
		set<int8_t> tmp;
		for(uint i = 0; i < 16; i++) tmp.emplace(i);
		for(auto const tmpval : vv){
			if(tmpval != -1) tmp.erase(tmpval);
		}
		if(tmp.size() > 1)
			cerr << "Error, more than one value is not determined, shouldn't happened" << endl;
		else{
			for(uint i = 0; i < 16; i++){
				if(vv[i] == -1)
					vv[i] = *(tmp.begin());
			}
		}


		v.emplace_back(move(vv));
	}
	return v;
}

vector<vector<uint8_t>> getGraphFromFile(uint const indexCentralDiGraph,
										 bool const withSSB){
	//Get the graph of index indexCentralDiGraph as an adjacency list
	//i.e. g[i] contains j iif there is an edge i->j

	vector<vector<uint8_t>> g(16);
	string filename = "./graphFiles/graphNoSSB"+to_string(indexCentralDiGraph);
	if(withSSB)
		filename = "./graphFiles/graph4444SSB0";
	ifstream infile(filename);
	string line;
	while(getline(infile, line)){
		istringstream ss(line);
		string val;
		vector<uint8_t> vv(2);
		uint i = 0;
		while(getline(ss, val, '>')){
			vv[i] = stoi(val);
			i++;
		}
		g[vv[0]].emplace_back(vv[1]);
	}
	return g;
}

void completePermutation(vector<int8_t> pp,
						 vector<vector<uint8_t>> const & g,
						 uint const nbPerm,
						 vector<vector<int8_t>> & listPerm,
						 bool const randomize){

	//Get the index of the first element that is not determined
	uint index = 0;
	while(index < 64){
		if(pp[index] == -1) 
			break;
		index++;
	}

	if(index == 64){ //Permutation is determined
		listPerm.emplace_back(pp);
	}
	else{
		uint i = index/4; //Index of the sbox
		// uint j = index%4; //Index of the bit in the sbox

		//Get the list of sboxes remaining to hit
		set<uint8_t> targetsbox(g[i].begin(), g[i].end());
		for(uint k = 0; k < 4; k++){
			if(pp[4*i+k] != -1)
				targetsbox.erase(pp[4*i+k]/4);
		}

	    std::random_device rd;
    	std::mt19937 rng(rd());
		vector<uint8_t> vsbox(targetsbox.begin(), targetsbox.end());
		if(randomize)
			shuffle(vsbox.begin(), vsbox.end(), rng);

		for(auto const isbox : vsbox){
			//Get the list of values in the target sbox that are not yet hit
			vector<uint8_t> targetvalues;
			for(uint k = 0; k < 4; k++){
				if(find(pp.begin(), pp.end(), 4*isbox+k) == pp.end())
					targetvalues.emplace_back(4*isbox+k);
			}
			for(auto const val : targetvalues){
				pp[index] = val;
				completePermutation(pp,g,nbPerm,listPerm);
				if(listPerm.size() >= nbPerm)
					break;
			}
			if(listPerm.size() >= nbPerm)
				break;
		}
	}
}

vector<vector<uint8_t>> genPermutationsFromIsomorphism(vector<int8_t> const & partialp,
													   vector<int8_t> const & iso,
													   vector<vector<uint8_t>> const & g,
													   uint const nbPerm,
													   bool const randomize){

	//Apply the isomoprhism to the partial permutation
	//if p[4i+j] = 4a+b then p'[4*iso[i]+j] = 4*iso[a]+b
	vector<int8_t> pp(64,-1);
	for(uint i = 0; i < 16; i++){
		for(uint j = 0; j < 4; j++){
			if(partialp[4*i+j] != -1){
				uint a = partialp[4*i+j]/4;
				uint b = partialp[4*i+j]%4;
				pp[4*iso[i]+j] = 4*iso[a] + b;
			}
		}
	}

	//Complete the permutation to generate up to #nbPerm permutations
	vector<vector<int8_t>> listPerm;
	listPerm.reserve(nbPerm);
	completePermutation(pp,g,nbPerm,listPerm,randomize);

	//Apply the inverse isomoprhism to get the final list of permutations
	//if p[4i+j] = 4a+b then p'[4*inviso[i]+j] = 4*inviso[a]
	vector<int8_t> inviso(16);
	for(uint i = 0; i < 16; i++)
		inviso[iso[i]] = i;

	vector<vector<uint8_t>> listFinalPerm(listPerm.size());
	for(uint iperm = 0; iperm < listPerm.size(); iperm++){
		vector<uint8_t> finalp(64);
		auto & pp_i = listPerm[iperm];

		for(uint i = 0; i < 16; i++){
			for(uint j = 0; j < 4; j++){
				uint a = pp_i[4*i+j]/4;
				uint b = pp_i[4*i+j]%4;
				finalp[4*inviso[i]+j] = 4*inviso[a] + b;
			}
		}
		listFinalPerm[iperm] = move(finalp);
	}
	return listFinalPerm;
}


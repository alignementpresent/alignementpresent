#include "aux_function.hpp"

using namespace std;

vector<vector<uint8_t>> getNewPermutations(uint const nbPermPerIso,
						uint const nbIso,
						bool const withSSB,
						bool const randomize,
						uint const MaxNbGoodGraph){
	//For each graph, check nbIso isomorphisms and build nbPermPerIso permutations for each isomoprhism
	//withSSB == true : same base graph as the original present permutation
	//withSSB == false : build permutations for each graph leading to no SSB
	//if randomizeIso == false, just pick the first isomorphism(s), and build the permutation deterministically
	//Otherwise pick the isomorphisms at random, and build the permutations in a random order

	std::random_device rd;
 	std::mt19937 rng(rd());

	//Generate the graph files
	// exec("sage genCentralDiGraphFiles.sage");

	//Original PRESENT permutation
	vector<uint8_t> p({0,16,32,48,1,17,33,49,2,18,34,50,3,19,35,51,4,20,36,52,5,21,37,53,6,22,38,54,7,23,39,55,8,24,40,56,9,25,41,57,10,26,42,58,11,27,43,59,12,28,44,60,13,29,45,61,14,30,46,62,15,31,47,63});
	//Sbox
	vector<uint8_t> s({0xc,0x5,0x6,0xb,0x9,0x0,0xa,0xd,0x3,0xe,0xf,0x8,0x4,0x7,0x1,0x2});
	auto ddt = computeDDT(s,4);

	//Generate the Sbox connectivity matrix for the partial permutations where values for the LSB at the input/output of the sboxes are not determined
	//Also generate the partial permutation, where undefined vales are replaced by -1
	vector<int8_t> partialp(64);
	vector<vector<uint8_t>> M(16, vector<uint8_t>(16,0));
	for(uint i = 0; i < 64; i++){
		if(i%4 != 0 && p[i]%4 != 0){
			M[i/4][p[i]/4] = 1;
			partialp[i] = p[i];
		}
		else
			partialp[i] = -1;
	}

	//Generate the graph file for M
	ofstream partialGraphFile("partialGraph");
	for(uint i = 0; i < 16; i++){
		for(uint j = 0; j < 16; j++){
			if(M[i][j] != 0)
				partialGraphFile << i << ">" << j << endl;
		}
	}
	partialGraphFile.close();

	vector<vector<uint8_t>> allPerms;

	if(withSSB){
		//There is only 1 DiGraph in this case
		auto allIso = getAllSubgraphIsomorphism(0,true);
		auto graph = getGraphFromFile(0,true);
		if(randomize){
			shuffle(allIso.begin(), allIso.end(), rng);
		}
		//Generate permutations for the first nbIso isomoprhisms
		for(uint indexIso = 0; indexIso < nbIso; indexIso++){
			auto listNewPerm = genPermutationsFromIsomorphism(partialp,allIso[indexIso],graph,nbPermPerIso,randomize);
			allPerms.insert(allPerms.end(), listNewPerm.begin(), listNewPerm.end());
		}
	}
	else{
		//There are 2027 digraphs to test
		uint nbGoodGraph = 0;
		for(uint iDigraph = 0; iDigraph < 2027; iDigraph++){
			auto allIso = getAllSubgraphIsomorphism(iDigraph);
			if(allIso.size() > 0){ //If this graph can generate satisfying permutations
				auto graph = getGraphFromFile(iDigraph);
				nbGoodGraph++;

				if(randomize){
					shuffle(allIso.begin(), allIso.end(), rng);
				}

				//Generate permutations for the first nbIso isomoprhisms
				for(uint indexIso = 0; indexIso < nbIso; indexIso++){
					auto listNewPerm = genPermutationsFromIsomorphism(partialp,allIso[indexIso],graph,nbPermPerIso,randomize);
					allPerms.insert(allPerms.end(), listNewPerm.begin(), listNewPerm.end());
				}
			}
			if(MaxNbGoodGraph > 0 && nbGoodGraph == MaxNbGoodGraph) break;
		}
	}

	return allPerms;
}

void runComputations(HistoType const htype,
					 uint const nbRound,
					 uint const nbSboxMax,
					 vector<vector<uint8_t>> const & allPerms,
					 string const & filename){

	vector<uint8_t> ogp({0,16,32,48,1,17,33,49,2,18,34,50,3,19,35,51,4,20,36,52,5,21,37,53,6,22,38,54,7,23,39,55,8,24,40,56,9,25,41,57,10,26,42,58,11,27,43,59,12,28,44,60,13,29,45,61,14,30,46,62,15,31,47,63}); //original perm
	vector<uint8_t> s({0xc,0x5,0x6,0xb,0x9,0x0,0xa,0xd,0x3,0xe,0xf,0x8,0x4,0x7,0x1,0x2});
	auto ddt = computeDDT(s,4);

	if(htype == HistoType::Trail){
		vector<map<double, uint64_t>> allHisto;
		map<double, uint64_t> histoOG;
		if(nbRound == 2){
			for(auto const & p : allPerms){
				PresentData PD(p,ddt);
				allHisto.emplace_back(PD.trailHistogram2R(nbSboxMax));
			}
			PresentData PD(ogp,ddt);
			histoOG = PD.trailHistogram2R(nbSboxMax);
		}
		if(nbRound == 3){
			for(auto const & p : allPerms){
				PresentData PD(p,ddt);
				allHisto.emplace_back(PD.trailHistogram3R(nbSboxMax));
			}
			PresentData PD(ogp,ddt);
			histoOG = PD.trailHistogram3R(nbSboxMax);
		}
		else if(nbRound == 4){
			for(auto const & p : allPerms){
				PresentData PD(p,ddt);
				allHisto.emplace_back(PD.trailHistogram4R(nbSboxMax));
			}
			PresentData PD(ogp,ddt);
			histoOG = PD.trailHistogram4R(nbSboxMax);
		}
		logHistogram(allHisto, histoOG, filename);
	}

	else if(htype == HistoType::BoxWeight){
		vector<map<uint64_t, uint64_t>> allHisto;
		map<uint64_t, uint64_t> histoOG;
		if(nbRound == 2){
			for(auto const & p : allPerms){
				PresentData PD(p,ddt);
				allHisto.emplace_back(PD.boxWeightHistogram2R(nbSboxMax));
			}
			PresentData PD(ogp,ddt);
			histoOG = PD.boxWeightHistogram2R(nbSboxMax);
		}
		if(nbRound == 3){
			for(auto const & p : allPerms){
				PresentData PD(p,ddt);
				allHisto.emplace_back(PD.boxWeightHistogram3R(nbSboxMax));
			}
			PresentData PD(ogp,ddt);
			histoOG = PD.boxWeightHistogram3R(nbSboxMax);
		}
		else if(nbRound == 4){
			for(auto const & p : allPerms){
				PresentData PD(p,ddt);
				allHisto.emplace_back(PD.boxWeightHistogram4R(nbSboxMax));
			}
			PresentData PD(ogp,ddt);
			histoOG = PD.boxWeightHistogram4R(nbSboxMax);
		}
		logHistogram(allHisto, histoOG, filename);
	}

	else if(htype == HistoType::CorePattern){
		vector<map<uint64_t, uint64_t>> allHisto;
		map<uint64_t, uint64_t> histoOG;
		if(nbRound == 2){
			for(auto const & p : allPerms){
				PresentData PD(p,ddt);
				allHisto.emplace_back(PD.corePatternHistogram2R(nbSboxMax));
			}
			PresentData PD(ogp,ddt);
			histoOG = PD.corePatternHistogram2R(nbSboxMax);
		}
		if(nbRound == 3){
			for(auto const & p : allPerms){
				PresentData PD(p,ddt);
				allHisto.emplace_back(PD.corePatternHistogram3R(nbSboxMax));
			}
			PresentData PD(ogp,ddt);
			histoOG = PD.corePatternHistogram3R(nbSboxMax);
		}
		else if(nbRound == 4){
			for(auto const & p : allPerms){
				PresentData PD(p,ddt);
				allHisto.emplace_back(PD.corePatternHistogram4R(nbSboxMax));
			}
			PresentData PD(ogp,ddt);
			histoOG = PD.corePatternHistogram4R(nbSboxMax);
		}
		logHistogram(allHisto, histoOG, filename);
	}

	else if(htype == HistoType::ClusterTrail){
		vector<map<double, uint64_t>> allHisto;
		map<double, uint64_t> histoOG;
		if(nbRound == 2){
			for(auto const & p : allPerms){
				PresentData PD(p,ddt);
				allHisto.emplace_back(PD.trailClusterHistogram2R(nbSboxMax));
			}
			PresentData PD(ogp,ddt);
			histoOG = PD.trailClusterHistogram2R(nbSboxMax);
		}
		if(nbRound == 3){
			for(auto const & p : allPerms){
				PresentData PD(p,ddt);
				allHisto.emplace_back(PD.trailClusterHistogram3R(nbSboxMax));
			}
			PresentData PD(ogp,ddt);
			histoOG = PD.trailClusterHistogram3R(nbSboxMax);
		}

		logCumulativeHistogram(allHisto, histoOG, filename);
	}
}

void timingNbSbox(HistoType const htype,
				  uint const nbRound,
				  uint const minSbox){

	vector<uint8_t> ogp({0,16,32,48,1,17,33,49,2,18,34,50,3,19,35,51,4,20,36,52,5,21,37,53,6,22,38,54,7,23,39,55,8,24,40,56,9,25,41,57,10,26,42,58,11,27,43,59,12,28,44,60,13,29,45,61,14,30,46,62,15,31,47,63}); //original perm
	vector<uint8_t> s({0xc,0x5,0x6,0xb,0x9,0x0,0xa,0xd,0x3,0xe,0xf,0x8,0x4,0x7,0x1,0x2});
	auto ddt = computeDDT(s,4);

	PresentData PD(ogp,ddt);

	for(uint nbSboxMax = minSbox; nbSboxMax < 25; nbSboxMax++){
		auto t1 = chrono::high_resolution_clock::now();

		if(htype == HistoType::Trail){
			if(nbRound == 2){
				cout << "Trails for 2 rounds with " << nbSboxMax << " max sboxes" << endl;
				auto histo = PD.trailHistogram2R(nbSboxMax);
				printHisto(histo);
			}
			else if(nbRound == 3){
				cout << "Trails for 3 rounds with " << nbSboxMax << " max sboxes" << endl;
				auto histo = PD.trailHistogram3R(nbSboxMax);
				printHisto(histo);
			}
			else if(nbRound == 4){
				cout << "Trails for 4 rounds with " << nbSboxMax << " max sboxes" << endl;
				auto histo = PD.trailHistogram4R(nbSboxMax);
				printHisto(histo);
			}
		}
		else if(htype == HistoType::BoxWeight){
			if(nbRound == 2){
				cout << "Box Weight for 2 rounds with " << nbSboxMax << " max sboxes" << endl;
				auto histo = PD.boxWeightHistogram2R(nbSboxMax);
				printHisto(histo);
			}
			else if(nbRound == 3){
				cout << "Box Weight for 3 rounds with " << nbSboxMax << " max sboxes" << endl;
				auto histo = PD.boxWeightHistogram3R(nbSboxMax);
				printHisto(histo);
			}
			else if(nbRound == 4){
				cout << "Box Weight for 4 rounds with " << nbSboxMax << " max sboxes" << endl;
				auto histo = PD.boxWeightHistogram4R(nbSboxMax);
				printHisto(histo);
			}
		}
		else if(htype == HistoType::CorePattern){
			if(nbRound == 2){
				cout << "Core Patterns for 2 rounds with " << nbSboxMax << " max sboxes" << endl;
				auto histo = PD.corePatternHistogram2R(nbSboxMax);
				printHisto(histo);
			}
			else if(nbRound == 3){
				cout << "Core Patterns for 3 rounds with " << nbSboxMax << " max sboxes" << endl;
				auto histo = PD.corePatternHistogram3R(nbSboxMax);
				printHisto(histo);
			}
			else if(nbRound == 4){
				cout << "Core Patterns for 4 rounds with " << nbSboxMax << " max sboxes" << endl;
				auto histo = PD.corePatternHistogram4R(nbSboxMax);
				printHisto(histo);
			}
		}
		else if(htype == HistoType::ClusterTrail){
			if(nbRound == 2){
				cout << "ClusterTrail for 2 rounds with " << nbSboxMax << " max sboxes" << endl;
				auto histo = PD.trailClusterHistogram2R(nbSboxMax);
				printHisto(histo);
			}
			else if(nbRound == 3){
				cout << "ClusterTrail for 3 rounds with " << nbSboxMax << " max sboxes" << endl;
				auto histo = PD.trailClusterHistogram3R(nbSboxMax);
				printHisto(histo);
			}
		}

		auto t2 = chrono::high_resolution_clock::now();
		cout << "Run time for nbSboxMax = " << nbSboxMax << " : " << chrono::duration_cast<chrono::seconds>(t2 - t1).count() << " s" << endl << endl;
	}
}

void genData(HistoType const htype,
			 uint const nbRound, 
			 uint const nbSboxMax, 
			 uint const nbPermPerIso,
			 uint const nbIso,
			 bool const withSSB,
			 bool const randomize,
			 uint const MaxNbGoodGraph){

	auto allPerms = getNewPermutations(nbPermPerIso, nbIso, withSSB, randomize, MaxNbGoodGraph);
	string filename = "./results/data";
	if(htype == HistoType::Trail)
		filename += "_Trail";
	else if(htype == HistoType::BoxWeight)
		filename += "_BoxWeight";
	else if(htype == HistoType::CorePattern)
		filename += "_CorePattern";
	else if(htype == HistoType::ClusterTrail)
		filename += "_ClusterTrail";
	filename += "_" + to_string(nbRound) + "r";
	filename += "_bound" + to_string(nbSboxMax);
	filename += "_" + to_string(nbPermPerIso) + "ppiso";
	filename += "_" + to_string(nbIso) + "iso";
	if(MaxNbGoodGraph == 0)
		filename += "_allGraph";
	else
		filename += "_" + to_string(MaxNbGoodGraph) + "Graph";
	if(withSSB)
		filename += "_SSB";
	else
		filename += "_noSSB";

	runComputations(htype, nbRound, nbSboxMax, allPerms, filename);
}

#ifndef H_PERMGEN
#define H_PERMGEN

#include <iostream>
#include <vector>
#include <cstdint>
#include <string>
#include <fstream>
#include <array>
#include <memory>
#include <sstream>
#include <set>
#include <map>
#include <algorithm>
#include <random>

typedef unsigned int uint;

std::string exec(const char* cmd);

std::vector<std::vector<int8_t>> getAllSubgraphIsomorphism(uint const indexCentralDiGraph,
														   bool const withSSB = false);
	//Get all possible isomorphisms from the specified CentralDigraph
	//There are 2027 such digraphs, indexed from 0 to 2026, the graph files are generated using the genCentralDiGraphFiles.sage script
	//if withSSB == true, instead just get the isomorphism from the (single) graph with SSB pattern 4444

std::vector<std::vector<uint8_t>> getGraphFromFile(uint const indexCentralDiGraph,
												   bool const withSSB = false);
	//Get the graph of index indexCentralDiGraph as an adjacency list
	//i.e. g[i] contains j iif there is an edge i->j
	//if witthSSB == true, instead just get the graph with SSB pattern 4444

void completePermutation(std::vector<int8_t> pp,
						 std::vector<std::vector<uint8_t>> const & g,
						 uint const nbPerm,
						 std::vector<std::vector<int8_t>> & listPerm,
						 bool const randomize = false);

std::vector<std::vector<uint8_t>> 
genPermutationsFromIsomorphism(std::vector<int8_t> const & partialp,
							   std::vector<int8_t> const & iso,
							   std::vector<std::vector<uint8_t>> const & g,
							   uint const nbPerm,
							   bool const randomize = false);




#endif
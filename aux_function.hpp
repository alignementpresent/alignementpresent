#ifndef H_AUXFUNCTION
#define H_AUXFUNCTION

#include <vector>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <cmath>
#include <map>
#include <unordered_map>
#include <iterator>
#include <random>
#include <chrono>
#include <set>
#include <fstream>

#include "PresentData.hpp"
#include "permGen.hpp"

enum class HistoType{
	Trail=0,
	BoxWeight,
	CorePattern,
	ClusterTrail
};
typedef unsigned int uint;

std::vector<std::vector<uint8_t>> 
getNewPermutations(uint const nbPermPerIso,
				   uint const nbIso,
				   bool const withSSB,
				   bool const randomize,
				   uint const MaxNbGoodGraph = 0);

void runComputations(HistoType const htype,
					 uint const nbRound, 
					 uint const nbSboxMax, 
					 std::vector<std::vector<uint8_t>> const & allPerms,
					 std::string const & filename);

void timingNbSbox(HistoType const htype,
				  uint const nbRound,
				  uint const minSbox);

void genData(HistoType const htype,
			 uint const nbRound, 
			 uint const nbSboxMax, 
			 uint const nbPermPerIso,
			 uint const nbIso,
			 bool const withSSB,
			 bool const randomize,
			 uint const MaxNbGoodGraph);

template <class T>
std::set<T> 
getAllWeight(std::vector<std::map<T, uint64_t>> const & allHisto,
			 std::map<T, uint64_t> const & histoOG){
	std::set<T> allWeight;
	for(auto const & histo : allHisto){
		for(auto const & p_ctr : histo)
			allWeight.emplace(p_ctr.first);
	}
	//Also checks the weight from the OG perm
	for(auto const & p_ctr : histoOG)
		allWeight.emplace(p_ctr.first);
	return allWeight;
}

template <class T>
std::map<T, uint64_t> 
smoothHistogram(std::map<T, uint64_t> const & histo,
				std::set<T> const & allWeight){
	std::map<T, uint64_t> smoothHisto;
	for(auto const w : allWeight){
		if(histo.find(w) != histo.end()){
			smoothHisto[w] = histo.at(w);
		}
		else{
			smoothHisto[w] = 0;
		}
	}
	return smoothHisto;
}

template <class T>
std::map<T, uint64_t>
cumuSumHisto(std::map<T, uint64_t> const & histo){
	std::map<T, uint64_t> cumuHisto;
	uint64_t sum = 0;
	for(auto const & x : histo){
		sum += x.second;
		cumuHisto[x.first] = sum;
	}
	return cumuHisto;
}

template <class T>
void printHisto(std::map<T, uint64_t> const & histo,
				std::ostream & os = std::cout){
	for(auto const & x : histo){
		os << "[" << x.first << "," << x.second << "],";
	}
	os << std::endl;
}

template <class T>
void logInFile(std::set<T> const & allWeight,
			   std::vector<std::map<T, uint64_t>> const & allHisto,
			   std::map<T, uint64_t> const & histoOG,
			   std::string const & filename){
	//Output
	std::ofstream outfile;
	outfile.open(filename, std::ios::out);
	//First write the weights
	for(auto const w : allWeight)
		outfile << w << ",";
	outfile << std::endl;
	//Then the histograms
	for(auto const & histo : allHisto){
		for(auto const & p_ctr : histo){
			outfile << p_ctr.second << ",";
		}
		outfile << std::endl;
	}
	outfile.close();

	//Output for the OG permutation
	std::ofstream outfileOG;
	outfileOG.open(filename+"_OG", std::ios::out);
	//First write the weights
	for(auto const w : allWeight)
		outfileOG << w << ",";
	outfileOG << std::endl;
	//Then the histograms
	for(auto const & p_ctr : histoOG){
		outfileOG << p_ctr.second << ",";
	}
	outfileOG << std::endl;
	outfileOG.close();
}

template <class T>
void logHistogram(std::vector<std::map<T, uint64_t>> const & allHisto,
				  std::map<T, uint64_t> const & histoOG,
				  std::string const & filename){
	//Output the data in two files filename and filename+"_OG"
	//The histograms are "smoothed" to fit a panda dataframe csv
	//Meaning if a given histogram has no value associated with a specific key, it's set to 0
	//This is very unlikely to happen for Trails, BoxWeight and CorePatterns, but just to be safe and have proper format output
	//For ClusterTrail, it's better to use the version that output the cumulative histogram

	//Get all weight
	auto allWeight = getAllWeight(allHisto, histoOG);
	//Smooth the histograms
	std::vector<std::map<T, uint64_t>> allHistoSmooth(allHisto.size());
	for(uint indexHisto = 0; indexHisto < allHisto.size(); indexHisto++)
		allHistoSmooth[indexHisto] = smoothHistogram(allHisto[indexHisto],allWeight);
	
	auto smoothOGHisto = smoothHistogram(histoOG,allWeight);

	logInFile(allWeight, allHistoSmooth, smoothOGHisto, filename);
}

template <class T>
void logCumulativeHistogram(std::vector<std::map<T, uint64_t>> const & allHisto,
							std::map<T, uint64_t> const & histoOG,
							std::string const & filename){
	//Output the data in two files filename and filename+"_OG"
	//This output the cumulative histograms

	//Get all weight
	auto allWeight = getAllWeight(allHisto, histoOG);
	//Smooth the histograms
	std::vector<std::map<T, uint64_t>> allHistoSmooth(allHisto.size());
	for(uint indexHisto = 0; indexHisto < allHisto.size(); indexHisto++)
		allHistoSmooth[indexHisto] = smoothHistogram(allHisto[indexHisto],allWeight);
	
	auto smoothOGHisto = smoothHistogram(histoOG,allWeight);

	//Get the cumulative histograms
	std::vector<std::map<double, uint64_t>> allHistoCumu(allHistoSmooth.size());
	for(uint indexHisto = 0; indexHisto < allHisto.size(); indexHisto++){
		allHistoCumu[indexHisto] = cumuSumHisto(allHistoSmooth[indexHisto]);
	}
	auto cumuOGHisto = cumuSumHisto(smoothOGHisto);

	logInFile(allWeight, allHistoCumu, cumuOGHisto, filename);
}



#endif
#ifndef H_PRESENTDATA
#define H_PRESENTDATA

#include <array>
#include <vector>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <map>
#include <unordered_map>
#include <omp.h>
#include <set>
#include <unordered_set>
#include <chrono>
#include <algorithm>
#include <random>

/*
	Truncated trails are weird for Present as you actually need to know the actual diffrential to propagate through the permutation.
	Thus we can see two definition for the histogram.
	The first one (boxWeightHistogram* functions) counts how many trails there are with a given number of active sboxes, disregarding the probability
	However it includes the ones resulting from sbox layer convolution, i.e.
	0x1 -S-> 0x2 -P-> 0x3 -S-> 0x4 and 0x5 -S-> 0x2 -P-> 0x3 -S-> 0x6
	count as *2* separate trails
	The second definition for the histogram (corePatternHistogram* functions) focus on trail cores, in which the two trails above would only count as one (for the trail core 0x2 -P-> 0x3).
*/

typedef unsigned int uint;

//Boost inspired hasher
struct pairHasher {
    std::size_t operator () (std::pair<uint64_t,uint64_t> const & k) const;
};
using mapStateHisto = std::unordered_map<uint64_t, std::map<double,uint64_t>>;
using mapPairStateHisto = std::unordered_map<std::pair<uint64_t,uint64_t>, std::map<double,uint64_t>, pairHasher>;

inline uint64_t getNibble(uint64_t const x, uint const i){ //Get the value of the i-th nibble of x
    return (x >> (4*i)) & 0xFULL;
}

inline uint64_t setNibble(uint64_t const x, uint const i, uint64_t const v){ //Set the i-th nibble of x to v, while keeping the other bits unchanged. If v > 16, force the use of the 4 LSB 
	return (x & ~(0xFULL << (4*i))) | ((v & 0xFULL) << (4*i));
}

class PresentData{

public:

	PresentData(std::vector<uint8_t> const & xp,
				std::vector<std::vector<double>> const & xddt);

	std::map<double, uint64_t> trailHistogram2R(unsigned int const nbSboxMax);
	std::map<double, uint64_t> trailHistogram3R(unsigned int const nbSboxMax);
	std::map<double, uint64_t> trailHistogram4R(uint const nbSboxMax);
	//Return a map m such that m[p] = the number of trails with weight p, i.e. proba 2^(-p) over 2/3/4 rounds
	//Limit the computations to trails with at most #nbSboxMax active sboxes

	std::map<uint64_t, uint64_t> boxWeightHistogram2R(uint const nbSboxMax);
	std::map<uint64_t, uint64_t> boxWeightHistogram3R(uint const nbSboxMax);
	std::map<uint64_t, uint64_t> boxWeightHistogram4R(uint const nbSboxMax);
	//Return a map m such that m[x] = the number of trails with x active sboxes over 2/3/4 rounds
	//Limit the computations to trails with at most #nbSboxMax active sboxes

	std::map<uint64_t, uint64_t> corePatternHistogram2R(uint const nbSboxMax);
	std::map<uint64_t, uint64_t> corePatternHistogram3R(uint const nbSboxMax);
	std::map<uint64_t, uint64_t> corePatternHistogram4R(uint const nbSboxMax);
	//Return a map m such that m[x] = the number of core patterns with x active sboxes over 2 rounds
	//Limit the computations to trails with at most #nbSboxMax active sboxes

	std::map<uint64_t, std::map<uint64_t, uint64_t>> clusterHistogram(uint const nbSboxMax);
	std::map<double, uint64_t> trailClusterHistogram2R(uint const nbSboxMax);
	std::map<double, uint64_t> trailClusterHistogram3R(uint const nbSboxMax);
	//Return the histogram for trail clustering
	//Limit the computations to trails with at most #nbSboxMax active sboxes

	uint64_t applyPerm(uint64_t x);
	std::map<double,uint64_t> convolutionRound(uint64_t const x,
											   double const boundProba);

private:

	std::vector<uint8_t> p;
	std::vector<uint8_t> invp;
	std::vector<std::vector<double>> ddt;
	std::array<uint64_t,16> sboxMasks;
	std::array<uint64_t,16> pSboxMasks;
	std::array<uint64_t,16> invpSboxMasks;
	std::vector<std::vector<uint16_t>> uint16WeightMap; //uint16WeightMap[x] contains all uint16_t of hamming weight x
	std::vector<std::map<double,uint64_t>> sboxTrails; //sboxTrails[x][p] contains the number of transitions with probability 2^-p != 1 for S
	std::vector<std::map<double,uint64_t>> sboxInvTrails; //sboxInvTrails[x][p] contains the number of transitions with probability 2^-p != 1 for S^-1
	std::vector<std::vector<uint64_t>> ddtPropagation;
	std::vector<std::vector<uint64_t>> ddtInvPropagation;
	std::vector<double> minProbaInvSbox;
	double minWeightSbox; //Minimum non zero weight in the DDT
	std::vector<std::vector<bool>> isValidSboxIOSubspace;
	//isValidSboxIOSubspace[bx][by] contains whether or not using bx/by as input/output basis for differentials over one sbox leads to at least one valid differential
	//e.g. bx = 0b0011 ; by = 0b0001
	//If there is at least one valid transition {0b0001, 0b0010, 0b0011} -> {0b0001} then [bx][by] = true, otherwise false
	std::vector<std::vector<double>> minWeightIOSubspace;
	//minWeightIOSubspace[bx][by] contains the minimum weight over all valid transitions made from using bx/by as basiss for the input/output differentials
	//Equal to "16" on impossible basis transitions


	//------------------------------
	//-- Generic useful functions --
	//------------------------------
	// uint64_t applyPerm(uint64_t x);
	//Return p(x), where x represents the state as 16 uint64_t (only the 4 LSB are relevant)
	//This is to avoid having to unnecessarily cast the state to a single uint64_t
	uint64_t applyInvPerm(uint64_t x);
	//Return invp(x), where p is applied ot the bits of x
	//index 0 is the LSB of x, i.e. x & 1.

	uint16_t getActivityPattern(uint64_t const x);
	//Return an uint16_t representing which Sboxes are active in x
	//i-th bit = 1 <-> i-th Sbox active

	unsigned int countActiveSbox(uint64_t const x);
	//Return the number of active Sboxes in x

	// std::map<double,uint64_t> convolutionRound(uint64_t const x,
	// 										   double const boundProba);
	//Return the histogram for all trails over a single Sbox Layer with input differential x
	//Ignore any trail with a weight over boundWeight
	std::map<double,uint64_t> convolutionInvRound(uint64_t const x,
											      double const boundProba);
	//Return the histogram for all trails over a single Sbox Layer with output differential x
	//Ignore any trail with a weight over boundWeight

	bool checkCompatibleBasis(uint64_t bx, uint64_t by);
	//Return true is bx and by are compatible basis, meaning that at least one valid differential can be generated from them

	double getMinWeightBasis(uint64_t bx, uint64_t by);
	//Return the minimal weight of a differential made from the basis pair bx/by
	//Doesn't check if the basis is compatible, should be done beforehand. Non-comptatible basis will add a factor of 16 to the weight.
	
	//-----------------------------------------------------
	//-- Internals for the trail histogram computations --
	//-----------------------------------------------------
	void exhaustTrails2R(uint64_t y1,
						 uint64_t x2,
						 uint16_t patterny1,
						 uint const nbSboxR1,
						 unsigned int const nbSboxMax,
						 double const boundProba,
						 double currentMinProbay1,
						 std::map<double,uint64_t> & histo);
	//Internal recursive function to search trails over 2 rounds

	void enumerateCorePatterns4R(std::array<uint64_t,4> const & bx2_ay1_by2_nbSbox,
								 uint64_t const ax3,
								 std::map<double, uint64_t> & histo,
								 uint const nbSboxMax,
								 double const boundWeight);
	void enumeratePartialCorePatterns3R(std::array<std::map<uint64_t, std::vector<uint64_t>>, 16> const & Ldiff,
										uint64_t const ay1,
										uint64_t const ax3,
										uint64_t const nbSboxR123,
										std::map<double, uint64_t> & histo,
										uint const currentIndex_x2,
										uint const currentIndex_y2,
										uint64_t x2,
										uint64_t y2,
										uint const nbSboxMax,
										double currentWeight,
										double const minWeightR1R3,
										double const boundWeight);
	void extendCore4r(uint64_t const x3,
					  uint64_t y3,
					  uint16_t patterny3,
					  uint64_t x4,
					  uint64_t const nbSboxR123,
					  uint64_t const nbSboxMax,
					  std::map<double, uint64_t> & histo,
					  std::map<double, uint64_t> const & histoR1,
					  double currentWeight,
					  double const boundWeight);

	//--------------------------------------------
	//-- Core Pattern and box weight histograms --
	//--------------------------------------------
	//Computations are extremely similar, so basically just the same functions with a flag to know what to compute
	std::map<uint64_t, uint64_t> corePatternHistogram2R_internal(uint const nbSboxMax,
																 bool const boxWeightHisto);
	std::map<uint64_t, uint64_t> corePatternHistogram3R_internal(uint const nbSboxMax,
																 bool const boxWeightHisto);
	std::map<uint64_t, uint64_t> corePatternHistogram4R_internal(uint const nbSboxMax,
																 bool const boxWeightHisto);
	// 2 rounds
	void corePatternExhaustTrails2R(uint64_t y1,
									uint64_t x2,
									uint16_t patterny1,
									uint const nbsboxR1,
									uint const nbSboxMax,
									std::map<uint64_t,uint64_t> & histo,
									bool const boxWeightHisto);

	// 3 rounds 
	void enumerateCorePatterns3R(std::array<uint64_t,4> const & bx2_ay1_by2_nbSbox,
								 uint64_t const ax3,
								 std::map<uint64_t, uint64_t> & histo,
								 bool const boxWeightHisto); //Core patterns and box weight
	void enumerateCorePatterns3R(std::array<uint64_t,4> const & bx2_ay1_by2_nbSbox,
								 uint64_t const ax3,
								 std::map<double, uint64_t> & histo,
								 double const boundWeight); //Overload for trails

	void enumerateCorePatterns3RFromCore(std::array<std::map<uint64_t, std::vector<uint64_t>>, 16> const & Ldiff,
										 uint64_t const ay1,
										 uint64_t const ax3,
										 uint64_t const nbSbox,
										 std::map<uint64_t, uint64_t> & histo,
										 uint const currentIndex_x2,
										 uint const currentIndex_y2,
										 uint64_t x2,
										 uint64_t y2,
										 bool const boxWeightHisto); //Core patterns and box weight
	void enumerateCorePatterns3RFromCore(std::array<std::map<uint64_t, std::vector<uint64_t>>, 16> const & Ldiff,
										 uint64_t const ay1,
										 uint64_t const ax3,
										 uint64_t const nbSbox,
										 std::map<double, uint64_t> & histo,
										 uint const currentIndex_x2,
										 uint const currentIndex_y2,
										 uint64_t x2,
										 uint64_t y2,
										 double currentWeight,
										 double const minWeightR1R3,
										 double const boundWeight); //Overload for trails

	// 4 rounds
	void enumerateCorePatterns4R(std::array<uint64_t,4> const & bx2_ay1_by2_nbSbox,
								 uint64_t const ax3,
								 std::map<uint64_t, uint64_t> & histo,
								 uint const nbSboxMax,
								 bool const boxWeightHisto);
	void enumeratePartialCorePatterns3R(std::array<std::map<uint64_t, std::vector<uint64_t>>, 16> const & Ldiff,
										uint64_t const ay1,
										uint64_t const ax3,
										uint64_t const nbSboxR123,
										std::map<uint64_t, uint64_t> & histo,
										uint const currentIndex_x2,
										uint const currentIndex_y2,
										uint64_t x2,
										uint64_t y2,
										uint const nbSboxMax,
										bool const boxWeightHisto);
	void extendCore4r(uint64_t const y1,
					  uint64_t const x3,
					  uint64_t y3,
					  uint16_t patterny3,
					  uint64_t x4,
					  uint64_t const nbSboxR123,
					  uint64_t const nbSboxMax,
					  std::map<uint64_t, uint64_t> & histo,
					  bool const boxWeightHisto);

	//------------------------
	//-- Cluster Histograms --
	//------------------------
	std::pair<uint64_t,std::unordered_set<uint64_t>> getClusterClassSize(uint16_t const patterna,
								 uint16_t const patternPa);
	void clusterExhaustTrails2R(uint64_t y1,
								uint64_t x2,
								uint16_t patterny1,
								uint const nbSboxR1,
								uint const nbSboxMax,
								uint16_t const basePy1,
								std::map<uint64_t, std::map<uint64_t, uint64_t>> & histo,
								std::unordered_map<uint16_t, uint64_t> & knownClusterSize,
								std::unordered_set<uint64_t> & knownElements);

	//----------------------
	//-- Trail Clustering --
	//----------------------

	// 2 rounds
	void exhaustTrailsForCluster2R(uint64_t y1,
								   uint64_t x2,
								   uint16_t patterny1,
								   uint const nbSboxR1,
								   uint const nbSboxMax,
								   mapPairStateHisto & allTrails_thread);
	void exhaustTrailsForCluster2RFromTrailCore(uint64_t x1,
												uint64_t const y1,
												uint64_t const x2,
												uint64_t y2,
												uint16_t patternx1,
												uint16_t patterny2,
												double currentProba,
												mapPairStateHisto & allTrails_thread);

	// 3 rounds
	void enumerateCorePatterns3R(std::array<uint64_t,4> const & bx2_ay1_by2_nbSbox,
								 uint64_t const ax3,
								 mapPairStateHisto & histo); //Overload for trail clustering
	void enumerateCorePatterns3RFromCore(std::array<std::map<uint64_t, std::vector<uint64_t>>, 16> const & Ldiff,
										 uint64_t const ay1,
										 uint64_t const ax3,
										 uint64_t const nbSbox,
										 mapPairStateHisto & histo,
										 uint const currentIndex_x2,
										 uint const currentIndex_y2,
										 uint64_t x2,
										 uint64_t y2,
										 double currentWeight); //Overload for trail clustering


};

std::vector<std::vector<double>> computeDDT(std::vector<uint8_t> const & s,
											unsigned int const n);
//Return the DDT of #s which is an #n-bit sbox


void binPrint(uint64_t const x,
			  bool const rightLSB = true,
			  bool const split = true);
void binPrint(std::array<uint8_t,16> const & x, 
			  bool const rightLSB = true, 
			  bool const split = true);
#endif
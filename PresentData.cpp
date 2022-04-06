#include "PresentData.hpp"

using namespace std;
typedef unsigned int uint;

std::size_t pairHasher::operator()(std::pair<uint64_t,uint64_t> const & a) const {
    std::size_t h = 0;
	h ^= std::hash<uint64_t>{}(a.first)  + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); 
	h ^= std::hash<uint64_t>{}(a.second)  + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); 
    return h;
}

// -- Constructor for data precomputation --
PresentData::PresentData(std::vector<uint8_t> const & xp,
						 std::vector<std::vector<double>> const & xddt) :
						 p(xp),
						 invp(xp.size()),
						 ddt(xddt),
						 uint16WeightMap(17),
						 sboxTrails(16),
						 sboxInvTrails(16),
						 ddtPropagation(xddt.size()),
						 ddtInvPropagation(xddt.size()),
						 minProbaInvSbox(xddt.size()),
						 minWeightSbox(16),
						 isValidSboxIOSubspace(16, vector<bool>(16,false)),
						 minWeightIOSubspace(16, vector<double>(16,16))
{

	for(uint i = 0; i < 64; i++)
		invp[p[i]] = i;

	for(uint i = 0; i < 16; i++){
		uint64_t x = (0xFULL << 4*i);
		sboxMasks[i] = x;
		pSboxMasks[i] = applyPerm(x);
		invpSboxMasks[i] = applyInvPerm(x);

	}
	
	uint16_t x = 0;
	do{
		uint16WeightMap[__builtin_popcount(x)].emplace_back(x);
		x++;
	}while(x != 0);


	for(uint dx = 0; dx < ddt.size(); dx++){
		auto & ddtPropagation_dx = ddtPropagation[dx];
		for(uint dy = 0; dy < ddt[dx].size(); dy++){
			if(ddt[dx][dy] !=-1){ //valid propagation
				ddtPropagation_dx.emplace_back(dy);
				ddtInvPropagation[dy].emplace_back(dx);

				sboxTrails[dx][ddt[dx][dy]]++;
				sboxInvTrails[dy][ddt[dx][dy]]++;
			}
		}
	}

	//minProbaInvSbox[dy] contains the minimal weight (-log2(proba)) of all trails dx -> dy
	for(uint dy = 0; dy < ddt.size(); dy++){
		double minProba = 16;
		for(uint dx = 0; dx < ddt.size(); dx++){
			if(ddt[dx][dy] > 0 && ddt[dx][dy] < minProba)
				minProba = ddt[dx][dy];
		}
		minProbaInvSbox[dy] = minProba;
	}

	for(auto const & ddt_dx : ddt){
		for(auto const & pdxdy : ddt_dx){
			if(pdxdy > 0 && pdxdy < minWeightSbox)
				minWeightSbox = pdxdy;
		}
	}

	//isValidSboxIOSubspace
	//Initialized to false everywhere
	//For bx=0 or by=0; only [0][0] is valid
	isValidSboxIOSubspace[0][0] = true;
	minWeightIOSubspace[0][0] = 0;
	for(uint dx = 1; dx < ddt.size(); dx++){
		for(uint dy = 1; dy < ddt[dx].size(); dy++){
			if(ddt[dx][dy] !=-1){ //valid propagation
				//dx is differential that can be generated from any basis b such that b & dx = dx
				//Same idea for dy
				for(uint bx = 1; bx < ddt.size(); bx++){
					if((bx & dx) == dx){
						for(uint by = 1; by < ddt.size(); by++){
							if((by & dy) == dy){
								isValidSboxIOSubspace[bx][by] = true;
								if(ddt[dx][dy] < minWeightIOSubspace[bx][by])
									minWeightIOSubspace[bx][by] = ddt[dx][dy];
							}
						}
					}
				}
			}
		}
	}
}

uint64_t PresentData::applyPerm(uint64_t x){
	//Return p(x), where p is applied ot the bits of x
	//index 0 is the LSB of x, i.e. x & 1.
	
	uint64_t y = 0;
	for(uint i = 0; i < 64; i++){
		y |= (x & 1ULL) << p[i];
		x >>= 1;
	}
	return y;
}

uint64_t PresentData::applyInvPerm(uint64_t x){
	//Return invp(x), where p is applied ot the bits of x
	//index 0 is the LSB of x, i.e. x & 1.
	
	uint64_t y = 0;
	for(uint i = 0; i < 64; i++){
		y |= (x & 1ULL) << invp[i];
		x >>= 1;
	}
	return y;
}

uint16_t PresentData::getActivityPattern(uint64_t const x){
	//Return an uint16_t representing which Sboxes are active in x
	//i-th bit = 1 <-> i-th Sbox active
	uint16_t pattern = 0;
	for(uint i = 0; i < 16; i++){
		if(x & sboxMasks[i])
			pattern |= (1 << i);
	}
	return pattern;
}

unsigned int PresentData::countActiveSbox(uint64_t const x){
	//Return the number of active Sboxes in x
	uint ctr = 0;
	for(uint i = 0; i < 16; i++){
		if(x & sboxMasks[i])
			ctr++;
	}
	return ctr;
}

map<double,uint64_t> PresentData::convolutionRound(uint64_t const x,
												   double const boundWeight){
	//Return the histogram for all trails over a single Sbox Layer with input differential x
	//Ignore any trail with a weight over boundWeight

	map<double,uint64_t> histo;
	histo[0] = 1;
	for(uint i = 0; i < 16; i++){
		map<double,uint64_t> tmp;
		for(auto const & p_ctr : sboxTrails[getNibble(x,i)]){
			for(auto const & p_ctr2 : histo){
				if(p_ctr.first+p_ctr2.first <= boundWeight)
					tmp[p_ctr.first+p_ctr2.first] += p_ctr.second * p_ctr2.second;
			}
		}
		histo = move(tmp);
	}

	return histo;
}

map<double,uint64_t> PresentData::convolutionInvRound(uint64_t const x,
												      double const boundWeight){
	//Return the histogram for all trails over a single Sbox Layer with output differential x
	//Ignore any trail with a weight over boundWeight

	map<double,uint64_t> histo;
	histo[0] = 1;
	for(uint i = 0; i < 16; i++){
		map<double,uint64_t> tmp;
		for(auto const & p_ctr : sboxInvTrails[getNibble(x,i)]){
			for(auto const & p_ctr2 : histo){
				if(p_ctr.first+p_ctr2.first <= boundWeight)
					tmp[p_ctr.first+p_ctr2.first] += p_ctr.second * p_ctr2.second;
			}
		}
		histo = move(tmp);
	}

	return histo;
}

bool PresentData::checkCompatibleBasis(uint64_t bx, uint64_t by){
	//Return true is bx and by are compatible basis, meaning that at least one valid differential can be generated from them
	for(uint i = 0; i < 16; i++){
		if(!isValidSboxIOSubspace[bx & 0xF][by & 0xF])
			return false;
		bx >>= 4;
		by >>= 4;
	}
	return true;
}

double PresentData::getMinWeightBasis(uint64_t bx, uint64_t by){
	//Return the minimal weight of a differential made from the basis pair bx/by
	//Doesn't check if the basis is compatible, should be done beforehand. Non-comptatible basis will add a factor of 16 to the weight.
	double minWeight = 0;
	for(uint i = 0; i < 16; i++){
		minWeight += minWeightIOSubspace[bx & 0xF][by & 0xF];
		bx >>= 4;
		by >>= 4;
	}
	return minWeight;
}

// -- Generic functions --
std::vector<std::vector<double>> computeDDT(std::vector<uint8_t> const & s,
											 uint const n){
//Return the DDT of #s which is an #n-bit sbox
//ddt[dx][dy] is -log2 of the probability if the probability is >0, -1 otherwise

	uint bound = (1 << n);
	//Compute the number of solution for each s[x]^s[x^dx] = dy
	vector<vector<double>> ddt(bound, vector<double>(bound,0));
	for(uint x = 0; x < bound; x++){
		for(uint dx = 0; dx < bound; dx++){
			ddt[dx][s[x]^s[x^dx]]++;
		}
	}

	//Compute -log2(prob)
	for(uint dx = 0; dx < bound; dx++){
		for(uint dy = 0;dy < bound;dy++){
			if(ddt[dx][dy] == 0) ddt[dx][dy] = -1;
			else if(ddt[dx][dy] == bound) ddt[dx][dy] = 0;
			else ddt[dx][dy] = -log2(ddt[dx][dy]/bound);
		}
	}

	return ddt;
}


void binPrint(uint64_t const x, bool const rightLSB, bool const split){
	if(rightLSB){
		for(int i = 63; i >= 0; i--){
			if(x & (1ULL << i)) cout << 1;
			else cout << 0;
			// else cout << "*";

			if(split && ((i)%4 == 0))
				cout << " ";
		}
		cout << endl;
	}
	else{
		for(uint i = 0; i < 64; i++){
			if(x & (1ULL << i)) cout << 1;
			else cout << 0;
			// else cout << "*";

			if(split && ((i+1)%4 == 0))
				cout << " ";
		}
		cout << endl;
	}
}

void binPrint(array<uint8_t,16> const & x, bool const rightLSB, bool const split){
	if(rightLSB){
		for(int i = 63; i >= 0; i--){
			if(x[i/4] & (1ULL << (i%4))) cout << 1;
			else cout << 0;
			// else cout << "*";

			if(split && ((i)%4 == 0))
				cout << " ";
		}
		cout << endl;
	}
	else{
		for(uint i = 0; i < 64; i++){
			if(x[i/4] & (1ULL << (i%4))) cout << 1;
			else cout << 0;
			// else cout << "*";

			if(split && ((i+1)%4 == 0))
				cout << " ";
		}
		cout << endl;
	}
}


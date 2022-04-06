#include "PresentData.hpp"

using namespace std;
typedef unsigned int uint;

std::map<uint64_t, uint64_t> PresentData::corePatternHistogram3R(uint const nbSboxMax){
	return corePatternHistogram3R_internal(nbSboxMax,false);
}
std::map<uint64_t, uint64_t> PresentData::boxWeightHistogram3R(uint const nbSboxMax){
	return corePatternHistogram3R_internal(nbSboxMax,true);
}
std::map<uint64_t, uint64_t> PresentData::corePatternHistogram3R_internal(uint const nbSboxMax,
																		  bool const boxWeightHisto){

	//x1 -S-> y1 -P-> x2 -S-> y2 -P-> x3 -S-> y3

	map<uint64_t, uint64_t> histo;

	//Generate all active patterns for y1 i.e. sbox i = 1111 if sbox i active
	//Compute x2 = p(y1)
	//Compute the corresponding active pattern patx2 (same as above)
	//Store (x2, patx2, ay1) in Tx2[i] for each i such that sbox i is active in x2
	//We need to keep ay1 as to be sure that we avoid duplicates later on, so that the actual core trail does follow ay1
	vector<vector<array<uint64_t,3>>> Tx2(16);

	//First round so at least ceil(wy1/4) sboxes in the second round + 1 for third
	for(uint wy1 = 1; wy1+ceil(double(wy1)/4)+1 <= nbSboxMax; wy1++){
		for(auto const & ay1 : uint16WeightMap[wy1]){
			//Compute x2 directly from the activity pattern of y1
			uint64_t x2 = 0;
			for(uint i = 0; i < 16; i++){
				if(ay1 & (1ULL << i))
					x2 |= pSboxMasks[i];
			}

			//Get patx2 and which sboxes are active
			vector<uint> activeSbox;
			activeSbox.reserve(16);
			uint64_t patx2 = 0;
			for(uint i = 0; i < 16; i++){
				if(x2 & sboxMasks[i]){
					patx2 |= sboxMasks[i];
					activeSbox.emplace_back(i);
				}
			}

			for(auto const & i : activeSbox)
				Tx2[i].push_back({x2,patx2,ay1});
		}
	}

	//Now generate all patterns for x3
	//Third round so :
	//- at least ceil(wx3/4) sboxes in the second round
	//- at least 1 sbox in the first round
	for(uint wx3 = 1; wx3+ceil(double(wx3)/4)+1 <= nbSboxMax; wx3++){
		cout << "Weight x3 : " << wx3 << endl;
		auto const & uint16WeightMap_wx3 = uint16WeightMap[wx3];
		uint64_t const bound = uint16WeightMap_wx3.size();

		#pragma omp parallel
		{
			map<uint64_t, uint64_t> histo_thread;
			#pragma omp for schedule(dynamic,1)
			for(uint64_t iax3 = 0; iax3 < bound; iax3++){
			// for(auto const & ax3 : uint16WeightMap[wx3]){
				auto ax3 = uint16WeightMap_wx3[iax3];
				//Compute y2 directly from the activity pattern of x3
				uint64_t y2 = 0;
				for(uint i = 0; i < 16; i++){
					if(ax3 & (1ULL << i))
						y2 |= invpSboxMasks[i];
				}

				//Get paty2 and which sboxes are active
				vector<uint> activeSbox;
				uint64_t paty2 = 0;
				for(uint i = 0; i < 16; i++){
					if(y2 & sboxMasks[i]){
						paty2 |= sboxMasks[i];
						activeSbox.emplace_back(i);
					}
				}

				//Compute which states for x2 *could* be compatible
				//This gives us a list of {bx2, ay1, by2} where 
				// - bx2 represent the basis of the vector space for x2
				// - ay1 represent the activity pattern that y1 must follow
				// - by2 represent the basis of the vector space for y2
				//Technically we also need ax3 for the activity pattern that x3 must follow, 
				// but we don't need to store it as it is local
				set<array<uint64_t,4>> compatibleStates;
				set<array<uint64_t,3>> checkedStates; //To avoid checking duplicates
				//x2 and y2 must have at least one active sbox in common
				for(auto const & isbox : activeSbox){
					for(auto const & x2_patx2_ay1 : Tx2[isbox]){
						if(checkedStates.find(x2_patx2_ay1) == checkedStates.end()){
							//Didn't check this state
							checkedStates.emplace(x2_patx2_ay1);

							//Compute the mask giving the common active sboxes
							uint64_t mask = paty2 & x2_patx2_ay1[1];

							//Check that the number of active sboxes in the middle round isn't too high
							if(__builtin_popcountll(getActivityPattern(mask)) <= nbSboxMax-2){

								//Compute the restricted basis for y2
								uint64_t by2 = y2 & mask;
								//Check that p(by2) still follows the activity pattern ax3
								//Doesn't mean that all induced differentials will follow ax3, 
								// but at least one (without considering sbox compatibility)
								if(getActivityPattern(applyPerm(by2)) == ax3){
									//Compute the restricted basis for x2
									uint64_t bx2 = x2_patx2_ay1[0] & mask;
									//Check that invp(bx2) still follows the activity pattern ay1
									//Doesn't mean that all induced differentials will follow ay1, 
									// but at least one (without considering sbox compatibility)
									uint64_t ay1 = x2_patx2_ay1[2];
									//Since we will enforce that ay1 and ax3 are verified, we already know the number of active sboxes in any trail generated from this
									//So check that it isn't too high
									//Also check that at least one differential is possible
									uint nbSbox = __builtin_popcountll(ay1) + countActiveSbox(mask) + __builtin_popcountll(ax3);
									if(nbSbox <= nbSboxMax &&
									getActivityPattern(applyInvPerm(bx2)) == ay1 &&
									checkCompatibleBasis(bx2,by2)){
										//bx2,ay1,by2,ax3 is valid so far
										compatibleStates.insert({bx2,ay1,by2,nbSbox});
									}
								}
							}
						}
					}
				}

				for(auto const & bx2_ay1_by2_nbSbox : compatibleStates){
					enumerateCorePatterns3R(bx2_ay1_by2_nbSbox, ax3, histo_thread,boxWeightHisto);
					// cout << "Finished enumerating all cores with this pattern" << endl << endl;
				}

			}

			//Merge in the global map
			#pragma omp critical
			{
				for(auto const & n_ctr : histo_thread){
					histo[n_ctr.first] += n_ctr.second;
				}
			}
		}
	}

	return histo;
}

void PresentData::enumerateCorePatterns3R(array<uint64_t,4> const & bx2_ay1_by2_nbSbox,
										  uint64_t const ax3,
										  std::map<uint64_t, uint64_t> & histo,
										  bool const boxWeightHisto){
	//x1 -S-> y1 -P-> x2 -S-> y2 -P-> x3 -S-> y3
	//Enumerate all core patterns such that :
	// - bx2 is a basis for the differential subspace over x2
	// - by2 is a basis for the differential subspace over y2
	// - The resulting difference in x2 must result in a difference compatible with activity pattern ay1 
	//   (i.e. getActivityPattern(invp(x2)) = ay1)
	// - Similarly, getActivityPattern(p(y2)) = ax3
	//And update the histogram accordingly
	//Note that we already know the number of active sboxes so no need to recompute it nor to check if it fits the bound

	uint64_t const bx2 = bx2_ay1_by2_nbSbox[0];
	uint64_t const ay1 = bx2_ay1_by2_nbSbox[1];
	uint64_t const by2 = bx2_ay1_by2_nbSbox[2];
	uint64_t const nbSbox = bx2_ay1_by2_nbSbox[3];

	//First, create a list Ldiff such that Ldiff[i] is a list of valid input/output differential pairs for Sbox i
	//This allows for early detection of incompatibilities
	//Each list would be of at most 2^8 elements so rather small
	//Morever, the stored differences will already be shifted to their position on the 64-bit state, such that
	//a difference over the full state can be obtained directly by xoring everything without additional shifts

	array<map<uint64_t, vector<uint64_t>>, 16> Ldiff;
	for(uint isbox = 0; isbox < 16; isbox++){
		// uint64_t bx2i = bx2 & sboxMasks[isbox]; //Basis for the input of i-th sbox
		// uint64_t by2i = by2 & sboxMasks[isbox]; //Basis for the input of i-th sbox
		auto & Ldiff_isbox = Ldiff[isbox];

		//Get the unit vectors (as uint64_t) in the basis for this sbox
		vector<uint64_t> vecbx2i;
		vector<uint64_t> vecby2i;
		uint64_t offset = 4*isbox;
		for(uint i = 0; i < 4; i++){
			uint64_t index = offset+i;
			if(bx2 & (1ULL << index))
				vecbx2i.emplace_back(1ULL << index);
			if(by2 & (1ULL << index))
				vecby2i.emplace_back(1ULL << index);
		}

		if((vecbx2i.size() == 0 && vecby2i.size() != 0) || (vecbx2i.size() != 0 && vecby2i.size() == 0)){
			//Should not happen with how the basis are generated but to be safe
			//In this case, incompatible basis, so no trails can be generated
			return;
		}
		if(vecbx2i.size() == 0 && vecby2i.size() == 0){
			//This means that the Sbox is inactive, so the only io difference pair is (0,0)
			Ldiff_isbox[0].emplace_back(0);
		}
		else{
			//Generate all possible io pairs for this sbox
			uint64_t in_bound = (1ULL << vecbx2i.size());
			uint64_t out_bound = (1ULL << vecby2i.size());
			for(uint64_t ix2 = 1; ix2 < in_bound; ix2++){
				//bit i of ix2 = 1 <=> current differential x2 uses the vector vecbx2i[i]
				uint64_t x2 = 0;
				for(uint i = 0; i < vecbx2i.size(); i++){
					if(ix2 & (1ULL << i))
						x2 |= vecbx2i[i];
				}
				//x2 is the shifted difference (i.e. align to the sbox)
				//get the value to check against the DDT
				uint64_t valx2 = (x2 >> offset);

				for(uint64_t iy2 = 1; iy2 < out_bound; iy2++){
					//bit i of iy2 = 1 <=> current differential y2 uses the vector vecby2i[i]
					uint64_t y2 = 0;
					for(uint i = 0; i < vecby2i.size(); i++){
						if(iy2 & (1ULL << i))
							y2 |= vecby2i[i];
					}
					//Again, y2 is the shifted value, get the actual value
					uint64_t valy2 = (y2 >> offset);

					if(ddt[valx2][valy2] != -1){//Valid differential
						Ldiff_isbox[x2].emplace_back(y2);
					}
				}
			}
		}

		//At this point all possible io pairs for this sbox have been generated
		//Check if there is at least one, if not, then no trail can be generated so stop here
		//Should not happen because of the previous checks, but to be safe
		if(Ldiff_isbox.size() == 0){
			// cout << "No possible transition on Sbox " << isbox << endl;
			return;
		}
	}
	//We have the list of all possible transitions for each sbox, and each sbox has at least one valid transition

	// //Time to go through all of them to see how many verify the activity patterns
	uint currentIndex_x2 = 0;
	uint currentIndex_y2 = 0;
	uint64_t x2 = 0;
	uint64_t y2 = 0;
	enumerateCorePatterns3RFromCore(Ldiff, ay1, ax3, nbSbox, histo, currentIndex_x2, currentIndex_y2, x2, y2,boxWeightHisto);

}

void PresentData::enumerateCorePatterns3RFromCore(array<map<uint64_t, vector<uint64_t>>, 16> const & Ldiff,
												  uint64_t const ay1,
												  uint64_t const ax3,
												  uint64_t const nbSbox,
												  std::map<uint64_t, uint64_t> & histo,
												  uint const currentIndex_x2,
												  uint const currentIndex_y2,
												  uint64_t x2,
												  uint64_t y2,
												  bool const boxWeightHisto){
	//Recursively build all possible trails from the transitions in Ldiff, 
	//and count how many verify the activity patterns ay1 and ax3
	//currentIndex_x2/y2 indicates the index of the sbox that needs to be guessed for the input/output, currentIndex = 16 marks the end of the recursion for the input/output (respectively)
	if(currentIndex_x2 == 16 && currentIndex_y2 == 16){
		//x2 and y2 are fully determined
		//Due to the workflow, it's already been checked that x2 verifies ay1 and y2 verifies ax3
		//So we have a valid core
		if(boxWeightHisto){
			//Only need to convolute the number of trails
			//Probability doesn't matter here, so arbitrary bound that should include any trail
			auto histoR1 = convolutionInvRound(applyInvPerm(x2),2*16*4);
			auto histoR3 = convolutionRound(applyPerm(y2),2*16*4);
			uint64_t nbTrailR1 = 0;
			uint64_t nbTrailR3 = 0;
			for(auto const & p_ctrR1 : histoR1)
				nbTrailR1 += p_ctrR1.second;
			for(auto const & p_ctrR3 : histoR3)
				nbTrailR3 += p_ctrR3.second;
			histo[nbSbox] += nbTrailR1*nbTrailR3;
		}
		else
			histo[nbSbox]++;
	}

	else{
		if(currentIndex_x2 < 16){ //Still need to guess the input x2

			uint64_t prev_x2 = x2;
			for(auto const & x2i_o : Ldiff[currentIndex_x2]){
				x2 = prev_x2 | x2i_o.first; //Already shifted

				//Either we will still have some guesses remaining so go forward
				// or we just finished guessing x2, check if it verifies ay1
				if(currentIndex_x2 < 15 || (currentIndex_x2 == 15 && getActivityPattern(applyInvPerm(x2)) == ay1)){
					enumerateCorePatterns3RFromCore(Ldiff, ay1, ax3, nbSbox, histo, currentIndex_x2+1, currentIndex_y2, x2, y2,boxWeightHisto);
				}
			}
		}

		else if(currentIndex_y2 < 16){ //Still need to guess the output y2

			//Get the input diff
			uint64_t const x2i = x2 & (0xFULL << (4*currentIndex_y2));
			uint64_t prev_y2 = y2;

			for(auto const y2i : Ldiff[currentIndex_y2].at(x2i)){ //.at for const access
				y2 = prev_y2 | y2i; //Already shifted
			
				//Either we will still have some guesses remaining so go forward
				// or we just finished guessing y2, check if it verifies ax3
				if(currentIndex_y2 < 15 || (currentIndex_y2 == 15 && getActivityPattern(applyPerm(y2)) == ax3)){
					enumerateCorePatterns3RFromCore(Ldiff, ay1, ax3, nbSbox, histo, currentIndex_x2, currentIndex_y2+1, x2, y2,boxWeightHisto);
				}
			}
		}
	}

}


map<double, uint64_t> PresentData::trailHistogram3R(unsigned int const nbSboxMax){

	//x1 -S-> y1 -P-> x2 -S-> y2 -P-> x3 -S-> y3

	map<double, uint64_t> histo;

	//Get the bound for the best probability for nbSboxMax active sboxes
	//This sets an upper bound on the relevant trails to exhaust
	double boundWeight = nbSboxMax*minWeightSbox;
	cout << "Bound on probability for 3 rounds: " << boundWeight << endl;

	//Generate all active patterns for y1 i.e. sbox i = 1111 if sbox i active
	//Compute x2 = p(y1)
	//Compute the corresponding active pattern patx2 (same as above)
	//Store (x2, patx2, ay1) in Tx2[i] for each i such that sbox i is active in x2
	//We need to keep ay1 as to be sure that we avoid duplicates later on, so that the actual core trail does follow ay1
	vector<vector<array<uint64_t,3>>> Tx2(16);

	//First round so at least ceil(wy1/4) sboxes in the second round + 1 for third
	for(uint wy1 = 1; wy1+ceil(double(wy1)/4)+1 <= nbSboxMax; wy1++){
		for(auto const & ay1 : uint16WeightMap[wy1]){
			//Compute x2 directly from the activity pattern of y1
			uint64_t x2 = 0;
			for(uint i = 0; i < 16; i++){
				if(ay1 & (1ULL << i))
					x2 |= pSboxMasks[i];
			}

			//Get patx2 and which sboxes are active
			vector<uint> activeSbox;
			activeSbox.reserve(16);
			uint64_t patx2 = 0;
			for(uint i = 0; i < 16; i++){
				if(x2 & sboxMasks[i]){
					patx2 |= sboxMasks[i];
					activeSbox.emplace_back(i);
				}
			}

			for(auto const & i : activeSbox)
				Tx2[i].push_back({x2,patx2,ay1});
		}
	}

	//Now generate all patterns for x3
	//Third round so :
	//- at least ceil(wx3/4) sboxes in the second round
	//- at least 1 sbox in the first round
	for(uint wx3 = 1; wx3+ceil(double(wx3)/4)+1 <= nbSboxMax; wx3++){
		cout << "Weight x3 : " << wx3 << endl;
		auto const & uint16WeightMap_wx3 = uint16WeightMap[wx3];
		uint64_t const bound = uint16WeightMap_wx3.size();

		#pragma omp parallel
		{
			map<double, uint64_t> histo_thread;
			#pragma omp for schedule(dynamic,1)
			for(uint64_t iax3 = 0; iax3 < bound; iax3++){
			// for(auto const & ax3 : uint16WeightMap[wx3]){
				auto ax3 = uint16WeightMap_wx3[iax3];
				//Compute y2 directly from the activity pattern of x3
				uint64_t y2 = 0;
				for(uint i = 0; i < 16; i++){
					if(ax3 & (1ULL << i))
						y2 |= invpSboxMasks[i];
				}

				//Get paty2 and which sboxes are active
				vector<uint> activeSbox;
				uint64_t paty2 = 0;
				for(uint i = 0; i < 16; i++){
					if(y2 & sboxMasks[i]){
						paty2 |= sboxMasks[i];
						activeSbox.emplace_back(i);
					}
				}

				//Compute which states for x2 *could* be compatible
				//This gives us a list of {bx2, ay1, by2} where 
				// - bx2 represent the basis of the vector space for x2
				// - ay1 represent the activity pattern that y1 must follow
				// - by2 represent the basis of the vector space for y2
				//Technically we also need ax3 for the activity pattern that x3 must follow, 
				// but we don't need to store it as it is local
				set<array<uint64_t,4>> compatibleStates;
				set<array<uint64_t,3>> checkedStates; //To avoid checking duplicates
				//x2 and y2 must have at least one active sbox in common
				for(auto const & isbox : activeSbox){
					for(auto const & x2_patx2_ay1 : Tx2[isbox]){
						if(checkedStates.find(x2_patx2_ay1) == checkedStates.end()){
							//Didn't check this state
							checkedStates.emplace(x2_patx2_ay1);

							//Compute the mask giving the common active sboxes
							uint64_t mask = paty2 & x2_patx2_ay1[1];

							//Check that the number of active sboxes in the middle round isn't too high
							if(__builtin_popcountll(getActivityPattern(mask)) <= nbSboxMax-2){

								//Compute the restricted basis for y2
								uint64_t by2 = y2 & mask;
								//Check that p(by2) still follows the activity pattern ax3
								//Doesn't mean that all induced differentials will follow ax3, 
								// but at least one (without considering sbox compatibility)
								if(getActivityPattern(applyPerm(by2)) == ax3){
									//Compute the restricted basis for x2
									uint64_t bx2 = x2_patx2_ay1[0] & mask;
									//Check that invp(bx2) still follows the activity pattern ay1
									//Doesn't mean that all induced differentials will follow ay1, 
									// but at least one (without considering sbox compatibility)
									uint64_t ay1 = x2_patx2_ay1[2];
									//Since we will enforce that ay1 and ax3 are verified, we already know the number of active sboxes in any trail generated from this
									//So check that it isn't too high
									//Also check that at least one differential is possible
									uint nbSbox = __builtin_popcountll(ay1) + countActiveSbox(mask) + __builtin_popcountll(ax3);
									if(nbSbox <= nbSboxMax &&
									getActivityPattern(applyInvPerm(bx2)) == ay1 &&
									checkCompatibleBasis(bx2,by2) &&
									getMinWeightBasis(bx2,by2) <= boundWeight){
										//bx2,ay1,by2,ax3 is valid so far
										compatibleStates.insert({bx2,ay1,by2,nbSbox});
									}
								}
							}
						}
					}
				}

				for(auto const & bx2_ay1_by2_nbSbox : compatibleStates){
					enumerateCorePatterns3R(bx2_ay1_by2_nbSbox, ax3, histo_thread, boundWeight);
					// cout << "Finished enumerating all cores with this pattern" << endl << endl;
				}

			}

			//Merge in the global map
			#pragma omp critical
			{
				for(auto const & n_ctr : histo_thread){
					histo[n_ctr.first] += n_ctr.second;
				}
			}
		}
	}

	return histo;
}

void PresentData::enumerateCorePatterns3R(array<uint64_t,4> const & bx2_ay1_by2_nbSbox,
										  uint64_t const ax3,
										  std::map<double, uint64_t> & histo,
										  double const boundWeight){
	//x1 -S-> y1 -P-> x2 -S-> y2 -P-> x3 -S-> y3
	//Enumerate all core patterns such that :
	// - bx2 is a basis for the differential subspace over x2
	// - by2 is a basis for the differential subspace over y2
	// - The resulting difference in x2 must result in a difference compatible with activity pattern ay1 
	//   (i.e. getActivityPattern(invp(x2)) = ay1)
	// - Similarly, getActivityPattern(p(y2)) = ax3
	//And update the histogram accordingly
	//Note that we already know the number of active sboxes so no need to recompute it nor to check if it fits the bound

	uint64_t const bx2 = bx2_ay1_by2_nbSbox[0];
	uint64_t const ay1 = bx2_ay1_by2_nbSbox[1];
	uint64_t const by2 = bx2_ay1_by2_nbSbox[2];
	uint64_t const nbSbox = bx2_ay1_by2_nbSbox[3];

	//First, create a list Ldiff such that Ldiff[i] is a list of valid input/output differential pairs for Sbox i
	//This allows for early detection of incompatibilities
	//Each list would be of at most 2^8 elements so rather small
	//Morever, the stored differences will already be shifted to their position on the 64-bit state, such that
	//a difference over the full state can be obtained directly by xoring everything without additional shifts

	array<map<uint64_t, vector<uint64_t>>, 16> Ldiff;
	for(uint isbox = 0; isbox < 16; isbox++){
		// uint64_t bx2i = bx2 & sboxMasks[isbox]; //Basis for the input of i-th sbox
		// uint64_t by2i = by2 & sboxMasks[isbox]; //Basis for the input of i-th sbox
		auto & Ldiff_isbox = Ldiff[isbox];

		//Get the unit vectors (as uint64_t) in the basis for this sbox
		vector<uint64_t> vecbx2i;
		vector<uint64_t> vecby2i;
		uint64_t offset = 4*isbox;
		for(uint i = 0; i < 4; i++){
			uint64_t index = offset+i;
			if(bx2 & (1ULL << index))
				vecbx2i.emplace_back(1ULL << index);
			if(by2 & (1ULL << index))
				vecby2i.emplace_back(1ULL << index);
		}

		if((vecbx2i.size() == 0 && vecby2i.size() != 0) || (vecbx2i.size() != 0 && vecby2i.size() == 0)){
			//Should not happen with how the basis are generated but to be safe
			//In this case, incompatible basis, so no trails can be generated
			return;
		}
		if(vecbx2i.size() == 0 && vecby2i.size() == 0){
			//This means that the Sbox is inactive, so the only io difference pair is (0,0)
			Ldiff_isbox[0].emplace_back(0);
		}
		else{
			//Generate all possible io pairs for this sbox
			uint64_t in_bound = (1ULL << vecbx2i.size());
			uint64_t out_bound = (1ULL << vecby2i.size());
			for(uint64_t ix2 = 1; ix2 < in_bound; ix2++){
				//bit i of ix2 = 1 <=> current differential x2 uses the vector vecbx2i[i]
				uint64_t x2 = 0;
				for(uint i = 0; i < vecbx2i.size(); i++){
					if(ix2 & (1ULL << i))
						x2 |= vecbx2i[i];
				}
				//x2 is the shifted difference (i.e. align to the sbox)
				//get the value to check against the DDT
				uint64_t valx2 = (x2 >> offset);

				for(uint64_t iy2 = 1; iy2 < out_bound; iy2++){
					//bit i of iy2 = 1 <=> current differential y2 uses the vector vecby2i[i]
					uint64_t y2 = 0;
					for(uint i = 0; i < vecby2i.size(); i++){
						if(iy2 & (1ULL << i))
							y2 |= vecby2i[i];
					}
					//Again, y2 is the shifted value, get the actual value
					uint64_t valy2 = (y2 >> offset);

					if(ddt[valx2][valy2] != -1){//Valid differential
						Ldiff_isbox[x2].emplace_back(y2);
					}
				}
			}
		}

		//At this point all possible io pairs for this sbox have been generated
		//Check if there is at least one, if not, then no trail can be generated so stop here
		//Should not happen because of the previous checks, but to be safe
		if(Ldiff_isbox.size() == 0){
			// cout << "No possible transition on Sbox " << isbox << endl;
			return;
		}
	}
	//We have the list of all possible transitions for each sbox, and each sbox has at least one valid transition

	// //Time to go through all of them to see how many verify the activity patterns
	uint currentIndex_x2 = 0;
	uint currentIndex_y2 = 0;
	uint64_t x2 = 0;
	uint64_t y2 = 0;
	double currentWeight = 0;
	//We know the number of active sboxes in R1 and R3, so we have a lower bound on the minimum weight for those rounds
	double minWeightR1R3 = __builtin_popcountll(ay1)*minWeightSbox + __builtin_popcountll(ax3)*minWeightSbox;

	enumerateCorePatterns3RFromCore(Ldiff, ay1, ax3, nbSbox, histo, currentIndex_x2, currentIndex_y2, x2, y2, currentWeight, minWeightR1R3, boundWeight);

}


void PresentData::enumerateCorePatterns3RFromCore(array<map<uint64_t, vector<uint64_t>>, 16> const & Ldiff,
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
												  double const boundWeight){
	//Recursively build all possible trails from the transitions in Ldiff, 
	//currentIndex_x2/y2 indicates the index of the sbox that needs to be guessed for the input/output, currentIndex = 16 marks the end of the recursion for the input/output (respectively)
	//Unlike for the boxWeight/corePattern enumeration, here we will build both input/output of a given sbox before moving to the next one
	//i.e. for boxWeight/corePattern, the indexes for x2/y2 was iterated over as x2 : 0 1 2 ... 15 16 then y2 : 0 1 2 ... 15 16
	//Here, we first do x2 : 0, y2 : 0 then x2 : 1, y2 : 1 etc.
	//Goal is to bound the weight of the trail early to cut off if necessary
	if(currentIndex_x2 == 16 && currentIndex_y2 == 16){
		//x2 and y2 are fully determined, and round 2 has weight currentWeight
		//Check if the activity pattern on the first and last round match
		auto histoR1 = convolutionInvRound(applyInvPerm(x2),boundWeight);
		auto histoR3 = convolutionRound(applyPerm(y2),boundWeight);
		for(auto const & p_ctrR1 : histoR1){
			double weightR1R2 = currentWeight+p_ctrR1.first;
			if(weightR1R2 <= boundWeight){
				for(auto const & p_ctrR3 : histoR3){
					double endWeight = weightR1R2+p_ctrR3.first;
					if(endWeight <= boundWeight)
						histo[endWeight] += p_ctrR1.second*p_ctrR3.second;
				}
			}
		}
	}

	else{
		if(currentIndex_x2 == currentIndex_y2){ //Guess on the input x2

			uint64_t prev_x2 = x2;
			for(auto const & x2i_o : Ldiff[currentIndex_x2]){
				x2 = prev_x2 | x2i_o.first; //Already shifted

				if(currentIndex_x2 < 15 || 
				  (currentIndex_x2 == 15 && getActivityPattern(applyInvPerm(x2)) == ay1)){
					enumerateCorePatterns3RFromCore(Ldiff, ay1, ax3, nbSbox, histo, currentIndex_x2+1, currentIndex_y2, x2, y2,currentWeight,minWeightR1R3,boundWeight);
				}
			}
		}
		else{ //Guess on the output y2

			//Get the input diff
			uint64_t offset = 4*currentIndex_y2;
			uint64_t const x2i = x2 & (0xFULL << offset);
			uint64_t x2isbox = x2i >> offset; //value of x2i between 0 and 15 to check ddt
			uint64_t prev_y2 = y2;
			double prev_currentWeight = currentWeight;

			for(auto const y2i : Ldiff[currentIndex_y2].at(x2i)){ //.at for const access
				y2 = prev_y2 | y2i; //Already shifted
				uint64_t y2isbox = y2i >> offset; //value of y2i between 0 and 15 to check ddt
				currentWeight = prev_currentWeight + ddt[x2isbox][y2isbox];

				if(((currentWeight + minWeightR1R3) <= boundWeight) && 
					(currentIndex_y2 < 15 || 
					(currentIndex_y2 == 15 && getActivityPattern(applyPerm(y2)) == ax3))){
					enumerateCorePatterns3RFromCore(Ldiff, ay1, ax3, nbSbox, histo, currentIndex_x2, currentIndex_y2+1, x2, y2,currentWeight,minWeightR1R3,boundWeight);
				}
			}
		}
	}
}

std::map<double, uint64_t> PresentData::trailClusterHistogram3R(uint const nbSboxMax){

	//x1 -S-> y1 -P-> x2 -S-> y2 -P-> x3 -S-> y3
	map<double, uint64_t> histo;

	//Generate all active patterns for y1 i.e. sbox i = 1111 if sbox i active
	//Compute x2 = p(y1)
	//Compute the corresponding active pattern patx2 (same as above)
	//Store (x2, patx2, ay1) in Tx2[i] for each i such that sbox i is active in x2
	//We need to keep ay1 as to be sure that we avoid duplicates later on, so that the actual core trail does follow ay1
	vector<vector<array<uint64_t,3>>> Tx2(16);

	//First round so at least ceil(wy1/4) sboxes in the second round + 1 for third
	for(uint wy1 = 1; wy1+ceil(double(wy1)/4)+1 <= nbSboxMax; wy1++){
		for(auto const & ay1 : uint16WeightMap[wy1]){
			//Compute x2 directly from the activity pattern of y1
			uint64_t x2 = 0;
			for(uint i = 0; i < 16; i++){
				if(ay1 & (1ULL << i))
					x2 |= pSboxMasks[i];
			}

			//Get patx2 and which sboxes are active
			vector<uint> activeSbox;
			activeSbox.reserve(16);
			uint64_t patx2 = 0;
			for(uint i = 0; i < 16; i++){
				if(x2 & sboxMasks[i]){
					patx2 |= sboxMasks[i];
					activeSbox.emplace_back(i);
				}
			}

			for(auto const & i : activeSbox)
				Tx2[i].push_back({x2,patx2,ay1});
		}
	}

	//Now generate all patterns for x3
	//Third round so :
	//- at least ceil(wx3/4) sboxes in the second round
	//- at least 1 sbox in the first round
	for(uint wx3 = 1; wx3+ceil(double(wx3)/4)+1 <= nbSboxMax; wx3++){
		cout << "Weight x3 : " << wx3 << endl;
		auto const & uint16WeightMap_wx3 = uint16WeightMap[wx3];
		uint64_t const bound = uint16WeightMap_wx3.size();

		#pragma omp parallel
		{
			map<double, uint64_t> histo_thread;
			#pragma omp for schedule(dynamic,1)
			for(uint64_t iax3 = 0; iax3 < bound; iax3++){
			// for(auto const & ax3 : uint16WeightMap[wx3]){
				auto ax3 = uint16WeightMap_wx3[iax3];
				//Compute y2 directly from the activity pattern of x3
				uint64_t y2 = 0;
				for(uint i = 0; i < 16; i++){
					if(ax3 & (1ULL << i))
						y2 |= invpSboxMasks[i];
				}

				//Get paty2 and which sboxes are active
				vector<uint> activeSbox;
				uint64_t paty2 = 0;
				for(uint i = 0; i < 16; i++){
					if(y2 & sboxMasks[i]){
						paty2 |= sboxMasks[i];
						activeSbox.emplace_back(i);
					}
				}

				//Compute which states for x2 *could* be compatible
				//This gives us a list of {bx2, ay1, by2} where 
				// - bx2 represent the basis of the vector space for x2
				// - ay1 represent the activity pattern that y1 must follow
				// - by2 represent the basis of the vector space for y2
				//Technically we also need ax3 for the activity pattern that x3 must follow, 
				// but we don't need to store it as it is local
				set<array<uint64_t,4>> compatibleStates;
				set<array<uint64_t,3>> checkedStates; //To avoid checking duplicates
				//x2 and y2 must have at least one active sbox in common
				for(auto const & isbox : activeSbox){
					for(auto const & x2_patx2_ay1 : Tx2[isbox]){
						if(checkedStates.find(x2_patx2_ay1) == checkedStates.end()){
							//Didn't check this state
							checkedStates.emplace(x2_patx2_ay1);

							//Compute the mask giving the common active sboxes
							uint64_t mask = paty2 & x2_patx2_ay1[1];

							//Check that the number of active sboxes in the middle round isn't too high
							if(__builtin_popcountll(getActivityPattern(mask)) <= nbSboxMax-2){

								//Compute the restricted basis for y2
								uint64_t by2 = y2 & mask;
								//Check that p(by2) still follows the activity pattern ax3
								//Doesn't mean that all induced differentials will follow ax3, 
								// but at least one (without considering sbox compatibility)
								if(getActivityPattern(applyPerm(by2)) == ax3){
									//Compute the restricted basis for x2
									uint64_t bx2 = x2_patx2_ay1[0] & mask;
									//Check that invp(bx2) still follows the activity pattern ay1
									//Doesn't mean that all induced differentials will follow ay1, 
									// but at least one (without considering sbox compatibility)
									uint64_t ay1 = x2_patx2_ay1[2];
									//Since we will enforce that ay1 and ax3 are verified, we already know the number of active sboxes in any trail generated from this
									//So check that it isn't too high
									//Also check that at least one differential is possible
									uint nbSbox = __builtin_popcountll(ay1) + countActiveSbox(mask) + __builtin_popcountll(ax3);
									if(nbSbox <= nbSboxMax &&
									getActivityPattern(applyInvPerm(bx2)) == ay1 &&
									checkCompatibleBasis(bx2,by2)){
										//bx2,ay1,by2,ax3 is valid so far
										compatibleStates.insert({bx2,ay1,by2,nbSbox});
									}
								}
							}
						}
					}
				}

				for(auto const & bx2_ay1_by2_nbSbox : compatibleStates){
					mapPairStateHisto allTrails_thread;
					//Next call will enumerate all trails that cluster together following the activity pattern of ay1 and ax3
					//Meaning that 2 separate calls enumerate trails that cannot belong in the same cluster, so safe to update the histogram after each call
					enumerateCorePatterns3R(bx2_ay1_by2_nbSbox, ax3, allTrails_thread);
					
					for(auto const & state_histo : allTrails_thread){
						double proba = 0;
						for(auto const & p_ctr : state_histo.second)
							proba += p_ctr.second * pow(2,-p_ctr.first);
						histo_thread[-log2(proba)]++;

						// cout << "dx : "; binPrint(state_histo.first.first);
						// cout << "dy : "; binPrint(state_histo.first.second);
						// cout << "proba : " << -log2(proba) << endl;
					}
					// getchar();
				}

			}

			//Merge in the global map
			#pragma omp critical
			{
				for(auto const & n_ctr : histo_thread){
					histo[n_ctr.first] += n_ctr.second;
				}
			}
		}
	}
	return histo;
}

void PresentData::enumerateCorePatterns3R(array<uint64_t,4> const & bx2_ay1_by2_nbSbox,
										  uint64_t const ax3,
										  mapPairStateHisto & histo){
	//x1 -S-> y1 -P-> x2 -S-> y2 -P-> x3 -S-> y3
	//Enumerate all core patterns such that :
	// - bx2 is a basis for the differential subspace over x2
	// - by2 is a basis for the differential subspace over y2
	// - The resulting difference in x2 must result in a difference compatible with activity pattern ay1 
	//   (i.e. getActivityPattern(invp(x2)) = ay1)
	// - Similarly, getActivityPattern(p(y2)) = ax3
	//And update the histogram accordingly
	//Note that we already know the number of active sboxes so no need to recompute it nor to check if it fits the bound

	uint64_t const bx2 = bx2_ay1_by2_nbSbox[0];
	uint64_t const ay1 = bx2_ay1_by2_nbSbox[1];
	uint64_t const by2 = bx2_ay1_by2_nbSbox[2];
	uint64_t const nbSbox = bx2_ay1_by2_nbSbox[3];

	//First, create a list Ldiff such that Ldiff[i] is a list of valid input/output differential pairs for Sbox i
	//This allows for early detection of incompatibilities
	//Each list would be of at most 2^8 elements so rather small
	//Morever, the stored differences will already be shifted to their position on the 64-bit state, such that
	//a difference over the full state can be obtained directly by xoring everything without additional shifts

	array<map<uint64_t, vector<uint64_t>>, 16> Ldiff;
	for(uint isbox = 0; isbox < 16; isbox++){
		// uint64_t bx2i = bx2 & sboxMasks[isbox]; //Basis for the input of i-th sbox
		// uint64_t by2i = by2 & sboxMasks[isbox]; //Basis for the input of i-th sbox
		auto & Ldiff_isbox = Ldiff[isbox];

		//Get the unit vectors (as uint64_t) in the basis for this sbox
		vector<uint64_t> vecbx2i;
		vector<uint64_t> vecby2i;
		uint64_t offset = 4*isbox;
		for(uint i = 0; i < 4; i++){
			uint64_t index = offset+i;
			if(bx2 & (1ULL << index))
				vecbx2i.emplace_back(1ULL << index);
			if(by2 & (1ULL << index))
				vecby2i.emplace_back(1ULL << index);
		}

		if((vecbx2i.size() == 0 && vecby2i.size() != 0) || (vecbx2i.size() != 0 && vecby2i.size() == 0)){
			//Should not happen with how the basis are generated but to be safe
			//In this case, incompatible basis, so no trails can be generated
			return;
		}
		if(vecbx2i.size() == 0 && vecby2i.size() == 0){
			//This means that the Sbox is inactive, so the only io difference pair is (0,0)
			Ldiff_isbox[0].emplace_back(0);
		}
		else{
			//Generate all possible io pairs for this sbox
			uint64_t in_bound = (1ULL << vecbx2i.size());
			uint64_t out_bound = (1ULL << vecby2i.size());
			for(uint64_t ix2 = 1; ix2 < in_bound; ix2++){
				//bit i of ix2 = 1 <=> current differential x2 uses the vector vecbx2i[i]
				uint64_t x2 = 0;
				for(uint i = 0; i < vecbx2i.size(); i++){
					if(ix2 & (1ULL << i))
						x2 |= vecbx2i[i];
				}
				//x2 is the shifted difference (i.e. align to the sbox)
				//get the value to check against the DDT
				uint64_t valx2 = (x2 >> offset);

				for(uint64_t iy2 = 1; iy2 < out_bound; iy2++){
					//bit i of iy2 = 1 <=> current differential y2 uses the vector vecby2i[i]
					uint64_t y2 = 0;
					for(uint i = 0; i < vecby2i.size(); i++){
						if(iy2 & (1ULL << i))
							y2 |= vecby2i[i];
					}
					//Again, y2 is the shifted value, get the actual value
					uint64_t valy2 = (y2 >> offset);

					if(ddt[valx2][valy2] != -1){//Valid differential
						Ldiff_isbox[x2].emplace_back(y2);
					}
				}
			}
		}

		//At this point all possible io pairs for this sbox have been generated
		//Check if there is at least one, if not, then no trail can be generated so stop here
		//Should not happen because of the previous checks, but to be safe
		if(Ldiff_isbox.size() == 0){
			// cout << "No possible transition on Sbox " << isbox << endl;
			return;
		}
	}
	//We have the list of all possible transitions for each sbox, and each sbox has at least one valid transition

	// //Time to go through all of them to see how many verify the activity patterns
	uint currentIndex_x2 = 0;
	uint currentIndex_y2 = 0;
	uint64_t x2 = 0;
	uint64_t y2 = 0;
	double currentWeight = 0;
	//We know the number of active sboxes in R1 and R3, so we have a lower bound on the minimum weight for those rounds

	enumerateCorePatterns3RFromCore(Ldiff, ay1, ax3, nbSbox, histo, currentIndex_x2, currentIndex_y2, x2, y2, currentWeight);
}

void PresentData::enumerateCorePatterns3RFromCore(array<map<uint64_t, vector<uint64_t>>, 16> const & Ldiff,
												  uint64_t const ay1,
												  uint64_t const ax3,
												  uint64_t const nbSbox,
												  mapPairStateHisto & histo,
												  uint const currentIndex_x2,
												  uint const currentIndex_y2,
												  uint64_t x2,
												  uint64_t y2,
												  double currentWeight){
	//Recursively build all possible trails from the transitions in Ldiff, 
	//currentIndex_x2/y2 indicates the index of the sbox that needs to be guessed for the input/output, currentIndex = 16 marks the end of the recursion for the input/output (respectively)
	//Unlike for the boxWeight/corePattern enumeration, here we will build both input/output of a given sbox before moving to the next one
	//i.e. for boxWeight/corePattern, the indexes for x2/y2 was iterated over as x2 : 0 1 2 ... 15 16 then y2 : 0 1 2 ... 15 16
	//Here, we first do x2 : 0, y2 : 0 then x2 : 1, y2 : 1 etc.
	//Goal is to bound the weight of the trail early to cut off if necessary
	if(currentIndex_x2 == 16 && currentIndex_y2 == 16){
		//x2 and y2 are fully determined, and round 2 has weight currentWeight
		//For trail clustering we need to actually enumerate all input/output differences
		//We can reuse the implementation for 2 rounds, just with a higher initial currentWeight (corrseponding here to the weight in the middle round)
		uint64_t x1 = 0;
		uint64_t y1 = applyInvPerm(x2);
		uint64_t x3 = applyPerm(y2);
		uint64_t y3 = 0;
		uint16_t patternx1 = getActivityPattern(y1); //Activity pattern of xi is the same as yi
		uint16_t patterny3 = getActivityPattern(x3);
		exhaustTrailsForCluster2RFromTrailCore(x1,y1,x3,y3,patternx1,patterny3,currentWeight,histo);
	}

	else{
		if(currentIndex_x2 == currentIndex_y2){ //Guess on the input x2

			uint64_t prev_x2 = x2;
			for(auto const & x2i_o : Ldiff[currentIndex_x2]){
				x2 = prev_x2 | x2i_o.first; //Already shifted

				if(currentIndex_x2 < 15 || 
				  (currentIndex_x2 == 15 && getActivityPattern(applyInvPerm(x2)) == ay1)){
					enumerateCorePatterns3RFromCore(Ldiff, ay1, ax3, nbSbox, histo, currentIndex_x2+1, currentIndex_y2, x2, y2,currentWeight);
				}
			}
		}
		else{ //Guess on the output y2

			//Get the input diff
			uint64_t offset = 4*currentIndex_y2;
			uint64_t const x2i = x2 & (0xFULL << offset);
			uint64_t x2isbox = x2i >> offset; //value of x2i between 0 and 15 to check ddt
			uint64_t prev_y2 = y2;
			double prev_currentWeight = currentWeight;

			for(auto const y2i : Ldiff[currentIndex_y2].at(x2i)){ //.at for const access
				y2 = prev_y2 | y2i; //Already shifted
				uint64_t y2isbox = y2i >> offset; //value of y2i between 0 and 15 to check ddt
				currentWeight = prev_currentWeight + ddt[x2isbox][y2isbox];

				if( (currentIndex_y2 < 15 || 
					(currentIndex_y2 == 15 && getActivityPattern(applyPerm(y2)) == ax3))){
					enumerateCorePatterns3RFromCore(Ldiff, ay1, ax3, nbSbox, histo, currentIndex_x2, currentIndex_y2+1, x2, y2,currentWeight);
				}
			}
		}
	}
}
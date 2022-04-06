#include "PresentData.hpp"

using namespace std;
typedef unsigned int uint;

std::map<uint64_t, uint64_t> PresentData::corePatternHistogram4R(uint const nbSboxMax){
	return corePatternHistogram4R_internal(nbSboxMax,false);
}
std::map<uint64_t, uint64_t> PresentData::boxWeightHistogram4R(uint const nbSboxMax){
	return corePatternHistogram4R_internal(nbSboxMax,true);
}

map<uint64_t, uint64_t> PresentData::corePatternHistogram4R_internal(uint const nbSboxMax,
																	 bool const boxWeightHisto){
	//More naive search for core patterns that's using a lot less memory than the previous implementation
	//Could be faster because a lot less memory management
	//Main idea is to generate core patterns for 3 rounds but extend them over 4 rounds

	//x1 -S-> y1 -P-> x2 -S-> y2 -P-> x3 -S-> y3 -P-> x4 -S-> y4

	map<uint64_t, uint64_t> histo;
	if(nbSboxMax <= 3) return histo;


	//Generate all active patterns for y1 i.e. sbox i = 1111 if sbox i active
	//Compute x2 = p(y1)
	//Compute the corresponding active pattern patx2 (same as above)
	//Store (x2, patx2, ay1) in Tx2[i] for each i such that sbox i is active in x2
	//We need to keep ay1 as to be sure that we avoid duplicates later on, so that the actual core trail does follow ay1
	vector<vector<array<uint64_t,3>>> Tx2(16);

	//First round so at least ceil(wy1/4) sboxes in the second round + 2 for third and fourth
	for(uint wy1 = 1; wy1+ceil(double(wy1)/4)+2 <= nbSboxMax; wy1++){
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
	//- at least ceil(wx3/4) sboxes in the fourth round
	//- at least 1 sbox in the first round
	for(uint wx3 = 1; wx3+2*ceil(double(wx3)/4)+1 <= nbSboxMax; wx3++){
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
							if(__builtin_popcountll(getActivityPattern(mask)) <= nbSboxMax-3){

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
									uint nbSboxR3 = __builtin_popcountll(ax3);
									uint nbSbox = __builtin_popcountll(ay1) + countActiveSbox(mask) + nbSboxR3;
									//Since we have nbSboxR3 sboxes in the third round, 
									//the 4th round will contain at least ceil(nbsboxR3/4) sboxes
									if(nbSbox+ceil(double(nbSboxR3)/4) <= nbSboxMax &&
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
					enumerateCorePatterns4R(bx2_ay1_by2_nbSbox, ax3, histo_thread, nbSboxMax, boxWeightHisto);
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

void PresentData::enumerateCorePatterns4R(array<uint64_t,4> const & bx2_ay1_by2_nbSbox,
										  uint64_t const ax3,
										  std::map<uint64_t, uint64_t> & histo,
										  uint const nbSboxMax,
										  bool const boxWeightHisto){

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
	enumeratePartialCorePatterns3R(Ldiff, ay1, ax3, nbSbox, histo, currentIndex_x2, currentIndex_y2, x2, y2, nbSboxMax, boxWeightHisto);

}

void PresentData::enumeratePartialCorePatterns3R(array<map<uint64_t, vector<uint64_t>>, 16> const & Ldiff,
												  uint64_t const ay1,
												  uint64_t const ax3,
												  uint64_t const nbSboxR123,
												  std::map<uint64_t, uint64_t> & histo,
												  uint const currentIndex_x2,
												  uint const currentIndex_y2,
												  uint64_t x2,
												  uint64_t y2,
												  uint const nbSboxMax,
												  bool const boxWeightHisto){
	//Recursively build all possible trails from the transitions in Ldiff, 
	//and count how many verify the activity patterns ay1 and ax3
	//currentIndex_x2/y2 indicates the index of the sbox that needs to be guessed for the input/output, currentIndex = 16 marks the end of the recursion for the input/output (respectively)
	if(currentIndex_x2 == 16 && currentIndex_y2 == 16){
		//x2 and y2 are fully determined
		//Due to the workflow, it's already been checked that x2 verifies ay1 and y2 verifies ax3
		//So we have a valid partial core over the first 3 rounds
		//Get x3
		auto y1 = applyInvPerm(x2);
		auto x3 = applyPerm(y2);
		auto patterny3 = getActivityPattern(x3); //Activity pattern of xi is the same as yi
		uint64_t y3 = 0;
		uint64_t x4 = 0;
		//Get through all possible resulting y3 to determine the core patterns
		extendCore4r(y1,x3,y3,patterny3,x4,nbSboxR123,nbSboxMax,histo,boxWeightHisto);
		
	}

	else{
		if(currentIndex_x2 < 16){ //Still need to guess the input x2

			uint64_t prev_x2 = x2;
			for(auto const & x2i_o : Ldiff[currentIndex_x2]){
				x2 = prev_x2 | x2i_o.first; //Already shifted

				//Either we will still have some guesses remaining so go forward
				// or we just finished guessing x2, check if it verifies ay1
				if(currentIndex_x2 < 15 || (currentIndex_x2 == 15 && getActivityPattern(applyInvPerm(x2)) == ay1)){
					enumeratePartialCorePatterns3R(Ldiff, ay1, ax3, nbSboxR123, histo, currentIndex_x2+1, currentIndex_y2, x2, y2,nbSboxMax,boxWeightHisto);
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
					enumeratePartialCorePatterns3R(Ldiff, ay1, ax3, nbSboxR123, histo, currentIndex_x2, currentIndex_y2+1, x2, y2,nbSboxMax,boxWeightHisto);
				}
			}
		}
	}
}

void PresentData::extendCore4r(uint64_t const y1,
							   uint64_t const x3,
							   uint64_t y3,
							   uint16_t patterny3,
							   uint64_t x4,
							   uint64_t const nbSboxR123,
							   uint64_t const nbSboxMax,
							   std::map<uint64_t, uint64_t> & histo,
							   bool const boxWeightHisto){
	if(patterny3 == 0){ //All guesses done, so we know both y3 and x4
		uint64_t totalNumberSbox = nbSboxR123 + countActiveSbox(x4);
		if(boxWeightHisto){
			//Only need to convolute the number of trails
			//Probability doesn't matter here, so arbitrary bound that should include any trail
			auto histoR1 = convolutionInvRound(y1,2*16*4);
			auto histoR4 = convolutionRound(x4,2*16*4);
			uint64_t nbTrailR1 = 0;
			uint64_t nbTrailR4 = 0;
			for(auto const & p_ctrR1 : histoR1)
				nbTrailR1 += p_ctrR1.second;
			for(auto const & p_ctrR4 : histoR4)
				nbTrailR4 += p_ctrR4.second;
			histo[totalNumberSbox] += nbTrailR1*nbTrailR4;
		}
		else{
			histo[totalNumberSbox]++;
		}
	}
	else{ //Still need to guess stuff
		//Get the index of the first stbox that needs to be guessed
		uint isbox = __builtin_ctz(patterny3);
		patterny3 ^= (1 << isbox); //Mark this sbox as being guessed
		uint offset = 4*isbox;

		uint64_t const x3i = getNibble(x3,isbox);
		auto const & ddtPropagation_x3i = ddtPropagation[x3i];

		//Backup value to restore after the guess
		uint64_t prev_x4 = x4;

		for(auto const val : ddtPropagation_x3i){
			x4 = prev_x4;

			y3 = setNibble(y3,isbox,val);
			//Apply the permutation on those 4 bits
			for(uint i = 0; i < 4; i++){
				if(val & (1ULL << i))
					x4 |= (1ULL << p[offset+i]);
			}

			//Check that we don't go over the bound
			if(nbSboxR123 + countActiveSbox(x4) <= nbSboxMax)
				extendCore4r(y1,x3,y3,patterny3,x4,nbSboxR123,nbSboxMax,histo,boxWeightHisto);
		}
	}
}



map<double, uint64_t> PresentData::trailHistogram4R(uint const nbSboxMax){
	//More naive search for trails that's using a lot less memory than the previous implementation
	//Could be faster because a lot less memory management
	//Main idea is to generate core patterns for 3 rounds but extend them over 4 rounds

	//x1 -S-> y1 -P-> x2 -S-> y2 -P-> x3 -S-> y3 -P-> x4 -S-> y4

	map<double, uint64_t> histo;
	if(nbSboxMax <= 3) return histo;

	double boundWeight = nbSboxMax*minWeightSbox;


	//Generate all active patterns for y1 i.e. sbox i = 1111 if sbox i active
	//Compute x2 = p(y1)
	//Compute the corresponding active pattern patx2 (same as above)
	//Store (x2, patx2, ay1) in Tx2[i] for each i such that sbox i is active in x2
	//We need to keep ay1 as to be sure that we avoid duplicates later on, so that the actual core trail does follow ay1
	vector<vector<array<uint64_t,3>>> Tx2(16);

	//First round so at least ceil(wy1/4) sboxes in the second round + 2 for third and fourth
	for(uint wy1 = 1; wy1+ceil(double(wy1)/4)+2 <= nbSboxMax; wy1++){
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
	//- at least ceil(wx3/4) sboxes in the fourth round
	//- at least 1 sbox in the first round
	for(uint wx3 = 1; wx3+2*ceil(double(wx3)/4)+1 <= nbSboxMax; wx3++){
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
							if(__builtin_popcountll(getActivityPattern(mask)) <= nbSboxMax-3){

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
									uint nbSboxR3 = __builtin_popcountll(ax3);
									uint nbSboxR1 = __builtin_popcountll(ay1);
									uint nbSbox = nbSboxR1 + countActiveSbox(mask) + nbSboxR3;
									//Since we have nbSboxR3 sboxes in the third round, 
									//the 4th round will contain at least ceil(nbsboxR3/4) sboxes
									if(nbSbox+ceil(double(nbSboxR3)/4) <= nbSboxMax &&
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
					enumerateCorePatterns4R(bx2_ay1_by2_nbSbox, ax3, histo_thread, nbSboxMax, boundWeight);
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

void PresentData::enumerateCorePatterns4R(array<uint64_t,4> const & bx2_ay1_by2_nbSbox,
										  uint64_t const ax3,
										  std::map<double, uint64_t> & histo,
										  uint const nbSboxMax,
										  double const boundWeight){

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
	enumeratePartialCorePatterns3R(Ldiff, ay1, ax3, nbSbox, histo, currentIndex_x2, currentIndex_y2, x2, y2, nbSboxMax,currentWeight, minWeightR1R3, boundWeight);

}

void PresentData::enumeratePartialCorePatterns3R(array<map<uint64_t, vector<uint64_t>>, 16> const & Ldiff,
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
												  double const boundWeight){
	//Recursively build all possible trails from the transitions in Ldiff, 
	//and count how many verify the activity patterns ay1 and ax3
	//currentIndex_x2/y2 indicates the index of the sbox that needs to be guessed for the input/output, currentIndex = 16 marks the end of the recursion for the input/output (respectively)
	if(currentIndex_x2 == 16 && currentIndex_y2 == 16){
		//x2 and y2 are fully determined
		//Due to the workflow, it's already been checked that x2 verifies ay1 and y2 verifies ax3
		//So we have a valid partial core over the first 3 rounds
		//We can already compute the convolution for the first round
		auto histoR1 = convolutionInvRound(applyInvPerm(x2),boundWeight);
		//Get x3
		auto x3 = applyPerm(y2);
		auto patterny3 = getActivityPattern(x3); //Activity pattern of xi is the same as yi
		uint64_t y3 = 0;
		uint64_t x4 = 0;
		//Get through all possible resulting y3 to determine the core patterns
		extendCore4r(x3,y3,patterny3,x4,nbSboxR123,nbSboxMax,histo,histoR1,currentWeight,boundWeight);
		
	}

	else{
		if(currentIndex_x2 == currentIndex_y2){ //Guess on the input x2

			uint64_t prev_x2 = x2;
			for(auto const & x2i_o : Ldiff[currentIndex_x2]){
				x2 = prev_x2 | x2i_o.first; //Already shifted

				if(currentIndex_x2 < 15 || 
				  (currentIndex_x2 == 15 && getActivityPattern(applyInvPerm(x2)) == ay1)){
					enumeratePartialCorePatterns3R(Ldiff, ay1, ax3, nbSboxR123, histo, currentIndex_x2+1, currentIndex_y2, x2, y2,nbSboxMax,currentWeight, minWeightR1R3, boundWeight);
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
					enumeratePartialCorePatterns3R(Ldiff, ay1, ax3, nbSboxR123, histo, currentIndex_x2, currentIndex_y2+1, x2, y2,nbSboxMax,currentWeight, minWeightR1R3, boundWeight);
				}
			}
		}
	}
}

void PresentData::extendCore4r(uint64_t const x3,
							   uint64_t y3,
							   uint16_t patterny3,
							   uint64_t x4,
							   uint64_t const nbSboxR123,
							   uint64_t const nbSboxMax,
							   std::map<double, uint64_t> & histo,
							   std::map<double, uint64_t> const & histoR1,
							   double currentWeight,
							   double const boundWeight){
	if(patterny3 == 0){ //All guesses done, so we know both y3 and x4
		//currentWeight holds the weight for round 2 + round 3
		//histoR1 has the convolution for the first round
		//just need to get the convolution for the last round and merge
		auto histoR4 = convolutionRound(x4,boundWeight);
		for(auto const p_ctrR1 : histoR1){
			double weightR123 = currentWeight+p_ctrR1.first;
			if(weightR123 <= boundWeight){
				for(auto const & p_ctrR4 : histoR4){
					double endWeight = weightR123+p_ctrR4.first;
					if(endWeight <= boundWeight)
						histo[endWeight] += p_ctrR1.second*p_ctrR4.second;
				}
			}
		}
	}
	else{ //Still need to guess stuff
		//Get the index of the first stbox that needs to be guessed
		uint isbox = __builtin_ctz(patterny3);
		patterny3 ^= (1 << isbox); //Mark this sbox as being guessed
		uint offset = 4*isbox;

		uint64_t const x3i = getNibble(x3,isbox);
		auto const & ddtPropagation_x3i = ddtPropagation[x3i];

		//Backup value to restore after the guess
		uint64_t prev_x4 = x4;
		double prev_currentWeight = currentWeight;

		for(auto const val : ddtPropagation_x3i){
			x4 = prev_x4;
			currentWeight = prev_currentWeight + ddt[x3i][val];

			y3 = setNibble(y3,isbox,val);
			//Apply the permutation on those 4 bits
			for(uint i = 0; i < 4; i++){
				if(val & (1ULL << i))
					x4 |= (1ULL << p[offset+i]);
			}

			//Check that we don't go over the bounds
			if((nbSboxR123 + countActiveSbox(x4) <= nbSboxMax) &&
			   (histoR1.begin()->first + currentWeight <= boundWeight)){ //Smallest weight in round 1 + weight on round 2 + weight on round 3 so far
				extendCore4r(x3,y3,patterny3,x4,nbSboxR123,nbSboxMax,histo,histoR1,currentWeight,boundWeight);
			}
		}
	}
}

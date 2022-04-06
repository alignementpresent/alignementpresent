#include "PresentData.hpp"

using namespace std;
typedef unsigned int uint;

map<double, uint64_t> PresentData::trailHistogram2R(uint const nbSboxMax){
	//Return a map m such that m[p] = the number of trails with weight p, i.e. proba 2^(-p) over 2 rounds
	//Limit the computations to trails with at most #nbSboxMax active sboxes

	
	// Since the linear layer is just a permutation, if there are n active sboxes in the first round, there are at least ceil(n/4) active sboxes in the second round
	// For example, assume that there are 14 active sboxes in the first round
	// The output weight is thus at least 14 (1 bit per sbox)
	// After the permutation, the weight remains the same, so 14 active bits leads to at least 4 active sboxes (n/4 = 3.5 -> ceil(n/4) = 4)
	
	//x1 -S-> y1 -P-> x2 -S-> y2

	//Bound for the maximum accepted weight of any single trail
	//If we would allow trails with a higher weight while keeping the same limit on the number of active sboxes, then we would have incomplete results, i.e. find X number of trais of weight W, but there are actually Y > X trails with weight W when we allow more active sboxes.
	//This gives *exact* histograms and also reduces the complexity
	double boundWeight = nbSboxMax*minWeightSbox;
	cout << "Bound on probability for 2 rounds : " << boundWeight << endl;

	map<double, uint64_t> histo;
	cout << "Computing for" << flush;
	for(uint nbSboxR1 = 1; nbSboxR1+ceil(double(nbSboxR1)/4) <= nbSboxMax; nbSboxR1++){
		cout << " " << nbSboxR1 << flush;
		auto const & uint16WeightMap_nbSboxR1 = uint16WeightMap[nbSboxR1];
		uint64_t bound = uint16WeightMap_nbSboxR1.size();

		#pragma omp parallel
		{	
			map<double, uint64_t> histo_thread;

			#pragma omp for schedule(dynamic,1)
			for(uint64_t i_patterny1 = 0; i_patterny1 < bound; i_patterny1++){
				auto const patterny1 = uint16WeightMap_nbSboxR1[i_patterny1];

				//Exhaust all trails over 2 rounds satisfying the sbox activity pattern on y1 given by #patterny1 and put the results in the histrogram
				//Only count trails that have at most #nbSboxMax active sboxes, and a weight of up to boundWeight

				//Initialize variables for the recursive calls
				uint64_t y1 = 0;
				uint64_t x2 = 0;
				double currentMinWeight_y1 = 0;
				exhaustTrails2R(y1,x2,patterny1,nbSboxR1,nbSboxMax,boundWeight,currentMinWeight_y1,histo_thread);
			}

			#pragma omp critical
			{
				for(auto const & p_ctr : histo_thread)
					histo[p_ctr.first] += p_ctr.second;
			}
		}
	}
	cout << endl;

	return histo;
}

void PresentData::exhaustTrails2R(uint64_t y1,
								  uint64_t x2,
								  uint16_t patterny1,
								  uint const nbSboxR1,
								  uint const nbSboxMax,
								  double const boundWeight,
								  double currentMinWeight_y1,
								  map<double,uint64_t> & histo){
	if(patterny1 == 0){ 
		//Means that all the values have been guessed for y1
		//Only need to convolute the number of trails
		auto histoR1 = convolutionInvRound(y1,boundWeight);
		auto histoR2 = convolutionRound(x2,boundWeight);
		for(auto const & p_ctrR1 : histoR1){
			for(auto const & p_ctrR2 : histoR2){
				if(p_ctrR1.first + p_ctrR2.first <= boundWeight)
					histo[p_ctrR1.first + p_ctrR2.first] += p_ctrR1.second * p_ctrR2.second;
			}
		}
	}
	else{ //Still need to guess stuff
		//Get the index of the first active Sbox that needs to be guessed
		uint isbox = __builtin_ctz(patterny1);
		patterny1 ^= (1 << isbox); //Mark this sbox as being guessed
		uint offset = 4*isbox;

		//Backup some values to restore them after a guess
		//Don't backup y1 as it's easy to just change on the fly
		uint64_t prev_x2 = x2;
		double prev_currentMinWeight_y1 = currentMinWeight_y1;

		//Go through all possible values for this sbox
		for(uint64_t val = 1; val < 16; val++){
			//Restore values before the guess
			x2 = prev_x2;
			currentMinWeight_y1 = prev_currentMinWeight_y1;

			y1 = setNibble(y1,isbox,val);
			currentMinWeight_y1 += minProbaInvSbox[val]; //update the best possible weight on y1 with the current guesses

			//Apply the permutation on those specific 4 bits for the next guess value
			for(uint i = 0; i < 4; i++){
				if(val & (1ULL << i))
					x2 |= (1ULL << p[offset+i]);
			}

			//Check if we don't exceed the number of active sboxes or the bound on the proba over the first round already
			if(nbSboxR1 + countActiveSbox(x2) <= nbSboxMax && currentMinWeight_y1 <= boundWeight){
				exhaustTrails2R(y1,x2,patterny1,nbSboxR1,nbSboxMax,boundWeight,currentMinWeight_y1,histo);
			}
		}
	}
}

std::map<uint64_t, uint64_t> PresentData::corePatternHistogram2R(uint const nbSboxMax){
	return corePatternHistogram2R_internal(nbSboxMax,false);
}
std::map<uint64_t, uint64_t> PresentData::boxWeightHistogram2R(uint const nbSboxMax){
	return corePatternHistogram2R_internal(nbSboxMax,true);
}

map<uint64_t, uint64_t> PresentData::corePatternHistogram2R_internal(uint const nbSboxMax,
																	 bool const boxWeightHisto){
	//Return a map m such that m[x] = the number of core patterns with x active sboxes
	//Limit the computations to trails with at most #nbSboxMax active sboxes

	map<uint64_t, uint64_t> histo;
	cout << "Computing for" << flush;
	for(uint nbSboxR1 = 1; nbSboxR1+ceil(double(nbSboxR1)/4) <= nbSboxMax; nbSboxR1++){
		cout << " " << nbSboxR1 << flush;
		auto const & uint16WeightMap_nbSboxR1 = uint16WeightMap[nbSboxR1];
		uint64_t bound = uint16WeightMap_nbSboxR1.size();

		#pragma omp parallel
		{	
			map<uint64_t, uint64_t> histo_thread;

			#pragma omp for schedule(dynamic,1)
			for(uint64_t i_patterny1 = 0; i_patterny1 < bound; i_patterny1++){
				auto const patterny1 = uint16WeightMap_nbSboxR1[i_patterny1];

				//Exhaust all trails over 2 rounds satisfying the sbox activity pattern on y1 given by #patterny1 and put the results in the histrogram
				//Only count trails that have at most #nbSboxMax active sboxes

				//Initialize variables for the recursive calls
				uint64_t y1 = 0;
				uint64_t x2 = 0;
				corePatternExhaustTrails2R(y1,x2,patterny1,nbSboxR1,nbSboxMax,histo_thread,boxWeightHisto);
			}

			#pragma omp critical
			{
				for(auto const & p_ctr : histo_thread)
					histo[p_ctr.first] += p_ctr.second;
			}
		}
	}
	cout << endl;

	return histo;
}

void PresentData::corePatternExhaustTrails2R(uint64_t y1,
											 uint64_t x2,
											 uint16_t patterny1,
											 uint const nbSboxR1,
											 uint const nbSboxMax,
											 map<uint64_t,uint64_t> & histo,
											 bool const boxWeightHisto){
	if(patterny1 == 0){ 
		//Means that all the values have been guessed for y1
		uint64_t totalNumberSbox = nbSboxR1 + countActiveSbox(x2);

		if(boxWeightHisto){
			//Only need to convolute the number of trails
			//Probability doesn't matter here, so arbitrary bound that should include any trail
			auto histoR1 = convolutionInvRound(y1,2*16*4);
			auto histoR2 = convolutionRound(x2,2*16*4);

			//Number of trails in the first and second round
			uint64_t nbTrailR1 = 0;
			uint64_t nbTrailR2 = 0;
			for(auto const & p_ctrR1 : histoR1)
				nbTrailR1 += p_ctrR1.second;
			for(auto const & p_ctrR2 : histoR2)
				nbTrailR2 += p_ctrR2.second;

			//Total number of trails for this value of y1
			histo[totalNumberSbox] += nbTrailR1*nbTrailR2;
		}
		else{
			histo[totalNumberSbox]++;
		}
		
	}
	else{ //Still need to guess stuff
		//Get the index of the first active Sbox that needs to be guessed
		uint isbox = __builtin_ctz(patterny1);
		patterny1 ^= (1 << isbox); //Mark this sbox as being guessed
		uint offset = 4*isbox;

		//Backup some values to restore them after a guess
		//Don't backup y1 as it's easy to just change on the fly
		uint64_t prev_x2 = x2;

		//Go through all possible values for this sbox
		for(uint8_t val = 1; val < 16; val++){
			x2 = prev_x2;

			y1 = setNibble(y1,isbox,val);
			//Apply the permutation on those specific 4 bits for the next guess value
			for(uint i = 0; i < 4; i++){
				if(val & (1ULL << i))
					x2 |= (1ULL << p[offset+i]);
			}

			//Check if we don't exceed the number of active sboxes
			if(nbSboxR1 + countActiveSbox(x2) <= nbSboxMax){
				corePatternExhaustTrails2R(y1,x2,patterny1,nbSboxR1,nbSboxMax,histo,boxWeightHisto);
			}
		}
	}
}

pair<uint64_t,unordered_set<uint64_t>> PresentData::getClusterClassSize(uint16_t const patterna, 
										  uint16_t const patternPa){
	//Return the size of the cluster class for a difference with pattern patterna
	// patternPa (= pattern of P(a)) is given as parameter as it will already be computed in practice

	// uint16_t patterna  = getActivityPattern(a);
	// uint16_t patternPa = getActivityPattern(Pa);

	//We are here considering a special case for the algorithm described in Section A.3 of the superbox paper
	//Rather than a generic linear mapping, here we only have a permutation
	//This means that both basis for P(V(a)) and V(P(a)) will only be unit vectors
	//As such, we don't need to use Zassenhaus' algorithm to compute the intersection of the two subspaces
	//Just need to compute the intersection of the two basis

	//Generate the set for a basis of P(V(a))
	set<uint64_t> PVa;
	uint16_t tmp_pattern = patterna;
	while(tmp_pattern != 0){
		uint isbox = __builtin_ctz(tmp_pattern);
		tmp_pattern ^= (1 << isbox);
		for(uint j = 0; j < 4; j++)
			PVa.emplace(1ULL << p[4*isbox+j]);
	}

	//Generate the set for V(P(a))
	set<uint64_t> VPa;
	tmp_pattern = patternPa;
	while(tmp_pattern != 0){
		uint isbox = __builtin_ctz(tmp_pattern);
		tmp_pattern ^= (1 << isbox);
		for(uint j = 0; j < 4; j++)
			VPa.emplace(1ULL << (4*isbox+j));
	}

	// auto basisInter = computeBasisIntersection(PVa,VPa);
	vector<uint64_t> basisInter;
	set_intersection(PVa.begin(), PVa.end(), VPa.begin(), VPa.end(),
                 std::back_inserter(basisInter));

	//Count how many elements v satisfy
	//v has the same activity pattern as Pa *and*
	//invP(v) has the same activity pattern as a

	uint64_t ctr = 0;
	uint64_t bound = (1ULL << basisInter.size());
	unordered_set<uint64_t> elements;
	for(uint64_t x = 0; x < bound; x++){
		//v = sum(basis[i] s.t. bit i of x = 1)
		uint64_t v = 0;
		for(uint i = 0; i < basisInter.size(); i++){
			if(x & (1ULL << i))
				v ^= basisInter[i];
		}

		//Check if v has the same activity pattern as P(a)
		uint16_t patternv = 0;
		for(uint i = 0; i < 16; i++){
			if(v & sboxMasks[i])
				patternv |= (1 << i);
		}
		if(patternv == patternPa){
			//Check if invP(v) has the same activity pattern as a
			uint64_t invPv = applyInvPerm(v);
			uint16_t pattern_invPv = 0;
			for(uint i = 0; i < 16; i++){
				if(invPv & sboxMasks[i])
					pattern_invPv |= (1 << i);
			}

			if(pattern_invPv == patterna){
				ctr++;
				elements.emplace(v);
			}
		}
	}

	return make_pair(ctr,elements);

}

map<uint64_t, map<uint64_t, uint64_t>> PresentData::clusterHistogram(uint const nbSboxMax){
	//Return the cluster histogram of the permutation as defined in Section 6.1 of the superbox paper
	//Limit the computations to trails with at most #nbSboxMax active sboxes

	map<uint64_t, map<uint64_t, uint64_t>> histo;
	cout << "Computing for" << flush;
	for(uint nbSboxR1 = 1; nbSboxR1+ceil(double(nbSboxR1)/4) <= nbSboxMax; nbSboxR1++){
		cout << " " << nbSboxR1 << flush;
		auto const & uint16WeightMap_nbSboxR1 = uint16WeightMap[nbSboxR1];
		uint64_t bound = uint16WeightMap_nbSboxR1.size();

		#pragma omp parallel
		{	
			map<uint64_t, map<uint64_t, uint64_t>> histo_thread;

			#pragma omp for schedule(dynamic,1)
			for(uint64_t i_patterny1 = 0; i_patterny1 < bound; i_patterny1++){
				auto const patterny1 = uint16WeightMap_nbSboxR1[i_patterny1];

				//Exhaust all trails over 2 rounds satisfying the sbox activity pattern on y1 given by #patterny1 and put the results in the histrogram
				//Only count trails that have at most #nbSboxMax active sboxes

				//Initialize variables for the recursive calls
				uint64_t y1 = 0;
				uint64_t x2 = 0;
				uint16_t basePy1 = patterny1;
				unordered_map<uint16_t, uint64_t> knownClusterSize;
				unordered_set<uint64_t> knownElements;
				clusterExhaustTrails2R(y1,x2,patterny1,nbSboxR1,nbSboxMax,basePy1,histo_thread,knownClusterSize,knownElements);
			}

			#pragma omp critical
			{
				for(auto const & w_h : histo_thread){
					for(auto const & p_ctr : w_h.second)
						histo[w_h.first][p_ctr.first] += p_ctr.second;
			}
			}
		}
	}
	cout << endl;

	return histo;
}

void PresentData::clusterExhaustTrails2R(uint64_t y1,
										 uint64_t x2,
										 uint16_t patterny1,
										 uint const nbSboxR1,
										 uint const nbSboxMax,
										 uint16_t const basePy1,
										 map<uint64_t, map<uint64_t, uint64_t>> & histo,
										 unordered_map<uint16_t, uint64_t> & knownClusterSize,
										 unordered_set<uint64_t> & knownElements){
	if(patterny1 == 0){ 
		if(knownElements.find(x2) == knownElements.end()){
			//Means that all the values have been guessed for y1
			uint64_t totalNumberSbox = nbSboxR1 + countActiveSbox(x2);

			// uint16_t py1 = getActivityPattern(y1); //= basePy1
			uint16_t px2 = getActivityPattern(x2);
			if(knownClusterSize.find(px2) == knownClusterSize.end()){
				auto clusterSize_elements = getClusterClassSize(basePy1,px2);
				knownElements.insert(clusterSize_elements.second.begin(), clusterSize_elements.second.end());
				knownClusterSize[px2] = clusterSize_elements.first;
				histo[totalNumberSbox][clusterSize_elements.first]++;
			}
			else{
				histo[totalNumberSbox][knownClusterSize[px2]]++;
			}
		}

		
		
	}
	else{ //Still need to guess stuff
		//Get the index of the first active Sbox that needs to be guessed
		uint isbox = __builtin_ctz(patterny1);
		patterny1 ^= (1 << isbox); //Mark this sbox as being guessed
		uint offset = 4*isbox;

		//Backup some values to restore them after a guess
		//Don't backup y1 as it's easy to just change on the fly
		uint64_t prev_x2 = x2;

		//Go through all possible values for this sbox
		for(uint8_t val = 1; val < 16; val++){
			x2 = prev_x2;
			
			y1 = setNibble(y1,isbox,val);
			//Apply the permutation on those specific 4 bits for the next guess value
			for(uint i = 0; i < 4; i++){
				if(val & (1ULL << i))
					x2 |= (1ULL << p[offset+i]);
			}

			//Check if we don't exceed the number of active sboxes or the bound on the proba over the first round already
			if(nbSboxR1 + countActiveSbox(x2) <= nbSboxMax){
				clusterExhaustTrails2R(y1,x2,patterny1,nbSboxR1,nbSboxMax,basePy1,histo,knownClusterSize,knownElements);
			}
		}
	}
}

//-- Trail Clustering --
map<double, uint64_t>
PresentData::trailClusterHistogram2R(uint const nbSboxMax){
	//Return the histogram for trail clustering
	//Limit the computations to trails with at most #nbSboxMax active sboxes
	//In the case of 2 rounds, these are exact results as two trails with a different number of active sboxes cannot cluster together

	//Returned maps
	map<double,uint64_t> retmap;

	cout << "Computing for" << flush;
	for(uint nbSboxR1 = 1; nbSboxR1+ceil(double(nbSboxR1)/4) <= nbSboxMax; nbSboxR1++){
		cout << " " << nbSboxR1 << flush;
		auto const & uint16WeightMap_nbSboxR1 = uint16WeightMap[nbSboxR1];
		uint64_t bound = uint16WeightMap_nbSboxR1.size();

		#pragma omp parallel
		{

			map<double, uint64_t> histo_thread;

			#pragma omp for schedule(dynamic,1)
			for(uint64_t i_patterny1 = 0; i_patterny1 < bound; i_patterny1++){
				auto const patterny1 = uint16WeightMap_nbSboxR1[i_patterny1];

				//Exhaust all trails over 2 rounds satisfying the sbox activity pattern on y1 given by #patterny1 and put the results in the histrogram
				//Only count trails that have at most #nbSboxMax active sboxes

				//Initialize variables for the recursive calls
				mapPairStateHisto allTrails_thread;
				allTrails_thread.reserve(1 << 19);
				uint64_t y1 = 0;
				uint64_t x2 = 0;
				exhaustTrailsForCluster2R(y1,x2,patterny1,nbSboxR1,nbSboxMax,allTrails_thread);

				for(auto const & state_histo : allTrails_thread){
					double proba = 0;
					for(auto const & p_ctr : state_histo.second)
						proba += p_ctr.second * pow(2,-p_ctr.first);
					histo_thread[-log2(proba)]++;
				}
			}

			#pragma omp critical 
			{
				for(auto const & p_ctr : histo_thread)
					retmap[p_ctr.first] += p_ctr.second;
			}
		}
	}
	cout << endl;

	return retmap;
}

void PresentData::exhaustTrailsForCluster2R(uint64_t y1,
								  			uint64_t x2,
								  			uint16_t patterny1,
								  			uint const nbSboxR1,
								  			uint const nbSboxMax,
								  			mapPairStateHisto & allTrails_thread){
	if(patterny1 == 0){ //Means that all the values have been guessed for y1, so time to exhaust the possible values for x1 and y2

		//Initialize the variables for the recursive calls
		uint64_t x1 = 0;
		uint64_t y2 = 0;
		uint16_t patternx1 = getActivityPattern(y1); //Activity pattern of xi is the same as yi
		uint16_t patterny2 = getActivityPattern(x2); //Activity pattern of xi is the same as yi
		double currentWeight = 0;
		exhaustTrailsForCluster2RFromTrailCore(x1,y1,x2,y2,patternx1,patterny2,currentWeight,allTrails_thread);

	}
	else{ //Still need to guess stuff
		//Get the index of the first active Sbox that needs to be guessed
		uint isbox = __builtin_ctz(patterny1);
		patterny1 ^= (1 << isbox); //Mark this sbox as being guessed
		uint offset = 4*isbox;

		//Backup some values to restore them after a guess
		//Don't backup y1 as it's easy to just change on the fly
		uint64_t prev_x2 = x2;

		//Go through all possible values for this sbox
		for(uint8_t val = 1; val < 16; val++){
			x2 = prev_x2;
			
			y1 = setNibble(y1,isbox,val);
			//Apply the permutation on those specific 4 bits for the next guess value
			for(uint i = 0; i < 4; i++){
				if(val & (1ULL << i))
					x2 |= (1ULL << p[offset+i]);
			}

			//Check if we don't exceed the number of active sboxes or the bound on the proba over the first round already
			if(nbSboxR1 + countActiveSbox(x2) <= nbSboxMax){
				exhaustTrailsForCluster2R(y1,x2,patterny1,nbSboxR1,nbSboxMax,allTrails_thread);
			}
		}
	}
}


void PresentData::exhaustTrailsForCluster2RFromTrailCore(uint64_t x1,
											   			 uint64_t const y1,
											   			 uint64_t const x2,
											   			 uint64_t y2,
											   			 uint16_t patternx1,
											   			 uint16_t patterny2,
											   			 double currentWeight,
											   			 mapPairStateHisto & allTrails_thread){

	if(patternx1 == 0 && patterny2 == 0){ //All guesses done, increase the proba counter
		allTrails_thread[make_pair(x1,y2)][currentWeight]++;
	}
	else{
		if(patternx1){ //Still have some guesses to do in x1
			//Get the index of the first active Sbox that needs to be guessed
			uint isbox = __builtin_ctz(patternx1);
			patternx1 ^= (1 << isbox); //Mark this sbox as being guessed

			//Go through all possible values for this sbox
			uint64_t const y1i = getNibble(y1,isbox);
			auto const & ddtInvPropagation_y1i = ddtInvPropagation[y1i];

			double prev_currentWeight = currentWeight;
			for(auto const val : ddtInvPropagation_y1i){
				currentWeight = prev_currentWeight;
				x1 = setNibble(x1,isbox,val);
				currentWeight += ddt[val][y1i];

				exhaustTrailsForCluster2RFromTrailCore(x1,y1,x2,y2,patternx1,patterny2,currentWeight,allTrails_thread);
			}
		}
		else{ //Still have guesses to do in y2
			//Get the index of the first active Sbox that needs to be guessed
			uint isbox = __builtin_ctz(patterny2);
			patterny2 ^= (1 << isbox); //Mark this sbox as being guessed

			//Go through all possible values for this sbox
			uint64_t const x2i = getNibble(x2,isbox);
			auto const & ddtPropagation_x2i = ddtPropagation[x2i];

			double prev_currentWeight = currentWeight;
			for(auto const val : ddtPropagation_x2i){
				currentWeight = prev_currentWeight;
				y2 = setNibble(y2,isbox,val);
				currentWeight += ddt[x2i][val];

				exhaustTrailsForCluster2RFromTrailCore(x1,y1,x2,y2,patternx1,patterny2,currentWeight,allTrails_thread);
			}
		}
	}
}

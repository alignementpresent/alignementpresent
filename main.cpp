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

#include "aux_function.hpp"

using namespace std;
typedef unsigned int uint;



//10 iso + 1 perm aligned for each case


int main(){

	omp_set_num_threads(1);

	// //Aligned Permutations
	// cout << "Trails 2 rounds" << endl;
	// genData(HistoType::Trail,2,12,1,10,true,true,1);
	// cout << "Trails 3 rounds" << endl;
	// genData(HistoType::Trail,3,16,1,10,true,true,1);
	// cout << "Trails 4 rounds" << endl;
	// genData(HistoType::Trail,4,16,1,10,true,true,1);

	// cout << "Core pattern 2 rounds" << endl;
	// genData(HistoType::CorePattern,2,15,1,10,true,true,1);
	// cout << "Core pattern 3 rounds" << endl;
	// genData(HistoType::CorePattern,3,17,1,10,true,true,1);
	// cout << "Core pattern 4 rounds" << endl;
	// genData(HistoType::CorePattern,4,18,1,10,true,true,1);

	cout << "Cluster 2 rounds" << endl;
	genData(HistoType::ClusterTrail,2,6,10,1,true,true,0);
	cout << "Cluster 3 rounds" << endl;
	genData(HistoType::ClusterTrail,3,8,10,1,true,true,0);
	
	
	// //2 rounds cluster
	// cout << "2 rounds cluster, 1 perm per iso, all graph" << endl;
	// genData(HistoType::ClusterTrail,2,6,1,1,false,true,0);
	// cout << "2 rounds cluster, 346 perm, 1 iso, 1 graph" << endl;
	// genData(HistoType::ClusterTrail,2,6,346,1,false,true,1);
	// //3 rounds cluster
	// cout << "3 rounds cluster, 1 perm per iso, all graph" << endl;
	// genData(HistoType::ClusterTrail,3,8,1,1,false,true,0);
	// cout << "3 rounds cluster, 346 perm, 1 iso, 1 graph" << endl;
	// genData(HistoType::ClusterTrail,3,8,346,1,false,true,1);

	// //2 rounds core pattern
	// cout << "2 rounds core pattern, 1 perm per iso, all graph" << endl;
	// genData(HistoType::CorePattern,2,15,1,1,false,true,0);
	// cout << "2 rounds core pattern, 346 perm, 1 iso, 1 graph" << endl;
	// genData(HistoType::CorePattern,2,15,346,1,false,true,1);
	// //3 rounds core pattern
	// cout << "3 rounds core pattern, 1 perm per iso, all graph" << endl;
	// genData(HistoType::CorePattern,3,17,1,1,false,true,0);
	// cout << "3 rounds core pattern, 346 perm, 1 iso, 1 graph" << endl;
	// genData(HistoType::CorePattern,3,17,346,1,false,true,1);
	// //4 rounds core pattern
	// cout << "4 rounds core pattern, 1 perm per iso, all graph" << endl;
	// genData(HistoType::CorePattern,4,9,1,1,false,true,0);
	// cout << "4 rounds core pattern, 346 perm, 1 iso, 1 graph" << endl;
	// genData(HistoType::CorePattern,4,9,346,1,false,true,1);

	// //2 rounds trail
	// cout << "2 rounds trail, 1 perm per iso, all graph" << endl;
	// genData(HistoType::Trail,2,12,1,1,false,true,0);
	// cout << "2 rounds trail, 346 perm, 1 iso, 1 graph" << endl;
	// genData(HistoType::Trail,2,12,346,1,false,true,1);
	// //3 rounds trail
	// cout << "3 rounds trail, 1 perm per iso, all graph" << endl;
	// genData(HistoType::Trail,3,16,1,1,false,true,0);
	// cout << "3 rounds trail, 346 perm, 1 iso, 1 graph" << endl;
	// genData(HistoType::Trail,3,16,346,1,false,true,1);
	// //4 rounds trail
	// cout << "4 rounds trail, 1 perm per iso, all graph" << endl;
	// genData(HistoType::Trail,4,10,1,1,false,true,0);
	// cout << "4 rounds trail, 346 perm, 1 iso, 1 graph" << endl;
	// genData(HistoType::Trail,4,10,346,1,false,true,1);



	// HistoType htype = HistoType::BoxWeight;
	// uint nbRound = 3; 
	// uint nbSboxMax = 7;
	// uint nbPermPerIso = 1;
	// uint nbIso = 1;
	// bool withSSB = false;
	// bool randomize = true;
	// uint MaxNbGoodGraph = 0;
	// genData(htype,nbRound,nbSboxMax,nbPermPerIso,nbIso,withSSB,randomize,MaxNbGoodGraph);




	// auto histo = PD.clusterHistogram(8);
	// for(auto const & w_h : histo){
	// 	cout << w_h.first << " : ";
	// 	uint64_t ctr = 0;
	// 	for(auto const & s_ctr : w_h.second){
	// 		cout << "(";
	// 		cout << s_ctr.second << " x " << s_ctr.first;
	// 		cout << ") ";
	// 		ctr += s_ctr.first * s_ctr.second;
	// 	}
	// 	cout << "= " << ctr << endl;
	// }

	// getchar();

}
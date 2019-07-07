#ifndef TOURGROUP_H
#define TOUPGROUP_H
#include "Auxiliary.h"
#include <vector>
#include <string>
using namespace::std;

class PoI;
class TourGroup {
	private:
		int tourist_num;
		int tourist_count;
		float Budget;
		float** interest;
		vector<int> path;
	public:
		TourGroup(int num, float restrict);
		~TourGroup();
		void TourAppend(vector<string> &field);
		void GetTour(const int id);
		float SatisfactionTour(int tour_id, int poi_id, PoI& tmpPoI);
		float SatisfactionPoI(int poi_id, PoI& tmpPoI);
		float VarianceCalculate(float*& tmp_satis);
		float SatisfactionSum(float*& tmp_satis);
		void SatisfactionPlus(float*& satis, float*& tmp_satis, int id, PoI& tmpPoI);
		void SatisfactionMinus(float*& satis, float*& tmp_satis, int id, PoI& tmpPoI);
		float SatisfactionCalculate(vector<int> &tmppth, PoI& tmpPoI);
		float PhiCalculate(float*& tmp_satis);
		float GreedyAlgorithm(int source, int terminal, PoI& tmpPoI);
		float BRBumaAlgorithm(int source, int terminal, PoI& tmpPoI);
		float BRBumaBRPlusPlusAlgorithm(int source, int terminal, PoI& tmpPoI);
		float BRplusBumaAlgorithm(int soruce, int terminal, PoI& tmpPoI);
		float BRplusplusBumaAlgorithm(int source, int terminal, PoI& tmpPoI);
		float EsAlgorithm(int source, int terminal, PoI& tmpPoI);
		friend void SetNumberTour(int num, TourGroup& tmpTG);
};

#endif // TOURGROUP_H


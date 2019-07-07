#ifndef GRAPH_H
#define GRAPH_H
#include "Auxiliary.h"
#include <string>
using namespace::std;

class PoI {
	private:
		int poi_num;
		int poi_count;
		float** cost;
		struct PoIFORM {
			string name;
			theme node_type;
			float popularity;
			float poi_cost;
			PoIFORM()
			{
				name = "";
				node_type = Empty;
				popularity = 0;
				poi_cost = 0;
			}
		};
		PoIFORM* PoIList;

	public:
		PoI(int num);
		~PoI();
		void PoIAppend(vector<string>& field);
		void DistanceAppend(vector <string>& field);
		float CostLeadin(const int source, const int terminal);
		float CostLeadinS(const int source, const int terminal);
		float CostLeadin(const int source, const int middle, const int terminal);
		friend class TourGroup;
		friend void SetNumberPoI(int num, PoI &tmpPoI);
};

#endif // GRAPH_H

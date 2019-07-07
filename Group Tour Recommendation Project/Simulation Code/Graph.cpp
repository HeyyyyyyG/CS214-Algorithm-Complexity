#include "Graph.h"
#include <stddef.h>

using namespace::std;
PoI::PoI(int num)
{
	poi_num = num;
	poi_count = 0;
	cost = new float* [poi_num];
	for (int i = 0; i < poi_num; ++i)
		cost[i] = new float[poi_num];
	PoIList = new PoIFORM[poi_num];
}

PoI::~PoI()
{
	for (int i = 0; i < poi_count; ++i)
		delete[] cost[i];
	delete[] cost;
	delete[] PoIList;
}

void PoI::PoIAppend(vector<string>& field)
{	
	PoIList[poi_count].node_type = theme(atoi(field[2].c_str()));
	PoIList[poi_count].poi_cost = atof(field[11].c_str()) + 0.5;
	PoIList[poi_count].popularity = atof(field[10].c_str());
	++poi_count;
}

void PoI::DistanceAppend(vector<string>& field)
{	
	int lmid = atoi(field[0].c_str());
	int rmid = atoi(field[1].c_str());
	if (lmid > 11) lmid -= 1;
	if (rmid > 11) rmid -= 1;
	cost[lmid][rmid] = atof(field[2].c_str());
}

float PoI::CostLeadin(const int source, const int terminal)
{
	return (cost[source][terminal]
		+ (PoIList[source].poi_cost 
			+ PoIList[terminal].poi_cost) / 2);
}

float PoI::CostLeadinS(const int source, const int terminal)
{
	return (cost[source][terminal] + PoIList[terminal].poi_cost);
}

float PoI::CostLeadin(const int source, const int middle, const int terminal)
{
	return (CostLeadin(source, middle) 
		+ CostLeadin(middle, terminal) 
		- CostLeadin(source, terminal));
}

 void SetNumberPoI(int num, PoI& tmpPoI)
{
	 tmpPoI.poi_num = num;
}


#include "TourGroup.h"
#include "Graph.h"
#include <stddef.h>
#include <iostream>
using namespace::std;

TourGroup::TourGroup(int num, float restrict)
{
	tourist_num = num;
	tourist_count = 0;
	Budget = restrict;
	interest = new float* [tourist_num];
	for (int i = 0; i < tourist_num; ++i)
		interest[i] = new float[THEME_NUM];
}

TourGroup::~TourGroup()
{
	for (int i = 0; i < tourist_num; ++i)
		delete[] interest[i];
	delete[] interest;
}

void TourGroup::TourAppend(vector<string>& field)
{
	for (int i = 0; i < THEME_NUM; ++i)
		interest[tourist_count][i] = atof(field[i + 1].c_str());
	++tourist_count;
}

void TourGroup::GetTour(const int id)
{	
	cout << id << '\t';
	for (int i = 0; i < THEME_NUM; ++i)
		cout << interest[id][i] << '\t';
	cout << endl;
}

float TourGroup::SatisfactionTour(int tour_id, int poi_id, PoI& tmpPoI)
{
	int poi_theme = tmpPoI.PoIList[poi_id].node_type;
	float poi_popularity = tmpPoI.PoIList[poi_id].popularity;

	return (b1 * interest[tour_id][poi_theme] + b2 * poi_popularity);
}

float TourGroup::SatisfactionPoI(const int id, PoI& originPoI)
{
	float sum = 0;
	int poi_theme = originPoI.PoIList[id].node_type;
	float poi_popularity = originPoI.PoIList[id].popularity;

	sum += b2 * tourist_num * poi_popularity;

	for (int i = 0; i < tourist_num; ++i)
		sum += b1 * interest[i][poi_theme];

	return sum;
}

float TourGroup::SatisfactionCalculate(vector<int> &tmppath, PoI &tmpPoI)
{
	float sum = 0;

	for (vector<int>::iterator iter = tmppath.begin(); iter != tmppath.end(); ++iter)
		sum += SatisfactionPoI(*iter, tmpPoI);

	return sum;
}

float TourGroup::SatisfactionSum(float*& tmp_satis)
{
	float sum = 0;
	
	for (int i = 0; i < tourist_num; ++i)
		sum += tmp_satis[i];

	return sum;
}

float TourGroup::VarianceCalculate(float*& tmp_satis)
{
	float average = 0, variance = 0;

	average = SatisfactionSum(tmp_satis) / tourist_num;
	for (int i = 0; i < tourist_num; ++i)
		variance += (tmp_satis[i] - average) * (tmp_satis[i] - average);

	return (variance / tourist_num);
	
}



float TourGroup::GreedyAlgorithm(int source, int terminal, PoI &tmpPoI)
{
	float satisfaction = 0, totalcost = 0, tmpsatis = 0, tmpcost = 0, ratio = 0, tmpratio = 0;
	int tempcost = 0, tmpPoIid = 0, tmpPoIid_2 = 0;
	vector<int>::iterator iter;
	bool flag_replace = false;
	bool *poi_status = new bool[tmpPoI.poi_num];

	for (int i = 0; i < tmpPoI.poi_num; ++i)
		poi_status[i] = false;

	path.push_back(source);
	path.push_back(terminal);
	poi_status[source] = true; 
	poi_status[terminal] = true;
	satisfaction += SatisfactionPoI(source, tmpPoI) + SatisfactionPoI(terminal, tmpPoI);
	totalcost += (tmpPoI.PoIList[source].poi_cost
		+ tmpPoI.PoIList[terminal].poi_cost) / 2
		+ tmpPoI.CostLeadin(source, terminal);
	while (totalcost < Budget) {
		iter = path.end() - 2;
		ratio = 0;
		flag_replace = false;
		for (int i = 0; i < tmpPoI.poi_num; ++i)
			if (!poi_status[i]) {
				tmpcost = totalcost + tmpPoI.CostLeadin(*iter, i, *(iter + 1));
				if (tmpcost > Budget)
					continue;
				tmpsatis = satisfaction + SatisfactionPoI(i, tmpPoI);
				tmpratio = tmpsatis / tmpcost;
				if (tmpratio > ratio) {
					ratio = tmpratio;
					tmpPoIid = i;
					flag_replace = true;
				}
			}
		if (flag_replace) {
			satisfaction += SatisfactionPoI(tmpPoIid, tmpPoI);
			totalcost += tmpPoI.CostLeadin(*iter, tmpPoIid, *(iter + 1));
			path.insert(iter + 1, tmpPoIid);
			poi_status[tmpPoIid] = true;
		}
		else break;
	}
	
	//cout << "BR cost: " << totalcost << endl;
	//cout << "BR satisfaction: " << satisfaction << endl;

	for (iter = path.begin() + 1; iter != path.end() - 1; ++iter) {
		ratio = satisfaction / totalcost;
		flag_replace = false;
		for (int i = 0; i < tmpPoI.poi_num; ++i)
			if (!poi_status[i]) {
				tmpcost = totalcost - tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
					+ tmpPoI.CostLeadin(*(iter - 1), i, *(iter + 1));
				if (tmpcost > Budget)
					continue;
				tmpsatis = satisfaction + 
					SatisfactionPoI(i, tmpPoI) - SatisfactionPoI(*iter, tmpPoI);
				tmpratio = tmpsatis / tmpcost;
				if (tmpratio > ratio) {
					ratio = tmpratio;
					tmpPoIid = i;
					flag_replace = true;
				}
			}
		if (flag_replace) {
			satisfaction += SatisfactionPoI(tmpPoIid, tmpPoI) - SatisfactionPoI(*iter, tmpPoI);
			totalcost -= tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
				+ tmpPoI.CostLeadin(*(iter - 1), tmpPoIid, *(iter + 1));
			poi_status[*iter] = false;
			poi_status[tmpPoIid] = true;
			*iter = tmpPoIid;
		}
	}

	//cout << "BR+ cost: " << totalcost << endl;
	//cout << "BR+ satisfaction: " << satisfaction << endl;

	for (iter = path.begin() + 1; iter != path.end() - 1; ++iter) {
		ratio = satisfaction / totalcost;
		flag_replace = false;
		for (int i = 0; i < tmpPoI.poi_num; ++i) {
			if (!poi_status[i]) {
				tempcost = totalcost - tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
					+ tmpPoI.CostLeadin(*(iter - 1), i, *(iter + 1));
				if (tempcost > Budget)
					continue;
				for (int j = 0; j < tmpPoI.poi_num; ++j) {
					if (j == i) continue;
					if (!poi_status[j]) {
						tmpcost = tempcost + 
							tmpPoI.CostLeadin(*(path.end() - 2), j, *(path.end() - 1));
						if (tmpcost > Budget)
							continue;
						tmpsatis = satisfaction + SatisfactionPoI(i, tmpPoI) 
							- SatisfactionPoI(*iter, tmpPoI) + SatisfactionPoI(j, tmpPoI);
						tmpratio = tmpsatis / tmpcost;
						if (tmpratio > ratio) {
							ratio = tmpratio;
							tmpPoIid = i;
							tmpPoIid_2 = j;
							flag_replace = true;
						}
					}
				}
			}
		}
		if (flag_replace) {
			satisfaction += SatisfactionPoI(tmpPoIid, tmpPoI) - 
				SatisfactionPoI(*iter, tmpPoI) + SatisfactionPoI(tmpPoIid_2, tmpPoI);
			totalcost =  totalcost - tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
				+ tmpPoI.CostLeadin(*(iter - 1), tmpPoIid, *(iter + 1))
				+ tmpPoI.CostLeadin(*(path.end() - 2), tmpPoIid_2, *(path.end() - 1));
			poi_status[*iter] = false;
			poi_status[tmpPoIid] = true;
			poi_status[tmpPoIid_2] = true;
			*iter = tmpPoIid;
			path.insert(path.end() - 1, tmpPoIid_2);
			iter = path.begin();
			while (*iter != tmpPoIid) ++iter;
		}
	}

	//cout << "BR++ cost: " << totalcost << endl;
	//cout << "BR++ satisfaction: " << satisfaction << endl;

	delete[] poi_status;

	return satisfaction;
}

float TourGroup::BRBumaAlgorithm(int source, int terminal, PoI &tmpPoI)
{	
	path.erase(path.begin(), path.end());
	float satisfaction = 0, totalcost = 0, tmpsatis = 0, tmpcost = 0;
	float ratio = 0, tmpratio = 0, tempsatis = 0, tempcost = 0;
	int tmpPoIid = 0, tmpPoIid_2 = 0;
	vector<int> tmppath;
	vector<int>::iterator iter;
	bool flag_replace = false;
	bool* poi_status = new bool[tmpPoI.poi_num];
	
	for (int i = 0; i < tourist_num; ++i) {
		memset(poi_status, 0, tmpPoI.poi_num);
		tmppath.erase(tmppath.begin(), tmppath.end());
		tmppath.push_back(source);
		tmppath.push_back(terminal);
		poi_status[source] = true;
		poi_status[terminal] = true;
		tmpsatis = 0;
		tmpcost = 0;
		tempsatis = 0;
		tempcost = 0;
		tempsatis += SatisfactionTour(i, source, tmpPoI) + SatisfactionTour(i, terminal, tmpPoI);
		tempcost += (tmpPoI.PoIList[source].poi_cost
			+ tmpPoI.PoIList[terminal].poi_cost) / 2
			+ tmpPoI.CostLeadin(source, terminal);
		while (tempcost < Budget) {
			iter = tmppath.end() - 2;
			ratio = 0;
			flag_replace = false;
			for (int j = 0; j < tmpPoI.poi_num; ++j)
				if (!poi_status[j]) {
					tmpcost = tempcost + tmpPoI.CostLeadin(*iter, j, *(iter + 1));
					if (tmpcost > Budget) continue;
					tmpsatis = tempsatis + SatisfactionTour(i, j, tmpPoI);
					tmpratio = tmpsatis / tmpcost;
					if (tmpratio > ratio) {
						ratio = tmpratio;
						tmpPoIid = j;
						flag_replace = true;
					}
				}
			if (flag_replace) {
				tempsatis += SatisfactionTour(i, tmpPoIid, tmpPoI);
				tempcost += tmpPoI.CostLeadin(*iter, tmpPoIid, *(iter + 1));
				tmppath.insert(iter + 1, tmpPoIid);
				poi_status[tmpPoIid] = true;
			}
			else
				break;
		}
		tempsatis = SatisfactionCalculate(tmppath, tmpPoI);
		if (tempsatis > satisfaction) {
			path.assign(tmppath.begin(), tmppath.end());
			satisfaction = tempsatis;
			totalcost = tempcost;
		}
	}

	delete[] poi_status;

	return satisfaction;
}

float TourGroup::BRBumaBRPlusPlusAlgorithm(int source, int terminal, PoI& tmpPoI)
{
	path.erase(path.begin(), path.end());
	float satisfaction = 0, totalcost = 0, tmpsatis = 0, tmpcost = 0;
	float ratio = 0, tmpratio = 0, tempsatis = 0, tempcost = 0;
	int tmpPoIid = 0, tmpPoIid_2 = 0;
	vector<int> tmppath;
	vector<int>::iterator iter;
	bool flag_replace = false;
	bool* poi_status = new bool[tmpPoI.poi_num];

	for (int i = 0; i < tourist_num; ++i) {
		memset(poi_status, 0, tmpPoI.poi_num);
		tmppath.erase(tmppath.begin(), tmppath.end());
		tmppath.push_back(source);
		tmppath.push_back(terminal);
		poi_status[source] = true;
		poi_status[terminal] = true;
		tmpsatis = 0;
		tmpcost = 0;
		tempsatis = 0;
		tempcost = 0;
		tempsatis += SatisfactionTour(i, source, tmpPoI) + SatisfactionTour(i, terminal, tmpPoI);
		tempcost += (tmpPoI.PoIList[source].poi_cost
			+ tmpPoI.PoIList[terminal].poi_cost) / 2
			+ tmpPoI.CostLeadin(source, terminal);
		while (tempcost < Budget) {
			iter = tmppath.end() - 2;
			ratio = 0;
			flag_replace = false;
			for (int j = 0; j < tmpPoI.poi_num; ++j)
				if (!poi_status[j]) {
					tmpcost = tempcost + tmpPoI.CostLeadin(*iter, j, *(iter + 1));
					if (tmpcost > Budget) continue;
					tmpsatis = tempsatis + SatisfactionTour(i, j, tmpPoI);
					tmpratio = tmpsatis / tmpcost;
					if (tmpratio > ratio) {
						ratio = tmpratio;
						tmpPoIid = j;
						flag_replace = true;
					}
				}
			if (flag_replace) {
				tempsatis += SatisfactionTour(i, tmpPoIid, tmpPoI);
				tempcost += tmpPoI.CostLeadin(*iter, tmpPoIid, *(iter + 1));
				tmppath.insert(iter + 1, tmpPoIid);
				poi_status[tmpPoIid] = true;
			}
			else
				break;
		}
		tempsatis = SatisfactionCalculate(tmppath, tmpPoI);
		if (tempsatis > satisfaction) {
			path.assign(tmppath.begin(), tmppath.end());
			satisfaction = tempsatis;
			totalcost = tempcost;
		}
	}

	for (iter = path.begin() + 1; iter != path.end() - 1; ++iter) {
		ratio = satisfaction / totalcost;
		flag_replace = false;
		for (int i = 0; i < tmpPoI.poi_num; ++i) {
			if (!poi_status[i]) {
				tempcost = totalcost - tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
					+ tmpPoI.CostLeadin(*(iter - 1), i, *(iter + 1));
				if (tempcost > Budget)
					continue;
				for (int j = 0; j < tmpPoI.poi_num; ++j) {
					if (j == i) continue;
					if (!poi_status[j]) {
						tmpcost = tempcost +
							tmpPoI.CostLeadin(*(path.end() - 2), j, *(path.end() - 1));
						if (tmpcost > Budget)
							continue;
						tmpsatis = satisfaction + SatisfactionPoI(i, tmpPoI)
							- SatisfactionPoI(*iter, tmpPoI) + SatisfactionPoI(j, tmpPoI);
						tmpratio = tmpsatis / tmpcost;
						if (tmpratio > ratio) {
							ratio = tmpratio;
							tmpPoIid = i;
							tmpPoIid_2 = j;
							flag_replace = true;
						}
					}
				}
			}
		}
		if (flag_replace) {
			satisfaction += SatisfactionPoI(tmpPoIid, tmpPoI) -
				SatisfactionPoI(*iter, tmpPoI) + SatisfactionPoI(tmpPoIid_2, tmpPoI);
			totalcost = totalcost - tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
				+ tmpPoI.CostLeadin(*(iter - 1), tmpPoIid, *(iter + 1))
				+ tmpPoI.CostLeadin(*(path.end() - 2), tmpPoIid_2, *(path.end() - 1));
			poi_status[*iter] = false;
			poi_status[tmpPoIid] = true;
			poi_status[tmpPoIid_2] = true;
			*iter = tmpPoIid;
			path.insert(path.end() - 1, tmpPoIid_2);
			iter = path.begin();
			while (*iter != tmpPoIid) ++iter;
		}
	}

	delete[] poi_status;

	return satisfaction;
}

float TourGroup::BRplusBumaAlgorithm(int source, int terminal, PoI& tmpPoI)
{
	path.erase(path.begin(), path.end());
	float satisfaction = 0, totalcost = 0, tmpsatis = 0, tmpcost = 0, ttmpcost = 0;
	float ratio = 0, tmpratio = 0, tempsatis = 0, tempcost = 0;
	int tmpPoIid = 0, tmpPoIid_2 = 0;
	vector<int> tmppath;
	vector<int>::iterator iter;
	bool flag_replace = false;
	bool* poi_status = new bool[tmpPoI.poi_num];

	for (int i = 0; i < tourist_num; ++i) {
		memset(poi_status, 0, tmpPoI.poi_num);
		tmppath.erase(tmppath.begin(), tmppath.end());
		tmppath.push_back(source);
		tmppath.push_back(terminal);
		poi_status[source] = true;
		poi_status[terminal] = true;
		tmpsatis = 0;
		tmpcost = 0;
		tempsatis = 0;
		tempcost = 0;
		tempsatis += SatisfactionTour(i, source, tmpPoI) + SatisfactionTour(i, terminal, tmpPoI);
		tempcost += (tmpPoI.PoIList[source].poi_cost
			+ tmpPoI.PoIList[terminal].poi_cost) / 2
			+ tmpPoI.CostLeadin(source, terminal);
		while (tempcost < Budget) {
			iter = tmppath.end() - 2;
			ratio = 0;
			flag_replace = false;
			for (int j = 0; j < tmpPoI.poi_num; ++j)
				if (!poi_status[j]) {
					tmpcost = tempcost + tmpPoI.CostLeadin(*iter, j, *(iter + 1));
					if (tmpcost > Budget) continue;
					tmpsatis = tempsatis + SatisfactionTour(i, j, tmpPoI);
					tmpratio = tmpsatis / tmpcost;
					if (tmpratio > ratio) {
						ratio = tmpratio;
						tmpPoIid = j;
						flag_replace = true;
					}
				}
			if (flag_replace) {
				tempsatis += SatisfactionTour(i, tmpPoIid, tmpPoI);
				tempcost += tmpPoI.CostLeadin(*iter, tmpPoIid, *(iter + 1));
				tmppath.insert(iter + 1, tmpPoIid);
				poi_status[tmpPoIid] = true;
			}
			else
				break;
		}

		for (iter = tmppath.begin() + 1; iter != tmppath.end() - 1; ++iter) {
			ratio = tempsatis / tempcost;
			flag_replace = false;
			for (int j = 0; j < tmpPoI.poi_num; ++j)
				if (!poi_status[j]) {
					tmpcost = tempcost - tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
						+ tmpPoI.CostLeadin(*(iter - 1), j, *(iter + 1));
					if (tmpcost > Budget)
						continue;
					tmpsatis = tempsatis +
						SatisfactionTour(i, j, tmpPoI) - SatisfactionTour(i, *iter, tmpPoI);
					tmpratio = tmpsatis / tmpcost;
					if (tmpratio > ratio) {
						ratio = tmpratio;
						tmpPoIid = j;
						flag_replace = true;
					}
				}
			if (flag_replace) {
				tempsatis += SatisfactionTour(i, tmpPoIid, tmpPoI) - SatisfactionTour(i, *iter, tmpPoI);
				tempcost -= tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
					+ tmpPoI.CostLeadin(*(iter - 1), tmpPoIid, *(iter + 1));
				poi_status[*iter] = false;
				poi_status[tmpPoIid] = true;
				*iter = tmpPoIid;
			}
		}

		tempsatis = SatisfactionCalculate(tmppath, tmpPoI);
		if (tempsatis > satisfaction) {
			path.assign(tmppath.begin(), tmppath.end());
			satisfaction = tempsatis;
			totalcost = tempcost;
		}
	}
	delete[] poi_status;

	return satisfaction;
}

float TourGroup::BRplusplusBumaAlgorithm(int source, int terminal, PoI& tmpPoI)
{
	path.erase(path.begin(), path.end());
	float satisfaction = 0, totalcost = 0, tmpsatis = 0, tmpcost = 0, ttmpcost = 0;
	float ratio = 0, tmpratio = 0, tempsatis = 0, tempcost = 0;
	int tmpPoIid = 0, tmpPoIid_2 = 0;
	vector<int> tmppath;
	vector<int>::iterator iter;
	bool flag_replace = false;
	bool* poi_status = new bool[tmpPoI.poi_num];

	for (int i = 0; i < tourist_num; ++i) {
		memset(poi_status, 0, tmpPoI.poi_num);
		tmppath.erase(tmppath.begin(), tmppath.end());
		tmppath.push_back(source);
		tmppath.push_back(terminal);
		poi_status[source] = true;
		poi_status[terminal] = true;
		tmpsatis = 0;
		tmpcost = 0;
		tempsatis = 0;
		tempcost = 0;
		tempsatis += SatisfactionTour(i, source, tmpPoI) + SatisfactionTour(i, terminal, tmpPoI);
		tempcost += (tmpPoI.PoIList[source].poi_cost
			+ tmpPoI.PoIList[terminal].poi_cost) / 2
			+ tmpPoI.CostLeadin(source, terminal);
		while (tempcost < Budget) {
			iter = tmppath.end() - 2;
			ratio = 0;
			flag_replace = false;
			for (int j = 0; j < tmpPoI.poi_num; ++j)
				if (!poi_status[j]) {
					tmpcost = tempcost + tmpPoI.CostLeadin(*iter, j, *(iter + 1));
					if (tmpcost > Budget) continue;
					tmpsatis = tempsatis + SatisfactionTour(i, j, tmpPoI);
					tmpratio = tmpsatis / tmpcost;
					if (tmpratio > ratio) {
						ratio = tmpratio;
						tmpPoIid = j;
						flag_replace = true;
					}
				}
			if (flag_replace) {
				tempsatis += SatisfactionTour(i, tmpPoIid, tmpPoI);
				tempcost += tmpPoI.CostLeadin(*iter, tmpPoIid, *(iter + 1));
				tmppath.insert(iter + 1, tmpPoIid);
				poi_status[tmpPoIid] = true;
			}
			else
				break;
		}

		for (iter = tmppath.begin() + 1; iter != tmppath.end() - 1; ++iter) {
			ratio = tempsatis / tempcost;
			flag_replace = false;
			for (int j = 0; j < tmpPoI.poi_num; ++j)
				if (!poi_status[j]) {
					tmpcost = tempcost - tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
						+ tmpPoI.CostLeadin(*(iter - 1), j, *(iter + 1));
					if (tmpcost > Budget)
						continue;
					tmpsatis = tempsatis +
						SatisfactionTour(i, j, tmpPoI) - SatisfactionTour(i, *iter, tmpPoI);
					tmpratio = tmpsatis / tmpcost;
					if (tmpratio > ratio) {
						ratio = tmpratio;
						tmpPoIid = j;
						flag_replace = true;
					}
				}
			if (flag_replace) {
				tempsatis += SatisfactionTour(i, tmpPoIid, tmpPoI) - SatisfactionTour(i, *iter, tmpPoI);
				tempcost -= tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
					+ tmpPoI.CostLeadin(*(iter - 1), tmpPoIid, *(iter + 1));
				poi_status[*iter] = false;
				poi_status[tmpPoIid] = true;
				*iter = tmpPoIid;
			}
		}

		for (iter = tmppath.begin() + 1; iter != tmppath.end() - 1; ++iter) {
			ratio = tempsatis / tempcost;
			flag_replace = false;
			for (int j = 0; j < tmpPoI.poi_num; ++j) {
				if (!poi_status[i]) {
					ttmpcost = tempcost - tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
						+ tmpPoI.CostLeadin(*(iter - 1), j, *(iter + 1));
					if (ttmpcost > Budget)
						continue;
					for (int k = 0; k < tmpPoI.poi_num; ++k) {
						if (k == j) continue;
						if (!poi_status[j]) {
							tmpcost = ttmpcost +
								tmpPoI.CostLeadin(*(path.end() - 2), k, *(path.end() - 1));
							if (tmpcost > Budget)
								continue;
							tmpsatis = tempsatis + SatisfactionTour(i, j, tmpPoI)
								- SatisfactionTour(i, *iter, tmpPoI) + SatisfactionTour(i, k, tmpPoI);
							tmpratio = tmpsatis / tmpcost;
							if (tmpratio > ratio) {
								ratio = tmpratio;
								tmpPoIid = i;
								tmpPoIid_2 = j;
								flag_replace = true;
							}
						}
					}
				}
			}
			if (flag_replace) {
				tempsatis += SatisfactionTour(i, tmpPoIid, tmpPoI) -
					SatisfactionTour(i, *iter, tmpPoI) + SatisfactionTour(i, tmpPoIid_2, tmpPoI);
				tempcost = tempcost - tmpPoI.CostLeadin(*(iter - 1), *iter, *(iter + 1))
					+ tmpPoI.CostLeadin(*(iter - 1), tmpPoIid, *(iter + 1))
					+ tmpPoI.CostLeadin(*(path.end() - 2), tmpPoIid_2, *(path.end() - 1));
				poi_status[*iter] = false;
				poi_status[tmpPoIid] = true;
				poi_status[tmpPoIid_2] = true;
				*iter = tmpPoIid;
				tmppath.insert(tmppath.end() - 1, tmpPoIid_2);
				iter = tmppath.begin();
				while (*iter != tmpPoIid) ++iter;
			}
		}
		tempsatis = SatisfactionCalculate(tmppath, tmpPoI);
		if (tempsatis > satisfaction) {
			path.assign(tmppath.begin(), tmppath.end());
			satisfaction = tempsatis;
			totalcost = tempcost;
		}
	}
	delete[] poi_status;

	return satisfaction;
}


float TourGroup::EsAlgorithm(int source, int terminal, PoI& tmpPoI)
{	
	vector<int> tmpstack;
	float totalcost = 0, tempcost = 0, tmpcost = 0, satisfaction = 0;
	int PoInum = tmpPoI.poi_num, i, j;
	float tmpvari = 0, variance = 10000;
	int tmpPoIid, lastPoIid;
	float* tour_satis = new float[tourist_num];
	float* tmp_satis = new float[tourist_num];
	bool* poi_instack = new bool[PoInum];
	bool** poi_visit = new bool*[PoInum];
	
	for (i = 0; i < PoInum; ++i)
		poi_visit[i] = new bool[PoInum];

	for (i = 0; i < tourist_num; ++i) {
		tmp_satis[i] = 0;
		tour_satis[i] = 0;
	}

	for (i = 0; i < PoInum; ++i)
		for (j = 0; j < PoInum; ++j)
			poi_visit[i][j] = 0;

	memset(poi_instack, 0, PoInum);

	tmpstack.push_back(source);
	tempcost = tmpPoI.PoIList[source].poi_cost;
	poi_instack[source] = true;

	cout << "tempcost: " << tempcost << endl;
	while (!tmpstack.empty()) {
		tmpPoIid = -1;
		lastPoIid = tmpstack.back();
		for (i = 0; i < PoInum; ++i) {
			if (poi_instack[i]) continue;
			if (poi_visit[lastPoIid][i]) continue;
			tmpcost = tempcost + tmpPoI.CostLeadinS(lastPoIid, i);
			if (tmpcost > Budget) continue;
			tmpPoIid = i;
			break;
		}
		if (tmpPoIid != -1) {
			if (tmpPoIid == terminal) {
				tmpcost = tempcost + tmpPoI.CostLeadinS(lastPoIid, terminal);
				tmpstack.push_back(terminal);
				for (i = 0; i < tourist_num; ++i)
					for (j = 0; j < tmpstack.size(); ++j)
						tmp_satis[i] += SatisfactionTour(i, tmpstack[j], tmpPoI);
				tmpvari = VarianceCalculate(tmp_satis);
				if (tmpvari < variance) {
					totalcost = tmpcost;
					satisfaction = SatisfactionSum(tmp_satis);
					variance = tmpvari;
					path.assign(tmpstack.begin(), tmpstack.end());
					cout << "EStmp: " << satisfaction << endl;
					cout << "EsCost: " << totalcost << endl;
				}
				tmpstack.pop_back();
				poi_visit[lastPoIid][terminal] = true;
			}
			else {
				tmpstack.push_back(tmpPoIid);
				poi_instack[tmpPoIid] = true;
				tempcost += tmpPoI.CostLeadinS(lastPoIid, tmpPoIid);
				lastPoIid = tmpPoIid;
				cout << lastPoIid << '\t' << tempcost << endl;
			}
		}
		else {
			tmpPoIid = lastPoIid;
			tmpstack.pop_back();
			memset(poi_visit[tmpPoIid], 0, PoInum);
			lastPoIid = tmpstack.back();
			poi_instack[tmpPoIid] = false;
			poi_visit[lastPoIid][tmpPoIid] = true;
			tempcost -= tmpPoI.CostLeadinS(lastPoIid, tmpPoIid);
		}
	}
		
	delete[] tour_satis;
	delete[] tmp_satis;
	delete[] poi_instack;
	for (i = 0; i < PoInum; ++i)
		delete[] poi_visit[i];
	delete[] poi_visit;

	return satisfaction;
}



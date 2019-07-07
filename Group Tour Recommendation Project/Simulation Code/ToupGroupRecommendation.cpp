#include <iostream>
#include <fstream>
#include <sstream>
#include "Graph.h"
#include "TourGroup.h"
#include <vector>
#include <chrono>
using namespace chrono;

int main()
{	
	int i = 0;
	TourGroup Tour_Exm(40, 6);
	PoI PoI_Exm(28);
	ifstream fin;
	fin.open("Vien//interest.csv", ios::in); 
	string line;

	while (getline(fin, line))
	{	
		if (i == 0) {
			++i;
			continue;
		}
		istringstream sin(line); 
		vector<string> fields;
		string field;
		while (getline(sin, field, ','))
			fields.push_back(field);
		Tour_Exm.TourAppend(fields);
		++i;
		if (i > 40) break;
	}
	cout << "Initalize TourGroup Complete!" << endl;
	fin.close();

	i = 0;
	fin.open("Vien//Vienpop.csv", ios::in);
	while (getline(fin, line))
	{
		if (i == 0) {
			++i;
			continue;
		}
		istringstream sin(line);
		vector<string> fields;
		string field;
		while (getline(sin, field, ','))
			fields.push_back(field);
		PoI_Exm.PoIAppend(fields);
	}
	cout << "Initialize PoIList Complete!" << endl;
	fin.close();

	i = 0;
	fin.open("Vien//time.csv", ios::in);
	while (getline(fin, line))
	{
		if (i == 0) {
			++i;
			continue;
		}
		istringstream sin(line);
		vector<string> fields;
		string field;
		while (getline(sin, field, ','))
			fields.push_back(field);
		PoI_Exm.DistanceAppend(fields);
	}
	cout << "Initialize cost Complete!" << endl;
	fin.close();

	for (int j = 1; j < 64; j = j * 2) {
		SetNumberTour(j, Tour_Exm);
		//SetNumberPoI(20, PoI_Exm);
		auto start = system_clock::now();
		cout << Tour_Exm.GreedyAlgorithm(0, 11, PoI_Exm) << endl;
		auto end = system_clock::now();
		auto duration = duration_cast<microseconds>(end - start);
		cout << "Greedy\n" << double(duration.count()) * microseconds::period::num / microseconds::period::den << "Ãë" << endl;

		start = system_clock::now();
		cout << Tour_Exm.BRplusplusBumaAlgorithm(0, 11, PoI_Exm) << endl;
		end = system_clock::now();
		duration = duration_cast<microseconds>(end - start);
		cout << "Power\n" << double(duration.count()) * microseconds::period::num / microseconds::period::den << "Ãë" << endl;

		start = system_clock::now();
		cout << Tour_Exm.EsAlgorithm(0, 11, PoI_Exm) << endl;
		end = system_clock::now();
		duration = duration_cast<microseconds>(end - start);
		cout << "ENumeration\n" << double(duration.count()) * microseconds::period::num / microseconds::period::den << "Ãë" << endl;
		cout << endl << endl;
	}

	return 0;
}

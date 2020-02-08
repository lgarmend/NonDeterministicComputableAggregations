// lgarmend@fdi.ucm.es
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <random> // For random functions
#include <ctime>
#include <windows.h> // For sleep
using namespace std;
//Population size
const int PopulationSize = 1000;
// Array (List) of population data
typedef double PopulationData[PopulationSize];
// Times to execute an aggregation to a Population Data
const int AggExecutionTimes = 100;
// Array (List) to store the ExecutionTimes outputs from computing an aggregation to a Population Data 
typedef double ExecutionsOutputData[AggExecutionTimes];
// Probability for samples to be used in aggregations
const double sampleProb = 0.1;

//Functions
void generateUniformPopulationData(PopulationData p);
void generateNormalPopulationData(PopulationData p);
void savePopulationData(PopulationData d);

//Aggregation: List<[0, 1]> -> [0, 1]
double agg(PopulationData d);
double aggSampleAVG(PopulationData d, double prob);
double aggSampleMax(PopulationData d, double prob);
double aggSampleMin(PopulationData d, double prob);
double generateNormal(double mean, double stdDev);
double aggNoisyAVG(PopulationData d, double mean, double stdDev);
double generateUniform(double min, double max);
double aggPureRandom(); 
double aggBoundedPureRandom(PopulationData d, double min, double max);

// Computing ExecutionTimes an Aggregation to a Population Data
// Input: Population Data  Output: ExecutionsOutputData
void generateAggExecutions(PopulationData p, ExecutionsOutputData o);
void generateSampleAVG(PopulationData p, ExecutionsOutputData o, double prob);
void generateSampleMax(PopulationData p, ExecutionsOutputData o, double prob); 
void generateSampleMin(PopulationData p, ExecutionsOutputData o, double prob);
void generateNoisyAVG(PopulationData p, ExecutionsOutputData o, double mean, double stdDev);
void generatePureRandom(PopulationData p, ExecutionsOutputData o);
void generateBoundedPureRandom(PopulationData p, ExecutionsOutputData o);

void saveAllExecutionsOutputData(PopulationData p, ExecutionsOutputData sampleAVG, 
	ExecutionsOutputData sampleMax, ExecutionsOutputData sampleMin, ExecutionsOutputData noisyAVG, 
	ExecutionsOutputData pureRandom, ExecutionsOutputData boundedPureRandom);
void bye();

// Main Program
int main()
{   //Storage Data Declaration
	PopulationData pd;
	ExecutionsOutputData sampleAVG;
	ExecutionsOutputData sampleMax;
	ExecutionsOutputData sampleMin;
	ExecutionsOutputData noisyAVG;
	ExecutionsOutputData pureRandom;
	ExecutionsOutputData boundedPureRandom;

	//Population Data generation
	generateNormalPopulationData(pd);
	//Save Population Data in file GeneratedOutputData.txt
	savePopulationData(pd);

	//Execute ExecutionTimes times the aggregations to the population data
	cout << "Number of Aggregation Executions: " << AggExecutionTimes << endl;
	cout << "Executing Agregations...... "  << endl;
	generateSampleAVG(pd, sampleAVG, sampleProb);
	generateSampleMax(pd, sampleMax, sampleProb);
	generateSampleMin(pd, sampleMin, sampleProb);
	generateNoisyAVG(pd, noisyAVG, 0, 1);
	generatePureRandom(pd, pureRandom);
	generateBoundedPureRandom(pd, boundedPureRandom);

	//Save the executions results in file ExecutionOutputData.txt
	saveAllExecutionsOutputData(pd,sampleAVG,sampleMax,sampleMin,noisyAVG,pureRandom, boundedPureRandom);
	bye();
	int n; cin >> n;
	return 0;
}

void generateUniformPopulationData(PopulationData p) {
	for (int i = 0; i < PopulationSize; i++) {
		 p[i] = generateUniform(0,1);
	}
}
void generateNormalPopulationData(PopulationData p) {
	for (int i = 0; i < PopulationSize; i++) {
		p[i] = -1;
		while (p[i]<0 || p[i]>1)
			p[i] = generateNormal(0.5, 0.1);
	}
}


//Aggregations: List<[0, 1]> -> [0, 1]
double agg(PopulationData d) {
	return 0;
}

double aggSampleAVG(PopulationData d, double prob) {
	int n = 0;
	double sum = 0;
	for (int i = 0; i < PopulationSize; i++) 
		if (generateUniform(0,1) < prob) {
			sum += d[i];
			n++;
		}
	return sum / n;
}
double aggSampleMax(PopulationData d, double prob) {
	double max = 0;
	for (int i = 0; i < PopulationSize; i++) 
		if (generateUniform(0, 1) < prob)
			if(max< d[i]) max= d[i];	
	return max;
}

double aggSampleMin(PopulationData d, double prob) {
	double min = 1;
	for (int i = 0; i < PopulationSize; i++) 
		if (generateUniform(0, 1) < prob)
			if (min> d[i]) min = d[i];
	return min;
}


double generateNormal(double mean, double stdDev) {
	std::random_device rd;
	std::mt19937 generator(rd());
	std::normal_distribution<double> n(mean, stdDev);
	return  n(generator);
}

double generateUniform(double min, double max) {
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<double> u(min, max);
	return  u(generator);
}


double aggNoisyAVG(PopulationData d, double mean, double stdDev) {
	double sum = 0;
	for (int i = 0; i < PopulationSize; i++)
			sum += d[i]+generateNormal(mean,stdDev);
	return sum / PopulationSize;
}

double aggPureRandom() {
	return generateUniform(0, 1);
}
double aggBoundedPureRandom(PopulationData d, double min, double max) {
	return generateUniform(min, max);
}

//General Execution aggregation agg on Population Data AggExecutionTimes Times
// Input: Population Data  Output: ExecutionsOutputData
void generateAggExecutions(PopulationData d, ExecutionsOutputData o) {
	for (int i = 0; i < AggExecutionTimes; i++) {
		o[i] = agg(d);
	}
}

void generateSampleAVG(PopulationData d, ExecutionsOutputData o, double prob) {
	for (int i = 0; i < AggExecutionTimes; i++) {
		o[i] = aggSampleAVG(d, prob);
	}
}
void generateSampleMax(PopulationData d, ExecutionsOutputData o, double prob) {
	for (int i = 0; i < AggExecutionTimes; i++) {
		o[i] = aggSampleMax(d, prob);
	}
}

void generateSampleMin(PopulationData d, ExecutionsOutputData o, double prob) {
	for (int i = 0; i < AggExecutionTimes; i++) {
		o[i] = aggSampleMin(d, prob);
	}
}
void generateNoisyAVG(PopulationData d, ExecutionsOutputData o, double mean,double stdDev) {
	for (int i = 0; i < AggExecutionTimes; i++) {
		o[i] =aggNoisyAVG(d, mean, stdDev);
	}
}
void generatePureRandom(PopulationData d, ExecutionsOutputData o) {
	for (int i = 0; i < AggExecutionTimes; i++) {
		o[i] = generateUniform(0, 1);
	}
}
void generateBoundedPureRandom(PopulationData d, ExecutionsOutputData o) {
	double min = aggSampleMin(d, 1);
	double max = aggSampleMax(d, 1);
	for (int i = 0; i < AggExecutionTimes; i++) {
		o[i] = generateUniform(min,max);
	}
}

void bye()
{
	cout << endl << "*** Thanks for testing Non Computable Aggregations ***" << endl;
	Sleep(2000);
}

//Saving data to files

void savePopulationData(PopulationData d) {
	ofstream file;
	file.open("GeneratedPopulationData.txt");
	if (file.is_open()) {
		cout << "Population Size: " <<PopulationSize<< endl;
		file << left << "PopulationData" << endl;
		file << fixed << setprecision(8);
		for (int i = 0; i < PopulationSize; i++) {
			file <<left<<setw(10)<< d[i] <<  endl;
		}
		file.close();
		cout << "File GeneratedPopulationData Saved" << endl;
	}
	else cout << "Can not open GeneratedPopulationData file. "<<endl;
}


void saveAllExecutionsOutputData(PopulationData p, ExecutionsOutputData sampleAVG, ExecutionsOutputData sampleMax, 
	ExecutionsOutputData sampleMin, ExecutionsOutputData noisyAVG, ExecutionsOutputData pureRandom, ExecutionsOutputData boundedPureRandom){
	ofstream file;
	file.open("GeneratedExecutionsLists.txt");
	if (file.is_open()) {
		file <<left<< setw(10) << "SampleAVG" << setw(10) << "SampleMax" << setw(10) << "SampleMin" <<
			setw(10) << "NoisyAVG" << setw(11) << "PureRandom" << setw(20) << "BoundedRandom" << endl;
		file << fixed << setprecision(6);
		for (int i = 0; i < AggExecutionTimes; i++) {
			file << setw(10) << sampleAVG[i] << setw(10) << sampleMax[i]<< setw(10) << sampleMin[i] <<
				setw(10) << noisyAVG[i] << setw(11) << pureRandom[i] << setw(20) << boundedPureRandom[i] << endl;
		}
		file.close();
		cout << "File GeneratedExecutionsLists Saved" << endl;
	}
	else cout << "Can not open GeneratedExecutionsLists file. "<<endl;
}

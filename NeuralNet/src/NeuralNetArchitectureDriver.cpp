//============================================================================
// Name        : NeuralNet.cpp
// Author      : Aaron Wilson
// Version     :
// Copyright   :
//


//Runing a few times with 2 layers and 5 nodes using 10 containers I've noticed that total nodes across all layers hovers around an
//average of 7 for the best architecture. The best attempt was 6 nodes:
//2 layers with Nodes 4:2 and an error of 0.000247646.
//
//There was one slightly better than that with 7 nodes, but it was less than 0.000005 better.
//============================================================================

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <sstream>
#include <time.h>
#include "NeuralNet.h"
using namespace std;
struct errorCalc{
	int numLayers;
	int *layerNodes;
	double error;
};
double learningTime, learningRate;
int maxHiddenLayers;
int maxNodes;
int K;
string data;
void readFile(char[]);
void printGlobal();
long getRounds();
void sort();
void mergeSort(errorCalc*,int, int);
void merge(errorCalc*, int, int, int);




int main(int argc, char *argv[]) {
	//time(NULL);
	//cout << argv[1] << endl;
	readFile(argv[1]);
	printGlobal();
	;
	//std::cout << "testing\n";
	dataSet *set = NeuralNet::readFile(data.c_str());
	set = NeuralNet::randomizeData(set);
	dataSet *splitSet = NeuralNet::splitData(set, K);
	std::cout << std::endl << std::endl;
	//dataSet *comb = NeuralNet::joinData(splitSet, K, 1);
	//NeuralNet::printData(comb);
	//NeuralNet *net;
	//net = new NeuralNet((comb), 2, layers, learningRate, learningTime);
	//net->printSolutions();
	//net->testPoints(splitSet + 1);
	//std::cout << net->getError(splitSet+1) << std::endl;
	//double errors[2];
	int numRounds = getRounds();
	errorCalc errors[numRounds];
	int numLayers = 1;
	int layerDiv = maxNodes;
	//int layerOffset = 1;
	cout << "Total Rounds: " <<getRounds() << endl;
	for (int round = 0; round < numRounds; round++){
		cout << "round " << round + 1 << endl;
		int *layers = new int[maxHiddenLayers];
		if (round / layerDiv == 1){
			numLayers += 1;
			layerDiv += (pow(maxNodes, numLayers));
		}

		//layerOffset = numLayers;
		layers[0] = (round) % maxNodes + 1;
		for (int x = 1; x < numLayers; x++){
			int section = round;
			for (int y = x; y > 0; y--){
				section -= pow(maxNodes, y);
			}
			cout << section << " / " << pow(maxNodes, x) << endl;
			layers[x] = section / (int)pow(maxNodes, x) % maxNodes + 1;

		}
		errors[round].layerNodes = layers;
		errors[round].numLayers = numLayers;
		double error = 0;
		cout << "Running Layers: " << numLayers << endl;
		for (int x = 0; x < numLayers; x++){
			cout << "Layer " << x << " nodes: " << layers[x] << endl;
		}
		//cout << "Trial data: ";
		for (int x = 0; x < K; x++){
			cout << "Trial data: " << x << endl;
			NeuralNet net(NeuralNet::joinData(splitSet, K, x), numLayers, layers, learningRate, learningTime);
			error = net.getError(splitSet + x);
		}
		//cout << endl;
		errors[round].error = error / K;
		cout << "Average error: " << errors[round].error << endl << endl;
	}

	mergeSort(errors, 0, numRounds - 1);

	cout << "Top: " << endl;
	int top = 5;
	if (top > numRounds) top = numRounds;
	for (int x = 0; x < top; x++){
		cout << "layers: " << errors[x].numLayers << endl;
		for (int layer = 0; layer < errors[x].numLayers; layer++){
			cout << "node " << layer << ": " << errors[x].layerNodes[layer] << endl;
		}
		cout << "Error: " << errors[x].error << endl << endl;
	}


	NeuralNet::deleteData(set);
	//delete net;


	std::cout << "Fin\n";

	return 0;
}
void sort(){

		//mergeSort(population, 0, populationSize - 1);

	}

void mergeSort(errorCalc* array, int low, int high){

		int mid = (low + high)/2;
		if (low < high){
			mergeSort(array, low, mid);
			mergeSort(array, mid + 1, high);
			merge(array, low, mid, high);
		}
	}

void merge(errorCalc* array, int low, int mid, int high){
	errorCalc left[mid - low + 2];
	for (int x = 0; x <= (mid - low); x++)
	{
		left[x] = array[x+low];


	}
	left[mid - low + 1].error = 1999999999999999999999.0;
	errorCalc right[high - mid + 1];
	for (int x = 0; x < (high - mid + 1); x++)
	{
		right[x] = array[x+mid+1];
	}
	right[high - mid].error = 19999999999999999.0;
	int i = 0;
	int j = 0;
	for (int x = low; x <= high; x++ ){
		if (abs(left[i].error) < abs(right[j].error)){
			array[x] = left[i];
			i++;
		}
		else{
			array[x] = right[j];
			j++;
		}
		//cout << array[x].numLayers<< endl;
		//cout << array[x].layerNodes[0]<< endl;
		//cout << array[x].error<<endl <<endl;

	}


}

long getRounds(){
	long rounds = 0;
	for (int x = 1; x <= maxHiddenLayers; x++){
		rounds += pow(maxNodes, x);
	}
	return rounds;

}


void readFile(char fileName[]){

	ifstream input;
	input.open(fileName);

	if(!input.is_open()){
		cout << "Problem opening file: exiting\n";

	}
	else{

		getline(input, data);
		input >> K;
		input >> maxHiddenLayers;
		input >> maxNodes;
		input >> learningRate;
		input >> learningTime;

		input.close();
	}

}

void printGlobal(){
	cout << "Data: " << data << endl
		 << "K: " << K << endl
		 << "Max Hidden Layers: " << maxHiddenLayers << endl
		 << "Max Hidden Nodes:" << maxNodes << endl;




	cout << "Learning Rate: " << learningRate <<endl
		 << "Time: " << learningTime <<endl;
}

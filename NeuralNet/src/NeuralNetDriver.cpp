//============================================================================
// Name        : NeuralNet.cpp
// Author      : Aaron Wilson
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <time.h>
#include "NeuralNet.h"
using namespace std;

double learningTime, learningRate;
int numHiddenLayers;
int* layers;
string data;
void readFile(char[]);
void printGlobal();


int main(int argc, char *argv[]) {
	//time(NULL);
	//cout << argv[1] << endl;
	readFile(argv[1]);
	printGlobal();

	NeuralNet *net;
	net = new NeuralNet(NeuralNet::readFile(data.c_str()), numHiddenLayers, layers, learningRate, learningTime);
	net->printSolutions();

	delete net;


	std::cout << "Fin\n";

	return 0;
}

void readFile(char fileName[]){

	ifstream input;
	input.open(fileName);

	if(!input.is_open()){
		cout << "Problem opening file: exiting\n";

	}
	else{

		getline(input, data);

		input >> numHiddenLayers;
		layers = new int[numHiddenLayers];

		for (int x = 0; x < numHiddenLayers; x++){
			input >> layers[x];
		}

		input >> learningRate;
		input >> learningTime;

		input.close();
	}

}

void printGlobal(){
	cout << "Data: " << data <<endl
		 << "Hidden Layers: " << numHiddenLayers << endl;

		for (int x = 0; x < numHiddenLayers; x++){
			cout << "layer " << x << ": " << layers[x] << endl;
		}

	cout << "Learning Rate: " << learningRate <<endl
		 << "Time: " << learningTime <<endl;
}

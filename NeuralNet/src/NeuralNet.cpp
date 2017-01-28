/*
 * NeuralNet.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: Aaron
 */

#include "NeuralNet.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <sstream>
#include <time.h>
#include <cmath>
//#include <random>


NeuralNet::NeuralNet(dataSet *workingSet, int numLayers, int* layers, double learningRate, double learningTime) {
	srand(time(NULL));
	//srand(11228374);
	this->data = workingSet;
	this->numLayers = numLayers;
	this->layers = layers;
	this->learningRate = learningRate;
	this->learningTime = learningTime;

	this->normilization = this->data->scale;
	this->inputNodes = this->data->inputs;
	this->outputNodes = this->data->outputs;
	//std::cout<< "test\n";
	this->neurons = new neuron*[numLayers+1];
	initLayers(this->neurons);

	initBiasWeight(this->neurons);
	//printNode(neurons[1][0], 0);
	//printNetwork();
	//std::cout << randomWeight();
	//solve(data->dataPoints, 1, 1);
	//printWeights();
	learn();
	//feedForward(data->dataPoints, 1,1);
	//printSolutions();
}
NeuralNet::~NeuralNet() {
	deleteData(data);
	deleteNet();


	// TODO Auto-generated destructor stub
}

void NeuralNet::deleteNet(){
	for (int x = 0; x < (numLayers); x++){
			for (int y = 0; y < layers[x]; y++){
				if (x == 0){
					delete neurons[x][y].biasWeight;
				}
				else{
					delete neurons[x][y].biasWeight;
				}//*/
			}
		}
		for (int y = 0; y < outputNodes; y++){
			//std::cout << "testing out: " << layers[numLayers - 1] << std::endl;
			delete neurons[numLayers][y].biasWeight;
		}

	delete neurons;
}
void NeuralNet::deleteData(dataSet* data){
	int entry = data->inputs + data->outputs;
	int total = data->points;
	total = total * entry;

	delete data->dataPoints;

	delete data;
}
void NeuralNet::printNode(neuron &node, int layer){
	std::cout << "node: \n";
	std::cout << "weights: ";
	int rounds;
	if (layer == 0) rounds = inputNodes;
	else rounds = layers[layer - 1];
	for (int x = 0; x < rounds; x++){
		std::cout << "(" << x << ")"<<node.biasWeight[x] << " ";
	}
	std::cout << std::endl;

}
void NeuralNet::printNetwork(){
	std::cout << "Network: " << numLayers << std::endl;
	for (int x = 0; x < (numLayers); x++){
		std::cout << "(" << layers[x] <<") ";
		for (int y = 0; y < layers[x]; y++){
			std::cout << neurons[x][y].out << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "output: ";
	for (int x = 0; x < outputNodes; x++){
		std::cout << neurons[numLayers][x].out << " ";
	}
	std::cout << std::endl;
}
double NeuralNet::randomWeight(){
	double randomD = (0.1 - -0.1) * ( (double)rand() / (double)RAND_MAX ) + -0.1;\
	//std::cout << randomD << std::endl;
	return randomD;
}
void NeuralNet::initBiasWeight(neuron **nodes){
	for (int y = 0; y < layers[0]; y++){
		nodes[0][y].biasWeight = new double[inputNodes];
		nodes[0][y].nodeWeight = randomWeight();
		nodes[0][y].out = 0;
	}
	for (int x = 1; x < numLayers; x++){
		for (int y = 0; y < layers[x]; y++){
			nodes[x][y].biasWeight = new double[layers[x - 1]];
			nodes[x][y].nodeWeight = randomWeight();
			nodes[x][y].out = 0;
		}
	}
	for (int y = 0; y < outputNodes; y++){
		nodes[numLayers][y].biasWeight = new double[layers[numLayers - 1]];
		nodes[numLayers][y].nodeWeight = randomWeight();
		nodes[numLayers][y].out = 0;
	}

	//*****FILL WEIGHTS*****
	for (int y = 0; y < layers[0]; y++){
		for (int z = 0; z < inputNodes; z++){
			nodes[0][y].biasWeight[z] = randomWeight();
		}
	}
	for (int x = 1; x < numLayers; x++){
		for (int y = 0; y < layers[x]; y++){
			for (int z = 0; z < layers[x - 1]; z++){
				nodes[x][y].biasWeight[z] = randomWeight();
			}
		}
	}
	for (int y = 0; y < outputNodes; y++){
		for (int z = 0; z < layers[numLayers - 1]; z++){
			nodes[numLayers][y].biasWeight[z] = randomWeight();
		}
	}


}
void NeuralNet::initLayers(neuron **nodes){
	//layers[0] = new neuron[this->inputNodes];
	for (int x = 0; x < (numLayers); x++){
		nodes[x] = new neuron[this->layers[x]];
	}
	nodes[numLayers] = new neuron[outputNodes];
}
void NeuralNet::feedForward(double *array, int inputs, int outputs){
	//double * answer = new double[inputs + outputs];
	//printWeights();
	//printPoint(array, 3, 1);
	for (int node = 0; node < layers[0]; node++){
		double sum = neurons[0][node].nodeWeight;
		for (int weight = 0; weight < inputs; weight++){
			sum += neurons[0][node].biasWeight[weight] * array[weight];
		}
		//std::cout << sum << std::endl;
		neurons[0][node].out = sigmoid(sum);
		//std::cout << "(0," << node << ") " << neurons[0][node].out << std::endl;
	}
	for (int layer = 1; layer < numLayers; layer++){
		for (int node = 0; node < layers[layer]; node++){
			double sum = neurons[layer][node].nodeWeight;
			for (int weight = 0; weight < layers[layer - 1]; weight++){
				sum += neurons[layer][node].biasWeight[weight] * neurons[layer - 1][weight].out;
			}
			//std::cout << sum << std::endl;
			neurons[layer][node].out = sigmoid(sum);
			//std::cout << "(" << layer << "," << node << ") " << neurons[layer][node].out << std::endl;
		}
	}
	for (int node = 0; node < outputs; node++){
			double sum = neurons[numLayers][node].nodeWeight;
			for (int weight = 0; weight < layers[numLayers - 1]; weight++){
				sum += neurons[numLayers][node].biasWeight[weight] * neurons[numLayers - 1][weight].out;
			}
			neurons[numLayers][node].out = sigmoid(sum);
			//std::cout << "(" << numLayers << "," << node << ") " << neurons[numLayers][node].out << std::endl;
		}




}
double *NeuralNet::solve(double *array, int inputs, int outputs){
		feedForward(array, inputs,outputs);
		double *answer = new double[inputs+outputs];
		for (int in = 0; in < inputs; in++){
			answer[in] = array[in];
		}
		for (int node = 0; node < (outputs); node++){
			//std::cout << neurons[numLayers][node].out << std:: endl;
			answer[node + inputs] = neurons[numLayers][node].out * normilization;
		}
		return answer;


}
double NeuralNet::sigmoid(double input){
	return 1/(1+pow(2.71828182845904523536, -input));
}
double NeuralNet::sigmoidPrime(double input){
	return input * (1 - input);
}
void NeuralNet::backProp(neuron *working, neuron *inputs, int workingNodes, int inputNodes){
	for (int node = 0; node < workingNodes; node++){
		double sum = 0;
		for (int weights = 0; weights < inputNodes; weights++){
			//std::cout << weights << ":" << inputs[weights].biasWeight[node] << std::endl;
			sum += inputs[weights].biasWeight[node] * inputs[weights].delta;

		}
		//std::cout << sum << std::endl;
		working[node].delta = sigmoidPrime(working[node].out) * sum;
		//std::cout << working[node].delta << std::endl;

	}
}
void NeuralNet::printPoint(double *point, int inputs, int outputs){
	for (int x = 0; x < (inputs + outputs); x++){
		std::cout << point[x] << ", ";
	}
	std::cout << std::endl;
}
void NeuralNet::learn(){
	time_t start = time(NULL);
	long long duration = 0;
	long long endTime = learningTime * 60;
	while (duration < endTime){
		for (int inputs = 0; inputs < (data->points * (data->inputs + data->outputs));inputs +=(inputNodes + outputNodes)){
			//std::cout << *(data->dataPoints + inputs) << std::endl;
			//printPoint((data->dataPoints + inputs));
			feedForward((data->dataPoints + inputs), inputNodes, outputNodes);
			for (int node = 0; node < outputNodes; node++){
				neurons[numLayers][node].delta = sigmoidPrime(neurons[numLayers][node].out) * (*(data->dataPoints + inputs + inputNodes)/normilization - neurons[numLayers][node].out);
				//std::cout << sigmoidPrime(neurons[numLayers][node].out) << ":" << neurons[numLayers][node].delta << std::endl;
			}
			backProp(neurons[numLayers - 1], neurons[numLayers], layers[numLayers - 1], outputNodes);
			for (int layer = (numLayers - 2); layer >= 0; layer--){
				backProp(neurons[layer], neurons[layer + 1], layers[layer],layers[layer + 1]);
			}
			updateWeights((data->dataPoints + inputs));
		}

		duration = difftime(time(0), start);
	}
}
void NeuralNet::updateWeights(double *input){
	//update first layer
	for (int node = 0; node < layers[0]; node++){
		//std::cout << "(" << 0 << "," << node << ") " << neurons[0][node].nodeWeight;
		neurons[0][node].nodeWeight += learningRate * neurons[0][node].delta;
		//std::cout << " -> " << neurons[0][node].nodeWeight << std::endl;
		for (int weight = 0; weight < inputNodes; weight++){
			//std::cout << "(" << 0 << "," << node << "," << weight << ") " << neurons[0][node].biasWeight[weight];
			neurons[0][node].biasWeight[weight] += (learningRate * *(input + weight) * neurons[0][node].delta);
			//std::cout << " -> " << neurons[0][node].biasWeight[weight] <<std::endl;

		}

	}
	//update middle layers
	for (int layer = 1; layer < numLayers; layer++){
		for (int node = 0; node < layers[layer]; node++){
			//std::cout << "(" << layer << "," << node << ") " << neurons[layer][node].nodeWeight;
			neurons[layer][node].nodeWeight += learningRate * neurons[layer][node].delta;
			//std::cout << " -> " << neurons[layer][node].nodeWeight << std::endl;
			for (int weight = 0; weight < layers[layer - 1]; weight++){
				//std::cout << "(" << layer << "," << node << "," << weight << ") " << neurons[layer][node].biasWeight[weight];
				neurons[layer][node].biasWeight[weight] += (learningRate * neurons[layer - 1][weight].out * neurons[layer][node].delta);
				//std::cout << " -> " << neurons[layer][node].biasWeight[weight] <<std::endl;

			}

		}
	}
	//update last layer
	for (int node = 0; node < outputNodes; node++){
		//std::cout << "(" << numLayers << "," << node << ") " << neurons[numLayers][node].nodeWeight;
		neurons[numLayers][node].nodeWeight += learningRate * neurons[numLayers][node].delta;
		//std::cout << " -> " << neurons[numLayers][node].nodeWeight << std::endl;
		for (int weight = 0; weight < layers[numLayers - 1]; weight++){
			//std::cout << "(" << numLayers << "," << node << "," << weight << ") " << neurons[numLayers][node].biasWeight[weight];
			neurons[numLayers][node].biasWeight[weight] += (learningRate * neurons[numLayers - 1][weight].out * neurons[numLayers][node].delta);
			//std::cout << " -> " << neurons[numLayers][node].biasWeight[weight] <<std::endl;
		}
	}
}
void NeuralNet::printWeights(){
	for (int node = 0; node < layers[0]; node++){
		std::cout << "(" << 0 << "," << node << ") " << neurons[0][node].nodeWeight<< std::endl;
		for (int weight = 0; weight < inputNodes; weight++){
			std::cout << "(" << 0 << "," << node << ") " << neurons[0][node].biasWeight[weight] << std::endl;
		}

	}
	for (int layer = 1; layer < numLayers; layer++){
		for (int node = 0; node < layers[layer]; node++){
			std::cout << "(" << layer << "," << node << ") " << neurons[layer][node].nodeWeight<< std::endl;
			for (int weight = 0; weight < layers[layer - 1]; weight++){
				std::cout << "(" << layer << "," << node << ") " << neurons[layer][node].biasWeight[weight] << std::endl;
			}
		}
	}
	for (int node = 0; node < outputNodes; node++){
		std::cout << "(" << numLayers << "," << node << ") " << neurons[numLayers][node].nodeWeight<< std::endl;
		for (int weight = 0; weight < layers[numLayers - 1]; weight++){
			std::cout << "(" << numLayers << "," << node << ") " << neurons[numLayers][node].biasWeight[weight] << std::endl;
		}
	}
}
void NeuralNet::printSolutions(){
	testPoints(data);
}

double NeuralNet::getError(dataSet *data){
	double error = 0;
	for (int point = 0; point < data->points; point++){
		double pointDiff = 0;
		double *test = solve((data->dataPoints + point * (data->inputs + data->outputs)), data->inputs, data->outputs);
		for (int answer = data->inputs; answer < (data->inputs + data->outputs); answer++){
			pointDiff += (data->dataPoints + point * (data->inputs + data->outputs))[answer]/normilization - test[answer]/normilization;
			//std::cout << pointDiff;
		}
		error += pointDiff;
		//std::cout << pointDiff << std::endl;
	}
	//std::cout << "testing\n";
	return error;
}
void NeuralNet::testPoints(dataSet *data){
	int pointWidth = data->inputs + data->outputs;
	for (int point = 0; point < data->points; point++){
		std::cout << "Inputs: ";
		for (int input = 0; input < data->inputs; input++){
			std::cout << *(data->dataPoints + point * pointWidth + input);
			if (input < (data->inputs - 1)){
				std::cout << ", ";
			}
		}
		std::cout << std::endl;
		std::cout << "output: ";
		double *answer = solve((data->dataPoints + (point * pointWidth)), data->inputs, data->outputs);
		for (int output = data->inputs; output < (data->inputs + data->outputs); output++){
			std::cout << answer[output];
				if (output < (data->inputs + data->outputs - 1)){
					std::cout << ", ";
				}
		}
		std::cout << std::endl;
		std::cout << "Expected: ";
		for (int output = data->inputs; output < (data->inputs + data->outputs); output++){
			std::cout << data->dataPoints[point * pointWidth + output];
			if (output < (data->inputs + data->outputs - 1)){
				std::cout << ", ";
			}
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
}
dataSet* NeuralNet::splitData(dataSet* data,int numContainers){
	int dataSize = data->points/numContainers;
	dataSet* containers = new dataSet[numContainers];

	for (int container = 0; container < numContainers; container ++){
		containers[container].inputs = data->inputs;
		containers[container].outputs = data->outputs;
		containers[container].scale = data->scale;
		containers[container].points = dataSize;
		//containers[container].dataPoints = new double[(dataSize * (data->inputs + data->outputs))];

		containers[container].dataPoints = (data->dataPoints + (container * (dataSize * (data->inputs + data->outputs))));

		//std::cout << "Container " << container << ":" << std::endl;
		//printData(containers + container);
	}

	return containers;
}
dataSet* NeuralNet::randomizeData(dataSet* data){
	srand(time(NULL));
	dataSet *randomOrder = new dataSet;
	randomOrder->inputs = data->inputs;
	randomOrder->outputs = data->outputs;
	randomOrder->points = data->points;
	randomOrder->scale = data->scale;
	randomOrder->dataPoints = new double[data->points * (data->inputs + data->outputs)];
	bool used[randomOrder->points] = {0};

	int out = rand() % randomOrder->points;
	for (int point = 0; point < data->points; point++){

		while (used[out] == true){
			out = rand() % randomOrder->points;
		}
		//std::cout << out << std::endl;
		for (int item = 0; item < (data->inputs + data->outputs); item++){
			randomOrder->dataPoints[out * (data->inputs + data->outputs) + item] = data->dataPoints[point * (data->inputs + data->outputs) + item];
		}

		used[out] = true;
	}
	//printData(randomOrder);
	deleteData(data);//*/
	return randomOrder;
}
dataSet* NeuralNet::joinData(dataSet* sets, int numSets ,int exclude){
	dataSet *combined = new dataSet;
	combined->points = sets[0].points * (numSets - 1);
	combined->inputs = sets[0].inputs;
	combined->outputs = sets[0].outputs;
	combined->scale = sets[0].scale;
	combined->dataPoints = new double[(combined->points * (combined->inputs + combined->outputs))];
	for (int combinddIndex = 0;  combinddIndex < (combined->points * (combined->inputs + combined->outputs));){
		for (int set = 0; set < numSets; set++){
			if (set != exclude){
				for (int setIndex = 0; setIndex < (sets[set].points * (sets[set].inputs + sets[set].outputs)); combinddIndex++, setIndex++){
									combined->dataPoints[combinddIndex] = sets[set].dataPoints[setIndex];
									//std::cout << combined->dataPoints[combinddIndex] << std::endl;

				}

			}
		}
	}
	//printData(combined);
	return combined;
}
void NeuralNet::printData(dataSet* data){
	for (int point = 0; point < data->points; point++){
		for (int inputs = 0; inputs < (data->inputs + data->outputs); inputs++){
			std::cout << data->dataPoints[(point * (data->inputs + data->outputs) + inputs)] << " , ";

		}
		std::cout << std::endl;
	}
}


dataSet* NeuralNet::readFile(const char* name){
	dataSet *set;
	set = new dataSet;
	int numLines = countLines(name);
	numLines -= 2;
	set->points = numLines;
	std::fstream input;
	input.open(name);
	std::string line;


	input >> set->scale;

	std::getline(input, line, ',');
	//std::cout << set->scale << std::endl;
	std::stringstream fileForm(line);
	fileForm >> set->inputs;
	//fileForm >> set->outputs;
	//set->inputs = atoi(line.c_str());
	input >> set->outputs;

	set->dataPoints = new double[numLines * (set->inputs + set->outputs)];
	//char delim[] = { ' ' , ',' , '\n'};

	for (int x = 0; x < (numLines * (set->inputs + set->outputs)); x= x+(set->inputs + set->outputs)){
		for (int y = 0; y < (set->inputs + set->outputs - 1); y++){
			std::getline(input, line, ',');
			std::stringstream ss(line);
			ss >> set->dataPoints[x+y];

		}
		input >> set->dataPoints[x + (set->inputs + set->outputs - 1)];
		//std::cout << "Adding to array: " << set->dataPoints[x]<< " at " << x << std::endl;

		//set->dataPoints[x] = atoi(line.c_str());
	}

	/*
	for (int x = 0 ; x < (numLines * (set->inputs + set->outputs)); x = x + 4){
		std::cout << set->dataPoints[x] <<  " , " << set->dataPoints[x+1] << " , " << set->dataPoints[x+2] << " = " << set->dataPoints[x+3] << std::endl;
	}//*/

	//std::cout << "Testing: " << numLines <<  " "<< set->inputs << " " << set->outputs << std::endl;


	return set;
}
int NeuralNet::countLines(const char* filename){
	int count = 0;
	std::ifstream file;
	file.open(filename);
	std::string line;

	if (file.is_open()){
		while(std::getline(file, line))
			count++;
	}
	file.close();

	return count;
}

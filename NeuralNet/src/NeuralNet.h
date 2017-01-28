/*
 * NeuralNet.h
 *
 *  Created on: Nov 15, 2016
 *      Author: Aaron
 */

#ifndef NEURALNET_H_
#define NEURALNET_H_

struct neuron{
	double *biasWeight;
	double nodeWeight;
	double delta;


	double out;
};

struct dataSet{
	double scale;
	int inputs;
	int outputs;
	int points;

	double *dataPoints;
};

class NeuralNet {
public:
	NeuralNet(dataSet *workingSet, int numLayers, int* layers, double learningRate, double learningTime);
	virtual ~NeuralNet();
	static dataSet* readFile(const char* name);
	static void deleteData(dataSet* data);
	static dataSet* splitData(dataSet* data,int containers);
	static dataSet* randomizeData(dataSet* data);
	static dataSet* joinData(dataSet* sets, int numSets, int exclude);
	static void printData(dataSet* data);
	void testPoints(dataSet* data);
	double getError(dataSet *data);

	void printSolutions();

private:
	dataSet *data;
	neuron *inputs;

	int numLayers;
	int* layers;
	double learningRate;
	double learningTime;

	neuron **neurons;

	int normilization;
	int inputNodes, outputNodes;


	static int countLines(const char* filename);
	void initLayers(neuron **layers);
	void initBiasWeight(neuron **nodes);
	double randomWeight();
	void printNode(neuron &node, int layer);
	void printNetwork();
	void deleteNet();
	void runNet();
	void feedForward(double *array, int inputs, int outputs);
	double *solve(double *array, int inputs, int outputs);
	double sigmoid(double input);
	double sigmoidPrime(double input);
	void backProp(neuron *working, neuron *inputs, int workingNodes, int inputNodes);
	void learn();
	void updateWeights(double *input);
	void printWeights();
	static void printPoint(double *point, int inputs, int outputs);



};

#endif /* NEURALNET_H_ */

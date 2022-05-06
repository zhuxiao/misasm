#include <iostream>
#include <string>

#include "Paras.h"
#include "util.h"

#include "Genome.h"

using namespace std;

int main(int argc, char **argv){

	Paras paras(argc, argv);

	// output parameters
	paras.outputParas();

	Genome genome(&paras);

	// estimate the parameters for noisy background
	genome.estimateSVSizeNum();
	paras.outputEstParas("After estimation:");

	genome.generateGenomeBlocks();
	genome.genomeIlluminaMisIdentify();

	return 0;
}




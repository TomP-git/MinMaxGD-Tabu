#include<iostream>
#include "input_algorithm.h"
#include "parser_bcmp_instance.h"
#include "HDAG.h"
#include "TABU.h"


int main(int argc, char *argv[]){

	unsigned input = input_manager(argc, argv);

	parse_instance(argc, argv);

	HDAG I; // Instance of the problem, a hierarchical graph
	string instance = argv[1];
	instance = "parsed_" + instance;
	I.read_instance(instance.c_str()); // read the instance

	// case 4 corresponds to the Tabu Search. This is the only
	//meta-heuristic currently implemented in the present code
	switch (input) {
	case 4:

		TABU T(I, argv);

		break;

	}




	return 0;
}




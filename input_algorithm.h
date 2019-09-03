#include <iostream>
#include <sstream>
#include <iomanip>
#include <type_traits>
#include <string.h>
#include <map>

#include <cstdint>

#ifndef INPUT_ALGORITHM_H_
#define INPUT_ALGORITHM_H_

using namespace std;

unsigned check_algorithm_name(int argc, char *argv[]) {

	char algorithm[20];
	strcpy(algorithm, argv[2]);

	if (strcmp(algorithm, "cplex") == 0 ) return 0;

	if (strcmp(algorithm, "grasp") == 0 ) return 1;

	if (strcmp(algorithm, "tabu") == 0 ) return 4;

	return 5;
}

bool correct_cplex_parameter(int argc, char *argv[]) {

	// - check the parameter [3], the k value
	int k = atoi(argv[3]);
	if (k < 0) {
		cerr << " k-value (arg[3]) must be an integer positive value" << endl;
		exit(0);
	}
	// - check the parameter [6], the time-limit
	int tl = atoi(argv[6]);
	if (tl < 10 || tl > 3600) {
		cerr << " time-limit (arg[6]) must be an integer positive value in the range [10, 3600]" << endl;
		exit(0);
	}

	return true;

}

bool correct_grasp_parameter(int argc, char *argv[]) {

	// - check the parameter[3], constructive phase
	char *cc = argv[3];
	if ((strcmp(cc, "C1") != 0) && (strcmp(cc, "c1") != 0)) {
		if ((strcmp(cc, "C2") != 0) && (strcmp(cc, "c2") != 0)) {
			if ((strcmp(cc, "C3") != 0) && (strcmp(cc, "c3") != 0)) {
				cerr << " Selected constructive (arg[3]) incorrect!" << endl;
				cerr << " - Possible choices {C1 or c1, C2 or c2, C3 or c3}" << endl;
				exit(0);
			}
		}
	}

	// - check the parameter [4], the alpha value
	double alpha = atof(argv[4]);
	if (alpha < 0 || alpha > 1) {
		cerr << " apha-value (arg[4]) must be a real nonnegative value in the range [0,1]" << endl;
		exit(0);
	}

	// - check the parameter [5], local_search
	char *ls = argv[5];
	if (strcmp(ls, "ls") != 0 ) {
		if (strcmp(ls, "ts") != 0 ) {
			if (strcmp(ls, "lsx") != 0 ) {
				if (strcmp(ls, "no") != 0 ) {
					cerr << " Selected local-search (arg[5]) incorrect!" << endl;
					cerr << " - Possible choices {ls, lsx, ts, no}" << endl;
					exit(0);
				}
			}
		}
	}
	// - check the parameter [6], max_iteration
	int maxit = atoi(argv[6]);
	if (maxit < 1 || maxit > 10000000) {
		cerr << " max-it (arg[6]) must be an integer positive value in the range [10, 10000000]" << endl;
		exit(0);
	}
	// - check the parameter [10], the time-limit
	int tl = atoi(argv[10]);
	if (tl < 10 || tl > 3600) {
		cerr << " time-limit (arg[10]) must be an integer positive value in the range [10, 3600]" << endl;
		exit(0);
	}
	// - check the parameters [11], [12] and [13]:
	// - [11]: type of path relinking
	// - [12]: dimension of the elite set
	// - [13]: threshold for DIVERSITY DEGREE
	if (argc > 11) {
		if (argc < 14) {
			cerr << " the parameters [11], [12], and [13] must be :" << endl;
			cerr << " - [11]: type of path relinking must be a char in {f(or F), b(or B), m(or M)} "
					"(f/F: forward, b/B: backward, m/M: mixed)" << endl;
			cerr << " - [12]: dimension of the elite set, an integer positive value in the range [5, 20]" << endl;
			cerr << " - [13]: threshold for access to the elite set, a real value in the range [0, 1]" << endl;
			exit(0);
		}
		else {
			char *t_PR = argv[11];
			int dim_ES = atoi(argv[12]);
			double th = atof(argv[13]);
			if (strcmp(t_PR, "f") != 0 && strcmp(t_PR, "F") != 0 &&
				strcmp(t_PR, "b") != 0 && strcmp(t_PR, "B") != 0 &&
				strcmp(t_PR, "m") != 0 && strcmp(t_PR, "M") != 0 ) {
				cerr << " - [11]: type of path relinking must be a char in {f(or F), b(or B), m(or M)} "
						"(f/F: forward, b/B: backward, m/M: mixed)" << endl;
				exit(0);
			}
			if (dim_ES < 5 || dim_ES > 20) {
				cerr << " - [12]: dimension of the elite set, an integer positive value in the range [5, 20]" << endl;
				exit(0);
			}
			if (th < 0 || th > 1) {
				cerr << " - [13]: threshold for access to the elite set, a real value in the range [0, 1]" << endl;
				exit(0);
			}
		}
	}

	return true;

}

bool correct_tabu_parameter(int argc, char *argv[]) {

	// - check the parameter [3], the k value
	int k = atoi(argv[3]);
	if (k < 0) {
		cerr << " the tenure (arg[3]) must be a non-negative value" << endl;
		exit(0);
	}
	// - check the parameter [4], the alpha value
	double alpha = atof(argv[4]);
	if (alpha < 0 || alpha > 1) {
		cerr << " alpha-value (arg[4]) must be a real nonnegative value in the range [0,1]" << endl;
		exit(0);
	}
	// - check the parameter [5], post_processing
	char *ls = argv[5];
		if (strcmp(ls, "yes") != 0 ) {
			if (strcmp(ls, "no") != 0 ) {
					cerr << " Pre-processing parameter (arg[5]) incorrect!" << endl;
					cerr << " - Possible choices {yes, no}" << endl;
					exit(0);
			}
		}

	// - check the parameter [6], max_iteration
	int maxit = atoi(argv[6]);
	if (maxit < 1 || maxit > 10000000) {
		cerr << " max-it (arg[6]) must be an integer positive value in the range [10, 10000000]" << endl;
		exit(0);
	}
	// - check the parameter [10], the time-limit
	int tl = atoi(argv[10]);
	if (tl < 10 || tl > 3600) {
		cerr << " time-limit (arg[10]) must be an integer positive value in the range [10, 3600]" << endl;
		exit(0);
	}

	unsigned Max_no_improv = atoi(argv[11]);
	if (Max_no_improv <=0 ){
		cerr << "The maximum number without improvement allowed should be a positive integer."<<endl;
		exit(0);
	}

	unsigned MAX_DIVERSE = atoi(argv[12]);
	if (MAX_DIVERSE < 0){
		cerr<< "The maximum number of diversification step performed each time should be a positive integer, or zero for no diversification."<< endl;
		exit(0);
	}

	return true;

}

unsigned input_manager(char argc, char *argv[]) {

	if ( ( strcmp(argv[1], "h") == 0) ||
		 ( strcmp(argv[1], "help") == 0) ||
		 ( strcmp(argv[1], "-h") == 0) ||
		 ( strcmp(argv[1], "-help") == 0) ||
		 ( strcmp(argv[1], "H") == 0) ||
		 ( strcmp(argv[1], "HELP") == 0) ||
		 ( strcmp(argv[1], "-HELP") == 0) ||
		 ( strcmp(argv[1], "-H") == 0) ) {

		cout << " For TABU algorithm: " << endl;
		cout << "  ./exec " << endl;
		cout << "  argv[1]: instance.txt" << endl;
		cout << "  argv[2]: algorithm={tabu}" << endl;
		cout << "  argv[3]: tenure " << endl;
		cout << "  argv[4]: alpha=[0:1] " << endl;
		cout << "  argv[5]: post-processing={yes, no} " << endl;
		cout << "  argv[6]: max-it={0,..,10000}" << endl;
		cout << "  argv[7]: instance (without .txt) " << endl;
		cout << "  argv[8]: out-file-summary=<string>" << endl;
		cout << "  argv[9]: out-file-complete=<string>" << endl;
		cout << " argv[10]: time-limit={0,..,3600}" << endl;
		cout << " argv[11]: max-it without improvement before diversification"  << endl;
		cout << " argv[12]: max diversification moves"<< endl;
		cout << " argv[13]: random seed" << endl;
		cout << " EXAMPLE: ./exec instance.txt tabu 5 1 no 10000 instance out-f1.txt out-f2.txt 60 20 10 1" << endl;
		exit(0);
	}

	unsigned algorithm = check_algorithm_name(argc, argv);

	switch(algorithm) {
	case(0):
	// "Execute cplex";
	{

		if (correct_cplex_parameter(argc, argv))
			return 0;

 	}
	break;
	case(1):
	{

		if (correct_grasp_parameter(argc, argv))
			return 1;

	}
	break;
	case(2):
	{

		if (correct_grasp_parameter(argc, argv))
			return 2;

	}
	break;
	case(3):
	{

		if (correct_grasp_parameter(argc, argv))
			return 3;

	}
	break;
	case(4):
	{
		if (correct_tabu_parameter(argc, argv))
			return 4;
	}
	break;
	}

	return 0;
}


#endif /* INPUT_ALGORITHM_H_ */

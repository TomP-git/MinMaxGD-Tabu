#include<vector>
#include<unordered_map>
#include<limits>
#include<chrono>
#include"HDAG.h"
#include<random>


#ifndef TABU_H_
#define TABU_H_

class TABU {
public:
	TABU( HDAG &I, char *argv[]);
	virtual ~TABU();

	//Computes the maximum number of crossings for arcs incident to non-TABU
	//nodes in level l. Updates max_NonTabu_v[l] accordingly
	void compute_maxNonTabu_inLev(HDAG &S,unsigned l);

	//Computes the probabilities used in the diversification phase
	void compute_SelectionProb_v();

	//Fills CL[l] according to the value of alpha and maxNonTabu_v[l]
	void fill_CandidateList(HDAG &S, unsigned l);

	//run the TABU search considering only the layer with the maximum
	void max_layer_TABU(HDAG &S);

	//run the TABU search considering only layers with an high number of crossings
	void high_cross_layers_TABU(HDAG &S);

	void run_max_TABU(HDAG &I);

	//run the TABU search
	void run_TABU(HDAG &I);

	//Sweep the levels to find improving moves among non TABU nodes
	void sweep_TABU(HDAG &S);

	void diversificate_search(HDAG &S, int seed);

	void reduce_total_crossing(HDAG &S);

	/*******************Allocation Methods*************************/

	void allocate_all_structures(HDAG &I);

	void allocate_TabuTrack(unsigned l);

	//allocates the space for level l of TabuTrack.
	//Initializes all the values to 0
	void allocate_TabuTrackLev(unsigned l, unsigned t);

	//allocates the space for maxNonTabu, and initializes all the values to 0
	void allocate_maxNonTabu_v(unsigned l);

	void allocate_CL(unsigned l);

	void allocate_SwapCount(unsigned l){
		SwapCount.resize(l);
	}

	void allocate_Prob_v(unsigned l){
		SelectionProb_v.resize(l);
	}

	void allocate_SwapCountLevel(unsigned l, unsigned nod){
		SwapCount[l].resize(nod);
	}

	void allocate_Prob_v_lev(unsigned l, unsigned nod){
		SelectionProb_v[l].resize(nod);
	}






	/*******************Setters & Getters**************************/

	std::vector<std::vector<unsigned> > &getTabuTrack() {
		return TabuTrack;
	}

	std::vector<std::vector<unsigned> > &getCL() {
		return CL;
	}

	unsigned get_tenure() const {
		return tenure;
	}

	unsigned get_r_Seed() const {
		return r_seed;
	}

	double get_alpha() const {
		return alpha;
	}

	void set_tenure(unsigned t) {
		tenure = t;
	}

	void set_r_seed(unsigned r){
		r_seed = r;
	}

	void set_currentIter(unsigned it){
		currentIter = it;
	}

	void set_alpha(double a){
		alpha = a;
	}

	void set_TimeLimit(int tl){
		TIME_LIMIT = tl;
	}

	void set_MAXITER(int IT){
		MAX_ITER = IT;
	}

	void set_post_proc(char* pp){
		post_proc = pp;
	}

	void set_Max_no_improv(unsigned IT){
		Max_no_improv = IT;
	}

	void set_MAX_DIVERSE(unsigned IT){
		MAX_DIVERSE = IT;
	}

	chrono::time_point<chrono::system_clock> &get_start_iteration() {return start_iteration;}
	chrono::time_point<chrono::system_clock> &get_end_iteration() {return end_iteration;}
	chrono::time_point<chrono::system_clock> &get_start_algorithm() {return start_algorithm;}

	void set_time_to_best(double x) {time_to_best = x;}

	void set_start_iteration(chrono::time_point<chrono::system_clock> x) {start_iteration = x;}
	void set_end_iteration(chrono::time_point<chrono::system_clock> x) {end_iteration = x;}
	void set_start_algorithm(chrono::time_point<chrono::system_clock> x) {start_algorithm = x;}

	double get_total_time_algorithm() {
		int elapsed =
				chrono::duration_cast<chrono::milliseconds>(end_iteration - start_algorithm).count();
		return (double)((double)elapsed / (double)1000);
	}


	void set_best_out_file(char *filename) { out_f_best = fopen(filename, "a"); };
	void set_complete_out_file(char *filename) { out_f_complete = fopen(filename, "a"); };
    void print_best_data_out_file() {
    	fprintf(out_f_best, "%0.2f %0.2f %0.2f %0.3f %0.3f\n",
    			alpha,
    			best_solution.update_cost_bottleneck(),
				best_solution.compute_crossing_sum_bottleneck(),
				get_total_time_algorithm(),
				time_to_best);
    }

     void print_complete_data_out_file(double actual_z, double current_t) {
    	fprintf(out_f_complete, "%0.2f %d %0.3f %0.3f\n",
    			alpha, tenure, actual_z, current_t);
    }
	void close_out_file() { fclose(out_f_best); fclose(out_f_complete); };

	/**************************************************************/

private:

	//Counter for current Iteration.
	unsigned currentIter;

	//random seed for the diversification phase
	unsigned r_seed;

	//Maximum number of Iterations.
	unsigned MAX_ITER;

	//switch on/off post-processing
	char* post_proc;

	//maximum number of iterations without improvement allowed
	unsigned Max_no_improv;

	//Maximum number of diversification step allowed each time
	unsigned MAX_DIVERSE;

	//keep track of improvements
	unsigned last_improv;

	//Number of iteration in which a solution holds the TABU status
	unsigned tenure;

	//Time limit in seconds for the algorithm
	int TIME_LIMIT;

	//Multiplier (0 <= alpha <= 1) used for the construction of the Candidate
	//List. If for a node i there exists an arc [i,j] whose number of crossings
	//is greater than or equal to alpha multiplied the max cross among the arcs
	//incident to non-TABU nodes in that specific layer, then i is in CL.
	//the CL is used to decide which node is elected for a move.
	double alpha;

	//Structure to track the TABU status of each node in each layer
	//TabuTrack[l][i] holds the last time in which node i of layer l has
	//changed position
	std::vector< std::vector < unsigned > > TabuTrack;

	//Structure to count the number of swaps for each node in each layer
	//SwapCount[l][i] holds the number of swaps that involved node with
	//id i in level l.
	std::vector<std::vector< unsigned > > SwapCount;

	//Structure to hold the probabilities used for random selection in the
	//diversification phase. SelectionProb_v[l][i] holds the probability
	//of selectiong the node with id i of level l
	std::vector<std::vector< double > > SelectionProb_v;

	//Structure to hold the position of the nodes in each layer that are
	//candidates for a move. A position p[i] is inserted in CL[l] iff
	//exists a node i (of layer l) such that there is an arc (i,j)
	//whose number of crossings is greater than or equal to alpha multiplied
	//the max cross among the arcs incident to non-TABU nodes in that specific layer, then i is in CL.
	std::vector< std::vector <unsigned > > CL;

	//Vector to hold the maximum number of crossing among arcs incident to nodes
	//not in the TABU-list
	std::vector< unsigned> maxNonTabu_v;


	//HDAG to hold the best solution found in the search process
	HDAG best_solution;
	//time needed to reach best solution found
	double time_to_best;

	//Time related members
	chrono::time_point<chrono::system_clock> start_algorithm;
	chrono::time_point<chrono::system_clock> start_iteration;
	chrono::time_point<chrono::system_clock> end_iteration;

	// file to store parameters about the best solution found
	FILE *out_f_best;

	// file to store all the parameters relative to all the iterations done
	FILE *out_f_complete;



};

#endif /* TABU_H_ */

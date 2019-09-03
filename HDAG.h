#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>


using namespace std;

#ifndef HDAG_H_
#define HDAG_H_

class HDAG {
public:
	HDAG();
	virtual ~HDAG();

	// HDAG represents a hierarchical graph, to maintain all the graph are used five different structure
	// 1) LEVELS: a 3D structure where the element LEVELS[l][i][k] establishes:
	//    - the l^th level of the hierarchical graph
	//    - the i^th node in the l^th level of the hierarchical graph
	//    - the k^th node in the FORWARD star of the i^th node in the l^th level of the hierarchical graph
	// 2) B_LEVELS: a 3D structure where the element B_LEVELS[l][i][k] establishes:
	//    - the l^th level of the hierarchical graph
	//    - the i^th node in the l^th level of the hierarchical graph
	//    - the k^th node in the BACKWARD star of the i^th node in the l^th level of the hierarchical graph
	// 3) IDs: a 2D structure where the element IDs[l][i] establishes:
	//    - the identifier of the node in LEVELS[l][i]
	// 4) Pos: a 2D structure where given a level l and an identifier id establishes:
	//    - the position of the node "id" in the level l, let's suppose that the
	//      position of "id" in l is given by "pos" then the node "id" in level l
	//      will be in the location LEVELS[l][pos]
	// 5) Os: a 2D structure where the element Os[l][i] establishes:
	//    - if LEVELS[l][i] = 1 then the node is an original node
	//    - if LEVELS[l][i] = 0 then the node is an incremental node.

	// --------------------->>>>>>>>>> ATTENTION <<<<<<<<<<<-------------------------
	// For each function used in this class when in the input are used "l" and "i"
	// to identify the element LEVEL[l][i] is ever referred to the original instance,
	// usually named "I", while with "S" is indicated the actual solution.
	// --------------------->>>>>>>>>> ATTENTION <<<<<<<<<<<-------------------------

	// Read the instance from a file
	void read_instance(string filename);

	// Return the number of level of the hierarchical graph
	unsigned getLevNumber();

	// Return the maximum degree in the hierarchical graph
	// taking into account both in and out degree
	unsigned getMaxDegree();

	// A 3D matrix where for each element [i][j][k]
	// - i represents the i^th level of the hierarchical graph
	// - j represents the j^th node in the i^th level of the hierarchical graph
	// - k represents the k^th node in the forward star of the j^th node in the i^th level
	// - Remember that LEVELS[i][j][k] is the identifier of the k^th node in th FS of the j^th node
	//   in the i^th level.
	vector< vector< vector< unsigned> > > &getLEVELS ();

	// A 3D matrix where for each element [i][j][k]
	// - i represents the i^th level of the hierarchical graph
	// - j represents the j^th node in the i^th level of the hierarchical graph
	// - k represents the k^th node in the backward star of the j^th node in the i^th level
	// - Remember that B_LEVELS[i][j][k] is the identifier of the k^th node in th BS of the j^th node
	//   in the i^th level
	vector< vector< vector< unsigned> > > &getB_LEVELS ();

	// A 3D matrix, of the same size of LEVELS, such that each element in [i][j][k]
	// represents the current number of crossings for the arc connecting the j^th node
	// of the i^th level, with the k^th node of its forward star.
	vector< vector< vector< unsigned> > > &getCurrentCostMatrix ();

	vector< unordered_map < unsigned, unsigned> > &getHowManyWCost ();

	vector< vector < bool> > &getAnyWCost ();

	vector< unsigned> &getTotalCostDistr ();

	void build_maxdegree_node_structure();

	void build_critical_edges_structure();

	void allocate_BottleneckCostStructures();

	// Allocate the space for the level of the hierarchical graph
	void allocateLEVELS (unsigned l);

	// Allocate the space for a number of node "n" in the level "l"
	void allocateLevel (unsigned l, unsigned n);

	void allocateCurrentCostMatrix (unsigned l);

	void allocateCCMatrixLevel (unsigned l, unsigned n);

	void allocateTotalCostDistr (unsigned h);

	void allocateHowManyWCost (unsigned h);

	void allocateAnyWCost (unsigned h);

	void allocateAnyWCostLevel (unsigned h, unsigned n);

	// Allocate the space for a number of level of the hierarchical graph
	void allocateIDs (unsigned l);

	// Allocate the space for a number of node "n" in the level "l"
	void allocateIDs (unsigned l, unsigned n);

	// Allocate the space for a number of level of the hierarchical graph
	void allocateOs (unsigned l);

	// Allocate the space for a number of node "n" in the level "l"
	void allocateOs (unsigned l, unsigned n);

	// Allocate the space for a number of level of the hierarchical graph
	void allocatePos (unsigned l);

	// Allocate the space for a number of node "n" in the level "l"
	void allocatePos (unsigned l, unsigned n);

	// Allocate space for a number of h of coefficient for the computation
	// of the score function
	void allocateScoreCVector(unsigned h);

	// Return a specific level of the HDAG
	vector <vector< unsigned > > &getLevel(unsigned i);

	// Return the Identifiers Matrix of the HDAG
	vector <vector< unsigned > > &getIDs();

	// Return the Original nodes Maxtrix of the HDAG
	vector <vector< unsigned > > &getOs();

	// Return the Positions Matrix of the HDAG
	vector <vector< unsigned > > &getPos();

	// Return the DEG Matrix of the HDAG
	vector <vector< unsigned > > &getDEG();

	// Return the vector which contains the pair of vertices whit max degree
	vector < pair<unsigned, unsigned > > &get_v_maxdeg();

	void build_MAX() {
		n_MAX.resize(LEVELS.size() + 1);
		for (unsigned l = 1; l < LEVELS.size() + 1; ++l) {
			n_MAX[l].first = 0; n_MAX[l].second = 0;
		}
	}

	void buildBlev();

	void recompute_n_MAX_in_level (unsigned l);

	//Returns the max number of crossing of and arc incident to node i of level l
	unsigned get_MaxCrossing_incidentTo(unsigned l, unsigned i);

	vector< pair<double, unsigned> > &get_n_MAX() { return n_MAX; };

	// Establish if the arcs (i,j) and (i_,j_) are crossings in the level l
	bool areCrossingEdge(unsigned i, unsigned j, unsigned i_, unsigned j_, unsigned l);

	// Establishes if the node i^th in the level l of the hierarchical graph is
	// an original node
	bool isOriginalNode(unsigned i, unsigned l);

	// Establishes if the arc (i,j) in the level l of the hierarchical graph
	// is made up only original node
	bool isOriginalLink(unsigned i, unsigned j, unsigned l);

	// Establishes if the node i^th in the level l of the hierarchical graph is
	// an incremental node
	bool isIncrementalNode(unsigned i, unsigned l);

	// Swaps the nodes in position p1 and in position p2 in the level l of the
	// hierarchical graph. The swap involves all the structure used to maintain the
	// hierarchical graph: LEVELS, IDs, Os, Pos. The swap involves bottleneck specific
	// structure CurrentCostMatrix as well.
	void swapPositions(unsigned p1, unsigned p2, unsigned l);


	// Compute the degree of the i^th node in the l^th level of
	// the hierarchical graph
	unsigned compute_degree(unsigned l, unsigned i);

	// Compute the barycenter for the node in LEVELS[l][i]
	unsigned compute_barycenter(vector<vector<unsigned>> &CL, HDAG &I, unsigned l, unsigned i);

	// Print in a file filename (.tex) using tikz the hierarchical graph specified by G
	void printHDAG(HDAG &G, string filename);

	// Return the number of original node for in the hierarchical graph G
	unsigned get_n_original_nodes(HDAG &G);

	// Return the number of total nodes in the hierarchical graph
	unsigned get_total_nodes() { return total_nodes; }

	/*** ARE USED ONLY IN THE READING PHASE OF THE INSTANCE ****/
	// Is used to set the number of total nodes
	void set_total_nodes(unsigned x) { total_nodes = x; };
	// Is used to increment the number of total nodes by one
	void increase_total_nodes(unsigned x) { total_nodes+=x; };
	/**********************************************************/
	// Utility function to make a copy of a hierarchical graph
	void copy_HDAG(HDAG &G);
	void copy_level(HDAG &G, unsigned l);
	void dealloc_HDAG();
	/*********************************************************/

	//Updates the value of the variable max_crossing according to
	//the TotalCostDistr vector
	void update_max_crossing();

	double get_cost();

	double get_cost_bottleneck();

	//returns the max crossing sweeping backwards TotalCostDistr.
	double update_cost_bottleneck();

	//computes the total number of crossings using the bottleneck structure
	//TotalCostDistr
	double compute_crossing_sum_bottleneck();

	//updates bottleneck structures after swapping positions p1 and p2 in level l
	void update_bottleneck_str_afterSwap(unsigned p1, unsigned p2, unsigned l);

	//set to zero bottleneck structures in order to properly recompute cost
	//from scratch with get_cost_bottleneck().
	void reset_bottleneck_structures();

	void set_cost(double x) { cost_changed = false; cost = x; };
	bool is_cost_changed() { return cost_changed; };
	void set_cost_changed(bool x) { cost_changed = x; };

	double get_total_MAX() { return total_MAX; };
	void set_total_MAX(double x) { total_MAX = x; };

	double get_crossings();

	const std::pair<unsigned, unsigned> &get_max_crossing();


	//Score computation methods
	void computeScoreCoefficient_v();
	void computeSolutionScore();
	const double getSolutionScore();
	std::vector<double> &getScoreCoefficient_v();


	double computeCrossingBound();
	void set_crossings(double x) { crossings = x; };
	double getCrossings_afterSwap(unsigned l, double cost);
	double getCrossings_beforeOneSwap(HDAG &I, unsigned l, unsigned i, unsigned i_, double cost);
	double getCrossings_afterOneSwap(HDAG &I, unsigned l, unsigned i, unsigned i_, double cost);
	double getCrossings_beforeOneShift(HDAG &I, unsigned l, unsigned i, unsigned i_, double cost);
	double getCrossings_afterOneShift(HDAG &I, unsigned l, unsigned i, unsigned i_, double cost);


	void initialize_crossingsForLevels(unsigned l) {
		crossingsForLevels.resize(l, 0);
		bestsCrossingsForLevels.resize(l, 0);
	};

	vector<double> &get_crossingsForLevels() { return crossingsForLevels; };
	vector<double> &get_bestsCrossingsForLevels() { return bestsCrossingsForLevels; };

	/** EXPERIMENTAL methods **/
	void get_all_crossings_edge(vector<unsigned> &CCV, unsigned &maxCCV);
	/**************************/

private:

	// Number of levels of the hierarchical graph
	unsigned _l;

	// Maximum degree indeg+outdeg between all the nodes in the instance
	unsigned maxdeg;


	// LEVELS: a 3D structure where the element LEVELS[l][i][k] establishes:
	// - the l^th level of the hierarchical graph
	// - the i^th node in the l^th level of the hierarchical graph
	// - the k^th node in the FORWARD star of the i^th node in the l^th level of the hierarchical graph
	vector< vector< vector< unsigned> > > LEVELS;

	// B_LEVELS: a 3D structure where the element B_LEVELS[l][i][k] establishes:
	// - the l^th level of the hierarchical graph
	// - the i^th node in the l^th level of the hierarchical graph
	// - the k^th node in the BACKWARD star of the i^th node in the l^th level of the hierarchical graph
	vector< vector< vector< unsigned> > > B_LEVELS;


	// CurrentCostMatrix: a 3D structure of the same size of LEVELS, CurrentCostMatrix[l][i][k]
	// indicates the value of the number of crossing of the arc connecting the i^th node in the l^th
	// level with the k^th node in its FORWARD star.
	vector< vector < vector< unsigned> > > CurrentCostMatrix;


	// HowManyWCost: a vector of unordered_maps, HowManyWCost[h] is a unordered map that has as
	// keys the levels, and as value the number of forward arcs starting from each level with cost h. If there
	// is not any arc with origin in the l^th level and cost h, then the l^th level is deleted from the corresponding
	// unordered_map
	vector< unordered_map < unsigned, unsigned> > HowManyWCost;


	//AnyWCost: a vector of vector of boolean values. Given an upper bound on the cost function UB,
	//AnyWCost is a structure of size O(UB*(LEVELS.size()-1)), AnyWCost[h][l], is true if there is at
	//least an arc connecting levels l and l+1 with h crossings, being false otherwise.
	vector< vector < bool> > AnyWCost;



	//TotalCostDistr: given an upper bound on the objective function value UB, TotalCostDistr
	// is a vector of UB+1 unsigned such that TotalCostDistr[h] indicates the total number of the
	// arcs in the HDAG with h crossing.
	vector< unsigned> TotalCostDistr;

	// IDs: a 2D structure where the element IDs[l][i] establishes:
	// - the identifier of the node in LEVELS[l][i]
	vector< vector< unsigned > > IDs;

	// Os: a 2D structure where the element Os[l][i] establishes:
	// - if LEVELS[l][i] = 1 then the node is an original node
	// - if LEVELS[l][i] = 0 then the node is an incremental node.
	vector< vector< unsigned > > Os; // index of original vertices 1=original - 0=not orignal;

	// Pos: a 2D structure where given a level l and an identifier id establishes:
	// - the position of the node "id" in the level l, let's suppose that the
	//   position of "id" in l is given by "pos" then the node "id" in level l
	//   will be in the location LEVELS[l][pos]
	vector< vector< unsigned > > Pos; // position of each node;

	// DEG: a 2D structure where:
	// - let "u" be the node IDs[l][i], then DEG[l][i] = |FS(u)| + |BS(u)|
	//   according to the instance I.
	vector< vector < unsigned > > DEG;

	// vector which contains all the nodes with maximum in_deg+out_deg degree to respect
	// the instance in input
	vector<pair<unsigned, unsigned>> v_maxdeg;

	// vector which specifies how many critical arcs are present in each layer
	vector< pair<double, unsigned> > n_MAX;

	// total nodes in the hierarchical graph
	unsigned total_nodes;

 	double cost;
 	bool cost_changed;
 	double total_MAX;


 	//structures to evaluate and hold the score of a solution
 	double solutionScore;
 	std::vector<double> scoreCoefficient_v;


 	//max_crossing is a pair of unsigned, whose first component indicates
 	//the maximum number of crossing for a single arc in the HDAG, while
 	//the second component refers to the number of arcs in the HDAG with
 	//that number of crossings.
 	std::pair <unsigned, unsigned> max_crossing;

 	double crossings;
	vector<double> crossingsForLevels;
	vector<double> bestsCrossingsForLevels;

};

#endif /* HDAG_H_ */

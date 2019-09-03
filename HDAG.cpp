#include "HDAG.h"
#include <math.h>
#include "string.h"

using namespace std;

HDAG::HDAG():_l(0), maxdeg(0), total_nodes(0), cost(INT_FAST64_MAX),
		cost_changed(true), total_MAX(0), solutionScore(INT_FAST64_MAX), crossings(INT_FAST64_MAX) {
	// TODO Auto-generated constructor stub

}


void HDAG::read_instance(string filename) {

	ifstream f;
	f.open(filename.c_str());

	try {

		string line;

		getline(f, line);

		line = line + "\n";

		std::istringstream iss(line);
		iss >> _l;

		LEVELS.resize(_l);
		B_LEVELS.resize(_l);
		CurrentCostMatrix.resize(_l);
		IDs.resize(_l);
		Os.resize(_l);
		Pos.resize(_l);
		DEG.resize(_l);

		getline(f, line);

		line = line + "\n";

		iss.clear();
		iss.str(line);
		for (unsigned i = 0; i < _l; ++i) {
			unsigned x;
			iss >> x;
			LEVELS[i].resize(x);
			B_LEVELS[i].resize(x);
			CurrentCostMatrix[i].resize(x);
			IDs[i].resize(x);
			Os[i].resize(x);
			Pos[i].resize(x);
			DEG[i].resize(x, 0);
			total_nodes = total_nodes + x;
		}


		unsigned n_lev = 0;
		for (; n_lev < _l - 1; ++n_lev) {

			for (unsigned i = 0; i < LEVELS[n_lev].size(); ++i) {

				getline(f, line);

				line = line + "\n";

				iss.clear();
				iss.str(line);
				unsigned o, u, v;
				iss >> o >> u;

				IDs[n_lev][i] = u;
				Os[n_lev][i] = o;
				Pos[n_lev][u] = i;

				while (1) {

					iss >> v;

					if (iss.eof() == 1)
						break;

					LEVELS[n_lev][i].push_back(v);
					CurrentCostMatrix[n_lev][i].push_back(0);

					DEG[n_lev][i]++;

				}


				unsigned deg = compute_degree(n_lev, i);
				if (deg > maxdeg)
					maxdeg = deg;

			}

		}


		for (unsigned i = 0; i < LEVELS[n_lev].size(); ++i) {

			getline(f, line);
			iss.clear();
			iss.str(line);
			unsigned o, u;
			iss >> o >> u;
			IDs[n_lev][i] = u;
			Os[n_lev][i] = o;
			Pos[n_lev][u] = i;

		}

		double UB = computeCrossingBound()+1;
		TotalCostDistr.resize(UB);
		HowManyWCost.resize(UB);
		AnyWCost.resize(UB);
		scoreCoefficient_v.resize(UB);

		for (unsigned i = 0;i < UB; i++){
			AnyWCost[i].resize(_l);
			for (unsigned lev = 0; lev < _l; lev++){
				AnyWCost[i][lev] = false;
			}
			TotalCostDistr[i]=0;
		}


	} catch (ifstream::failure &e) {

		cerr << "Failure when reading configuration-file: " << e.what() << endl;
		exit (0);

	}



	// BUILD BACKWARD_STAR
	for (unsigned l = 0; l < LEVELS.size() - 1; ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			unsigned id = IDs[l][i];
			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
				unsigned id_ = LEVELS[l][i][j];
				unsigned pos_id_ = Pos[l + 1][id_];

				B_LEVELS[l + 1][pos_id_].push_back(id);
				DEG[l + 1][pos_id_]++;

			}
		}
	}
}


void HDAG::buildBlev() {

	B_LEVELS.resize(LEVELS.size());
	DEG.resize(LEVELS.size());


	for(unsigned l=0; l<LEVELS.size();++l){

		B_LEVELS[l].resize(LEVELS[l].size());
		DEG[l].resize(LEVELS[l].size());
	}

	// BUILD BACKWARD_STAR
	for (unsigned l = 0; l < LEVELS.size() - 1; ++l) {

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			unsigned id = IDs[l][i];

			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
				unsigned id_ = LEVELS[l][i][j];
				unsigned pos_id_ = Pos[l + 1][id_];

				B_LEVELS[l + 1][pos_id_].push_back(id);
				DEG[l + 1][pos_id_]++;

			}
		}
	}
}

HDAG::~HDAG() {
	// TODO Auto-generated destructor stub
}

void HDAG::allocate_BottleneckCostStructures(){

	unsigned ln = this->getLevNumber();


	allocateCurrentCostMatrix(ln);

	for(unsigned l = 0; l < ln; l++){
		allocateCCMatrixLevel(l,LEVELS[l].size());
		for(unsigned i = 0; i < LEVELS[l].size(); i++){
			CurrentCostMatrix[l][i].resize(LEVELS[l][i].size());
		}
	}



	double UB = computeCrossingBound();
	TotalCostDistr.resize(UB);
	HowManyWCost.resize(UB);
	AnyWCost.resize(UB);
	scoreCoefficient_v.resize(UB);

	for (unsigned i = 0;i < UB; i++){
		AnyWCost[i].resize(ln);
		for (unsigned lev = 0; lev < ln; lev++){
			AnyWCost[i][lev] = false;
		}
		TotalCostDistr[i]=0;
	}
}

void HDAG::allocateLEVELS(unsigned l) {

	this->LEVELS.resize(l);
	this->_l = l;

}

void HDAG::allocateLevel(unsigned l, unsigned n) {

	this->LEVELS[l].resize(n);

}


void HDAG::allocateCurrentCostMatrix (unsigned l){

	this->CurrentCostMatrix.resize(l);
}

void HDAG::allocateCCMatrixLevel(unsigned l, unsigned n) {

	this->CurrentCostMatrix[l].resize(n);

}

void HDAG::allocateTotalCostDistr (unsigned h){
	this->TotalCostDistr.resize(h);
}

void HDAG::allocateHowManyWCost (unsigned h){
	this->HowManyWCost.resize(h);
}

void HDAG::allocateAnyWCost (unsigned h){
	this->AnyWCost.resize(h);
}

void HDAG::allocateAnyWCostLevel(unsigned h, unsigned n){

	this->AnyWCost[h].resize(n);
}


void HDAG::allocateIDs(unsigned l) {

	this->IDs.resize(l);

}

void HDAG::allocateIDs(unsigned l, unsigned n) {

	this->IDs[l].resize(n, this->LEVELS[l].size());

}

void HDAG::allocateOs(unsigned l) {

	this->Os.resize(l);

}

void HDAG::allocateOs(unsigned l, unsigned n) {

	this->Os[l].resize(n, this->LEVELS[l].size());

}

void HDAG::allocatePos(unsigned l) {

	this->Pos.resize(l);

}

void HDAG::allocatePos(unsigned l, unsigned n) {

	this->Pos[l].resize(n, this->LEVELS[l].size());

}

void HDAG::allocateScoreCVector(unsigned h) {

	this->scoreCoefficient_v.resize(h);
}

std::vector< std::vector< std::vector< unsigned> > > &HDAG::getLEVELS() {

	return this->LEVELS;

}


std::vector< std::vector< std::vector< unsigned> > > &HDAG::getB_LEVELS() {

	return this->B_LEVELS;

}


std::vector< std::vector < std::vector< unsigned> > > &HDAG::getCurrentCostMatrix(){

	return this-> CurrentCostMatrix;
}

std::vector< std::unordered_map< unsigned, unsigned > > &HDAG::getHowManyWCost(){

	return this-> HowManyWCost;
}

std::vector< std::vector< bool> > &HDAG::getAnyWCost(){

	return this->AnyWCost;
}

std::vector< unsigned> &HDAG::getTotalCostDistr(){
	return this-> TotalCostDistr;
}

unsigned HDAG::getLevNumber() {

	return this->LEVELS.size();

}

unsigned HDAG::getMaxDegree() {

	return this->maxdeg;

}

std::vector<std::vector<unsigned> > &HDAG::getLevel(unsigned i) {

	return this->LEVELS[i];

}

std::vector<std::vector<unsigned> > &HDAG::getIDs() {

	return this->IDs;

}

std::vector<std::vector<unsigned> > &HDAG::getOs() {

	return this->Os;

}

std::vector<std::vector<unsigned> > &HDAG::getPos() {

	return this->Pos;

}

std::vector<std::vector<unsigned> > &HDAG::getDEG() {

	return this->DEG;

}

vector<std::pair<unsigned, unsigned>> &HDAG::get_v_maxdeg() {

	return this->v_maxdeg;

}


void HDAG::build_maxdegree_node_structure() {

	auto &LEVELS = this->getLEVELS();
	auto &B_LEVELS = this->getB_LEVELS();
	for (unsigned l = 0; l < LEVELS.size(); ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			unsigned in_deg = LEVELS[l][i].size();
			unsigned out_deg = B_LEVELS[l][i].size();
			unsigned tot_deg = in_deg + out_deg;
			if (tot_deg > maxdeg)
				maxdeg = tot_deg;
		}
	}


	for (unsigned l = 0; l < LEVELS.size(); ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			unsigned in_deg = LEVELS[l][i].size();
			unsigned out_deg = B_LEVELS[l][i].size();
			unsigned tot_deg = in_deg + out_deg;

			if (tot_deg == maxdeg) {

				pair<unsigned, unsigned> p;
				p.first = l;
				p.second = i;
				v_maxdeg.push_back(p);

			}
		}
	}


	return;

}


bool HDAG::areCrossingEdge(unsigned i, unsigned j, unsigned i_, unsigned j_,
		unsigned l) {


	if((i==i_) || (j==j_)) return false;


	if ((this->Pos[l][i] < this->Pos[l][i_]) && (this->Pos[l + 1][j] > this->Pos[l + 1][j_]))
		return true;

	if ((this->Pos[l][i] > this->Pos[l][i_]) && (this->Pos[l + 1][j] < this->Pos[l + 1][j_]))
		return true;

	return false;

}

bool HDAG::isOriginalLink(unsigned i, unsigned j, unsigned l) {

	if (this->Os[l][this->Pos[l][i]] == 1 && this->Os[l + 1][this->Pos[l + 1][j]] == 1)
		return true;

	return false;

}

bool HDAG::isOriginalNode(unsigned i, unsigned l) {

	if (this->Os[l][i] == 1)
		return true;

	return false;

}

bool HDAG::isIncrementalNode(unsigned i, unsigned l) {

	if (!isOriginalNode(i, l))
		return true;

	return false;

}

void HDAG::swapPositions(unsigned p1, unsigned p2, unsigned l) {


	this->Pos[l][this->IDs[l][p1]] = p2;
	this->Pos[l][this->IDs[l][p2]] = p1;


	this->LEVELS[l][p1].swap(this->LEVELS[l][p2]);
	this->CurrentCostMatrix[l][p1].swap(this->CurrentCostMatrix[l][p2]);



	unsigned tmp = this->Os[l][p1];
	this->Os[l][p1] = this->Os[l][p2];
	this->Os[l][p2] = tmp;



	tmp = this->IDs[l][p1];
	this->IDs[l][p1] = this->IDs[l][p2];
	this->IDs[l][p2] = tmp;



	return;

}



unsigned HDAG::compute_degree(unsigned l, unsigned i) {

	auto &LEVELS = this->getLEVELS();
	unsigned deg = 0;

	if (l > 0) {

		for (unsigned i_ = 0; i_ < LEVELS[l - 1].size(); ++i_) {
			for (unsigned j_ = 0; j_ < LEVELS[l - 1][i_].size(); ++j_) {
				if 	(LEVELS[l - 1][i_][j_] == IDs[l][i]) {
					deg++;
					break;
				}
			}
		}
	}

	deg += LEVELS[l][i].size();

	return deg;

}

unsigned HDAG::compute_barycenter(vector<vector<unsigned>> &CL, HDAG &I, unsigned l, unsigned i) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	auto &Pos_S = this->getPos();

	double BS_bc = 0;
	double FS_bc = 0;
	unsigned BS_elem = 0;
	unsigned FS_elem = 0;

	// compute the barycenter in backward sense
	for (unsigned j = 0; j < B_LEVELS[l][i].size(); ++j) {
		unsigned id_j = B_LEVELS[l][i][j];
		unsigned original_pos_id_j = Pos[l - 1][id_j];
		// if the node is already in the new solution
		if (CL[l - 1][original_pos_id_j] == 0) {
			unsigned new_pos_id_j = Pos_S[l - 1][id_j];
			BS_bc = BS_bc + new_pos_id_j;
			BS_elem++;
		}
	}

	// compute the barycenter in forward sense
	for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
		unsigned id_j = LEVELS[l][i][j];
		unsigned original_pos_id_j = Pos[l + 1][id_j];
		// if the node is already in the new solution
		if (CL[l + 1][original_pos_id_j] == 0) {
			unsigned new_pos_id_j = Pos_S[l + 1][id_j];
			FS_bc = FS_bc + new_pos_id_j;
			FS_elem++;
		}
	}



	if (BS_elem != 0)
		BS_bc = BS_bc / BS_elem;
	if (FS_elem != 0)
		FS_bc = FS_bc / FS_elem;

	if (BS_elem == 0 && FS_elem == 0) return 0;


	//rounding to the nearest integer
	unsigned bc = (round)((BS_bc + FS_bc) / 2);


	return bc;
}

unsigned HDAG::get_n_original_nodes(HDAG &G) {

	auto &LEVELS = G.getLEVELS();
	unsigned n = 0;

	for (unsigned l = 0; l < LEVELS.size(); ++l)
		for (unsigned i = 0; i < LEVELS[l].size(); ++i)
			if (G.isOriginalNode(i, l))
				n++;

	return n;

}

void HDAG::printHDAG(HDAG &G, string filename) {

	ofstream f(filename.c_str());

	f << "\\documentclass{standalone}" << '\n';
	f << "\\usepackage{tikz}" << '\n';
	f << "\\usepackage{amsmath}" << '\n';
	f << "\\usepackage{amsfonts}" << '\n';
	f << "\\usepackage{comment}" << '\n';
	f << "\\begin{document}" << '\n';
	f << "\\begin{tikzpicture}[thick]" << '\n';

	// Print the nodes
	unsigned index_l = 0;

	unsigned max_dim_level = 0;


	f << "% VERTICES" << '\n';

	for (unsigned i = 0; i < G.getLevNumber(); ++i) {

		int index_i = 0;

		f << "% Node of level " << i << '\n';

		if (G.getLevel(i).size() > max_dim_level)
			max_dim_level = G.getLevel(i).size();

		for (unsigned j = 0; j < G.getLevel(i).size(); ++j) {

			if (G.Os[i][j] == 1) {

				if (G.IDs[i][j] < 10) {

					f << "\\node [draw, circle, inner sep = 0.85] (" << G.IDs[i][j] << "-" << i
						<< ") at (" << index_l << ", " << index_i << ") {\\tiny "
						<< G.IDs[i][j] << "};" << '\n';

				}
				else {

					f << "\\node [draw, circle, inner sep = 0] (" << G.IDs[i][j] << "-" << i
						<< ") at (" << index_l << ", " << index_i << ") {\\tiny "
						<< G.IDs[i][j] << "};" << '\n';

				}



			} else {

				if (G.IDs[i][j] == LEVELS[i].size()) {

					f << "\\node [draw, circle, dashed, inner sep = 0] (" << G.IDs[i][j] << "-"
							<< i << ") at (" << index_l << ", " << index_i << ") {\\tiny x };" << '\n';

				}
				else {
					f << "\\node [draw, circle, dashed, inner sep = 0] (" << G.IDs[i][j] << "-"
						<< i << ") at (" << index_l << ", " << index_i << ") {\\tiny "
						<< G.IDs[i][j] << "};" << '\n';
				}
			}

			index_i++;

		}

		f << '\n';

		index_l = index_l + 3;

	}

	f << "% EDGES" << '\n';

	//Print edges
	for (unsigned i = 0; i < G.getLevNumber() - 1; ++i) {

		f << "% Connection between levels " << i << " and " << i + 1 << '\n';

		for (unsigned j = 0; j < G.getLevel(i).size(); ++j) {

			for (unsigned k = 0; k < G.getLevel(i)[j].size(); ++k) {

				// Find link (u, v)
				unsigned u = G.IDs[i][j];
				unsigned v = G.LEVELS[i][j][k];

				if (isOriginalLink(u, v, i)) {

					f << "\\draw (" << G.IDs[i][j] << "-" << i << ") -- ("
							<< G.getLevel(i)[j][k] << "-" << i + 1 << ");"
							<< '\n';

				} else {

					f << "\\draw[dashed] (" << G.IDs[i][j] << "-" << i
							<< ") -- (" << G.getLevel(i)[j][k] << "-" << i + 1
							<< ");" << '\n';

				}

			}

		}

	}


	f << "\\end{tikzpicture}" << '\n';
	f << "\\end{document}" << '\n';

	f.close();

	return;

}

void HDAG::copy_HDAG(HDAG &G) {

	//TODO: qui bisogna copiare anche le strutture bottleneck?
	auto &LEVELS = G.getLEVELS();
	auto &IDs = G.getIDs();
	auto &Os = G.getOs();
	auto &CurrentCostMatrix = G.getCurrentCostMatrix();

	auto &IDs_new = this->getIDs();
	auto &Pos_new = this->getPos();
	auto &Os_new = this->getOs();

	unsigned ln, vn;

	ln = G.getLevNumber();

	this->allocateLEVELS(ln);
	this->allocateIDs(ln);
	this->allocateOs(ln);
	this->allocatePos(ln);
	this->allocateCurrentCostMatrix(ln);

	for (unsigned l = 0; l < ln; ++l) {

		vn = LEVELS[l].size();

		this->allocateLevel(l, vn);
		this->allocateIDs(l, vn);
		this->allocateOs(l, vn);
		this->allocatePos(l, vn);
		this->allocateCCMatrixLevel(l, vn);


		for (unsigned i = 0; i < vn; ++i) {

			unsigned id_u = IDs[l][i];
			IDs_new[l][i] = id_u;
			Pos_new[l][id_u] = i;
			Os_new[l][i] = Os[l][i];

			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

				this->getLEVELS()[l][i].push_back(LEVELS[l][i][j]);
				this->getCurrentCostMatrix()[l][i].push_back(CurrentCostMatrix[l][i][j]);

			}

		}

	}



	auto &TotalCostDistribution = G.getTotalCostDistr();
	auto &HMWCost = G.getHowManyWCost();
	auto &AnyWCost = G.getAnyWCost();

	auto UB = TotalCostDistribution.size();

	auto &TotalCostDistr_new = this->TotalCostDistr;
	auto &HMWCost_new = this->HowManyWCost;
	auto &AnyWCost_new = this->AnyWCost;

	this->allocateTotalCostDistr(UB);
	this->allocateAnyWCost(UB);
	this->allocateHowManyWCost(UB);

	for(unsigned h = 0; h < UB ; h++){
		TotalCostDistr_new[h] = TotalCostDistribution[h];
		this->allocateAnyWCostLevel(h,ln);
		HMWCost_new[h]=HMWCost[h];

		for(unsigned l = 0; l < ln; l++){
			AnyWCost_new[h][l]=AnyWCost[h][l];

		}


	}

	return;

}

void HDAG::copy_level(HDAG &G, unsigned l) {

	auto &LEVELS_target = G.getLEVELS();
	auto &IDs_target = G.getIDs();

	auto &Pos_source = this->getPos();

	for (unsigned i = 0; i < LEVELS_target[l].size(); ++i) {

		unsigned id_v_in_target = IDs_target[l][i];
		unsigned pos_v_in_target = i;

		unsigned id_v_in_source = id_v_in_target;
		unsigned pos_v_in_source = Pos_source[l][id_v_in_source];

		if (pos_v_in_target != pos_v_in_source) {

			this->swapPositions(pos_v_in_source, pos_v_in_target, l);

		}

	}

	return;

}

void HDAG::dealloc_HDAG() {

	for (unsigned l = 0; l < this->getLEVELS().size(); ++l) {
		for (unsigned i = 0; i < this->getLEVELS()[l].size(); ++i) {
			while (this->getLEVELS()[l][i].size() > 0)
				this->getLEVELS()[l][i].pop_back();
			delete &this->LEVELS[l][i];
		}
		delete &this->LEVELS[l];
	}
	delete this;

}

void HDAG::get_all_crossings_edge(vector<unsigned> &CCV, unsigned &maxCCV) {

	double new_cost = 0;

	for (unsigned l = 0; l < LEVELS.size() - 1; ++l) {

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

				double actual_cost = 0;

				for (unsigned i_ = 0; i_ < LEVELS[l].size(); ++i_) {

					for (unsigned j_ = 0; j_ < LEVELS[l][i_].size(); ++j_) {

						// arc (i, j)
						unsigned id_i = IDs[l][i];
						unsigned id_j = LEVELS[l][i][j];

						// arc (i_, j_)
						unsigned id_i_ = IDs[l][i_];
						unsigned id_j_ = LEVELS[l][i_][j_];

						if (areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {

							actual_cost++;

						}

					}

				}

				CCV[actual_cost]++;

				if (actual_cost > new_cost) {
					new_cost = actual_cost;
				}

			}

		}

	}

	cost = new_cost;
	maxCCV = cost;
	return;

}

double HDAG::get_cost() {

	if (!this->is_cost_changed()) return this->cost;

	double new_cost = 0;
	total_MAX = 0;

	for (unsigned l = 0; l < LEVELS.size() - 1; ++l) {

		double MAX_in_level = 0;
		unsigned n_MAX_in_level = 0;

		n_MAX[l + 1].first = 0;
		n_MAX[l + 1].second = 0;

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

				double actual_cost = 0;

				for (unsigned i_ = 0; i_ < LEVELS[l].size(); ++i_) {

					for (unsigned j_ = 0; j_ < LEVELS[l][i_].size(); ++j_) {

						// arc (i, j)
						unsigned id_i = IDs[l][i];
						unsigned id_j = LEVELS[l][i][j];

						// arc (i_, j_)
						unsigned id_i_ = IDs[l][i_];
						unsigned id_j_ = LEVELS[l][i_][j_];

						if (areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {

							actual_cost++;

						}

					}

				}

				if (actual_cost > new_cost) {
					new_cost = actual_cost;
				}

				if (MAX_in_level < actual_cost) {
					MAX_in_level = actual_cost;
					n_MAX_in_level = 1;
				}
				else if (MAX_in_level == actual_cost) {
					n_MAX_in_level++;
				}

			}

		}

		n_MAX[l + 1].first = MAX_in_level;
		n_MAX[l + 1].second = n_MAX_in_level;

	}

	cost_changed = false;
	cost = new_cost;

	total_MAX = 0;
	for (unsigned l = 1; l < LEVELS.size() + 1; ++l)
		if (n_MAX[l].first == cost)
			total_MAX = total_MAX + n_MAX[l].second;



	return cost;

}


double HDAG::get_cost_bottleneck() {

	if (!this->is_cost_changed()) return this->cost;

	double new_cost = 0;
	total_MAX = 0;

	for (unsigned l = 0; l < LEVELS.size() - 1; ++l) {

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

				double actual_cost = 0;

				for (unsigned i_ = 0; i_ < LEVELS[l].size(); ++i_) {

					for (unsigned j_ = 0; j_ < LEVELS[l][i_].size(); ++j_) {

						// arc (i, j)
						unsigned id_i = IDs[l][i];
						unsigned id_j = LEVELS[l][i][j];

						// arc (i_, j_)
						unsigned id_i_ = IDs[l][i_];
						unsigned id_j_ = LEVELS[l][i_][j_];

						if (areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {

							actual_cost++;

						}

					}

				}

				CurrentCostMatrix[l][i][j]=actual_cost;
				TotalCostDistr[actual_cost]++;

				if(AnyWCost[actual_cost][l]==true){
					//update map
					HowManyWCost[actual_cost][l]++;
				}

				else{
					//insert level in map, set AnyWCost to true
					HowManyWCost[actual_cost].insert({l,1});
					AnyWCost[actual_cost][l] = true;
				}


				if (actual_cost > new_cost) {
					new_cost = actual_cost;
				}

			}

		}

	}

	cost_changed = false;
	cost = new_cost;


	return cost;

}

double HDAG::update_cost_bottleneck(){

	unsigned i = TotalCostDistr.size();
	bool found = false;

	while(found == false && i > 0){
		i--;

		if(TotalCostDistr[i]>0){
			found = true;
		}
	}
	return i;
}

double HDAG::compute_crossing_sum_bottleneck(){
	double total_crossing_sum = 0;

		for(int h = 1; h <= this->update_cost_bottleneck(); h++){
			total_crossing_sum = total_crossing_sum + (h*TotalCostDistr[h]);
		}

		total_crossing_sum = total_crossing_sum/2;
		return total_crossing_sum;
}

void HDAG::reset_bottleneck_structures(){

	for (auto &x : CurrentCostMatrix){
		for(auto &y : x){
			std::fill(y.begin(), y.end(), 0);
		}
	}

	for(unsigned costAddress = 0 ; costAddress < TotalCostDistr.size(); costAddress++){

		TotalCostDistr[costAddress] = 0;

		std::fill(AnyWCost[costAddress].begin(), AnyWCost[costAddress].end(),false);

		HowManyWCost[costAddress].clear();
	}

	this -> set_cost_changed(true);

}


unsigned HDAG::get_MaxCrossing_incidentTo(unsigned l, unsigned i){

	unsigned MaxCrossingIncident = 0;
	const unsigned id_i = IDs[l][i];



	if(l < this->getLevNumber()){
		for (unsigned j = 0; j < CurrentCostMatrix[l][i].size(); j++){
			if(CurrentCostMatrix[l][i][j] > MaxCrossingIncident){
				MaxCrossingIncident = CurrentCostMatrix[l][i][j];
			}
		}
	}




	if(l > 0){
		unsigned id_k = 0;
		unsigned pos_k = 0 ;

		//to hold the position of i in  the outgoing star of node k of level l-1
		unsigned relative_pos_i = 0;



		for (unsigned k = 0; k < this->getB_LEVELS()[l][id_i].size(); k++){


			id_k = getB_LEVELS()[l][id_i][k];
			pos_k = Pos[l-1][id_k];


			for(unsigned j_= 0; j_ < this->getLEVELS()[l-1][pos_k].size();j_++){


				if(LEVELS[l-1][pos_k][j_]== id_i){
					relative_pos_i = j_;
					break;
				}
			}


			if(CurrentCostMatrix[l-1][pos_k][relative_pos_i] > MaxCrossingIncident){
				MaxCrossingIncident = CurrentCostMatrix[l-1][pos_k][relative_pos_i];
			}
		}

	}


	return MaxCrossingIncident;

}


void HDAG::recompute_n_MAX_in_level (unsigned l) {

	if (n_MAX[l + 1].first == cost)
		total_MAX -= n_MAX[l + 1].second;

	n_MAX[l + 1].first = 0;
	n_MAX[l + 1].second = 0;

	double MAX_in_level = 0;
	unsigned n_MAX_in_level = 0;

	for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

		for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

			double actual_cost = 0;

			for (unsigned i_ = 0; i_ < LEVELS[l].size(); ++i_) {

				for (unsigned j_ = 0; j_ < LEVELS[l][i_].size(); ++j_) {

					// arc (i, j)
					unsigned id_i = IDs[l][i];
					unsigned id_j = LEVELS[l][i][j];

					// arc (i_, j_)
					unsigned id_i_ = IDs[l][i_];
					unsigned id_j_ = LEVELS[l][i_][j_];

					if (areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {

						actual_cost++;

					}

				}

			}

			if (MAX_in_level < actual_cost) {
				MAX_in_level = actual_cost;
				n_MAX_in_level = 1;
			}
			else if (MAX_in_level == actual_cost) {
				n_MAX_in_level++;
			}

		}

	}

	n_MAX[l + 1].first = MAX_in_level;
	n_MAX[l + 1].second = n_MAX_in_level;

	if (n_MAX[l + 1].first == cost)
		total_MAX += n_MAX[l + 1].second;
	else if (n_MAX[l + 1].first > cost) {
		total_MAX = n_MAX[l + 1].second;
	}

	return;

}

double HDAG::get_crossings() {

	const auto &LEVELS = this->getLEVELS();
	const auto &IDs = this->getIDs();

	double new_crossings = 0;

	for (unsigned l = 0; l < LEVELS.size() - 1; ++l) {

		for (unsigned i = 0; i < LEVELS[l].size() - 1; ++i) {

			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

				for (unsigned i_ = i + 1; i_ < LEVELS[l].size(); ++i_) {

					for (unsigned j_ = 0; j_ < LEVELS[l][i_].size(); ++j_) {

						// arc (i, j)
						unsigned id_i = IDs[l][i];
						unsigned id_j = LEVELS[l][i][j];

						// arc (i_, j_)
						unsigned id_i_ = IDs[l][i_];
						unsigned id_j_ = LEVELS[l][i_][j_];

						if (this->areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {

							new_crossings++;

						}

					}

				}

			}

		}

	}

	crossings = new_crossings;
	return crossings;

}

double HDAG::computeCrossingBound(){

	auto LEVELS = this->getLEVELS();

		unsigned maxedge_for_layer = 0;
		for (unsigned l = 0; l < LEVELS.size(); ++l) {
			unsigned actual_number_edge = 0;
			for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
				actual_number_edge += LEVELS[l][i].size();
			}
			if (actual_number_edge > maxedge_for_layer)
				maxedge_for_layer = actual_number_edge;

		}

		//maxedge_for_layer--;
		return maxedge_for_layer;
}


void HDAG::update_max_crossing(){

	int i = TotalCostDistr.size()-1;

	while(true){

		if(TotalCostDistr[i] > 0){
			max_crossing.first = i;
			max_crossing.second = TotalCostDistr[i];
			break;
		}

		i++;
	}

}

const std::pair <unsigned, unsigned> &HDAG::get_max_crossing(){

	return this-> max_crossing;
}

const double HDAG::getSolutionScore(){

	return this->solutionScore;

}

void HDAG::computeScoreCoefficient_v(){

	double product = 1;
	double factor = computeCrossingBound()+1;


	for(unsigned i = 0; i < scoreCoefficient_v.size(); i++){

		factor +=i;
		product	= product * factor;
		scoreCoefficient_v[i] = product;
	}

}

void HDAG::computeSolutionScore(){

	double currentMaxCross = update_cost_bottleneck();

	solutionScore = currentMaxCross;

	for(unsigned i = 0; i <= currentMaxCross; i++){
		solutionScore += double(TotalCostDistr[currentMaxCross - i])
				/ (scoreCoefficient_v[i]);
	}

}

std::vector<double> &HDAG::getScoreCoefficient_v(){

	return this->scoreCoefficient_v;
}

double HDAG::getCrossings_afterSwap(unsigned l, double cost) {

	double new_cost = cost;

	if (l > 0) {

		l = l - 1;
 	    new_cost = new_cost - crossingsForLevels[l];
 	   crossingsForLevels[l] = 0;
		for (unsigned i = 0; i < LEVELS[l].size() - 1; ++i) {
			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
				for (unsigned i_ = i + 1; i_ < LEVELS[l].size(); ++i_) {
					for (unsigned j_ = 0; j_ < LEVELS[l][i_].size(); ++j_) {
						// arc (i, j)
						unsigned id_i = IDs[l][i];
						unsigned id_j = LEVELS[l][i][j];
						// arc (i_, j_)
						unsigned id_i_ = IDs[l][i_];
						unsigned id_j_ = LEVELS[l][i_][j_];
						if (areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {
							crossingsForLevels[l]++;
						}
					}
				}
			}
		}
		new_cost = new_cost + crossingsForLevels[l];
		l = l + 1;

	}

	if (l < LEVELS.size() - 1) {

	 	new_cost = new_cost - crossingsForLevels[l];
	 	crossingsForLevels[l] = 0;
		for (unsigned i = 0; i < LEVELS[l].size() - 1; ++i) {
			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
				for (unsigned i_ = i + 1; i_ < LEVELS[l].size(); ++i_) {
					for (unsigned j_ = 0; j_ < LEVELS[l][i_].size(); ++j_) {
						// arc (i, j)
						unsigned id_i = IDs[l][i];
						unsigned id_j = LEVELS[l][i][j];
						// arc (i_, j_)
						unsigned id_i_ = IDs[l][i_];
						unsigned id_j_ = LEVELS[l][i_][j_];
						if (areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {
							crossingsForLevels[l]++;
						}
					}
				}
			}
		}
		new_cost = new_cost + crossingsForLevels[l];
	}

	return new_cost;

}


void HDAG::update_bottleneck_str_afterSwap(unsigned p1, unsigned p2, unsigned l){

	//establishing ranges of positions affected
	//by the swap
	unsigned min_affected = min(p1,p2);
	unsigned max_affected = max(p1,p2);


	//if l is not the last layer, then the outgoing arcs
	//have to be updated
	if(l < this->getLevNumber()){

		for (unsigned affected = min_affected; affected <= max_affected;
				affected++) {
			for(unsigned j = 0; j <this->LEVELS[l][affected].size(); j++){

				unsigned previousCost = CurrentCostMatrix[l][affected][j];


				unsigned currentArcCost = 0;

				//subtract one from TotalCostDistr istogram in the
				//position relative to previous cost of the arc
				this->TotalCostDistr[previousCost]--;

				//update the number of arcs starting from l
				//with cost previousCost
				this->HowManyWCost[previousCost][l]--;

				if(this->HowManyWCost[previousCost][l] == 0){
					this->AnyWCost[previousCost][l] = false;
					this->HowManyWCost[previousCost].erase(l);
				}


				//compute new cost for the arc [affected, j]

				for(unsigned i_ = 0; i_ <LEVELS[l].size(); i_++){
					for(unsigned j_ = 0; j_ < LEVELS[l][i_].size(); j_++){


						//spot the arc(id_affected, id_j)
						unsigned id_affected = IDs[l][affected];
						unsigned id_j = LEVELS[l][affected][j];


						//spot the arc (id_i, id_j_)
						unsigned id_i_= IDs[l][i_];
						unsigned id_j_ = LEVELS[l][i_][j_];

						if(areCrossingEdge(id_affected,id_j,id_i_,id_j_,l)){
							currentArcCost++;
								}
						}
				}


				CurrentCostMatrix[l][affected][j] = currentArcCost;
				TotalCostDistr[currentArcCost]++;


				if(AnyWCost[currentArcCost][l] == false){
					AnyWCost[currentArcCost][l] = true;
					HowManyWCost[currentArcCost].insert({l,1});
				}
				else{
					HowManyWCost[currentArcCost][l]++;
				}

			}
		}
	}


	//if l is not the first layer, then the ingoing arcs
	//have to be updated
	if (l > 0){
		l--;
		for(unsigned pos = 0; pos < LEVELS[l].size(); pos++){
			for(unsigned j = 0; j < LEVELS[l][pos].size(); j++){

				unsigned previousCost = CurrentCostMatrix[l][pos][j];

				unsigned currentArcCost = 0;



				//subtract one from TotalCostDistr istogram in the
				//position relative to previous cost of the arc
				this->TotalCostDistr[previousCost]--;

				//update the number of arcs starting from l
				//with cost previousCost
				this->HowManyWCost[previousCost][l]--;

				if(this->HowManyWCost[previousCost][l] == 0){
					this->AnyWCost[previousCost][l] = false;
					this->HowManyWCost[previousCost].erase(l);
				}



				//compute new cost for arc (id_pos, id_j)

				for(unsigned i_ = 0; i_ < LEVELS[l].size(); i_++){
					for(unsigned j_= 0; j_ < LEVELS[l][i_].size(); j_++){

						//spot the arc (id_pos, id_j)
						unsigned id_pos = IDs[l][pos];
						unsigned id_j = LEVELS[l][pos][j];

						//spot the arc (id_i_, id_j_)
						unsigned id_i_ = IDs[l][i_];
						unsigned id_j_ = LEVELS[l][i_][j_];

						if(areCrossingEdge(id_pos, id_j,id_i_,id_j_,l)){
							currentArcCost++;
						}

					}
				}

				CurrentCostMatrix[l][pos][j] = currentArcCost;
				TotalCostDistr[currentArcCost]++;

				if(AnyWCost[currentArcCost][l] == false){
					AnyWCost[currentArcCost][l] = true;
					HowManyWCost[currentArcCost].insert({l,1});
				}
				else{
					HowManyWCost[currentArcCost][l]++;
				}
			}
		}
	}
}

double HDAG::getCrossings_beforeOneSwap(HDAG &I, unsigned l, unsigned i, unsigned i_, double cost) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	auto &LEVELS_S = this->getLEVELS();
	auto &IDs_S = this->getIDs();

	unsigned x_min = min(i, i_);
	unsigned x_max = max(i, i_);

	if (l > 0) {
		for (unsigned x = x_min; x <= x_max - 1; ++x) {
			unsigned id_x = IDs_S[l][x];
			unsigned original_pos_x = Pos[l][id_x];
			for (unsigned j = 0; j < B_LEVELS[l][original_pos_x].size(); ++j) {
				unsigned id_j = B_LEVELS[l][original_pos_x][j];
				for (unsigned x_ = x + 1; x_ <= x_max; ++x_) {
					unsigned id_x_ = IDs_S[l][x_];
					unsigned original_pos_x_ = Pos[l][id_x_];
					for (unsigned j_ = 0;
							j_ < B_LEVELS[l][original_pos_x_].size(); ++j_) {
						unsigned id_j_ = B_LEVELS[l][original_pos_x_][j_];
						if (id_j != id_j_) {
							if (this->areCrossingEdge(id_j, id_x, id_j_, id_x_,
									l - 1)) {

								cost--;
							}
						}
					}
				}
			}
		}
	}
	if (l < LEVELS.size() - 1) {
		for (unsigned x = x_min; x <= x_max - 1; ++x) {
			unsigned id_x = IDs_S[l][x];
			for (unsigned j = 0; j < LEVELS_S[l][x].size(); ++j) {
				// arc (x,j)
				unsigned id_j = LEVELS_S[l][x][j];
				for (unsigned x_ = x + 1; x_ <= x_max; ++x_) {
					unsigned id_x_ = IDs_S[l][x_];
					for (unsigned j_ = 0; j_ < LEVELS_S[l][x_].size(); ++j_) {
						// arc (x_, j_)
						unsigned id_j_ = LEVELS_S[l][x_][j_];
						if (id_j != id_j_) {

							if (this->areCrossingEdge(id_x, id_j, id_x_, id_j_, l)) {

								cost--;
							}
						}
					}
				}
			}
		}
	}

	return cost;

}

double HDAG::getCrossings_afterOneSwap(HDAG &I, unsigned l, unsigned i, unsigned i_, double cost) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	auto &LEVELS_S = this->getLEVELS();
	auto &IDs_S = this->getIDs();

	unsigned x_min = min(i, i_);
	unsigned x_max = max(i, i_);

	if (l > 0) {
		for (unsigned x = x_min; x <= x_max - 1; ++x) {
			unsigned id_x = IDs_S[l][x];
			unsigned original_pos_x = Pos[l][id_x];
			for (unsigned j = 0; j < B_LEVELS[l][original_pos_x].size(); ++j) {
				unsigned id_j = B_LEVELS[l][original_pos_x][j];
				for (unsigned x_ = x + 1; x_ <= x_max; ++x_) {
					unsigned id_x_ = IDs_S[l][x_];
					unsigned original_pos_x_ = Pos[l][id_x_];
					for (unsigned j_ = 0;
							j_ < B_LEVELS[l][original_pos_x_].size(); ++j_) {
						unsigned id_j_ = B_LEVELS[l][original_pos_x_][j_];
						if (id_j != id_j_) {
							if (this->areCrossingEdge(id_j, id_x, id_j_, id_x_,
									l - 1)) {
								cost++;
							}
						}
					}
				}
			}
		}
	}
	if (l < LEVELS.size() - 1) {
		for (unsigned x = x_min; x <= x_max - 1; ++x) {
			unsigned id_x = IDs_S[l][x];
			for (unsigned j = 0; j < LEVELS_S[l][x].size(); ++j) {
				// arc (x,j)
				unsigned id_j = LEVELS_S[l][x][j];
				for (unsigned x_ = x + 1; x_ <= x_max; ++x_) {
					unsigned id_x_ = IDs_S[l][x_];
					for (unsigned j_ = 0; j_ < LEVELS_S[l][x_].size(); ++j_) {
						// arc (x_, j_)
						unsigned id_j_ = LEVELS_S[l][x_][j_];
						if (id_j != id_j_) {
							if (this->areCrossingEdge(id_x, id_j, id_x_, id_j_, l)) {

								cost++;
							}
						}
					}
				}
			}
		}
	}
	return cost;

}

double HDAG::getCrossings_beforeOneShift(HDAG &I, unsigned l,
		unsigned i, unsigned i_, double cost) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	auto &LEVELS_S = this->getLEVELS();
	auto &IDs_S = this->getIDs();

	if (l > 0) {
		unsigned id_i = IDs_S[l][i];
		unsigned id_i_ = IDs_S[l][i_];
		unsigned original_pos_i = Pos[l][id_i];
		unsigned original_pos_i_ = Pos[l][id_i_];
		for (unsigned j = 0; j < B_LEVELS[l][original_pos_i].size(); ++j) {
			for (unsigned j_ = 0; j_ < B_LEVELS[l][original_pos_i_].size();
					++j_) {
				unsigned id_j = B_LEVELS[l][original_pos_i][j];
				unsigned id_j_ = B_LEVELS[l][original_pos_i_][j_];

				if (id_j != id_j_) {
					if (this->areCrossingEdge(id_j, id_i, id_j_, id_i_, l - 1)) {

						cost = cost - 1;

					}
				}
			}
		}
	}
	if (l < LEVELS.size() - 1) {

		for (unsigned j = 0; j < LEVELS_S[l][i].size(); ++j) {
			for (unsigned j_ = 0; j_ < LEVELS_S[l][i_].size(); ++j_) {
				// arc (i,j)
				unsigned id_i = IDs_S[l][i];
				unsigned id_j = LEVELS_S[l][i][j];

				// arc (i_,j_)
				unsigned id_i_ = IDs_S[l][i_];
				unsigned id_j_ = LEVELS_S[l][i_][j_];

				if (id_j != id_j_) {
					if (this->areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {
						cost = cost - 1;

					}
				}
			}
		}
	}

	return cost;

}

double HDAG::getCrossings_afterOneShift(HDAG &I, unsigned l,
		unsigned i, unsigned i_, double cost) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	auto &LEVELS_S = this->getLEVELS();
	auto &IDs_S = this->getIDs();

	if (l > 0) {
		unsigned id_i = IDs_S[l][i];
		unsigned id_i_ = IDs_S[l][i_];
		unsigned original_pos_i = Pos[l][id_i];
		unsigned original_pos_i_ = Pos[l][id_i_];

		for (unsigned j = 0; j < B_LEVELS[l][original_pos_i].size(); ++j) {
			for (unsigned j_ = 0; j_ < B_LEVELS[l][original_pos_i_].size();
					++j_) {
				unsigned id_j = B_LEVELS[l][original_pos_i][j];
				unsigned id_j_ = B_LEVELS[l][original_pos_i_][j_];

				if (id_j != id_j_) {
					if (this->areCrossingEdge(id_j, id_i, id_j_, id_i_, l - 1)) {
						cost = cost + 1;

					}
				}
			}
		}
	}
	if (l < LEVELS.size() - 1) {

		for (unsigned j = 0; j < LEVELS_S[l][i].size(); ++j) {
			for (unsigned j_ = 0; j_ < LEVELS_S[l][i_].size(); ++j_) {
				// arc (i,j)
				unsigned id_i = IDs_S[l][i];
				unsigned id_j = LEVELS_S[l][i][j];

				// arc (i_,j_)
				unsigned id_i_ = IDs_S[l][i_];
				unsigned id_j_ = LEVELS_S[l][i_][j_];

				if (id_j != id_j_) {
					if (this->areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {
						cost = cost + 1;

					}
				}
			}
		}
	}

	return cost;

}

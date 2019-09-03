#include "TABU.h"
#include "string.h"

TABU::TABU(HDAG &I, char *argv[]) :
		currentIter(0), r_seed(0), MAX_ITER(0), Max_no_improv(
				0), MAX_DIVERSE(0), last_improv(0), tenure(0), alpha(0) {


	set_tenure(atoi(argv[3]));
	set_currentIter(tenure + 1);
	set_alpha(atof(argv[4]));
	set_post_proc(argv[5]);
	set_MAXITER(atoi(argv[6])+tenure);
	set_start_algorithm(chrono::system_clock::now());
	set_best_out_file(argv[8]);
	set_complete_out_file(argv[9]);
	set_TimeLimit(atoi(argv[10]));
	set_Max_no_improv(atoi(argv[11]));
	set_MAX_DIVERSE(atoi(argv[12]));
	set_r_seed(atoi(argv[13]));

	allocate_all_structures(I);

	I.get_cost_bottleneck();
	I.computeScoreCoefficient_v();

	this->run_TABU(I);


	if(strcmp(post_proc, "yes") == 0){
		std::cout<<"Starting post-processing..."<<std::endl;
		reduce_total_crossing(best_solution);
		std::cout<<"The max-cross is: "<<best_solution.update_cost_bottleneck()<<std::endl;
	}

}

void TABU::run_max_TABU(HDAG &I){

	best_solution.copy_HDAG(I);
	best_solution.set_cost(INT64_MAX);
	HDAG &S = I;

	while(currentIter < MAX_ITER){

		for(unsigned l = 0; l < TabuTrack.size(); l++){
			compute_maxNonTabu_inLev(S,l);
			fill_CandidateList(S,l);
			}
		max_layer_TABU(S);
		currentIter++;


		int currentTimer = chrono::duration_cast<chrono::milliseconds>(
						chrono::system_clock::now() - start_algorithm).count();

				if((double)((double)currentTimer / (double)1000) >= TIME_LIMIT){
					break;
				}
	}

	set_end_iteration(chrono::system_clock::now());
	print_best_data_out_file();
	close_out_file();

	std::cout<<best_solution.update_cost_bottleneck()<<std::endl;

}

void TABU::run_TABU(HDAG &I){

	best_solution.copy_HDAG(I);
	best_solution.set_cost(INT64_MAX);
	HDAG &S = I;

	unsigned random_seed = r_seed;
	std::default_random_engine globalGenerator(random_seed);
	std::uniform_int_distribution<int> distribution(0,INT32_MAX);

	while(currentIter < MAX_ITER){

		last_improv++;

		for(unsigned l = 0; l < TabuTrack.size(); l++){

			compute_maxNonTabu_inLev(S,l);
			fill_CandidateList(S,l);
			}

		//possible change of setup here with
		// max_layer_TABU(S) and high_cross_layers_TABU(S)
		sweep_TABU(S);



		currentIter++;

		int currentTimer = chrono::duration_cast<chrono::milliseconds>(
				chrono::system_clock::now() - start_algorithm).count();

		if((double)((double)currentTimer / (double)1000) >= TIME_LIMIT){
			break;
		}


		if(last_improv > Max_no_improv){
			int rand_int_seed = distribution(globalGenerator);
			diversificate_search(S, rand_int_seed);
		}

	}



	set_end_iteration(chrono::system_clock::now());
	print_best_data_out_file();
	close_out_file();
	std::cout<<"Best "<<best_solution.update_cost_bottleneck()<<std::endl;


}

void TABU::max_layer_TABU(HDAG &S){

	S.get_cost_bottleneck();
	S.computeScoreCoefficient_v();
	S.computeSolutionScore();

	double currentBestScore = S.getSolutionScore();
	bool move_found = false;

	//node position to swap for best move in each level
	unsigned best_m_p;

	//If no improving move is found: pick the least worsening one to
	//excape from local minima
	unsigned least_worsening_p1 = 0;
	unsigned least_worsening_p2 = 0;

	bool perform_worsening = true;


	double least_worsening_score;

	unsigned currentCost = S.update_cost_bottleneck();

	auto mapElement = S.getHowManyWCost()[currentCost].begin();

	//select level with max cost
	unsigned l = mapElement->first;

	least_worsening_score = std::numeric_limits<double>::max();


	for (unsigned pos : CL[l]) {

		for (unsigned i = 0; i < S.getLEVELS()[l].size(); i++) {

			if (currentIter - TabuTrack[l][i] > tenure) {
				//Try a move
				S.swapPositions(pos, i, l);
				//S.reset_bottleneck_structures();
				//S.get_cost_bottleneck();
				S.update_bottleneck_str_afterSwap(pos,i,l);

				S.computeSolutionScore();

				if (S.getSolutionScore() < currentBestScore) {

					best_m_p = i;
					currentBestScore = S.getSolutionScore();
					move_found = true;

				}

				else if ((S.getSolutionScore() < least_worsening_score)
						&& (pos != i)) {
					least_worsening_p1 = pos;
					least_worsening_p2 = i;
					least_worsening_score = S.getSolutionScore();
				}

				//undo the move
				S.swapPositions(pos, i, l);
				S.update_bottleneck_str_afterSwap(pos,i,l);
				S.computeSolutionScore();
			}
		}


		if (move_found) {
			S.swapPositions(pos, best_m_p, l);
			S.update_bottleneck_str_afterSwap(pos, best_m_p,l);
			TabuTrack[l][pos] = currentIter;
			TabuTrack[l][best_m_p] = currentIter;
			move_found = false;
			perform_worsening = false;
		}

		if(S.update_cost_bottleneck() < best_solution.update_cost_bottleneck()){

			set_end_iteration(chrono::system_clock::now());

			int elapsed =
					chrono::duration_cast<chrono::milliseconds>(
							get_end_iteration()
									- get_start_algorithm()).count();
			double current = (double) ((double) elapsed
					/ (double) 1000);

			best_solution = S;
			print_complete_data_out_file(S.update_cost_bottleneck(), current);
			set_time_to_best(current);


		}

	}

	if(perform_worsening){
		S.swapPositions(least_worsening_p1, least_worsening_p2, l);
		S.update_bottleneck_str_afterSwap(least_worsening_p1,
				least_worsening_p2, l);
		TabuTrack[l][least_worsening_p1] = currentIter;
		TabuTrack[l][least_worsening_p2] = currentIter;
	}

}



void TABU::high_cross_layers_TABU(HDAG &S) {

	S.get_cost_bottleneck();
	S.computeScoreCoefficient_v();
	S.computeSolutionScore();
	double currentBestScore = S.getSolutionScore();
	bool move_found;


	//node position to swap for best move in each level
	unsigned best_m_p;




	bool perform_worsening;

	bool enough_nodes;

	double least_worsening_score;



	for (unsigned l = 0; l < S.getLevNumber(); l++){

		int currentMaxCross = S.get_cost_bottleneck();
		double perc = 0.3;
		int costToCheck = currentMaxCross;
		while(!S.getAnyWCost()[costToCheck][l]){
			costToCheck--;
			if(costToCheck == 0) break;
		}

		if(costToCheck <= perc*currentMaxCross) continue;

		move_found = false;
		perform_worsening = true;
		enough_nodes = false;

		least_worsening_score = std::numeric_limits<double>::max();
		//If no improving move is found: pick the least worsening one to
		//excape from local minima
		unsigned least_worsening_p1 = 0;
		unsigned least_worsening_p2 = 0;

		for (unsigned pos : CL[l]) {

			for (unsigned i = 0; i < S.getLEVELS()[l].size(); i++) {

				if (currentIter - TabuTrack[l][i] > tenure) {

					enough_nodes=true;

					//Try a move
					S.swapPositions(pos, i, l);
					S.update_bottleneck_str_afterSwap(pos,i,l);

					S.computeSolutionScore();

					if (S.getSolutionScore() < currentBestScore) {

						best_m_p = i;
						currentBestScore = S.getSolutionScore();
						move_found = true;

					}

					else if ((S.getSolutionScore() < least_worsening_score)
							&& (pos != i)) {
						least_worsening_p1 = pos;
						least_worsening_p2 = i;
						least_worsening_score = S.getSolutionScore();
					}

					//undo the move
					S.swapPositions(pos, i, l);
					S.update_bottleneck_str_afterSwap(pos,i,l);
					S.computeSolutionScore();
				}
			}



			if (move_found) {
				S.swapPositions(pos, best_m_p, l);
				S.update_bottleneck_str_afterSwap(pos, best_m_p,l);
				TabuTrack[l][pos] = currentIter;
				TabuTrack[l][best_m_p] = currentIter;
				move_found = false;
				perform_worsening = false;

				auto CurrentID = S.getIDs();

				SwapCount[l][CurrentID[l][pos]]++;
				SwapCount[l][CurrentID[l][best_m_p]]++;

			}

			if(S.update_cost_bottleneck() < best_solution.update_cost_bottleneck()){

				set_end_iteration(chrono::system_clock::now());

				int elapsed =
						chrono::duration_cast<chrono::milliseconds>(
								get_end_iteration()
										- get_start_algorithm()).count();
				double current = (double) ((double) elapsed
						/ (double) 1000);

				best_solution = S;
				print_complete_data_out_file(S.update_cost_bottleneck(), current);
				set_time_to_best(current);

				last_improv = 0;
			}



		}

		if(perform_worsening && enough_nodes){

			S.swapPositions(least_worsening_p1, least_worsening_p2, l);
			S.update_bottleneck_str_afterSwap(least_worsening_p1,
					least_worsening_p2, l);
			TabuTrack[l][least_worsening_p1] = currentIter;
			TabuTrack[l][least_worsening_p2] = currentIter;

			auto CurrentID = S.getIDs();

			SwapCount[l][CurrentID[l][least_worsening_p1]]++;
			SwapCount[l][CurrentID[l][least_worsening_p2]]++;
		}

	}


}

void TABU::sweep_TABU(HDAG &S) {

	S.get_cost_bottleneck();
	S.computeScoreCoefficient_v();
	S.computeSolutionScore();
	double currentBestScore = S.getSolutionScore();
	bool move_found;


	//node position to swap for best move in each level
	unsigned best_m_p;


	bool perform_worsening;

	bool enough_nodes;

	double least_worsening_score;



	for (unsigned l = 0; l < S.getLevNumber(); l++){

		move_found = false;
		perform_worsening = true;
		enough_nodes = false;

		least_worsening_score = std::numeric_limits<double>::max();
		//If no improving move is found: pick the least worsening one to
		//excape from local minima
		unsigned least_worsening_p1 = 0;
		unsigned least_worsening_p2 = 0;

		for (unsigned pos : CL[l]) {


			for (unsigned i = 0; i < S.getLEVELS()[l].size(); i++) {

				if (currentIter - TabuTrack[l][i] > tenure) {

					enough_nodes=true;

					//Try a move
					S.swapPositions(pos, i, l);

					S.update_bottleneck_str_afterSwap(pos,i,l);

					S.computeSolutionScore();

					if (S.getSolutionScore() < currentBestScore) {

						best_m_p = i;
						currentBestScore = S.getSolutionScore();
						move_found = true;

					}

					else if ((S.getSolutionScore() < least_worsening_score)
							&& (pos != i)) {
						least_worsening_p1 = pos;
						least_worsening_p2 = i;
						least_worsening_score = S.getSolutionScore();
						}

					//undo the move
					S.swapPositions(pos, i, l);
					S.update_bottleneck_str_afterSwap(pos,i,l);
					S.computeSolutionScore();
				}
			}



			if (move_found) {
				S.swapPositions(pos, best_m_p, l);
				S.update_bottleneck_str_afterSwap(pos, best_m_p,l);
				TabuTrack[l][pos] = currentIter;
				TabuTrack[l][best_m_p] = currentIter;
				move_found = false;
				perform_worsening = false;

				auto CurrentID = S.getIDs();

				SwapCount[l][CurrentID[l][pos]]++;
				SwapCount[l][CurrentID[l][best_m_p]]++;

			}

			if(S.update_cost_bottleneck() < best_solution.update_cost_bottleneck()){

				set_end_iteration(chrono::system_clock::now());

				int elapsed =
						chrono::duration_cast<chrono::milliseconds>(
								get_end_iteration()
										- get_start_algorithm()).count();
				double current = (double) ((double) elapsed
						/ (double) 1000);
				best_solution = S;
				print_complete_data_out_file(S.update_cost_bottleneck(), current);
				set_time_to_best(current);

				last_improv = 0;
			}



		}

		if(perform_worsening && enough_nodes){

			S.swapPositions(least_worsening_p1, least_worsening_p2, l);
			S.update_bottleneck_str_afterSwap(least_worsening_p1,
					least_worsening_p2, l);
			TabuTrack[l][least_worsening_p1] = currentIter;
			TabuTrack[l][least_worsening_p2] = currentIter;

			auto CurrentID = S.getIDs();

			SwapCount[l][CurrentID[l][least_worsening_p1]]++;
			SwapCount[l][CurrentID[l][least_worsening_p2]]++;
		}

	}


}

TABU::~TABU() {
	// TODO Auto-generated destructor stub
}


void TABU::allocate_all_structures(HDAG &I){

	unsigned ln = I.getLevNumber();
	unsigned nel;
	allocate_TabuTrack(ln);
	allocate_SwapCount(ln);
	allocate_Prob_v(ln);

	for(unsigned l = 0; l < ln; l++){

		nel = I.getLEVELS()[l].size();
		allocate_TabuTrackLev(l, nel);
		allocate_SwapCountLevel(l,nel);
		allocate_Prob_v_lev(l,nel);

	}

	allocate_CL(ln);
	allocate_maxNonTabu_v(ln);


}

void TABU::allocate_TabuTrack(unsigned l){

	this->TabuTrack.resize(l);

}

void TABU::allocate_TabuTrackLev(unsigned l, unsigned t){

	this->TabuTrack[l].resize(t,0);

}

void TABU::allocate_maxNonTabu_v(unsigned l){

	this->maxNonTabu_v.resize(l,0);

}


void TABU::allocate_CL(unsigned l){

	this->CL.resize(l);

}

void TABU::compute_maxNonTabu_inLev(HDAG &S, unsigned l) {

	unsigned currentMax = 0;
	unsigned maxForNode = 0;
	unsigned pos_i = 0;

	for (unsigned i = 0; i < TabuTrack[l].size(); i++) {


		//check if Tabu or not
		if (currentIter - TabuTrack[l][i] > tenure) {

			pos_i = S.getPos()[l][i];

			//get max crossing among arcs incident to i
			maxForNode = S.get_MaxCrossing_incidentTo(l, pos_i);

			//Update the maximum in the level if needed
			if (maxForNode > currentMax) {
				currentMax = maxForNode;
				maxNonTabu_v[l] = maxForNode;
			}
		}
	}
}


void TABU::fill_CandidateList(HDAG &S, unsigned l){


	CL[l].clear();


	for(unsigned i = 0; i < TabuTrack[l].size(); i++){
		//check if Tabu or not
		if (currentIter - TabuTrack[l][i] > tenure) {

			//if the maximum crossing for arcs incident
			//to i is big enough, insert it in CL
			if(S.get_MaxCrossing_incidentTo(l,i) >= alpha*maxNonTabu_v[l]){
				CL[l].push_back(i);
			}

		}
	}
}

void TABU::compute_SelectionProb_v(){

	for(unsigned l = 0 ; l < SwapCount.size(); l++){
		unsigned totalSwap = 0;
		for(unsigned nod = 0; nod < SwapCount[l].size(); nod++){
			totalSwap = totalSwap + SwapCount[l][nod];
		}



		totalSwap = totalSwap/2;

		//denominator used to normalize the quantities TotalSwap-swap_count
		unsigned normalizingDen = 0;

		for(unsigned id = 0; id < SelectionProb_v[l].size(); id++){
			SelectionProb_v[l][id] = double(totalSwap - SwapCount[l][id]);
			normalizingDen+=SelectionProb_v[l][id];
		}

		for(unsigned id = 0; id < SelectionProb_v[l].size(); id++){
			SelectionProb_v[l][id] = double(totalSwap - SwapCount[l][id])/double(normalizingDen);
		}
	}

}

void TABU::diversificate_search(HDAG &S, int random_seed) {

	unsigned current_diverse = 0;
	bool improve_found = false;


	std::default_random_engine generator(random_seed);

	//while the number of diversification is less than the max allowed
	while (current_diverse < MAX_DIVERSE) {

		//compute selection prob
		compute_SelectionProb_v();

		//go through the layers and make a random move according to selection prob
		for (unsigned l = 0; l < S.getLEVELS().size(); l++) {

			/***************************/
			//If the level has a single node the diversification skips the level
			  if(S.getLEVELS()[l].size()<=2) continue;
			/***************************/
			std::discrete_distribution<int> current_distr(
					SelectionProb_v[l].begin(), SelectionProb_v[l].end());

			unsigned first_id = current_distr(generator);
			unsigned second_id = current_distr(generator);

			unsigned escapeCount=0;

			while (second_id == first_id) {
				second_id = current_distr(generator);
				escapeCount++;
				if(escapeCount++ > 10)break;

			}

			unsigned pos_first = S.getPos()[l][first_id];
			unsigned pos_second = S.getPos()[l][second_id];


			S.swapPositions(pos_first, pos_second, l);
			S.update_bottleneck_str_afterSwap(pos_first, pos_second, l);
			TabuTrack[l][pos_first] = currentIter;
			TabuTrack[l][pos_second] = currentIter;

			SwapCount[l][first_id]++;
			SwapCount[l][second_id]++;


			//compute the o.f. value of the new solution, if better break.
			if (S.update_cost_bottleneck()
					< best_solution.update_cost_bottleneck()) {

				int elapsed = chrono::duration_cast<chrono::milliseconds>(
						get_end_iteration() - get_start_algorithm()).count();
				double current = (double) ((double) elapsed / (double) 1000);
				std::cout
						<< "Ho improvato durante la diversificazione! new cost: "
						<< S.update_cost_bottleneck() << std::endl;
				best_solution = S;
				print_complete_data_out_file(S.update_cost_bottleneck(),
						current);
				set_time_to_best(current);

				improve_found = true;
				break;
			}

		}
		current_diverse++;

		if (improve_found) {
			break;
		}
	}

	last_improv = 0;


}

void TABU::reduce_total_crossing(HDAG &S){


	unsigned currentMaxCrossing = S.update_cost_bottleneck();

	unsigned bestTotalCross = 0;

	for(unsigned h = 0; h < S.getTotalCostDistr().size();h++){
		bestTotalCross += h*S.getTotalCostDistr()[h];
	}

	bestTotalCross= bestTotalCross/2;

	std::cout<<"The initial total crossing is: "<<bestTotalCross<<std::endl;

	unsigned noImprov = 0;
	unsigned maxNoImpr = 100;


	while(noImprov < maxNoImpr){

		noImprov++;

		for(unsigned lev = 0; lev < S.getLEVELS().size(); lev++){
			for(unsigned i = 0; i < S.getLEVELS()[lev].size(); i++){

				unsigned bestNode = 0;
				bool moveFound = false;

				//for(unsigned i_ = i + 1; i_ < S.getLEVELS()[lev].size(); i_++){
					for(unsigned i_ = 0; i_ < S.getLEVELS()[lev].size(); i_++){

					//try a move
					S.swapPositions(i,i_,lev);
					S.update_bottleneck_str_afterSwap(i,i_,lev);


					//if the move worsen the current main objective, undo the move
					if(S.update_cost_bottleneck() > currentMaxCrossing){
						S.swapPositions(i,i_,lev);
						S.update_bottleneck_str_afterSwap(i,i_,lev);
					}


					else{
						//to store the total crossings of the HDAG after the move
						unsigned totalInConsideration = 0;

						//compute the total after the move
						for(unsigned h = 0; h < S.getTotalCostDistr().size();h++){
							totalInConsideration += h*S.getTotalCostDistr()[h];
						}
						totalInConsideration = totalInConsideration /2;

						//if the move improves the total crossing, store it in memory
						if(totalInConsideration <= bestTotalCross){
							bestNode = i_;
							bestTotalCross = totalInConsideration;
							moveFound = true;
						}
						//undo the move
						S.swapPositions(i,i_,lev);
						S.update_bottleneck_str_afterSwap(i,i_,lev);



					}

				}

				if(moveFound){
					std::cout<<"The total crossing has been improved!"<<std::endl;
					std::cout<<"New total crossing: "<<bestTotalCross<<std::endl;
					S.swapPositions(i, bestNode,lev);
					S.update_bottleneck_str_afterSwap(i, bestNode, lev);
					//noImprov = 0;
				}

			}
		}

	}

}






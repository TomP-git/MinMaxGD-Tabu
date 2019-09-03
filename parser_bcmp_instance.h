#include <string>
#include <iostream>
#include <vector>

using namespace std;

pair<unsigned, unsigned> rescale(int i, vector<vector<vector< unsigned> > > &LEVELS) {

	pair<unsigned, unsigned> tmp;

	for (unsigned l = 0; l < LEVELS.size(); ++l) {
		if (i - (int)LEVELS[l].size() <= 0) {
			tmp.first = l;
			break;
		}
		i = i - (int)LEVELS[l].size();
	}

	tmp.second = i-1;

	return tmp;

}

void parse_instance(int argc, char *argv[]) {

	char *instance_to_parse = argv[1];

	string parsed_instance = argv[1];
	parsed_instance = "parsed_"+ parsed_instance;


	FILE *F_in = fopen(instance_to_parse, "r");

    if (F_in == NULL) {
        fprintf(stderr, "error opening: %s\n", instance_to_parse);
        return;
    }

	FILE *F_out = fopen(parsed_instance.c_str(), "w");


	int n_nodes, n_edges, n_layers;


	// read the number of nodes, number of edges and number of layers
	fscanf(F_in, "%d %d %d\n", &n_nodes, &n_edges, &n_layers);


	// Structure HDAG
	vector<vector<vector< unsigned> > > LEVELS; // levels of the HDAG
	// Allocate the space for the layers
	LEVELS.resize(n_layers);

	// Read all the dimensions for the layers
	for (int l = 0; l < n_layers; ++l) {
		int dim;
		fscanf(F_in, "%d ", &dim);
		// Allocate dimension for each layer
		LEVELS[l].resize(dim);

	}

	fscanf(F_in, "\n");

	// Build the forward star for each layers
	while (!feof(F_in)) {
		unsigned id_i;
		unsigned id_j;
		fscanf(F_in, "%d %d\n", &id_i, &id_j);
		pair<unsigned, unsigned> tmp_id_i = rescale(id_i, LEVELS);
		pair<unsigned, unsigned> tmp_id_j = rescale(id_j, LEVELS);
		LEVELS[tmp_id_i.first][tmp_id_i.second].push_back(tmp_id_j.second);
	}

	fclose(F_in);

	fprintf(F_out, "%d\n", n_layers);

	for (unsigned l = 0; l < LEVELS.size(); ++l) {
		fprintf(F_out, "%d ", (int)LEVELS[l].size());
	}
	fprintf(F_out, "\n");

	for (unsigned l = 0; l < LEVELS.size(); ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			fprintf(F_out, "1 %d", i);
			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
				fprintf(F_out, " %d", LEVELS[l][i][j]);
			}
			fprintf(F_out, "\n");
		}
	}

	fclose(F_out);

	return;

}

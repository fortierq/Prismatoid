#ifndef PRISMATOID_H
#define PRISMATOID_H

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/layered/DfsAcyclicSubgraph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/planarity/PlanarizationLayout.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/planarity/FastPlanarSubgraph.h>
#include <ogdf/orthogonal/OrthoLayout.h>
#include <ogdf/planarity/EmbedderMinDepthMaxFaceLayers.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/graphalg/MinimumCut.h>

#include "set.h"
#include <assert.h>
#include <queue>
#include <z3.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <time.h>
#include <algorithm>

class Prismatoid_Z3
{
private:
	//std::vector<int> nFacet;
	int dim;
	Z3_context ctx;
	std::vector<std::vector<Z3_func_decl> > faces; // k-faces (k=0 or 1)
	std::vector<std::vector<int> > adj; // adj[i] contain the indices of the vertices adjacent to faces[1][i] (they may be > 2)
	std::vector<std::vector<int> > adj_e; // adj_e[i] contain the indices of the edges adjacent to faces[0][i] 

	void get_facets(int var, std::vector<std::string> &F);
	void get_facets(std::string &var, std::vector<std::string> &F);
	// \returns true if the vertex v is adjacent to the edge e
	bool isAdj(std::string &v, std::string &e);

	std::string decl_to_str_var(Z3_func_decl &d);
	int get_val(Z3_model &m, Z3_func_decl &d);
	bool isSubset(std::string &s1, std::string &s2);
	/*! \returns 1 (2) if vertex belongs to G+ (G-), 0 otherwise */
	int belong_to(Z3_func_decl &face);
	ogdf::String get_color(int vertex);
	int get_nFacets(std::string &var);
	int get_dim(std::string &var);
	std::string inter(std::string &s1, std::string &s2);

	void make_graph(ogdf::Graph &G, Z3_model &m);
	void make_basis(ogdf::Graph &G, Z3_model &m, int nBasis);
	void write_graph(ogdf::Graph &G, std::string file_name);
	void write_smt(std::string file_name, Z3_model &m);

	int check_enforce_kconnectivity(ogdf::Graph &G);
	bool check_enforce_convexity(ogdf::Graph &G, bool enforce);
	void ensure_convexity(int s, int t, std::vector<std::string> &F_inter, std::vector< std::vector<ogdf::node> > &cut, 
		ogdf::Graph &G);
	bool bfs_in_face(std::vector<ogdf::node> &v, std::vector< std::vector<ogdf::node> > &cut, ogdf::Graph &G, bool enforce);
	
	Z3_ast toBool(Z3_ast a, bool b = true);
	bool bit_vector;

public:
	Prismatoid_Z3(std::string file_name, int dim, bool bit_vector);
	void run(int min_connectivity, bool find_all, bool enforce);
};

#endif

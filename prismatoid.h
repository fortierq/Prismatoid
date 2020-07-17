#ifndef PRISMATOID_H
#define PRISMATOID_H

#include "set.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <time.h>
#include <algorithm>
#include <cmath>

class Prismatoid 
{
private:
	std::string toVariable(std::vector<Set> &s);
	std::string toBoolVar(std::vector<Set> &s, bool b = true);
	std::string neg(std::string s);
	std::string formula(std::string s);

	int dim, maxFacets, nGeo; // nGeo == 2
	std::vector<int> nFacet;
	std::string file_name;
	bool verbose;
	bool use_lib2;
	bool bit_vector; 

public:
	Prismatoid(bool bit_vector, int dim, std::vector<int> nFacet, std::string file_name, bool verbose);

	void set_options();
	void get_output();
	void declare_face(int d);
	void assert_no_dstep();
	void assert_0_or_1(int d);
	void assert_degree(); 
	/*! \param k dimension of face on which vertices belong */
	void assert_vertex_on_face(int k); 
	void assert_vertex_on_facet(int geo); // for G+ and G- 

	static std::string setToString(Set &s, int k);
	//template<typename type> std::string toString(type x);

	// For tests:
	void assert_number_vertices(int geo, int d, int n);
	void assert_cube(int geo, bool only_bases);
};

#endif

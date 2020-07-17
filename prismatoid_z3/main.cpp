#include "Prismatoid_Z3.h"
#include <iostream>
#include <string> 

using namespace std;

int toInt(string a)
{
        istringstream iss(a);
        int tmp;
        iss>>tmp;
        return tmp;
}

int main(int argc, char *argv[])
{
	int dim = 3;
	if(argc >= 2)
		dim = toInt(argv[1]);
	vector<int> nFacet;
	Prismatoid_Z3 p("test.smt", dim, true);
	p.run(dim+1, false, true);
}

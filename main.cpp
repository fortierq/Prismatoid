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
	Prismatoid_Z3 p("../prismatoid/Release/test.smt", dim, true);
	p.run(dim + 1, false, true);
}

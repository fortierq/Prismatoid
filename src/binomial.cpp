#include "binomial.h"

#include <vector>

using namespace std;

static vector<unsigned int> lut;
static unsigned int lut_n = 0; // invariant: lut is filled in for all n < lut_n

static void fill_lut(unsigned int n)
{
	if (n < lut_n)
		return;

	lut.resize((n+1)*(n+2)/2);
	while(n >= lut_n) {
		// fill lut in for n == lut_n
		if (lut_n == 0) {
			lut[0] = 1;
		} else {
			unsigned int ofs_prev = (lut_n-1)*lut_n/2;
			lut[ofs_prev+lut_n] = 1;
			lut[ofs_prev+2*lut_n] = 1;
			for(unsigned int k = 1; k < lut_n; ++k)
				lut[ofs_prev+lut_n+k] = lut[ofs_prev+k-1] + lut[ofs_prev+k];
		}

		lut_n++;
	}
}

unsigned int binomial(int n, int k)
{
	if (n < 0 || k < 0 || k > n)
		return 0;

	unsigned int ofs = n*(n+1)/2 + k;

	if (ofs >= lut.size())
		fill_lut(n);

	return lut[ofs];
}

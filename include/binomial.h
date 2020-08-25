#ifndef guard_binomial_h
#define guard_binomial_h

/**
 * Return the binomial coefficient "n over k".
 *
 * \note This is an implementation using look-up tables, so
 * effectively O(1) runtime.
 */
unsigned int binomial(int n, int k);

#endif

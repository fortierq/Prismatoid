#ifndef guard_set_h
#define guard_set_h

#include "binomial.h"
#include <cassert>
#include <cstdlib>

//#include "binomial.h"

/**
 * Efficient manipulation of subsets of the set {0,1,..,max_element},
 * where \ref max_element is at least the number of bits of an \ref int.
 *
 * The philosophy of this class is that every set is an immutable object.
 * That is, all operations on sets generate new sets instead of modifying
 * existing instances.
 *
 * \note This is completely untested for sizeof(int) != 4.
 */
class Set {
	unsigned int bits;

public:
	enum { max_element = 31 };

	/**
	 * The default constructor creates the empty set.
	 */
	Set() : bits(0) {}

	/**
	 * Returns the set that contains the element x iff the bit (1 << x)
	 * is set in \ref _bits
	 */
	static Set from_bits(unsigned int _bits) {
		Set s;
		s.bits = _bits;
		return s;
	}

	/**
	 * Returns the set that contains only the element \ref x
	 */
	static Set singleton(unsigned int x) {
		assert(x <= max_element);
		return from_bits(1 << x);
	}

	/**
	 * Returns the set {0,1,..,n-1}.
	 */
	static Set first_n(unsigned int n) {
		assert(n <= max_element+1);
		if (n <= max_element)
			return from_bits((1 << n)-1);
		else
			return from_bits(0xffffffff);
	}

	// Basic set operations
	Set operator&(const Set other) const {
		return from_bits(bits & other.bits);
	}
	Set operator|(const Set other) const {
		return from_bits(bits | other.bits);
	}
	Set operator+(unsigned int x) const {
		assert(x <= max_element);
		return from_bits(bits | (1 << x));
	}
	Set operator-(unsigned int x) const {
		assert(x <= max_element);
		return from_bits(bits & ~(1 << x));
	}

	/**
	 * \return the number of elements in the set.
	 */
	unsigned int cardinality() const {
		unsigned int c = (bits & 0x55555555) + ((bits & 0xAAAAAAAA)>>1);
		c = (c & 0x33333333) + ((c & 0xCCCCCCCC)>>2);
		c = (c & 0x0F0F0F0F) + ((c & 0xF0F0F0F0)>>4);
		c = (c & 0x00FF00FF) + ((c & 0xFF00FF00)>>8);
		c = (c & 0x0000FFFF) + ((c & 0xFFFF0000)>>16);
		return c;
	}

	/**
	 * \return \c true iff the set is empty.
	 */
	bool empty() const {
		return bits == 0;
	}

	/**
	 * \return \c true iff the set contains the element x.
	 */
	bool contains(unsigned int x) const {
		assert(x <= max_element);
		return bits & (1 << x);
	}

	/**
	 * \return \c true iff the set is contained in the \ref other set.
	 */
	bool subset(const Set& other) const {
		return (bits & other.bits) == bits;
	}

	/**
	 * Calculate and return the 0-based index of this set in the colex order
	 * of all sets of this cardinality.
	 *
	 * Recall that e.g. the colex-order of 3-sets is
	 * 012, 013, 023, 123, 014, 024, 124, 034, 134, 234, 015, ...
	 *
	 * \note Running time is O(max_element).
	 */
	unsigned int colex_index() const {
		unsigned int idx = 0;
		unsigned int card = 0;
		for(unsigned int b = 0; b <= max_element; ++b) {
			if (contains(b)) {
				card++;
				idx = binomial(b,card) + idx;
			}
		}
		return idx;
	}

	/**
	 * Return the next set in colex order.
	 *
	 * Typical usage is to loop through all k-element subsets of {0..n-1}
	 * as follows:

	 * for(Set S = Set::first_n(k); !S.contains(n); S = S.colex_next()) {
	 *	...
	 * }
	 *
	 * \note Running time is O(max_element). We should be able to get
	 * amortized O(1) with a generator approach.
	 */
	Set colex_next() const {
		unsigned int card = 0;
		for(unsigned int b = 0; b < max_element; ++b) {
			if (contains(b)) {
				if (!contains(b+1)) {
					return from_bits((bits & ~((1 << (b+1))-1)) | (1 << (b+1)) | ((1 << card)-1));
				}
				card++;
			}
		}
		abort(); // exception?
	}

	/**
	 * Shift the set as if the given elements from the groundsets have been
	 * removed.
	 *
	 * Example:
	 * 347.remove_groundset(058) == 246
	 * 125.remove_groundset(014) == 02
	 */
	Set remove_groundset(const Set groundset) const {
		Set s(*this);
		for(int b = max_element; b >= 0; --b) {
			if (groundset.contains(b)) {
				unsigned int mask = (1<<b) - 1;
				s.bits = (s.bits & mask) | ((s.bits & ~(mask|(1<<b)))>>1);
			}
		}
		return s;
	}
};

#endif

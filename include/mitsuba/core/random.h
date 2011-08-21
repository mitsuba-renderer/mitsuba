/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__RANDOM_H)
#define __RANDOM_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/cobject.h>

/* 
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)  
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and 
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and 
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
 * \ingroup libcore
*/

/* Period parameters */  
#define MT_N 312
#define MT_M 156
#define MT_MATRIX_A 0xB5026F5AA96619E9ULL /* constant vector a */
#define MT_UPPER_MASK 0xFFFFFFFF80000000ULL /* most significant 33 bits */
#define MT_LOWER_MASK 0x7FFFFFFFULL /* least significant 31 bits */

MTS_NAMESPACE_BEGIN

/**
 * \brief %Random number generator based on Mersenne Twister
 * by Takuji Nishimura and Makoto Matsumoto.
 * \ingroup libpython
 */
class MTS_EXPORT_CORE Random : public SerializableObject {
public:
	/**
	 * \brief Construct a new seeded random generator.
	 *
	 * Uses the default seed on Windows and '/dev/urandom' 
	 * on OSX and Linux.
	 */
	Random();

	/**
	 * \brief Construct a random generator with a custom seed
	 */
	Random(uint64_t seed);

	/// Construct a new random generator seeded from a pre-existing one
	Random(Random *random);

	/// Unserialize a random generator
	Random(Stream *stream, InstanceManager *manager);

	/// Copy the state from another random generator
	void set(Random *random);

	/// Seed the random generator with a single 64bit value
	void seed(uint64_t value = 5489ULL);

	/// Seed the random generator from another random generator
	void seed(Random *random);

	/**
	 * \brief Seed the random generator from an array
	 * \remark This function is currently not exposed 
	 * by the Python bindings
	 */
	void seed(uint64_t *values, uint64_t length);

	/// Return an integer on the [0, 2^63-1]-interval 
	uint64_t nextULong();

	/// Return an integer on the [0, n)-interval 
	uint32_t nextUInt(uint32_t n);

	/// Return an integer on the [0, n)-interval 
	size_t nextSize(size_t n);

	/// Return a floating point value on the [0, 1) interval
	Float nextFloat();

	/**
	 * \brief Draw a uniformly distributed permutation and permute the 
	 * given STL container.
	 *
	 * See Knuth, TAoCP Vol. 2 (3rd 3d), Section 3.4.2.
	 *
	 * \remark This function is currently not exposed 
	 * by the Python bindings
	 */
	template <typename Iterator> void shuffle(Iterator it1, Iterator it2) {
		for (Iterator it = it2 - 1; it > it1; --it) 
			std::iter_swap(it, it1 + nextSize((size_t) (it-it1)));
	}

	/// Serialize a random generator to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Random() { }
private:
	uint64_t mt[MT_N]; /* the array for the state vector  */
	int mti;
};


MTS_NAMESPACE_END

#endif /* __RANDOM_H */

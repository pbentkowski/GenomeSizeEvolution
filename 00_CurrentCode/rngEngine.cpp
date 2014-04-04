// 		rngEngine.cpp

// 		Copyright:
// 		(C) 1997 Makoto Matsumoto and Takuji Nishimura with additions from
// 		Topher Cooper and Marc Rieffel in July-Aug. 1997
// 		(C) 1998 Brian Gough
// 		(C) 2011 Piotr Bentkowski <bentkowski.piotr@googlemail.com>
// 		
// 		This program is free software; you can redistribute it and/or modify
// 		it under the terms of the GNU General Public License as published by
// 		the Free Software Foundation; either version 2 of the License, or
// 		(at your option) any later version.
// 		
// 		This program is distributed in the hope that it will be useful,
// 		but WITHOUT ANY WARRANTY; without even the implied warranty of
// 		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// 		GNU General Public License for more details.
//
// 		You should have received a copy of the GNU General Public License
// 		along with this program; if not, write to the Free Software
// 		Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// 		MA 02110-1301, USA.
//
// 		This is Mersenne Twister random number generator singleton 
// 		implementation based on the C code of the GNU Scientific Library
// 		version GSL-1.14 . The seeding procedure matches the 2002 release
// 		of MT19937 .
//
// 		The original Makoto Matsumoto's code included the comment: "When 
// 		you use this, send an email to: matumoto@math.keio.ac.jp with an
// 		appropriate reference to your work".
//
// 		Makoto Matsumoto has a web page with more information about the
// 		generator, http://www.math.keio.ac.jp/~matumoto/emt.html. 
//
// 		The paper below has details of the algorithm.
//
// 		From: Makoto Matsumoto and Takuji Nishimura, "Mersenne Twister: A
// 		623-dimensionally equidistributerd uniform pseudorandom number
// 		generator". ACM Transactions on Modeling and Computer Simulation,
// 		Vol. 8, No. 1 (Jan. 1998), Pages 3-30
//
// 		You can obtain the paper directly from Makoto Matsumoto's web page.
//
// 		The period of this generator is 2^{19937} - 1.

#include "rngEngine.h"

rngEngine* rngEngine::s_pInstance = NULL;

unsigned long rngEngine::mt[N_rng];
int rngEngine::mti = N_rng + 1;

/** 
 * @brief Core method. Constructor of the RNG.
 */
rngEngine::rngEngine() {
}

/** 
 * @brief Core method. Destructor of the RNG.
 */
rngEngine::~rngEngine() {
}

/**
 * @brief Core method. Gets instance of the RNG.
 */
rngEngine* rngEngine::getInstance() {
    if (NULL == s_pInstance) {
        s_pInstance = new rngEngine();
    }
    return s_pInstance;
}

/**
 * @brief Core method. Releases instance of the RNG.
 */
void rngEngine::release() {
    if (NULL != s_pInstance) {
        delete s_pInstance;
        s_pInstance = NULL;
    }
}

/** 
 * @brief Seeding the Marsenne Twister random number generator.
 * 
 * Allows the user to seed the Marsenne Twister random number generator.
 * @param - the seed (unsigned long int)
 */
void rngEngine::rngMTset(unsigned long int s) {
    if (s == 0) {
        // the default seed is 4357
        s = 4357;
    }
    mt[0] = s & 0xffffffffUL;
    for (mti = 1; mti < N_rng; mti++) {
        // See Knuth's "Art of Computer Programming" Vol. 2, 3rd Ed. p.106
        // for multiplier.
        mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
        mt[mti] &= 0xffffffffUL;
    }
}

/** * @brief Marsenne Twister random number generator's guts
 * 
 * @param void *vstate
 * @return - pseudo-randomly generated integer
 */
unsigned long rngEngine::rngMTget() {
    unsigned long y;
#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)
    if (mti >= N_rng) {
        // generate N words at one time
        int kk;
        // if rngMTset() has not been called, a default initial seed is used
        if (mti == (N_rng + 1)) {
            rngMTset(5489UL);
        }
        for (kk = 0; kk < (N_rng - M_rng); kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + M_rng ] ^ (y >> 1) ^ MAGIC(y);
        }
        for (; kk < N_rng - 1; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (M_rng - N_rng)] ^ (y >> 1) ^ MAGIC(y);
        }
        y = (mt[N_rng - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N_rng - 1] = mt[M_rng - 1] ^ (y >> 1) ^ MAGIC(y);

        mti = 0;
    }
    // Tempering 
    y = mt[mti++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/** * @brief Generates a random double precision number
 * 
 * Generates a double precision number on [0,10)-real-interval from
 * Marsenne Twister random number generator
 * @param void *vstate
 * @return pseudo-random double precision number
 */
double rngEngine::rngMTgetDauble() {
    // divided by 2^32
    return rngMTget() / 4294967296.0;
}

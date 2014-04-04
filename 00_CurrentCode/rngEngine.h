// 		rngEngine.h

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


#ifndef _RNG_ENGINE_
#define _RNG_ENGINE_

#include <stdio.h>

#define N_rng (624) // Period parameters for Marsenne Twister
#define M_rng (397)
#define	UPPER_MASK (0x80000000UL) // most significant w-r bits
#define	LOWER_MASK (0x7fffffffUL) // least significant r bits 
#define	A (0x9908b0dfUL) // constant vector A

/** 
 * @class rngEngine
 *
 * @brief A singleton class which runs the Marsenne Twister random number
 * generator based on the C code of the GNU Scientific Library 
 * implementation of the MT RNG algorithm (version: GSL-1.14).
 */
class rngEngine {
private:
    rngEngine();
    virtual ~rngEngine();

public:
    static rngEngine* getInstance();
    void release();
    void rngMTset(unsigned long int s);
    double rngMTgetDauble();

private:
    unsigned long int rngMTget();
    static rngEngine* s_pInstance;
    // array for the state vector
    static unsigned long mt[N_rng];
    static int mti;
};

#endif //_RNG_ENGINE_

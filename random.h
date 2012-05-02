/* random.h

RandLib -- Wrapper classes for several RNG algorithms.
This file is part of the OpenQubit project.

Copyright (C) 1998-1999 OpenQubit.org
Yan Pritzker <yan@pritzker.ws>

Please see the CREDITS file for contributors.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

//! author="Joe Nelson" 
//! lib="RandLib"

template<class T>
class RandGenerator
//: Interface to which all RNG classes in RandLib must conform
{
public:
	//: Seeds the randomizer with a beginning value
	// In some children of RandGenerator, this function
	// may be meaningless. Also, seedVal2 need not be
	// used by all children.
	virtual void Seed(const unsigned int seedVal, const unsigned int seedVal2) = 0;

	//: Gets a rand between minVal and maxVal inclusive
	// PRECONDITIONS:  								<BR>
	// (minVal <= maxVal), 							<BR>
	// (minVal >= LowestRand), 					<BR>	
	// (maxVal <= HighestRand),	 				<BR>
	// (RandGenerationRange >= maxVal-minVal)	<BR>
	virtual T GetRandBetween(const T minVal, const T maxVal)
	{
		assert(minVal <= maxVal);
		assert(minVal >= LowestRand());
		assert(maxVal <= HighestRand());
		assert(RandGenerationRange() >= (long double)maxVal - (long double)minVal);		

		if(minVal == maxVal)
			return minVal;

		return minVal + GetRand(maxVal-minVal);
	}

	//: Highest random number that will be generated
	virtual const T HighestRand() const = 0;

	//: Lowest random number that will be generated
	virtual const T LowestRand() const = 0;

	//: The range (HighestRand - LowestRand) that can be generated
	virtual const T RandGenerationRange() const = 0;
protected:
	virtual T GetRand(const long double maxVal) = 0;
};

class IntStdRandGenerator : public RandGenerator<int>
//: Integer RNG encapsulating standard functions like rand()
// For details on methods, see the documentation
// for its interface: RandGenerator
{
public:
	IntStdRandGenerator(const unsigned int seedVal = 0, const unsigned int seedVal2 = 0)
	{ Seed(seedVal, 0); }
	virtual void Seed(const unsigned int seedVal, const unsigned int seedVal2)
	{
		if(seedVal == 0)
			srand(time(NULL));
		else
			srand(seedVal);
	}
	virtual const int HighestRand() const
	{ return RAND_MAX; }
	virtual const int LowestRand() const
	{ return -RAND_MAX; }
	virtual const int RandGenerationRange() const
	{ return RAND_MAX; }
protected:
	virtual int GetRand(const long double maxVal)
	{
		return static_cast<int>(maxVal*rand() / RAND_MAX);
	}
};

class DblUniformRandGenerator : public RandGenerator<double>
//: Double-precision RNG using "uniform" algorithm
// Encapsulates code compiled and "slightly" modified by Paul Bourk
// http://www.mhri.edu.au/~pdb/unixsoftware/random/
// For details on methods, see the documentation
// for its interface: RandGenerator
{
public:
	DblUniformRandGenerator(const unsigned int seedVal = 0, const unsigned int seedVal2 = 0)
	{ Seed(seedVal, seedVal2); }

	// PRECONDITIONS:  
	// (seedVal <= 31328), 
	// (seedVal2 <= 30081)
	virtual void Seed(unsigned int seedVal, unsigned int seedVal2)
	{
		double s,t;
		int ii,i,j,k,l,jj,m;

		/*
		Handle the seed range errors
		First random number seed must be between 0 and 31328
		Second seed must have a value between 0 and 30081
		*/
		assert(seedVal <= 31328);    // changed to assertions
		assert(seedVal2 <= 30081);

		if(seedVal == 0 || seedVal2 == 0) //added by Joe
		{
			seedVal  = time(NULL) % 31328;
			seedVal2 = time(NULL) % 30081;
		}

		i = (seedVal / 177) % 177 + 2;
		j = (seedVal % 177)       + 2;
		k = (seedVal2 / 169) % 178 + 1;
		l = (seedVal2 % 169);

		for (ii=0; ii<97; ii++) {
			s = 0.0;
			t = 0.5;
			for (jj=0; jj<24; jj++) {
				m = (((i * j) % 179) * k) % 179;
				i = j;
				j = k;
				k = m;
				l = (53 * l + 1) % 169;
				if (((l * m % 64)) >= 32)
					s += t;
				t *= 0.5;
			}
			u[ii] = s;
		}

		c    = 362436.0 / 16777216.0;
		cd   = 7654321.0 / 16777216.0;
		cm   = 16777213.0 / 16777216.0;
		i97  = 97;
		j97  = 33;
	}

	virtual const double HighestRand() const
	{ return DBL_MAX; }
	virtual const double LowestRand() const
	{ return -DBL_MAX; }
	virtual const double RandGenerationRange() const
	{ return DBL_MAX; }	
		
protected:
	double u[97],c,cd,cm;
	int i97,j97;

	virtual double GetRand(const long double maxVal)
	{ return maxVal * RandomUniform(); }

	double RandomUniform()
	{
		double uni;

		uni = u[i97-1] - u[j97-1];
		if (uni <= 0.0)
			uni++;
		u[i97-1] = uni;
		i97--;
		if (i97 == 0)
			i97 = 97;
		j97--;
		if (j97 == 0)
      		j97 = 97;
		c -= cd;
		if (c < 0.0)
			c += cm;
		uni -= c;
		if (uni < 0.0)
			uni++;

		return(uni);
	}
};

///////////////////////////////////////////////////
// This is a wrapper class for the
// DblUniformRandGenerator
class IntUniformRandGenerator : public RandGenerator<int>
//: Integer RNG using "uniform" algorithm
// Encapsulates code compiled and "slightly" modified by Paul Bourk
// http://www.mhri.edu.au/~pdb/unixsoftware/random/
// For details on methods, see the documentation
// for its interface: RandGenerator
{
public:
	IntUniformRandGenerator(const unsigned int seedVal = 0, const unsigned int seedVal2 = 0)
	{
		_inner = new DblUniformRandGenerator(seedVal, seedVal2);
	}	
	virtual ~IntUniformRandGenerator()
	{ delete _inner; }

	virtual void Seed(unsigned int seedVal, unsigned int seedVal2)
	{ _inner->Seed(seedVal, seedVal2); }

	virtual const int HighestRand() const
	{ return INT_MAX; }
	virtual const int LowestRand() const
	{ return -INT_MAX; }
	virtual const int RandGenerationRange() const
	{ return INT_MAX; }	
			
protected:
	RandGenerator<double> *_inner;

	virtual int GetRand(const long double maxVal)
	{ return static_cast<int>(_inner->GetRandBetween(0, maxVal)); }
};

#endif

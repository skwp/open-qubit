/* utility.cc

Implementation of various utility algorithms.
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

#include "utility.h"

bool IsBitSet(int n, unsigned short i)
//: Check if a bit is set in a number
   { return (n >> i) & 1; }

char *dtob(unsigned long value, unsigned short pad = 0)
{
	int numBits = count_bits(value), index = 0;
	int length = (pad > numBits) ? pad : numBits;
	char *result = new char[length+1];
	unsigned int mask = 1 << (numBits - 1);

	if(pad > numBits)
	{
		for(index = 0; index < pad-numBits; index++)
			result[index] = '0';
		result[index] = '\0';
	}

	// if pad>numBits, then index will be greater than zero at
	// this point (in fact, it will be positioned right after
	// the pad of zeroes), but if not, then index will have been
	// initialized to zero at its declaration
	for( ; mask && index < length; index++, mask >>= 1)
	{
		result[index] = (value & mask) ? '1' : '0';
	}
	result[index] = '\0';

   return result;
}

inline int count_bits(unsigned long value)
//: Count number of bits in a number
{
    int nbits = 0;
    while (value) {
        nbits++;
        value >>= 1;
    }

    // If we "need" 0 bits, return 1 so we at least print a 0.
    return nbits ? nbits : 1;
}


unsigned long CreateMask(int bits[])
{
	int bit;
	unsigned long mask;
	while(bit = *bits++) mask ^= (1 << bit);
	return mask;
}

int GCD(int a, int b)
// Greatest common divisor of a,b [Euclidean algorithm]
{
  while (a % b !=0)
  {
    int d = a % b;
    a = b;  
    b = d;
  }
  return(b);
}

int PeriodExtract(int v, int M, int domain)

// the function extracts period guess from FFT result
// this is continued fraction expansion taken from
// quant-ph/9809016

// v is the result from FFT measurement
// M is factorized number (it sets limit of the expansion)
// domain is number of states used in FFT

{
  int a0,a1,a2,p0,p1,p2,q0,q1,q2;
  double e0,e1,e2;


  if (v!=0) // if the period guess is 0, we will get nothing with it
  {
    int divisor = GCD(v,domain);
    v /= divisor;
    domain /= divisor;

    if (domain >= M)    // if the fraction cannot be cancelled 
         // or after cancelation the denominator is >=M
         // so it cannot be the true period and
                        // the result is not precise so use the recurense
                        // of continued fraction expansion
                   // otherwise this is the correct result
    {

// starting values
      a0 = int((double(v)/double(domain)));
      D("a0 = %d\n",a0);
      e0 = fabs(double(v)/double(domain) - a0);
      D("e0 = %f\n",e0);
      a1 = int((1/e0));
      D("a1 = %d\n",a1);
      e1 = fabs(1/e0 - a1);
      D("e1 = %f\n",e1);
      p0 = a0;
      p1 = a1*a0 + 1;
      q0 = 1;
      q1 = a1;
      D("p1 = %d\n",p1);
      D("q1 = %d\n",q1);
      q2=0;
// recurense starts
      while ((e1>1/domain) && (q2<M) ) //the first condition in order to
                                       //prevent overflows
      {
        a2 = int((1/e1));
        p2 = a2*p1 + p0;
        q2 = a2*q1 + q0;
        e2 = fabs(1/e1 - a2);
        e1 = e2;
        q0 = q1;
        p0 = p1;
        q1 = q2;
        p1 = p2;
        D("p2 = %d\n",p2);
        D("q2 = %d\n",q2);
      }
      if (q1==q2) // value from q1 was moved to q0
      {
        q1=q0;
        p1=p0;
      }
      D("q1 = %d\n",q1);
    }
   else  //exact value
    {
      p1 = v;
      q1 = domain;
    }
  return (q1/GCD(p1,q1)); //the fraction is cancelled here
  }
  else
    return(0);
}

int Reverse(int num, int nbits)
// reverse bits in num
// nbits is number of bits 

{
  int result=0;

  for (int i=0;i<nbits;i++)
    if (IsBitSet(num,i))
      result+=1<<(nbits-1-i);

  return(result);
}

// Russian Peasant modular exponentiation
int modexp(int x, int y, int m) {
        int xx = 1;
        int p = x % m;
        while(y > 0) {
                if(y & 1) xx = (xx*p) % m;
                p = p*p % m;
                y >>= 1;
        }
        return xx;
}

/*
bool IsNotPrime(int n)
// the function uses Fermat's theorem to check if the number is not prime
// For each prime x^(r-1) mod r == 1, if x<>r

{
  if (n!=3) //it does't work for 3
  {
    int result = modexp(3,n-1,n); //x=3 is good since most of non-primes
                                //gives result!=1, see D. E. Knuth,
				//"The Art of Computer Programming",
				//Vol.2 "Seminumerical algorithms",
				//Second ed., Addison-Wesley
    if (result==1)
      return(0);
    else
      return(1);
  }
  else
    return(0);
} 
*/

bool IsPrime(int n) {
	//this uses code borrowed from Bernhard Oemer
	int i;
	if (n<=1) return false;
	for (i=2; i<=floor(sqrt(n)); i++)
		if ((n%i)==0) return false;
	return true;
}

bool IsPrimePower(int n) {
	//code by Bernhard Oemer
	int i;
	int f;
	i=2;
	while (i<=floor(sqrt(n)) && f==0) {
		if((n%i)==0) f=i;
		i++;
	}

	for(i=2; i<=floor(log(n)/log(f)); ++i) 
		if(pow(f,i)==n) return true;

	return false;
}

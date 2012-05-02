/* main.cc

Main testing code.
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
#include <string.h>
#include <iostream.h>
#include "quantum"
#include <vector>

int Count(QState &q)
{
   static int count;
   int temp=0;
   for(int i=0; i<q.Outcomes(); ++i)
      if(abs(real(q[i])) > ROUND_ERR ||
         abs(imag(q[i])) > ROUND_ERR)
         temp++;

   count = (temp > count) ? temp : count;
	return count;
}

int main() {

	//some operators we will be using     
	ModExp MX;
	FFT	fft;

	int x;
	int M;
	char diag;	
	
	cout  << "OpenQubit version 0.2.0, Copyright (C) 1999 OpenQubit.org\n"
			<< "The OpenQubit library comes with ABSOLUTELY NO WARRANTY; "
			<< "for details\nread COPYING. This is free software, and you "
			<< "are welcome to redistribute\nit under certain conditions. "
			<< "Please see the file COPYING for details.\n\n";

	printf("Shor's algorithm for factoring numbers\n");
	printf("Would you like array usage diagnostics? (y/n) ");
	scanf("%s",&diag);
	printf("Enter number to factorize\n");
	scanf("%d",&M);
	
	if (M % 2 !=1)
	{
	  printf("The number is even. Factors found\n");
	  printf("%d = 2 * %d\n",M,M/2);
	  exit(0);
	} 

	if (IsPrime(M)) 
	{
	   printf("The number is prime. You cannot factor a prime number.\n");
	   exit(0);
	}

	//this doesn't work!
	if (IsPrimePower(M))
	{
		printf("The number is a prime power. It cannot be factored.\n");
		exit(0);
	}

	//usually this is taken randomly
	printf("Enter a number from 1..%d \n",M-1);
	scanf("%d",&x);
      
	int factor = GCD(M,x);

	if (factor!=1 && factor!=M) 
	{
		printf("Factor found since GCD(%d,%d)=%d\n",M,x,factor);
		printf("%d = %d * %d\n",M,factor,M/factor);
		exit(0);
	}

	// first part size
	int first=count_bits(M*M); //continued fraction expansion 
										//needs enough bits to represent 
	                           //M^2
	int firstsize=1<<first;

	// total register size
	int bits=first + count_bits(M);//values after modular exponentiation are < M
	int size= 1<<bits;

	QState *qureg;
	std::vector<Complex> StateCoefs(size);
	// /equal superposition of all the states in register 1
	// and |0..0> in register 2
	for (int k=0; k<firstsize;k++)
		StateCoefs[k]=1/sqrt(firstsize); 

	qureg = new QState(bits,StateCoefs);
	 
	//benchmarking tool
 	Count(*qureg);

	printf("Preparing equal superposition in the first register\n");
	//	qureg->Print();
 
	// Act with modular expenentiation function
	MX(*qureg,x,M,first);
	Count(*qureg);
	printf("Modular exponentiation\n");
	//        qureg->Print();
	
	//Uncomment if you want to try measuring the second register
	//It'd probably be working too

	for (int k=first;k<bits;k++)
		Measure(*qureg,k);
              
	//	qureg->Print();	

	// Fourier tranform in first register
	printf("Fourier transformation of the first register\n");
	fft(*qureg, first);
	Count(*qureg);
	//        qureg->Print();

	int result = 0;

	printf("Measurement\n");
	result=Measure(*qureg);
	printf("Measured state is:\n");
	qureg->Print();
	Count(*qureg);

	// extracting result of the first register from the total
	result = result % firstsize;

	// reverse bits since FFT gives bits reversed
	// Procedure ReverseBits doesn't work
	result=Reverse(result,first);

	printf("The result is %d\n",result);
	printf("Fourier domain is %d\n",firstsize);
	// Use continued fraction expansion to extract the period

	int period=PeriodExtract(result,M,firstsize);
	printf("Extracting the period using continued fraction expansion\n");
	printf("Period guess is: %d\n",period);

	if (period !=0 && modexp(x,period,M) == 1) //we check if it is true period
		printf("Period guess is probably correct\n");
	else
	{
		printf("Period guess is incorrect.\n");
		printf("Try this algorithm again with same staring numer\n");
		exit(0);
	}

	if (period % 2 == 0) // period is even. That's OK
	{
		factor=GCD(modexp(x,period/2,M)+1,M);
		if (factor!=1 && factor!=M) 
		{
			printf("Factors found!\n");
			printf("%d = %d * %d\n",M,factor,M/factor);
		}
		else
			if (factor==1)
				printf("Procedure failed due to bad period guess\n");
			else
	    	{
				printf("Procedure failed since %d^%d mod %d == -1\n",x,period/2,M);
				printf("Try again with another number\n");
				exit(0);
			}
		}
	else
		printf("Procedure failed; period is odd\n");
	
	int nonzerocoefs = Count(*qureg);
	int zerocoefs = size-nonzerocoefs;
	if (diag == 'y' || diag == 'Y') {
		printf("\nDiagnostics for factoring the number %d...\n", M);
		printf("Used %d qubits to factor this number.\n", bits);
		printf("Total number of elements in coef array: \t %d\n", 1 << bits);
		printf("Maximum non-zero coefs in array: \t\t %d\n", nonzerocoefs);
		printf("Total size (in bytes) of the array: \t\t %d\n", 
					sizeof(Complex) * size);
		printf("Bytes of array used by non-zero coefs: \t %d\n",
					nonzerocoefs*sizeof(Complex));
		printf("Percent of array wasted by zero coefs: \t %d bytes (%f%%)\n",
				zerocoefs*sizeof(Complex), zerocoefs*100/(double)size);				
	}
	return 0;
}

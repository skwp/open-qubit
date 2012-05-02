/* utility.h

Definition of several utility algorithms.
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

#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <stdio.h>
#include <math.h>
#include "debug.h"

int PeriodExtract(int v, int M, int domain);
int GCD(int a, int b);
int Reverse(int num, int nbits);

int count_bits(unsigned long value); 
char *dtob(unsigned long value, unsigned short pad = 0);

unsigned long CreateMask(int bits[]);
bool IsBitSet(int n, unsigned short i);
int modexp(int x, int y, int m);
//bool IsNotPrime(int n);

bool IsPrime(int n);
bool IsPrimePower(int n);

#endif 


/* iomanip.cc

Stream operators for printing QStates.
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
#include <iostream.h>
#include <iomanip>
#include "qstate.h"

ostream& ket(ostream& os, char* ket_val)
{
  return os << " |" << ket_val << '>';
}

omanip<char*> ket(char* ket_val)
{
  return omanip<char*> (ket, ket_val);
} 

ostream& coeff(ostream& os, Complex c) 
{
  if( abs(imag(c)) > ROUND_ERR )
      return os << '(' << setiosflags(ios::showpoint | ios::fixed)
                << setprecision(6) << setw(8) << real(c) << ", "
                << imag(c) << ')';

  		return os << setiosflags(ios::showpoint | ios::fixed)
            << setprecision(6) << setw(8) << real(c);
}

omanip<Complex> coeff(Complex c)
{
  return omanip<Complex> (coeff, c);
}

//added by Yan Pritzker 
bool ImagOrReal(const Complex &c)
{
	return (abs(real(c)) > ROUND_ERR || abs(imag(c)) > ROUND_ERR);	
}

bool ImagAndReal(const Complex &c) {
	return (abs(real(c)) > ROUND_ERR && abs(imag(c)) > ROUND_ERR);
}

ostream& operator<< (ostream& out, const QState &q)
{	
	int pad=count_bits(q._nStates-1);
	
	/* added by Yan Pritzker...to make sure + is printed even
		if there is more than one 0 coefficient between states */
	int nonzero=0; //number of nonzero probability states (for plus sign)
	for( int i=0; i < q._nStates; i++ )
		if(ImagOrReal(q._qArray[i])) nonzero++;
 
 
  for(int i=0; i < q._nStates-1; i++) {
    if(!ImagOrReal(q._qArray[i]))
      continue; 

    out << coeff(q._qArray[i]) << ket(dtob(i,pad));

	if(--nonzero) cout << " + ";
  }

  if(!ImagOrReal(q._qArray[q._nStates-1])) {
    out << endl; 
    return out; 
  }

  out << coeff(q._qArray[q._nStates-1]) 
       << ket(dtob(q._nStates-1,pad)) << endl;
	return out;

}

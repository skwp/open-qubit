/* qstate.cc

Implementation of Data Model of a Quantum State/Register
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
#include "qstate.h"

/* collapse code contributed by Rafal Podeszwa */

void QState::SetState(const std::vector<Complex> &c)
{
	double totalprob=0;
	
	for(int i=0; i<_nStates; i++)
		{ totalprob+= norm(c[i]); } 

	D("Sum of probabilities: %1.25f\n",totalprob);
	assert(1-ROUND_ERR <= totalprob && totalprob <= 1+ROUND_ERR );

	for(int j=0; j<_nStates; j++) 
		_qArray[j] = c[j];
}

int QState::_Collapse()
//collapse entire register
//this is an implementation of "Bernhard's Collapse"
//as suggested by Peter Belkner
{
	double x=0.0;
	double rnd = RNG->GetRandBetween(0,norm(*this));

	RandGenerator<int> *INTRNG = new IntUniformRandGenerator;
	unsigned long start = INTRNG->GetRandBetween(0,_nStates);
	unsigned long i=start;
	D("Normalized amplitudes: %f\n",norm(*this));
	
	D("Collapsing register. Got rnd %1.5f\n", rnd);
	
	while ((x += norm(_qArray[i])) < rnd) {
		if(++i == _nStates)
			i=0;
		
		if (start==i) {
			//this should never happen
			cerr << "Unexpected error #7891\n";
			i = (0 == i) ? _nStates - 1 : i - 1;
			break;
		}
	}
	
	int result=i;
	D("Set register state to %d\n", result);

	/* added by Yan Pritzker -- set all other coefs to 0 since
		qubit is collapsed; set collapsed state to prob 1 */

	for(int i=0; i < _nStates; i++) _qArray[i]=0;
	_qArray[result]=1;

	return result;
}

int QState::_Collapse(int index)
//single qubit collapse
{
	//index is 0..nQubits-1, 0 being the most significant qubit,
	//e.g. in |0001>, the 1 is the least significant qubit and
	//its position is nQubits-1
	assert(0 <= index && index < _nQubits);	

	double p0,p1;	//probabilities of measuring this bit as 0 and 1
	p0=p1=0.0;
	
	for (int i=0; i < _nStates; i++) {
		if (IsBitSet(i,index))
			p1 += norm(_qArray[i]);
		else
			p0 += norm(_qArray[i]);
		}

	D("Got probabilities...p0=%1.3f, p1=%1.3f\n", p0, p1);
   D("Sum of probabilities: %1.5f; Normalized amplitudes: %1.5f\n",
		 p0+p1, norm(*this));
	
	//sum of normalized amplitudes and sum of added probabilities
	//should equal..if not we messed up somewhere in the addition/normalization
	D("norm(*this) + ROUND_ERR = %1.25f\n", norm(*this) + ROUND_ERR);
	D("norm(*this) - ROUND_ERR = %1.25f\n", norm(*this) - ROUND_ERR);
	D("p0+p1 should be between those and it is %1.10f\n", p0+p1);
	assert( 	p0+p1 <= norm(*this) + ROUND_ERR &&
				p0+p1 >= norm(*this) - ROUND_ERR);

	double rnd=RNG->GetRandBetween(0,norm(*this));
	bool on  = p0 < rnd;
	
	for (int i=0;i < _nStates; i++) 
		if (IsBitSet(i,index))
   		if(!on) _qArray[i] = 0.0;
			else _qArray[i] /= sqrt(p1);
		else
			if(on) _qArray[i] = 0.0;
			else _qArray[i] /= sqrt(p0);
	
	D("Set bit state to %d\n", on);
	return on;
}

int QState::_CollapseSet(unsigned long bits)
{
	int i; 
	int state=0;

	for (i = 0; i < _nQubits; i++)
		if ((1 << i) & bits) 
			state += _Collapse(i);

	return state;
}

void QState::PrintSTD() const
{
	static const char* withimag		= "(%1.6f,%1.6f) |%s>";
	static const char* withimag_last	= "(%1.6f,%1.6f) |%s>\n";
	static const char* noimag			= "%1.6f |%s>";
	static const char* noimag_last	= "%1.6f |%s>\n";
   
	int nonzero=0;	//number of nonzero probability states (for plus sign)
	for( int i=0; i < _nStates; i++ ) 
		if(real(_qArray[i]) || imag(_qArray[i])) nonzero++;

	int pad=count_bits(_nStates-1);

	for(int i=0; i < _nStates-1; i++) {
      if(!imag(_qArray[i]) && !real(_qArray[i])) { continue; }
		if(imag(_qArray[i])) {
         printf(	withimag,
               	real(_qArray[i]),
               	imag(_qArray[i]),
               	dtob(i,pad));
      } else {
         printf(noimag, real(_qArray[i]), dtob(i,pad));
   	}
		if(--nonzero) { printf(" + "); } 
	}
 
   if(!imag(_qArray[_nStates-1]) && !real(_qArray[_nStates-1])) 
		{ printf("\n"); return; }
	if(imag(_qArray[_nStates-1])) {
		 printf ( withimag_last,
               real(_qArray[_nStates-1]),
               imag(_qArray[_nStates-1]),
               dtob(_nStates-1,pad));
   } else {
      printf ( noimag_last,
               real(_qArray[_nStates-1]),
               dtob(_nStates-1,pad));
	}
	printf("\n");
}

void QState::Dump(char filename[]) const
{
	FILE *FH;
	if ((FH=fopen(filename,"w"))==NULL)
		cerr << "ERROR: could not open file " << filename << endl;

	fprintf(FH,"QSTATE SIZE %d\n", _nStates);

	for(int i=0; i < _nStates; i++) {
		if(real(_qArray[i]) || imag(_qArray[i]))
			fprintf(FH, "%+1.17f \t %+1.17f \t |0x%X>\n",
					 real(_qArray[i]),
					 imag(_qArray[i]),
					 i);
	}
	fclose(FH);
}

void QState::Read(char filename[])
{
	FILE *FH;
	
	if((FH=fopen(filename,"r"))==NULL)
		cerr << "ERROR: could not open file " << filename << endl;

	double real,imag;
	int index;
	int size;

	fscanf(FH, "QSTATE SIZE %d\n", &size);
	assert(this->_nStates==size);
	
	while(FH)
	{
		int r=fscanf(FH,"%+1.17f \t %+1.17f \t |0x%X>\n", 
					 	 &real, &imag, &index);

		if(r==0) break;

		D("Scanned %f %f %d\n",real,imag,index);
		_qArray[index] = Complex(real,imag);
	}

	fclose(FH);
	D("Error: %f",norm(*this));
	assert(1-ROUND_ERR <= norm(*this) && norm(*this) <= 1+ROUND_ERR);
}

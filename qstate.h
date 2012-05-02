/* qstate.h

Data Model of a Quantum State/Register
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

#ifndef _QSTATE_H_
#define _QSTATE_H_
#include <map>
#include <vector>
#include <stdio.h>
#include <iostream.h>

#include "utility.h"					//for decimal to binary string conversions
#include "complex.h"					//complex numbers
#include "debug.h"					//debugging stuff
#include "random.h"					//random number generator

#define _RNGT_ double						//what type of generator
#define _RNG_  DblUniformRandGenerator	//specifically what type
 
//i think an error in the 12th decimal place is reasonable
//to accept as a rounding error, and not something important
//e.g. when adding norm(sqrt(1/2)) + norm(sqrt(1/2)) the answer 
//deviates from 1 by about 2e-16

//: Acceptable rounding error when normalizing.
static const double ROUND_ERR = 1e-12;

//! author = "Yan Pritzker, Peter Belkner, Rafal Podeszwa, Chris Dawson"
//! lib = "Quantum State [OpenQubit Core]"

class QState 
//: Model of a quantum state/register
{
private:
	std::vector<Complex> _qArray;		//: array of complex amplitudes
	RandGenerator<_RNGT_> *RNG;		//: random number generator
	
	int _nQubits;							//: number of qubits
	int _nStates;							//: number of states = 2^nQubits

	int _Collapse();						//: collapse entire register
	int _Collapse(int);					//: collapse a certain qubit
	int _CollapseSet(unsigned long);//: collapse a set of bits

	//: creates a state with no coefficients
	void _Clear()
		{ for(int i=0; i<_nStates; i++) _qArray[i]=0; }
	
	//: initialize array of complex amplitudes
	void _init(const std::vector<Complex> &c)
		{ RNG = new _RNG_; SetState(c); }

	//: default initializer (set up default coefs 1/sqrt(nStates))
	void _default_init()
		{ Reset(*this); RNG = new _RNG_; }
		
public:
	
	//: Default constructor. (Create one qubit)
	QState()
		: _nQubits(1), _nStates(1), _qArray(1)
		{ _qArray[0]=sqrt(1.0/2.0); }

	//: Coefficients not specified. Set up Walsh-Hadamard state (1/sqrt(size)).
	QState(int size)
		: _nQubits(size), _nStates(1 << size), _qArray(1 << size)
		{ assert(size>=1); _default_init(); }

	//: Constructor with initialization to complex amplitude array.
	QState(int size, const std::vector<Complex> &c) 
		: _nQubits(size), _nStates(1 << size), _qArray(1 << size)
		{ assert(size>=1); _init(c); }

	//: Default destructor.
	~QState() {};	

	//: Set state to specified coefficients.
	void SetState(const std::vector<Complex> &c);
	
	//: Display outcomes and coefficients [does not disturb state]
	void PrintSTD() const;

	//: New version of above. (uses streams) [DOES NOT COLLAPSE]
	void Print() const
		{ cout << *this; }

	//: Same as above.
	friend ostream& operator<<(ostream&, const QState&);

	//: Dump state to file. (Oemer-like Format) [does not disturb state]
	void Dump(char filename[]) const;
	
	//: Read state from file. [convenience function]
	void Read(char filename[]);	

	//: Returns total number of outcomes.
	int  Outcomes() const { return _nStates; }	

	//: Returns total number of bits.
	int  Qubits() const { return _nQubits; }

	//: Reset to base state |00...0>.
	friend void Reset(QState &q) 
		{ q._Clear();
		  q._qArray[0] = Complex(1); }

	//: Destructive register measure.
	// Note that the collapsing functions do not return another
	// state, as suggested by Peter Belkner, but rather destroy
	// the *this state. This is more physically correct because
	// when we look at a quantum state in real life, we don't 
	// get 'another' state, but rather the state we are looking
	// at collapses. After the collapse all the coefficients are
	// set to 0 except for that of the state that the register
	// collapsed to.
	friend long Measure(QState &q)
		{ return q._Collapse(); }

	//: Destructive bit measure.
	friend short Measure(QState &q, int i)
		{ return q._Collapse(i); }

	//: Destructive measure of a set of bits.
	friend int MeasureSet(QState &q, unsigned long bits)
		{ q._CollapseSet(bits); }

	//: Sum of normalized amplitudes. 
	friend double norm(const QState &q)
		{ 	double n=0.0;
			for(int i=0; i<q._nStates; i++)
				n += norm(q._qArray[i]);
			return n;
		}

	//: Implicit conversion operator causes destructive measure.
	operator long()
		{ return Measure(*this); }

	//: Access to coefficients.
	Complex& operator[] (int index)
		{ return _qArray[index]; }

};

#endif

/* qop.h

Several quantum gates. (one bit and multi-controlled)

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

//! author="Rafal Podeszwa, Christopher Dawson, Peter Belkner, Yan Pritzker"
//! lib="Quantum Operators and Gates"

/*

OK, Here's an explanation of the mess you find below :)

There are two base classes: SingleBit and Controlled. They represent
one-bit gates and multi-controlled one-bit gates (multi-bit gates),
respectively. 

Every operator is derived from these bases. In order
to enable any operator to be controlled and one-bit,
I implemented something I call  'templated inheritance'
which allow for example to specify a negation (NOT) like so:

	typedef opNOT<SingleBit>			Not;

In fact, you will find this declaration at the end of this
file. Its controlled counterpart is declared like so:

	typedef opNOT<Controlled>			CNot;

This syntax allows you to write an operator once and immediately
have a controlled counterpart for it.

Each operator also has a method called Param() which allows
you to re-specify the parameters for this gate. This allows
the reusal of a certain gate. Here's an example:

	RotQubit Rx(PI/2);	//gate which rotates qubit by PI/2
	Rx(mystate,0);			//rotate the first qubit by PI/2
	Rx.Param(PI/3);		//now it will perform PI/3 rotations
	Rx(mystate,1);			//the second bit is rotated by PI/3

Multi-controlled means that several bits can be used
as controlling bits. This also means that the controlling bit 
argument to the controlled functions is going to be a mask, not just
the number of the bit. 

e.g.  if you wanted to CNot the third bit using bits 1 and 3 to control
		the mask will be 101 == 5
		the controlled bit is bit 2 (third bit):
	
		CNot(myqstate, 5, 2); 

Also found below is an interesting class called DoAllBits
for lack of a better name. This also is a templated class
taking as a type, any of the operators that are defined.
This class allows any operator to operate on all bits in 
a state. So the Walsh-Hadamard operator is defined like:

	typedef opHadamard<SingleBit>		Hadamard;
	typedef DoAllBits<Hadamard>		WalshHadamard;

That's the only application I've found for it but i'm sure
there are more. It just saves time..you don't have to write
loops.

NOTE: FFT and SPhaseShift have been implemented separately
        as their application algorithms don't conform to either
        of the base classes.

	--Yan Pritzker <yan@pritzker.ws>

WARNING: The following code may be strange, unusual, perhaps 
			even offensive. Coder discretion is strongly advised.
*/

#ifndef _QOP_H_
#define _QOP_H_

#include <math.h>
#include "utility.h"
#include "qstate.h"

class SingleBit
//: Base for one-bit gates.
{
public:
	virtual void operator() (QState &q, int bit);

	void SetMatrix(const Complex &a00, const Complex &a01,
						const Complex &a10, const Complex &a11)
		{ _a00 = a00; _a01 = a01; _a10 = a10; _a11 = a11; }	

protected:

	//: Constructor to create gate matrix (Identity by default)
	SingleBit(const Complex &a00=1, const Complex &a01=0,
		  	    const Complex &a10=0, const Complex &a11=1)
			{ SetMatrix(a00,a01,a10,a11); }

private:
	
	Complex _a00, _a01, _a10, _a11;

};

class Controlled
//: Base for multi-bit (multi-controlled one bit) gates.
{
public:
  
	//: Operator for application of gate to a QState.
	virtual void operator() (QState &q, int mask, int bit);

	//: This version allows bitmask to be specified in an array of int.
	virtual void operator() (QState &q, int bits[], int bit)
		{ 	unsigned long mask=CreateMask(bits);
			operator()(q,mask,bit); }

        //: Allows reusal of a defined gate by changing its matrix.
	void SetMatrix(const Complex &a00, const Complex &a01,
						const Complex &a10, const Complex &a11)
		{ _a00 = a00; _a01 = a01; _a10 = a10; _a11 = a11; }	

protected:
	
	//: Constructor to create gate matrix (Identity Matrix by default)
	Controlled(const Complex &a00=1, const Complex &a01=0,
 				  const Complex &a10=0, const Complex &a11=1)
		{ SetMatrix(a00,a01,a10,a11); }

private:
	
	Complex _a00, _a01, _a10, _a11;

};

template <class BaseClassT>
class opUnitary : public BaseClassT
//: General Unitary operator.
{
public:
	void Param(double alpha, double beta, double delta, double theta)
	{                
		BaseClassT::SetMatrix
                (Complex(cos(delta+alpha/2+beta/2)*cos(theta/2),
                         sin(delta+alpha/2+beta/2)*cos(theta/2)),
						
                Complex(cos(delta+alpha/2-beta/2)*sin(theta/2),
                        sin(delta+alpha/2-beta/2)*sin(theta/2)),

                Complex(-cos(delta-alpha/2+beta/2)*sin(theta/2),
                        -sin(delta-alpha/2+beta/2)*sin(theta/2)),

                Complex(cos(delta-alpha/2-beta/2)*cos(theta/2),
                        sin(delta-alpha/2-beta/2)*cos(theta/2)));
 	}

	opUnitary(double a=0, double b=0, double d=0, double t=0)
		{ Param(a,b,d,t); }
};

template <class BaseClassT>
class opRotQubit : public BaseClassT 
//: Rotation on y-axis
{
public:
	//: Set the parameters needed for the gate to function.
	void Param(double theta) {
		BaseClassT::SetMatrix( 
			 cos(theta/2),  sin(theta/2),
			-sin(theta/2),  cos(theta/2)
		);
	}
	
	opRotQubit(double theta=0) { Param(theta); }
}; 

template <class BaseClassT>
class opRotPhase : public BaseClassT
//: Rotation on z-axis
{
public:
	void Param(double alpha) {
		BaseClassT::SetMatrix(
			Complex(cos(alpha/2), sin(alpha/2)), 0,
			0, Complex(cos(alpha/2), -sin(alpha/2))
		);
	}

	opRotPhase(double alpha=0) { Param(alpha); }
};

template <class BaseClassT>
class opPhaseShift : public BaseClassT
//: Scalar multiplication on z-axis
{
public:
	Param(double delta) {
		BaseClassT::SetMatrix(
			Complex(cos(delta), sin(delta)),0,
			0, Complex(cos(delta), sin(delta))
		);
	}

	opPhaseShift(double delta=0) { Param(delta); }
};

template <class BaseClassT>
class opHadamard : public BaseClassT
//: Hadamard operator [|0> -> |0> + |1> and |1> -> |0> - |1>]
{
public:
	opHadamard() : BaseClassT(M_SQRT1_2, M_SQRT1_2,
									  M_SQRT1_2, -M_SQRT1_2) {};
};

template <class BaseClassT>
class opNOT : public BaseClassT
//: SingleBit == Pauli negation, Controlled == CNot
{
public:
	opNOT() : BaseClassT(0,1,1,0) {};
};

/*** Below are some gates that are different enough that they are implemented
	  without implementing from the base classes ***/

class opFFT 
//: Fast Fourier Transform
{
public:
	opFFT() {};
	void operator() (QState &q, int numbits=-1);
};

class opSPhaseShift 
//: Shor's phase shift operator (double-controlled phase shift);
{
// the matrix used for controlled single bit operation is
// (1             0             )
// (0   exp(i*Pi*delta/2^(k-j)) )

private:
	//opPhaseShift<Controlled> CPS;
	//opRotPhase<Controlled> Rz;
	opUnitary<Controlled> CU;
public:
	opSPhaseShift() {};
	
	void operator() (QState &q, int j, int k)
	{
		assert(j < k);    //see Shor's paper
		double delta = M_PI/(1 << (k-j));
		unsigned long mask;
		
		mask = (1 << j);
		
		//call general unitary controlled gate as suggested by Rafal
		//instead of the four lines following. (works faster)
		CU.Param(delta,0,-delta/2,0);
		CU(q,mask,k);
		//CPS.Param(-delta/2);
		//CPS(q,mask,k);
		//Rz.Param(delta);
		//Rz(q,mask,k); 
	}
};

class ModExp
//: Modular Exponentiation
{
public:
	
	ModExp() {};
	void operator() (QState &q, int a, int n, int b) {
		//results go here first
		vector<Complex> _qArrayTmp(q.Outcomes()); 

		for(int i = 0; i < q.Outcomes(); i++)
			{
				if (q[i] != 0)
					_qArrayTmp[i + (modexp(a,i,n) << b)] = q[i];
			}

		for (int k = 0; k < q.Outcomes(); k++)
			q[k] = _qArrayTmp[k]; //copy resulta
	}
};

template <class OperatorType>
class DoAllBits
//: This class allows any operator to work on all bits of a state.
// This will only work for single bit operators requiring no parameters
// (i.e. this will not work with RotQubit which requires an angle theta)
{
private:
	OperatorType op;

public:
        void operator() (QState &q) {
                for(int i=0; i < q.Qubits(); ++i) 
                        op(q,i);
	}
};

/* the controlled or single-bit counterparts of several of the functions
bellow are probably not needed but they are there just in case  */

//: Unitary Operator
typedef opUnitary<SingleBit> 		Unitary;		

//: Controlled Unitary Operator
typedef opUnitary<Controlled>		CUnitary;

//: Qubit Rotation (Ry)
typedef opRotQubit<SingleBit>		RotQubit;

//: Controlled Qubit Rotation (Ry)
typedef opRotQubit<Controlled>	CRotQubit;

//: Phase Rotation (Rz)
typedef opRotPhase<SingleBit>		RotPhase;

//: Controlled Phase Rotation
typedef opRotPhase<Controlled>	CRotPhase;

//: Phase Shift (scalar multiplication of z axis)
typedef opPhaseShift<SingleBit>	PhaseShift;

//: Controlled Phase Shift
typedef opPhaseShift<Controlled>	CPhaseShift;

//: Shor Phase Shift
typedef opSPhaseShift				SPhaseShift;

//: Hadamard Operator
typedef opHadamard<SingleBit>		Hadamard;

//: Controlled Hadamard Operator
typedef opHadamard<Controlled>	CHadamard;

//: Negation Operator
typedef opNOT<SingleBit>			Not;

//: Controlled Negation (CNot - a.k.a XOR)
typedef opNOT<Controlled>			CNot;

//: Fast Fourier Transform
typedef opFFT							FFT;

//: Hadamard on all bits == WalshHadamard
typedef DoAllBits<Hadamard> 		WalshHadamard;

#endif

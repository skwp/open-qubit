#include "qop.h"

void SingleBit::operator() (QState &q, int bit=0)
{
	std::vector<Complex> _qArrayTmp(q.Outcomes());  //results will first go here
	int maski  = 1 << bit;  // set bit mask
	int k;
  
	for (k = 0; k < q.Outcomes(); k++)
		if (k & maski)  		// bit set so we use second row of the matrix
			_qArrayTmp[k] = _a10 * q[k ^ maski] + _a11 * q[k];
		else
			_qArrayTmp[k] = _a00 * q[k] + _a01 * q[k ^ maski];
       
		for (k = 0; k < q.Outcomes(); k++)
			q[k] = _qArrayTmp[k]; //copy result
}

void Controlled::operator() (QState &q, int mask, int bit=0)
//multi-controlled (via mask) operator
//suggested by Rafal Podeszwa
{
	std::vector<Complex> _qArrayTmp(q.Outcomes());
	int maski = 1 << bit;
	D("Controlling: %d \t Controlled: %d\n",mask, maski);
	D("In common: %d\n", mask & maski);
	assert((mask & maski) == 0); //can't control controlling bit
	int k;

	for(k=0; k < q.Outcomes(); k++) {
		if (k & mask)           //if mask is set
			if (k & maski)  //if controlled bit set
				_qArrayTmp[k] = _a10 * q[k ^ maski] + _a11 * q[k];
			else
				_qArrayTmp[k] = _a00 * q[k] + _a10 * q[k ^ maski];

		if (k & mask)
			q[k] = _qArrayTmp[k];
	}
}

void opFFT::operator() (QState &q, int numbits=-1)
{ 
	if (numbits == -1) numbits = q.Qubits();
   int j,k; 
   assert(numbits>=2); 	//need at least 2 qubits
  
	//operators we will be needing
   opSPhaseShift S;      //conditional phase shift
   Hadamard H;           //hadamard operator
	
   for (j = numbits-1; j>=0; j--) 
   { 
      for (k = numbits-1; k>j; k--)
      {
    		S(q,j,k);		
		D("S(%d,%d)",j,k);
      }
      D("H(%d)",j);
      H(q,j);
	
   } 
}

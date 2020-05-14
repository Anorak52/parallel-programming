#include <iostream>
#include <tbb/tbb.h>
#include "Fox.h"

using namespace tbb;
using namespace std;

int main(int argc, char** argv)
{
	double* A;
	double* B;
	double* CSeq;
	double* CFox;

	double StartFoxAlg = 0, FinFoxAlg = 0, StartSeqAlg = 0, FinSeqAlg = 0; //time

	int size = 16; //matrix size
	int procNum = 4; //num of proc
	int BlockSize = int(sqrt(procNum));  //size of block
	
	//
	task_scheduler_init init(procNum);
	//
	if (BlockSize * BlockSize == procNum && size % BlockSize == 0) {


		//initialization
		A = new double[size * size];
		B = new double[size * size];
		CSeq = new double[size * size];
		CFox = new double[size * size];
		MatrixCreate(A, B, CFox, CSeq, size);
		//

		//print on screen
		if (size < 30) {
			cout << "Matrix A" << endl;
			PrintMatrix(A, size);
			cout << "Matrix B" << endl;
			PrintMatrix(B, size);
		}
		//

		//FoxAlgPar
		tick_count StartFoxAlg = tick_count::now();
		if (procNum != 1) {
			
			Fox(A, B, CFox, size, procNum);
			
		}
		tick_count FinFoxAlg = tick_count::now();
		//

		//FoxAlgSeq
		///task_scheduler_init init(1);
		tick_count StartSeqAlg = tick_count::now();
		SequentialAlgorithm(A, B, CSeq, size);
		tick_count FinSeqAlg = tick_count::now();
		//
		

		//print on screen
		if (size < 30) {
			cout << "Parallel Alg " << endl;
			PrintMatrix(CFox, size);
		}
		if (size < 30) {
			cout << "Seq Alg " << endl;
			PrintMatrix(CSeq, size);
		}
		cout << "===========================================" << endl;
		cout << "Time of parallel alg is " << (FinFoxAlg - StartFoxAlg).seconds() << endl;
		cout << "Time of seq alg is " << (FinSeqAlg - StartSeqAlg).seconds() << endl;
		cout << "Acceleration = " << (FinSeqAlg - StartSeqAlg).seconds() / (FinFoxAlg - StartFoxAlg).seconds() << endl;
		//
	}
	else
	{
		cout << "Error" << endl;
	}
}
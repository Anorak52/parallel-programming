#include <iostream>
#include <omp.h>
#include <ctime>
#include <random>

using namespace std;

void PrintMatrix(double* matrix, int N) {

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
		{
			cout.width(3);
			cout << matrix[i * N + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void MatrixCreate(double* A, double* B, double* C, double *CSeq, int size) {
	for (int i = 0; i < size * size; ++i) {
		A[i] = rand() % 10 + 1;
		B[i] = rand() % 10 + 1;
		C[i] = 0;
		CSeq[i] = 0;
	}
}

bool matComparison(double* A, double* B, int size)
{
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
		{
			if (A[i * size + j] != B[i * size + j])
				return false;
		}
	}
	return true;
}


void Multiplication(double* Ablock, double* Bblock, double* Cblock, int _blockSize)
{
	double temp;
	for (int i = 0; i < _blockSize; i++)
		for (int j = 0; j < _blockSize; j++) {
			temp = 0;
			for (int k = 0; k < _blockSize; k++)
				temp += Ablock[i * _blockSize + k] * Bblock[k * _blockSize + j];
			Cblock[i * _blockSize + j] += static_cast<int>(temp * 100) / 100;
		}
}

void SequentialAlgorithm(double* A, double* B, double* C, int size) {
	for (int i = 0; i < size * size; ++i)
		C[i] = 0;
	Multiplication(A, B, C, size);
}

void Fox(double* A, double* B, double* C, int size, int procNum) {
	int Grid = int(sqrt(procNum));
	int BlockSize = size / Grid;
	omp_set_num_threads(procNum);
#pragma omp parallel
	{
		int numThreads = omp_get_thread_num();
		int RowIndex = numThreads / Grid;
		int ColIndex = numThreads % Grid;
		for (int iter = 0; iter < Grid; iter++)
		{
			for (int i = RowIndex * BlockSize; i < (RowIndex + 1) * BlockSize; i++)
			{
				for (int j = ColIndex * BlockSize; j < (ColIndex + 1) * BlockSize; j++)
				{
					for (int k = iter * BlockSize; k < (iter + 1) * BlockSize; k++)
					{
						C[i * size + j] += A[i * size + k] * B[k * size + j];
					}
				}
			}
		}
	} //pragma parallel
}

int main(int argc, char* argv[])
{
	double* A;
	double* B;
	double* CSeq;
	double* CFox;
	int size, procNum;
	double StartFoxAlg = 0, FinFoxAlg = 0, StartSeqAlg = 0, FinSeqAlg = 0; //time

	cout << "input the size and then the procnun: " << endl;

	if (argc > 1) {
		size = atoi(argv[1]);
		procNum = atoi(argv[2]);
	}
	else {
		cout << "Input the size of matrix and then procnum: " << endl;
		cin >> size;
		cin >> procNum;
	}

	int BlockSize = int(sqrt(procNum));

	if (BlockSize * BlockSize == procNum && size % BlockSize == 0) {	


		//инициализация
		A = new double[size * size];
		B = new double[size * size];
		CSeq = new double[size * size];
		CFox = new double[size * size];
		MatrixCreate(A, B, CFox, CSeq,size);
		//
	
		//вывод на экран
		if (size < 30) {
			cout << "Matrix A" << endl;
			PrintMatrix(A, size);
			cout << "Matrix B" << endl;
			PrintMatrix(B, size);
		}
		//

		//FoxAlgPar
		if (procNum != 1) {
			StartFoxAlg = omp_get_wtime(); //omp_get_wtime ();
			Fox(A, B, CFox, size, procNum);
			FinFoxAlg = omp_get_wtime();
		}
		//

		//FoxAlgSeq
		StartSeqAlg = omp_get_wtime(); //omp_get_wtime ();
		SequentialAlgorithm(A, B, CSeq, size);
		FinSeqAlg = omp_get_wtime();
		//

		//вывод на экран
		if (size < 30) {
			cout << "Parallel Alg " << endl;
			PrintMatrix(CFox, size);
		}
		if (size < 30) {
				cout << "Seq Alg " << endl;
			PrintMatrix(CSeq, size);
		}

		if (matComparison(CFox, CSeq, size))
			cout << "Matrix identical\n" << endl;

		cout << "time of parallel alg is " << (FinFoxAlg - StartFoxAlg) << endl;
		cout << "time of seq alg is " << (FinSeqAlg - StartSeqAlg) << endl;
		//

		delete[] A, B, CSeq, CFox;
	}
	else
	{
		cout << "Error" << endl;
	}

}
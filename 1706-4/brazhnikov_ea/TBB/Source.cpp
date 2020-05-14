#include <iostream>
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include <tbb/tbb.h>
//#include "FoxAlg.h"

using namespace tbb;
using namespace std;

void PrintMatrix(double* matrix, int N) {

	int M = N;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++)
		{
			cout << matrix[i * M + j] << " ";
			for (int k = 0; k < (5 - log10((double)(matrix[i * M + j]))); k++)
				cout << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void MatrixCreate(double* A, double* B, double* C, double* CSeq, int size) {
	for (int i = 0; i < size * size; ++i) {
		A[i] = rand() % 10 + 1;
		B[i] = rand() % 10 + 1;
		C[i] = 0;
		CSeq[i] = 0;
	}
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

class ClassFoxAlg {
private:
	static int threadCounter;
	static mutex myMutex;
	const double* A;
	const double* B;
	double* C;
	int size;
	int procNum;
public:
	static void SetToZerothreadCounter() {
		threadCounter = 0;
	}

	ClassFoxAlg (const double* A_, const double* B_, double* C_, int size_, int procNum_) : A(A_), B(B_), C(C_), size(size_), procNum(procNum_) {}

	void operator()(const tbb::blocked_range<int>& iter) const {
		int Grid = int(sqrt(procNum));
		int BlockSize = size / Grid;
		{

			tbb::mutex::scoped_lock lock;
			lock.acquire(myMutex);
			int RowIndex = threadCounter / Grid;
			int ColIndex = threadCounter % Grid;
			threadCounter++;
			lock.release();

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
		}
	}
};


tbb::mutex ClassFoxAlg::myMutex;
int ClassFoxAlg::threadCounter = 0;
void Fox(const double *A, const double *B, double* C, int size, int procNum) {
	task_scheduler_init init(procNum);

	ClassFoxAlg::SetToZerothreadCounter();

	ClassFoxAlg processing(A, B, C, size, procNum);

	parallel_for(blocked_range<int>(0, procNum, 1), processing);
}

int main(int argc, char** argv)
{
	double* A;
	double* B;
	double* CSeq;
	double* CFox;

	double StartFoxAlg = 0, FinFoxAlg = 0, StartSeqAlg = 0, FinSeqAlg = 0; //время

	int size = 16; //размер матрицы 
	int procNum = 4; //количество процессов
	int BlockSize = int(sqrt(procNum));  //размер сетки 

	if (BlockSize * BlockSize == procNum && size % BlockSize == 0) {


		//инициализация
		A = new double[size * size];
		B = new double[size * size];
		CSeq = new double[size * size];
		CFox = new double[size * size];
		MatrixCreate(A, B, CFox, CSeq, size);
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
			
			Fox(A, B, CFox, size, procNum);
			
		}
		//

		//FoxAlgSeq
		
		SequentialAlgorithm(A, B, CSeq, size);
		
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
		cout << "time of parallel alg is " << (FinFoxAlg - StartFoxAlg) << endl;
		cout << "time of seq alg is " << (FinSeqAlg - StartSeqAlg) << endl;
		//
	}
	else
	{
		cout << "Error" << endl;
	}
}
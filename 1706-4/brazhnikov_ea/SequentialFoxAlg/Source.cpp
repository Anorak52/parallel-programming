#include <iostream>
#include <random>
#include <chrono> //замер времени

using namespace std;

void PrintMatrix(double* matrix, int N) {

	int M = N;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++)
		{
			cout.width(3);
			cout << matrix[i * M + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void MatrixCreate(double* A, double* B, int size) {
	for (int i = 0; i < size * size; ++i) {
		A[i] = rand() % 10 + 1;
		B[i] = rand() % 10 + 1;
	}
}

void SequentialAlgorithm(double* A, double* B, double* C, int size) {
	for (int i = 0; i < size * size; ++i)
		C[i] = 0;

	double temp;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			temp = 0;
			for (int k = 0; k < size; k++)
				temp += A[i * size + k] * B[k * size + j];
			C[i * size + j] = temp;
		}
	}
}

void blockMatrixMult(double* A, double* B, double* C, int size)
{
	int block_size = sqrt(size);

	for (int i = 0; i < size * size; ++i)
		C[i] = 0;

	double temp;

	for (int i = 0; i < size; i += block_size) 
	{ 
		for (int j = 0; j < size; j += block_size)
		{
			for (int k = 0; k < size; k += block_size) 
			{
				for (int a = i; a < min(size, i + block_size); a++) 
				{
					for (int b = j; b < min(size, j + block_size); b++) 
					{
						temp = 0;
						for (int c = k; c < min(size, k + block_size); c++) 
						{
							temp += A[a * size + c] * B[c * size + b];
						}
						C[a * size + b] += temp;
					}
				}
			}
		}
	}
}

int main()
{
	int size = 10;
	double* A;
	double* B;
	double* CSeq;
	double* CBlock;

	A = new double[size * size];
	B = new double[size * size];
	CSeq = new double[size * size];
	CBlock = new double[size * size];

	MatrixCreate(A, B, size);

	if (size < 30) {
		cout << "Matrix A" << endl;
		cout << endl;
		PrintMatrix(A, size);
		cout << "Matrix B" << endl;
		cout << endl;
		PrintMatrix(B, size);
		cout << "---------------------" << endl;
		cout << endl;
	}

	auto begin = chrono::steady_clock::now();
	SequentialAlgorithm(A, B, CSeq, size);
	auto end = std::chrono::steady_clock::now();

	auto elapsed_ms = chrono::duration_cast<chrono::milliseconds>(end - begin);

	cout << endl;
	cout << "SeqAlgTime: " << elapsed_ms.count() << endl;

	auto beginFox = chrono::steady_clock::now();
	blockMatrixMult(A, B, CBlock, size);
	auto endFox = std::chrono::steady_clock::now();
	auto elapsed_msFox = chrono::duration_cast<chrono::milliseconds>(endFox - beginFox);

	cout << endl;
	cout << "Block multiplication time:" << elapsed_msFox.count() << endl;

	if (size < 30) {
		cout << endl;
		cout << "Seq multiplication:" << endl;
		PrintMatrix(CSeq, size); 
		cout << endl;
		cout << "Block matrix multiplication:" << endl;
		PrintMatrix(CBlock, size);
	}

	delete A, B, CSeq, CBlock;
	
}

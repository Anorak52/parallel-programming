#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <stdexcept>


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

void MatrixCreate(double* A, double* B, int size) {
	for (int i = 0; i < size * size; ++i) {
		A[i] = rand() % 10 + 1;
		B[i] = rand() % 10 + 1;
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

int main()
{
	int size = 10;
	double tmp;
	double* A = &tmp;
	double* B = &tmp;
	double* CSeq = &tmp;

	A = new double[size * size];
	B = new double[size * size];
	CSeq = new double[size * size];

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
		
	SequentialAlgorithm(A, B, CSeq, size);
	PrintMatrix(CSeq, size);

}



#pragma once

void PrintMatrix(double* matrix, int N);
void MatrixCreate(double* A, double* B, double* C, double* CSeq, int size);
void Multiplication(double* Ablock, double* Bblock, double* Cblock, int _blockSize);
void SequentialAlgorithm(double* A, double* B, double* C, int size);
void Fox(const double* A, const double* B, double* C, int size, int procNum);
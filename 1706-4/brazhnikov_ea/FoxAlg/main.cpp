#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <mpi.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <stdexcept>


using namespace std;

int cartSize;
int blockSize;
MPI_Comm comm_row;
MPI_Comm comm_col;
int coords[2];



void PrintMatrix(double *matrix, int N) {
	
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

void Scatter(double* Matrix, double* block, int size) {
	double* rowbuff = new double[blockSize * size];
	if (coords[1] == 0) {
		MPI_Scatter(Matrix, blockSize * size, MPI_DOUBLE, rowbuff, blockSize * size, MPI_DOUBLE, 0, comm_col);
	}
	for (int i = 0; i < blockSize; i++) {
		MPI_Scatter(&rowbuff[i * size], blockSize, MPI_DOUBLE, &(block[i * blockSize]), blockSize, MPI_DOUBLE, 0, comm_row);
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

void FoxAlg(double* partA, double* AMatrixblock, double* partB, double* partC)
{
	MPI_Status status;
	for (int z = 0; z < cartSize; ++z) {
		//рассылка блоков матрицы A по строкам процессорной решетки
		int main = (coords[0] + z) % cartSize; //вед проц ведущ рассылку
		if (coords[1] == main) {
			for (int i = 0; i < blockSize * blockSize; ++i)
				partA[i] = AMatrixblock[i]; //AMB используется для начального распр данных
		}
		MPI_Bcast(partA, blockSize * blockSize, MPI_DOUBLE, main, comm_row);


		Multiplication(partA, partB, partC, blockSize);

		int dest = coords[0] - 1;
		if (coords[0] == 0) dest = cartSize - 1;
		int source = coords[0] + 1;
		if (coords[0] == cartSize - 1) source = 0;
		MPI_Sendrecv_replace(partB, blockSize * blockSize, MPI_DOUBLE, dest, 0,
			source, 0, comm_col, &status);
	}
}

void SequentialAlgorithm(double* A, double* B, double* C, int size) {
	for (int i = 0; i < size * size; ++i)
		C[i] = 0;
	Multiplication(A, B, C, size);
}

void Fox(double* A, double* B, double* C, int size) {
	int procNum;
	int rank;

	MPI_Comm comm_cart;
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (procNum == 1) {
		SequentialAlgorithm(A, B, C, size);
	}
	else {
		cartSize = sqrt(procNum);
		blockSize = size / cartSize;
		int dimSize[2] = { cartSize, cartSize }; //Количество процессов в каждом измерении решетки
		int periods[2] = { 0, 0 };  // =1 для каждого измерения, являющегося периодическим 

		//инициализация памяти
		int numElemBlock = blockSize * blockSize;
		double* Ablock = new double[numElemBlock];
		double* Bblock = new double[numElemBlock];
		double* Cblock = new double[numElemBlock];
		double* AMatrixblock = new double[numElemBlock];
		double* rowbuff = new double[size * blockSize];

		for (int i = 0; i < numElemBlock; i++)  //обнуление матриц Цэ
			Cblock[i] = 0;
		//

		//создание решетки
		MPI_Cart_create(MPI_COMM_WORLD, 2, dimSize, periods, true, &comm_cart);
		MPI_Cart_coords(comm_cart, rank, 2, coords);  // Определение координат процесса в решетке 
		
		int subdims[2];   // =1 для каждого измерения, оставляемого в подрешетке

		// Создание коммуникаторов для строк процессной решетки
		subdims[0] = false;
		subdims[1] = true;
		MPI_Cart_sub(comm_cart, subdims, &comm_row);
		
		// Создание коммуникаторов для столбцов процессной решетки
		subdims[0] = true;
		subdims[1] = false;
		MPI_Cart_sub(comm_cart, subdims, &comm_col);
		//

		//распределение блоков между процессами
		Scatter(A, AMatrixblock, size);
		Scatter(B, Bblock, size);
		//

		//выполнение параллельного фокса(перемножения)
		FoxAlg(Ablock, AMatrixblock, Bblock, Cblock);
		//


		//сбор матрицы на ведущем процессе
		for (int i = 0; i < blockSize; i++) {
			MPI_Gather(&Cblock[i * blockSize], blockSize, MPI_DOUBLE,
				&rowbuff[i * size], blockSize, MPI_DOUBLE, 0, comm_row);
		}
		if (coords[1] == 0) {
			MPI_Gather(rowbuff, blockSize * size, MPI_DOUBLE, C, blockSize * size,
				MPI_DOUBLE, 0, comm_col);
		}
		//
	}
}

int main(int argc, char *argv[])
{
	int x;
	int size = 8;
	double tmp;
	double* A = &tmp;
	double* B = &tmp;
	double* CSeq = &tmp;
	double* CFox = &tmp;
	int rank, procNum;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);

	int BlockSize = sqrt(procNum);
	if (BlockSize * BlockSize == procNum && size % BlockSize == 0) {
		if (rank == 0) 
		{
			A = new double[size * size];
			B = new double[size * size];
			CSeq = new double[size * size];
			CFox = new double[size * size];

			MatrixCreate(A, B, size);

			cout << "Matrix A" << endl;
			PrintMatrix(A, size);
			cout << "Matrix B" << endl;
			PrintMatrix(B, size);
		}
		//ВСТАВИТЬ ВРЕМЯ
		Fox(A, B, CFox, size);

		if (rank == 0) 
		{
			SequentialAlgorithm(A, B, CSeq, size);
			cout << "Parallel Alg " << endl;
			PrintMatrix(CFox, size);

			cout << "Seq Alg " << endl;
			PrintMatrix(CSeq, size);
		}
	}
	cin >> x;
}



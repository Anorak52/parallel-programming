#define MSMPI_NO_DEPRECATE_20
#include <iostream>
#include "mpi.h"
#include <typeinfo>
using namespace std;

int rProc(int r, int root, int procNum)
{
	return (r - 1 + root) % procNum;
}


void treeGather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
	int ProcNum, ProcRank;
	MPI_Comm_rank(comm, &ProcRank);
	MPI_Comm_size(comm, &ProcNum);

	MPI_Aint sendExtent;
	MPI_Aint recvExtent;

	MPI_Type_extent(sendtype, &sendExtent);
	MPI_Type_extent(recvtype, &recvExtent);

	int r = ((ProcRank - root) + ProcNum) % ProcNum + 1;
	int recvProc = rProc(2 * r, root, ProcNum);

	double *temp = new double[7];


	if (2 * r > ProcNum) //sendятся процессы с 0 по 3
	{
		MPI_Send(sendbuf, sendcount, sendtype, rProc(r / 2, root, ProcNum), 0, comm);

	}

	if (2 * r <= ProcNum && ProcRank != root)  //принимают процессы 5,6
	{
		for (int i = recvProc; i <= recvProc + 1; i++) {

			MPI_Recv((char*)recvbuf + i * recvcount * recvExtent, recvcount, recvtype, i, 0, comm, MPI_STATUS_IGNORE);
		}

	}

	if (ProcRank == root)
	{
		for (int j = 0; j < sendcount; j++)
		{
			memcpy((char*)recvbuf + root * recvcount * recvExtent, sendbuf, sendExtent*sendcount);
		}
		for (int i = recvProc; i <= recvProc + 1; i++) 
		{
			for (int j = 0; j < (ProcNum - 1) / 2; j++) 
			{
				int a;
				MPI_Recv(&a, recvcount, recvtype, i, 0, comm, MPI_STATUS_IGNORE);
				MPI_Recv((char*)recvbuf + a * recvcount * recvExtent, recvcount, recvtype, i, 0, comm, MPI_STATUS_IGNORE);
			}
		}
	}

	if (2*r <= ProcNum && ProcRank != root) //по идее сендят 5,6
	{
		for (int i = recvProc; i <= recvProc + 1; i++) 
		{
			MPI_Send(&i, 1, MPI_INT, rProc(r / 2, root, ProcNum), 0, comm);
			MPI_Send((char*)recvbuf + i * recvcount * recvExtent, sendcount, sendtype, rProc(r / 2, root, ProcNum), 0, comm);
		}

			MPI_Send(&ProcRank, 1, MPI_INT, rProc(r / 2, root, ProcNum), 0, comm);
			MPI_Send(sendbuf, sendcount, sendtype, rProc(r / 2, root, ProcNum), 0, comm);
		
	}
}

const int root = 4;

int n = 10;

int main()
{
	setlocale(LC_ALL, "rus");
	int ProcNum, ProcRank;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int part, x;
	double StartMy = 0, FinMy = 0, TimeMy = 0;
	double StartMPI = 0, FinMPI = 0, TimeMPI = 0;

	//int* recv = new int[7];
	//int data;
	//int k = 0;
	
	double* recv = new double [7];
	double* recv2 = new double [7];
	double data;
	double k = 0.3;
	part = 1;
	
	data = ProcRank + k;
	cout << "Rank " << ProcRank << " Send " << data <<endl;

	StartMy = MPI_Wtime();
	treeGather(&data, part, MPI_DOUBLE, recv, part, MPI_DOUBLE, root, MPI_COMM_WORLD);
	FinMy = MPI_Wtime();

	if (ProcRank == root)
	{
		cout << endl;
		for (int i = 0; i < ProcNum; i++)
			cout <<"Recv treeGather " << recv[i] << endl;
	}

	StartMPI = MPI_Wtime();
	MPI_Gather(&data, part, MPI_DOUBLE, recv2, part, MPI_DOUBLE, root, MPI_COMM_WORLD);
	FinMPI = MPI_Wtime();

	if (ProcRank == root)
	{
		cout << endl;
		for (int i = 0; i < ProcNum; i++)
			cout << "Recv MPI Gather " << recv2[i] << endl;
	}

	if (ProcRank == root) {
		cout << "\n" << endl;
		TimeMy = (FinMy - StartMy);
		cout << "Time of my gather is " << TimeMy << endl;
		TimeMPI = (FinMPI - StartMPI)*10;
		cout << "Time of MPI gather is " << TimeMPI << endl;
	}


	cin >> x;

	MPI_Finalize();
	return 0;
}
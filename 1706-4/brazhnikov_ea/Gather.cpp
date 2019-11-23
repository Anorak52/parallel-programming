#define MSMPI_NO_DEPRECATE_20
#include <iostream>
#include <mpi.h>
#include <time.h>
#include <cstdlib> 
#include <ctime>
#include <conio.h>
#include <math.h>

using namespace std;

int rProc(int r, int root, int procSize)
{
	return (r - 1 + root) % procSize;
}

void treeGather(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
	void* recvbuf, int recvcount, MPI_Datatype recvtype,
	int root, MPI_Comm comm)
{
	int ProcSize, ProcRank;
	MPI_Comm_rank(comm, &ProcRank);
	MPI_Comm_size(comm, &ProcSize);
	int x;
	MPI_Aint sendExtent;
	MPI_Aint recvExtent;
	MPI_Type_extent(sendtype, &sendExtent);
	MPI_Type_extent(recvtype, &recvExtent);

	int r = ((ProcRank - root) + ProcSize) % ProcSize + 1;
	int recvProc = rProc(2 * r, root, ProcSize);

	for (int i = 0; i < 2; i++)
	{
		if (2 * r <= ProcSize) {
			for (int j = recvProc; j <= recvProc + 1; j++)
			{
				MPI_Recv((char*)sendbuf + j * sendcount*recvExtent, recvcount, recvtype, j, 0, comm, MPI_STATUS_IGNORE);
			}
		}
	}
	
	if (ProcRank != root)
	{
		MPI_Send(sendbuf, sendcount, sendtype, rProc(r / 2, root, ProcSize), 0, comm);
	}
}

const int root = 4;

int n = 10;

int plus = 0;
int main(int argc, char *argv[])
{
	setlocale(LC_ALL, "rus");
	int ProcSize, ProcRank;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	int part;
	part = 7 / ProcSize;
	int x;
	double* num = new double[7];
	double data;
	double k = 0.3;

	data = k + ProcRank;

	treeGather(&data, part, MPI_DOUBLE, num, part, MPI_DOUBLE, root, MPI_COMM_WORLD);
	
	if (ProcRank == root)
	{
		for (int i = 0; i < 7; i++)
		{
			cout << "HEY " << num[i] << endl;
		}
	}
	MPI_Finalize();
	return 0;
}
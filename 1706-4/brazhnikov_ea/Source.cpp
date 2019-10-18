#include <iostream>
#include <mpi.h>
#include <time.h>
#include <cstdlib> 
#include <ctime>
#include <conio.h>

using namespace std;


void PrintMatrix(int *matrix, int N, int M) {
	cout << "Our matrix is:" << endl;
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


int* CreateMatrix(int N, int M) {
	int *matrix;
	matrix = new int[N*M];
	return matrix;
}

int* CreateVector(int N) {
	int *result;
	result = new int[N];
	return result;
}

void Rand(int *matrix, int *res, int N, int M) {
	srand(time(0));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++)
		{
			matrix[i * M + j] = rand() % 1000 + 1;
		}
	}
}

//============================================================================================


int main(int argc, char *argv[])
{
	int *matrix, *result, *result_fin;
	int lines, columns, tempres;
	int proc_rank, proc_size;
	double StartSeqAlg = 0, FinSeqAlg = 0, StartParAlg = 0, FinParAlg = 0, TimeSeqAlg = 0, TimeParAlg = 0;
	int x, n = 1, realcol = 0, countcolumn = 1;

	MPI_Datatype columntype;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_size); //всего процессов
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank); //возвр ранг процесса


	//============================================================================
	if (proc_rank == 0) {


		cout << "Input the num of line" << endl;
		cin >> lines;
		cout << "Input the num of column" << endl;
		cin >> columns;
		
		if ((columns % (proc_size - 1)) != 0)
		{
			if ((columns - 1) % (proc_size - 1) != 0)
			{
				n = (columns - 2) / (proc_size - 1); //узнаем сколько столбцов в вектор класть 
				realcol = 2;
			}
			else
			{
				n = (columns - 1) / (proc_size - 1); //узнаем сколько столбцов в вектор класть 
				realcol = 1; //на случай если попали сюда, ставим ограничение для фора отправки, послед столбец считается здесь же
			}
		}

		else if ((columns % (proc_size - 1) == 0) && (columns > proc_size - 1))
		{
			n = columns / (proc_size - 1);
		}

		matrix = CreateMatrix(lines, columns);
		result = CreateVector(lines);
		result_fin = CreateVector(lines);
		Rand(matrix, result, lines, columns);

		if (columns < 15)
			PrintMatrix(matrix, lines, columns);

		MPI_Type_vector(lines, n, lines, MPI_INT, &columntype);
		MPI_Type_commit(&columntype);

		StartParAlg = MPI_Wtime();

		int j = 1;
		for (int i = 0; i < columns - realcol; i += n) {
			MPI_Send(&n, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
			MPI_Send(&lines, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
			MPI_Send(&matrix[i], 1, columntype, j, 0, MPI_COMM_WORLD);
			if (j < proc_size)
			{
				j++;
			}
		} 

		int tmp = 0;
		for (int i = 1; i < proc_size; i++) {
			MPI_Recv(result, n, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
				for (int j = 0; j < n; j++)
				{
					if (result[j] > 0) {
						result_fin[tmp] = result[j];
						tmp++;
					}
				}
		}
		FinParAlg = MPI_Wtime();

		for (int j = 1; j < realcol+1; j++)
		{
			int tmp1 = lines - j;
				if (n > 1)
				{
					result_fin[lines - j] = TMP_MAX;
						for (int i = 0; i < lines; i++)
							if (matrix[i * lines + tmp1] < result_fin[tmp1])
								result_fin[tmp1] = matrix[i * lines + tmp1];
				}
		}


		cout << "Result of parallel alg is " << endl;
		for (int i = 0; i < columns; i++)
			cout << "|"<< result_fin[i] << endl;
		
		
		//Линейный алг---------------------------------------------------------------------

		StartSeqAlg = MPI_Wtime();
		for (int j = 0; j < lines; j++) //столбец, мы его фиксируем 
		{
			result[j] = INT_MAX;
			for (int i = 0; i < lines; i++)
				if (matrix[i * lines + j] < result[j])
				{
					result[j] = matrix[i * lines + j];
				}
		}
		FinSeqAlg = MPI_Wtime();

		//Вывод--------------------------------------------------------------------------------

		cout << "\n" << endl;
		cout << "Result of Lin alg is " << endl;
		for (int i = 0; i < columns; i++)
			cout << "||" << result[i] << endl;

		cout << "\n" << endl;
		TimeParAlg = (FinParAlg - StartParAlg);
		cout << "Time of par alg is " << TimeParAlg << endl;

		TimeSeqAlg = (FinSeqAlg - StartSeqAlg) * 1000;
		cout << "Time of seq alg is " << TimeSeqAlg << endl;
		cout << "\n" << endl;

	}
	//-------------------------------------------------------------------------------
	


	else
	{
		MPI_Recv(&n, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
		MPI_Recv(&lines, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
		matrix = CreateVector(lines*n);
		result = CreateVector(lines);
		MPI_Recv(&matrix[0], lines*n, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);

		for (int j = 0; j < n; j++)
			result[j] = INT_MAX;
		
			for (int j = 0; j < n; j++) {
				result[j] = INT_MAX;
				for (int i = j; i < lines*n; i += n)
				{
					if (matrix[i] < result[j])
					{
						result[j] = matrix[i];
					}
				}
			}
	
		/*процесс отсылает результат 0 процессу*/
		MPI_Send(result, n, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}


	MPI_Finalize();
	_getch();
	return 0;
}
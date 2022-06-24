#include <mpi.h>
#include<windows.h>
#include<fstream>
#include<iostream>
#include<cmath>
#include<algorithm>
#include<omp.h>
using namespace std;
int const MM = 2048;//读入文件大小
double arr[MM + 10][MM + 10];//矩阵大小
double a[MM + 10][MM + 10];
void m_read() {
	ifstream fin;
	int n = MM;
	fin.open("input2048.dat");
	for (int i = 1;i <= n;i++) {
		for (int j = 1;j <= n;j++) {
			fin >> a[i][j];
		}
	}
}
void m_init(int n) {
	for (int i = 1;i <= n;i++) {
		for (int j = 1;j <= n;j++) {
			arr[i][j] = a[i][j];
		}
	}
}
void m_print(int n) {
	ofstream fout;
	fout.open("output2048.dat");
	for (int i = 1;i <= n;i++) {
		for (int j = 1;j <= n;j++) {
			fout << arr[i][j] << ' ';
		}
		fout << endl;
	}
}
void test(int n) {
	int myid, numprocs;
	long long head, tail, freq;
	m_init(n);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	if (myid == 0) {
		QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
		QueryPerformanceCounter((LARGE_INTEGER*)&head);
	}
	int low_bd, high_bd;
	if (myid != numprocs - 1) {
		low_bd = myid * (n - n % numprocs) / numprocs;
		high_bd = low_bd + (n - n % numprocs) / numprocs;
	}
	else {
		low_bd = (numprocs - 1) * (n - n % numprocs) / numprocs;
		high_bd = n;
	}
		for (int k = 1;k<= n;k++) {
			if (myid == k) {
				for (int j = k + 1;j <= n;j++) {
					arr[k][j] = arr[k][j] / arr[k][k];
				}
				arr[k][k] = 1.0;
			}
			MPI_Bcast(&arr[k][k], n - k + 1, MPI_DOUBLE, myid, MPI_COMM_WORLD);
			
			int low_bd2 = max(low_bd, k + 1);
			for (int i = low_bd2;i <high_bd;i++) {
				for (int j = k + 1;j <= n;j++) {
					arr[i][j] -= arr[i][k] * arr[k][j];
				}
				arr[i][k] = 0;
			}
		}

		MPI_Request request;
		MPI_Status status;
		if (myid != 0) {
			for (int i = low_bd;i <=high_bd;i ++) {
				MPI_Send(&arr[i][1], n, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
			}
		}
		if (myid == 0) {
			for (int i = 1;i <= n;i++) {
				MPI_Recv(&arr[i][1], n, MPI_DOUBLE, MPI_ANY_SOURCE, i, MPI_COMM_WORLD,&status);
				//MPI_Wait(&request, &status);				}
			}
		}
	if (myid == 0) {
		QueryPerformanceCounter((LARGE_INTEGER*)&tail);
		cout << n << '\t' << (tail - head) * 1000.0 / freq << "ms" << endl;
	}
}
int main(int argc, char* argv[]) {
	int provided;
	//MPI_Init(&argc, &argv);
	MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
	if (provided < MPI_THREAD_MULTIPLE)
		MPI_Abort(MPI_COMM_WORLD, 1);

	m_read();
	MPI_Barrier(MPI_COMM_WORLD);
	int t = 4;
	int thread_count = 2;

	for (int i = 1;i <= 9;i++) {
		t *= 2;
		test(t);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}



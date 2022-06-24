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
	fout.open("2048.dat");
	for (int i = 1;i <= n;i++) {
		for (int j = 1;j <= n;j++) {
			fout << arr[i][j] << ' ';
		}
		fout << endl;
	}
}
void test(int n) {
	int myid, numprocs;
	int  namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	long long head, tail, freq;
	m_init(n);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);
	if (myid == 0) {
		QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
		QueryPerformanceCounter((LARGE_INTEGER*)&head);
	}
	/*int thread_id = omp_get_thread_num();
	int nthreads = omp_get_num_threads();
	cout << "threadid: "<<thread_id<<"and processid "<<myid <<" of numprocs "<<numprocs<< endl;*/
	int thread_count = 4;
#pragma omp parallel num_threads(thread_count) 
	{
		for (int k0 = 1;k0 <= n;k0 = k0 + numprocs) {
			int lim = min(n - k0 + 1, numprocs);
			for (int k1 = 0;k1 < lim;k1++) {
				int k, start;
					k = k0 + k1;
					start = k0 + myid;
					if (start <= k) start += numprocs;
#pragma omp single
				{
					{
						if (myid == k1) {

							for (int j = k + 1;j <= n;j++) {
								arr[k][j] = arr[k][j] / arr[k][k];
							}
							arr[k][k] = 1.0;

						}
						MPI_Bcast(&arr[k][k], n - k + 1, MPI_DOUBLE, k1, MPI_COMM_WORLD);
					}
				}
#pragma omp for
				
					for (int i = start;i <= n;i += numprocs) {
						for (int j = k + 1;j <= n;j++) {
							arr[i][j] -= arr[i][k] * arr[k][j];
						}
						arr[i][k] = 0;
					}
				

			}
		}

		MPI_Request request;
		MPI_Status status;
		if (myid != 0) {
#pragma omp for
			for (int i = 1 + myid;i <= n;i += numprocs) {
				MPI_Isend(&arr[i][1], n, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, &request);
			}
		}
		if (myid == 0) {
#pragma omp for
			for (int i = 1;i <= n;i++) {
				if (i % 4 != 1) {
					MPI_Irecv(&arr[i][1], n, MPI_DOUBLE, (i - 1) % 4, i, MPI_COMM_WORLD, &request);
					MPI_Wait(&request, &status);
				}
			}
		}
	}
	if (myid == 0) {
		QueryPerformanceCounter((LARGE_INTEGER*)&tail);
		cout << n << '\t' << (tail - head) * 1000.0 / freq << "ms" << endl;
		m_print(n);
	}
	
	//int thread_id = omp_get_thread_num();
	//int nthreads = omp_get_num_threads();
	//printf("Hello,The World!from thread number %d (on %d)for the MPI process number %d (on %d)\n", thread_id, nthreads, myid, numprocs);
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
	/*m_read();
	m_init(8);
	m_print(8);*/
	return 0;
}



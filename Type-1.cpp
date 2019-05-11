#include<bits/stdc++.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <mpi.h>

using namespace std;


void multiply(int **A, int **B, int **C, int n, int rootp){
	cout<<"Printing Input Matrix"<<endl;
	for( int i = 0 ; i < n ; i++ ){
		for( int j = 0 ; j < n ; j++ ){
			cout<<A[i][j]<<" ";
		}
		cout<<endl;
	}
	for( int i = 0 ; i < n ; i++ ){
		for( int j = 0 ; j < n ; j++ ){
			cout<<B[i][j]<<" ";
		}
		cout<<endl;
	}
}

int main(int argc, char* argv[]){
	srand(time(NULL));
	int n = atoi(argv[1]);
	int p = atoi(argv[2]);
	int rootp = sqrt(p);
	n = (1<<n);
	int *A = new int[n];
	int *B = new int[n];
	int *C = new int[n];
	for( int i = 0 ; i < n ; i++ ){
		A[i] = new int[n];
		B[i] = new int[n];
		C[i] = new int[n]; 
	}
	for( int i = 0 ; i < n; i++ ){
		for( int j = 0 ; j < n ; j++ ){
			A[i][j] = rand()%100;
			B[i][j] = rand()%100;
		}
	}
	MPI_Init(&argc,&argv);
	multiply(A,B,C,n,rootp);
	MPI_Finalize();
	return 0;
}
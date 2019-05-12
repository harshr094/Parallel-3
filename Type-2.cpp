#include<bits/stdc++.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <mpi.h>

#define synchronized(m) \
    for(std::unique_lock<std::recursive_mutex> lk(m); lk; lk.unlock())
using namespace std;

recursive_mutex m_mutex;

void setC(int **C, int x, int y, int len){
	for(int i = x; i< x+ len; i++){
		for(int j = y; j < y + len; j++){
			C[i][j]=0;
		}
	}
}


int* getSlice(int **array, int x, int y, int len){
	int *ret = new int[len*len];
	int count = 0;
	for(int i = 0 ; i < len; i++){
		for(int j = 0 ; j < len; j++){
			ret[count++] = array[i+x][j+y];
		}
	}
	return ret;
}
int getProcessorId(int x, int y, int rootp){
	return rootp*x+y;
}


void print(int **array, int n){
	for(int i = 0 ; i < n ;i++){
		for(int j = 0 ; j < n ;j++){
			cout<<array[i][j]<<" ";
		}
		cout<<endl;
	}
}

int** get_2d(int* array,int n){
	int** ret = new int*[n];
	for(int i = 0 ; i < n ; i++)
		ret[i] = new int[n];
	int count=0;
	for(int i = 0 ; i < n ; i++)
	for(int j = 0; j < n ; j++)
		ret[i][j] = array[count++];
	return ret;
}

void addFinal(int** C, int* data, int n){
	int** A = get_2d(data,n);
	for(int i = 0 ; i < n ; i++)
	for(int j = 0 ; j < n ; j++)
		C[i][j]+=A[i][j];
}

void multiplyMatrix(int** C, int* a, int* b, int n){
	int**  A = get_2d(a,n);
	int**  B = get_2d(b,n);
	for(int i = 0 ; i < n ; i++ )
	for(int j = 0 ; j < n ; j++ )
		C[i][j]=0;
	for(int i = 0 ; i < n ; i++){
		for(int j = 0 ; j < n ; j++){
			for(int k = 0 ; k < n ; k++){
				C[i][j]+=A[i][k]*B[k][j];
			}
		}
	}
}

void copyOriginal(int **C,int **copy, int x, int y, int size){
		cilk_for(int i = 0 ; i < size; i++)
		for(int j = 0 ; j < size; j++)
			C[i+x][j+y]+=copy[i][j];
	
}

void verify(int** A,int **B, int** C, int n){
	int** correct = new int*[n];
	cout<<"Verifying Result"<<endl;
	for(int i = 0 ; i < n ; i++){
		correct[i] = new int[n];
		for(int j = 0 ; j < n ;j++)
			correct[i][j]=0;
	}
	for(int i = 0 ; i < n ;i++){
		for(int j = 0 ; j < n ;j++){
			for(int k = 0 ; k < n ; k++){
				correct[i][j]+=A[i][k]*B[k][j];
			}
		}
	}
	bool incorrect = false;
	for(int i = 0 ; i < n ; i++){
		for( int j = 0 ; j < n ; j++){
			if(correct[i][j]!=C[i][j])
				incorrect = true;
		}
	} 
	if(incorrect)
		cout<<"WRONG ANSWER"<<endl;
	else
		cout<<"ACCEPTED"<<endl;
}

void multiply(int **A, int **B, int **C, int n, int rootp){
	int process_id;
	MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
	if(process_id >= rootp*rootp)
		return;
	double t1,t2;
	if(process_id==0)
	 	t1 = MPI_Wtime();
	int x = process_id/rootp;
	int y = process_id%rootp;
	//cout<<"Running Processor "<<process_id<<"~"<<x<<"~"<<y<<endl;
	int size = n/rootp;
	setC(C,x*size,y*size,size);
        //cout<<"C Matrix set 0"<<endl;	
	int* sendA = getSlice(A,x*size,y*size,size);
	int* sendB = getSlice(B,x*size,y*size,size);
	int** ccopy = new int*[size];
	for(int i = 0 ; i < size; i++)
		ccopy[i] = new int[size];
	for( int l = 1; l <= rootp; l++){
		int k = (y + l -1)%rootp;
		if(k==i){
			MPI_Bcast(&sendB[0],size*size,MPI_INT,k,MPI_COMM_WORLD);
		}
		multiplyMatrix(ccopy,sendA,sendB,size);
		copyOriginal(C,ccopy,x*size,y*size,size);
		if(l<rootp){
			int left_processor = getProcessorId(x,(rootp+y-1)%rootp,rootp);
			MPI_Isend(&sendA[0],size*size,MPI_INT,left_processor,0,MPI_COMM_WORLD,&sendreq[0]);
	        int* dummyA = new int[size*size];
			left_processor = getProcessorId(x,(rootp+y+1)%rootp,rootp);
			MPI_Irecv(&dummyA[0],size*size,MPI_INT,left_processor,0,MPI_COMM_WORLD,&recvreq[0]);
			MPI_Wait(&recvreq[0],&stat);
			sendA = dummyA;
	 	}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// if(process_id!=0){
	// 	int* send = getSlice(C,0,0,n);
	// 	MPI_Isend(&send[0],n*n,MPI_INT,0,0,MPI_COMM_WORLD,&sendreq[0]);
	// }
	// else{
	// 	t2 = MPI_Wtime();
	// 	int* data = new int[n*n];
	// 	MPI_Request reqst[rootp];
	//  	for(int i = 1; i < rootp*rootp; i++){
	// 		MPI_Irecv(&data[0],n*n,MPI_INT,i,0,MPI_COMM_WORLD,&reqst[i]);
	// 		MPI_Wait(&reqst[i],&stat);
	// 		addFinal(C,data,n);
	// 	}
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(process_id==0){
	// 	printf("Time Taken %.6lf\n",t2-t1);
		//cout<<"Multiplication Done"<<endl;
        	//cout<<"Print A"<<endl;
       	 	//print(A,n);
        	//cout<<"Print B"<<endl;
        	//print(B,n);
        	//cout<<"Print C"<<endl;
        	//print(C,n);
        	//verify(A,B,C,n);
	//}
}

int main(int argc, char* argv[]){
	srand(time(NULL));
	int n = atoi(argv[1]);
	int l = atoi(argv[2]);
	MPI_Init(&argc,&argv);
	int p;
	MPI_Comm_size(MPI_COMM_WORLD,&p); 
	int rootp = (1<<l);
	n = (1<<n);
	int **A = new int*[n];
	int **B = new int*[n];
	int **C = new int*[n];
	for( int i = 0 ; i < n ; i++ ){
		A[i] = new int[n];
		B[i] = new int[n];
		C[i] = new int[n]; 
	}
	for( int i = 0 ; i < n; i++ ){
		for( int j = 0 ; j < n ; j++ ){
			A[i][j] = -10+rand()%20;
			B[i][j] = -10+rand()%10;
			C[i][j]=0;
		}
	}
	multiply(A,B,C,n,rootp);
	MPI_Finalize();
	cout<<"Multiplication Done"<<endl;
    cout<<"Print A"<<endl;
    print(A,n);
    cout<<"Print B"<<endl;
    print(B,n);
    cout<<"Print C"<<endl;
    print(C,n);
    verify(A,B,C,n);
	return 0;
}

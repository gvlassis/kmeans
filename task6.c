#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

#define N 100000
#define K 100
#define D 1000
#define THRESHOLD 0.000001
#define MIN -1
#define MAX 1

void vectorsInitialization(double ** vectors);
void centroidsInitialization(double ** vectors, double ** centroids);
double assignment(double ** vectors, int * clustersMapping, double ** centroids);
void update(double ** vectors, int * clustersMapping, double ** centroids);
double calculateDistanceSquared(double * vector, double * centroid);
void setVector(double * vector,double value);
void accumulateSum(double * A,double * B);
void accumulateQuotient(double * vector,double value);
double randBetween(double min, double max);
double printVector(double * vector);
int * intArray1DInHeap(int dim1);
double * doubleArray1DInHeap(int dim1);
double ** doubleArray2DInHeap(int dim1, int dim2);

#pragma GCC optimize("O3","unroll-loops","omit-frame-pointer","inline","unsafe-math-optimizations")
#pragma GCC option("arch=native","tune=native","no-zero-upper")
int main(){
	double ** vectors=doubleArray2DInHeap(N,D);	//Represents our observations. Contains podoubleers to double arrays of D lenth. Each one of the D doubles is a component of the corresponding vector. e.g: vectors[0][2]: Component 2 of vector 1.
	int * clustersMapping=intArray1DInHeap(N);	//clusterMapping[i] holds the cluster in which vectors[i] is currently included.
	double ** centroids=doubleArray2DInHeap(K,D);	//centroids[i] is the current centroid of cluster i.
	double lossFunction=DBL_MAX;	//Initialized to DOUBLE_MAX so the condition for the while loop will certainly hold for the first iteration
	double lossFunctionNew;
	int iteration;
	
	time(NULL);
	srand(time(NULL));	//Initialize seed for subsequent calls of randBetween
	
	//Core logic=======================================================================================
	vectorsInitialization(vectors);
	centroidsInitialization(vectors,centroids);
	
	iteration=0;
	while( lossFunction-(lossFunctionNew=assignment(vectors,clustersMapping,centroids)) > THRESHOLD){
		lossFunction=lossFunctionNew;
		update(vectors,clustersMapping,centroids);
		printf("Iteration %d, lossFunction=%.8f \n",iteration,lossFunction);
		iteration++;
	}
	//=================================================================================================
	
	return 0;	
}

void vectorsInitialization(double ** vectors){
	int a,b;
	double * vector;	//If vector is declared as an Array, its value cannot be modified, since it is presumed to be constant. This is why it is declared as a podoubleer.
	
	for(a=0; a<N; a++){
		vector=vectors[a];
		for(b=0; b<D; b++){
			vector[b]=randBetween(MIN,MAX);
		}
	}
	
}

void centroidsInitialization(double ** vectors, double ** centroids){
	int a,b;
	double * centroid;
	int chosenInitialVector=rand()%N;	//A number between 0 and N-1 is chosen at random. Then the centroids are initialized as follows: centroid[0]=vectors[number], centroid[1]=vectors[number+1],... .This way no two centroids can be similar SUPPOSING VECTORS ARE DIFFERENT.
	int vectorToAssign;
	
	for(a=0; a<K; a++){
		centroid=centroids[a];
		vectorToAssign=(chosenInitialVector+a)%N;
		for(b=0; b<D; b++){
			centroid[b]=vectors[vectorToAssign][b];
		}
	}
	
}

//Assign each vector to the closest centroid
double assignment(double ** vectors, int * clustersMapping, double ** centroids){
	int a,b;
	double * vector;
	double distance, nearestDistance;
	int nearestCentroid;
	double lossFunction=0;
	
	for(a=0; a<N; a++){	//For each vector
		vector=vectors[a];
		//Start by matching vector with cluster 0
		nearestDistance=calculateDistanceSquared(vector,centroids[0]);
		nearestCentroid=0;
		for(b=1; b<K; b++){	//For each of the remaining centroids
			distance=calculateDistanceSquared(vector,centroids[b]);
			if(distance<nearestDistance){
				nearestDistance=distance;
				nearestCentroid=b;
			}
		}
		clustersMapping[a]=nearestCentroid;
		lossFunction+=nearestDistance;
	}
	
	return lossFunction;
	
}

//Update the centroids according to the newly formed clusters
void update(double ** vectors, int * clustersMapping, double ** centroids){
	int a,b;
	double * centroid;
	int members;
	
	for(a=0; a<K; a++){	//For each centroid/cluster
		centroid=centroids[a];
		setVector(centroid,0);
		members=0;
		
		for(b=0; b<N; b++){	//Search all the vectors to see which are included in the cluster
			if(clustersMapping[b]==a){
				accumulateSum(centroid,vectors[b]);
				members++;
			}
		}
		
		accumulateQuotient(centroid,members);
	}
	
}

double calculateDistanceSquared(double * vector, double * centroid){
	register int a;
	register double distance=0;
	register double dif;
	
	for(a=0; a<D; a++){
		dif=vector[a]-centroid[a];
		distance+=dif*dif;
	}
	
	return distance;
}

//Set all vector components to the same value
void setVector(double * vector,double value){
	int a;
	
	for(a=0; a<D; a++){
		vector[a]=value;
	}
}

void accumulateSum(double * A,double * B){
	int a=0;
	
	for(a=0; a<D; a++){
		A[a]=A[a]+B[a];
	}
	
}

void accumulateQuotient(double * vector,double value){
	int a;
	
	for(a=0; a<D; a++){
		vector[a]=vector[a]/value;
	}
	
}

double randBetween(double min, double max){
	return min+((double) rand()/RAND_MAX)*(max-min);	//(rand()/RAND_MAX) returns a random number between 0 and 1. Obviously, the random numbers we get this way do not cover the full range of double values.
}

double printVector(double * vector){
	int a;
	
	printf("(");
	for(a=0; a<D; a++){
		printf("%.8f ",vector[a]);
	}
	printf(")\n");
	
}

//Returns a pointer which points to a int 1D Array in the Heap.
//The caller is responsible for freeing the memory when it is done with it.
int * intArray1DInHeap(int dim1){
	int * pointer;
	
	pointer=(int *) malloc(dim1*sizeof(int));
	
	return pointer;
}

//Returns a pointer which points to a double 1D Array in the Heap.
//The caller is responsible for freeing the memory when it is done with it.
double * doubleArray1DInHeap(int dim1){
	double * pointer;
	
	pointer=(double *) malloc(dim1*sizeof(double));
	
	return pointer;
}

//Returns a pointer which points to a double 2D Array in the Heap.
//The caller is responsible for freeing the memory when it is done with it.
double ** doubleArray2DInHeap(int dim1, int dim2){
	int a;
	double ** pointer;
	
	pointer=(double**) malloc(dim1*sizeof(double *));	//Allocate space in the heap for dim1 pointers to double
	for(a=0; a<dim1; a++){
		pointer[a]=doubleArray1DInHeap(dim2);	//Allocate space in the heap for dim2 double values for each of the dim1 pointers
	}
	
	return pointer;
}
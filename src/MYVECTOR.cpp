#include <iostream>
#include "MYVECTOR.hpp"

myVector::myVector(int sz){
    size = sz;
    data = new double[size];

    for(int i=0; i<size; i++)
	data[i] = 0.0;
}

myVector::~myVector(){
    size = 0;
    delete[] data;
    data = NULL;
}


void myVector::Print(){
    
    for(int i=0; i<size; i++){
	printf(" Vector: %d %f\n", i, data[i]);
    }
}


myMatrix::myMatrix(int rows, int cols){

    nrows = rows;
    ncols = cols;

    data    = new double*[nrows];       //allocate double* indexing array
    data[0] = new double [nrows*ncols]; //Allocate as a contiguous block
    for(int i=1; i<nrows; i++)
	data[i] = new double [ncols];

    //Initialize
    for(int i=0; i<nrows; i++)
	for(int j=0; j<ncols; j++)
	 data[i][j] = 0;

}


void myMatrix::Print(){
    printf(" Matrix %d %d\n", nrows, ncols);
    for(int i=0; i<nrows; i++){
	for(int j=0; j<ncols; j++){
	    printf(" Matrix %d %d %f\n", i, j, data[i][j]);
	}
    }
}


void myMatrix::Initialize(double **v){
    
    for(int i=0; i<nrows; i++){
	for(int j=0; j<ncols; j++){
	    v[i][j]= 22.1;
	}
    }
}

myMatrix::~myMatrix(){
    
    delete[] data[0];
    delete[] data;
    data = NULL;
}

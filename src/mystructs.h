#ifndef MYSTRUCTS_H
#define MYSTRUCTS_H

/*
 * This file contains all the structs definitions
 * S. Marras, April 2014
 */
typedef struct {
    double *coords;
    double *weights;
    size_t size;
} st_lgl;

typedef struct  {
    double legendre;
    double dlegendre;
    double q;         //q  = legendre(p+1)  - legendre(p-1)
    double dq;        //dq = dlegendre(p+1) - dlegendre(p-1)
} st_legendre;

typedef struct {
    double *data;
    size_t size;
} st_vector;
 

typedef struct Vector {
    double *data;    // no place for x and n in readable code :-)
    size_t size;
} Vector;

#endif

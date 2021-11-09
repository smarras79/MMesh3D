#ifndef MESH_H
#define MESH_H
/*
 * This file contains all the structs for mesh generation
 * S. Marras, Nov 2021
 */
typedef struct {    
    /* 
     * 1D elements of order 'p' used to then build the 3D elements via tensor-product
     * ksi=eta=zeta EXCEPT for the case when the polynomial order 'p'
     * is different in each element direction 
     */
    double *ksi;     //ksi  coords of lgl=p+1 points in the logical elements ksi  = [-1,1]
    double *eta;     //eta  coords of lgl=p+1 points in the logical elements eta  = [-1,1]
    double *zeta;    //zeta coords of lgl=p+1 points in the logical elements zeta = [-1,1]

    //LGL or GL weights
    double *weights;
    
    /*
     * # lgl points in each element direction.
     *  WARNING: These may be different only for conforming cartesian grids
     */
    size_t size;
    size_t size_ksi;  
    size_t size_eta;
    size_t size_zeta;
    
} st_lgl;


typedef struct {
    /*
     * Algorightm 126 of Kopriva (adapted)
     */

    double **coords;   //coords[npoints][nsd]
    int    **conn;     //conn[nelem][ngl^3]
    int     *elements; //Array of elements

    int    **ElCornerNodes; //element connectivity of corner nodes: ElCornerNodes[nelem][2^{nsd}]
    int    **ElEdges;       //element connectivity of edges:        ElEdges[nelem][]
    size_t nelem;

    
    
} st_mesh;


/*typedef struct {
    /*
     * Algorithm 144 (struct version only.)
     * THe procedures are defined separately for lack of Classes in C
     *
    int *head;
    int *tail;
    int *current;

} st_LinkedList;
*/

typedef struct {

    int    id;         //element id: 1...,nelem
    int    *conn;      //conn[nelem][ngl^3]
    int     *element;  //nodes[nelm][4] --> global id's of the corner nodes of each element
    
    size_t npoints;
    size_t nelem;
    
} st_element;


//struct st_Record {
typedef struct {
    /*
     * Algorithm 143
     * THe procedures are defined separately for lack of Classes in C
     */

    int listData;
    struct st_LinkedList *next;
    
} st_Record;

//Procedures on linked lists:
int  Construct(void);
int  Add(st_Record *newRecord, int data);
void PrintList(st_Record *newRecord);


#endif

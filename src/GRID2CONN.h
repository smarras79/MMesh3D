#ifndef GRID2CONN_H
#define GRID2CONN_H

int GRID2CONN(unsigned int strting_node_number, unsigned int strting_elem_number, 
			  int nnodesx, int nnodesy, int **CONN, int *ELTYPE, int iblock);

int GRID2CONNtri(unsigned int strting_node_number, unsigned int strting_elem_number, 
                 int nnodesx, int nnodesy, int nelem, int **CONN, int *ELTYPE, int iblock);

int GRID2CONNhex(unsigned int strting_node_number, unsigned int strting_elem_number, 
		 int nnodesx, int nnodesy, int nnodesz, int **CONN, int *ELTYPE, int iblock);

int GRID2CONNwedge(unsigned int strting_node_number, unsigned int strting_elem_number, 
		   int nnodesx, int nnodesy, int nnodesz, int **CONN, int *ELTYPE, int iblock);

#endif

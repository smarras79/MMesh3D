#ifndef BUILD_CONN_H
#define BUILD_CONN_H

int BUILD_CONN(void);

int CGNS_ORDERING(int **CONN, int nelem);

int BUILD_EDGES(int **CONN, int nelem);

int ADD_HIGH_ORDER_NODES(void);

#endif

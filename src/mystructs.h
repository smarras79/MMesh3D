/*
 * This file contains all the structs definitions
 * S. Marras, April 2014
 */


/* GRID_INFO_T */
typedef struct
{
  MPI_Comm comm;       /* Communicator for entire grid */
  MPI_Comm row_com;    /* Communicator for my row      */
  MPI_Comm col_comm;   /* Communicator for my col      */
  
  int      nprocs;     /* Total number of processes    */
  int      q;          /* Order of the grid            */
  int      my_row;     /* My row number                */
  int      my_col;     /* My column number             */
  int      my_rank;    /* My rank in the grid comm     */

} GRID_INFO_T;


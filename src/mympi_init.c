/*
 * This file contains all the necessary functions to initialize
 * what is needed by MPI 
 *
 * Simone Marras, NPS April 2014
 */
#include "myinclude.h"

MPI_Status MYMPI_INIT(int argc, char *argv[])
{
  int my_rank;
  int size;            /* number of processes  */
  int source;          /* rank of sender       */
  int dest;            /* rank of receiver     */
  int tag = 0;         /* tag for messages     */
  char message[100];   /* storage for messages */
  
  MPI_Status status;   /* return status for receive */
  MPI_Request request; /* request for Isend */
  

  /* Get process rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
  /* Get number of processes */
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  //printf(" RANK: %d %d \n", my_rank, size);
  
  if(my_rank != 0)
    {
      sprintf(message, "Greetings from process %d ! \n", my_rank);
      dest = 0;
      
      /* Use strlen+1 so that  '\0' gets transmitted */
      MPI_Isend(message, strlen(message)+1, MPI_CHAR, dest, 999, MPI_COMM_WORLD, &request);
    }
  else
    {
      for(source = 1; source < size; source++)
	{
	  MPI_Recv(message, 100, MPI_CHAR, source, 999, MPI_COMM_WORLD, &status);
	  printf(" %s \n", message);
	}
    }
  
  return status;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "MESH.h"
#include "PRINT.h"

#define P4EST_CHILDREN 4
#define P4_TO_P8 1

/*
 * Procedure of Algorithm 145
 */
int Construct(void) {

    st_Record *head    = NULL;
    st_Record *tail    = NULL;
    st_Record *current = NULL;

    head    = (st_Record *) malloc(sizeof(st_Record));
    if (head == NULL) {
	PRINT_ERROR("HEAD struct not allocated");
	return 1;
    }
    
    tail    = (st_Record *) malloc(sizeof(st_Record));
    if (tail == NULL) {
	PRINT_ERROR(" TAIL struct not allocated");
	return 1;
    }
    
    current = (st_Record *) malloc(sizeof(st_Record));
    if (current == NULL) {
	PRINT_ERROR("CURRENT struct not allocated");
	return 1;
    }
    
    return 0;
}


int Add(st_Record *newRecord, int data) {

    st_Record *head    = NULL;
    st_Record *tail    = NULL;
    
    //Allocate new record to be added to the list:
    newRecord = (st_Record *) malloc(sizeof(st_Record));
    if (newRecord == NULL) {
	PRINT_ERROR("newRecord not allocated");
	return 1;
    }
    newRecord->next = NULL;

    if(tail == NULL) {
	head = newRecord;
	tail = newRecord;
    } else {
	tail->next = newRecord;
	tail = newRecord;
    }
    newRecord = tail;
    newRecord->next = NULL;
    newRecord->listData = data;
    //printf(" -- %d) \n",newRecord->listData);

    return 0;
}

/*
  struct st_Record MoveToNext(struct st_Record *current) {   
       	
  current = current->next;

  return current;
  }*/


void PrintList(st_Record *newRecord) {
    
    int error_flag;

    st_Record *head;
    st_Record *current;
   
    current = newRecord;
    /*    
	  head->listData = 12;
	  head->next = NULL;
    */
    
    printf("\n[ ");
   	
    //start from the beginning
    while(current != NULL) {
	printf("(%d) \n",newRecord->listData);
	//MoveToNext();
	current = current->next;
    }
	
    printf(" ]");
}

/*
char* READ_GMSH(FILE *stream)
{

    char *line = (char*) malloc(1024 * sizeof(char *)), *linep = line;
    size_t lenmax = 1024, len = lenmax;
    int c;

    if (line == NULL)
	return NULL;

    for (;;)
	{
	    c = fgetc(stream);
	    c = toupper(c);
	    if (c == EOF && linep == line)
		{
		    free(linep);
		    return NULL;
		}

	    if (--len == 0)
		{
		    len = lenmax;
		    lenmax *= 2;

		    char *linen = realloc(linep,lenmax*sizeof(char*));
		    //char *linen = P4EST_REALLOC(linep, char, lenmax);

		    if (linen == NULL)
			{
			    free(linep);
			    return NULL;
			}
		    line = linen + (line - linep);
		    linep = linen;
		}
	    if ((*line++ = c) == '\n')
		break;
	}
    *line = '\0';
    return linep;
    }*/

/* This function reads the mesh file and fills bc info*
static int read_inp_stream(FILE *stream,int *num_bc,
			   int *num_elem,
			   int *bc_physical_type,
			   int *bc_el_label,
			   int *bc_to_vertex,
			   int *vc_el_type,
			   int *vc_el_label)
*
int read_inp_stream(FILE *stream, int *num_bc, int *num_elem)
{
    int reading_elements = 0;
    int reading_elsets = 0;
    int reading_elsets_volume = 0;
    int reading_elements_volume = 0;
    int wrong_element = 0;
    int lines_read = 0, lines_free = 0;
    char *line;
    int num_physical_el = 0;
    int num_elements = 0;
    int num_elsets = 0;
    int num_elsets_volume = 0;
    int num_elements_volume = 0;
    //int fill_bc_data = (bc_physical_type != NULL && bc_to_vertex != NULL);
    int boundary_condition = 0;
    int *physical_type;
    int *element_label;
    int *vc_el_label_list;
    int ibc = 0;
    int ivc = 0;

     P4EST_ASSERT((bc_physical_type == NULL && bc_to_vertex == NULL) ||
		 (bc_physical_type != NULL && bc_to_vertex != NULL));

    if (fill_bc_data)
	{
	    physical_type = malloc(*num_bc * sizeof(int));
	    element_label = malloc(*num_bc * sizeof(int));
	    vc_el_label_list = malloc(*num_elem * sizeof(int));
	}
    
    for (;;)
	{
	    line = (stream);

	    if (line == NULL)
		{
		    break;
		}

	    ++lines_read;

	    // check for control line 
	    if (line[0] == '*')
		{
		    reading_elements = reading_elsets = reading_elsets_volume =
			reading_elements_volume = wrong_element = 0;
		    if (strstr(line, "*ELEMENT"))
			{

			    if (
#ifdef P4_TO_P8
				strstr(line, "TYPE=C2D4") || strstr(line, "TYPE=CPS4") ||
				strstr(line, "TYPE=S4") ||
				strstr(line, "TYPE=CPS3") // quad elements, mostly CPS4
#else
				strstr(line, "TYPE=T3D2")
#endif
				)
				{
				    if (strstr(line, "TYPE=CPS3")) // we don't want triangles
					wrong_element = 1;
				    reading_elements = 1;
				    ++lines_free;
				    free(line);
				    continue;
				}
			    else if (strstr(line, "TYPE=C3D8")) // just counting 3D elements
				{
				    reading_elements_volume = 1;
				    ++lines_free;
				    free(line);
				    continue;
				}
			}
		    else if (strstr(line, "*ELSET") &&
			     strstr(line, ":BC_")) // boundary conditions (faces)
			{

			    // find which bc is prescribed
			    boundary_condition = which_bc(line);
			    P4EST_ASSERT(boundary_condition > 0);

			    reading_elsets = 1;
			    ++lines_free;
			    free(line);
			    continue;
			}
		    else if (strstr(line, "*ELSET") &&
			     strstr(line, ":VC_")) // volume condition (sponge)
			{
			    // find which vc is prescribed
			    boundary_condition = which_bc(line);
			    P4EST_ASSERT(boundary_condition > 0);
			    reading_elsets_volume = 1;
			    ++lines_free;
			    free(line);
			    continue;
			}
		    
		}

	    if (reading_elements)
		{
		    if (fill_bc_data)
			{
			    long long int v[P4EST_CHILDREN / 2];
			    long long int n;
			    long long int el_num;
			    int retval;

			    if (num_elements >= *num_bc)
				{
				    PRINT_ERROR("Encountered element that will not fit into"
						 " bc_to_vertex array. More elements than expected.\n");
				    free(line);
				    return 1;
				}

			    if (wrong_element)
				{
				    retval = sscanf(line, "%lld, %lld, %lld, %lld", &el_num, &v[0], &v[1],
						    &v[2]);

				    if (retval != 4)
					{
					    PRINT_ERROR("Premature end of file");
					    free(line);
					    return 1;
					}

				    bc_el_label[num_elements] = el_num;
				    for (n = 0; n < 3; ++n)
					bc_to_vertex[P4EST_CHILDREN / 2 * num_elements + n] = v[n] - 1;
				    bc_to_vertex[P4EST_CHILDREN / 2 * num_elements + 3] = -1;
				}
			    else
				{
				    retval = sscanf(line, "%lld, %lld, %lld"
#ifdef P4_TO_P8
						    ", %lld, %lld"
#endif
						    ,
						    &el_num, &v[0], &v[1]
#ifdef P4_TO_P8
						    ,
						    &v[2], &v[3]
#endif
						    );

				    if (retval != P4EST_CHILDREN / 2 + 1)
					{
					    PRINT_ERROR("Premature end of file");
					    free(line);
					    return 1;
					}

				    bc_el_label[num_elements] = el_num;
				    for (n = 0; n < P4EST_CHILDREN / 2; ++n)
					bc_to_vertex[P4EST_CHILDREN / 2 * num_elements + n] = v[n] - 1;
				}
				}
		    ++num_elements;
		}
	    else if (reading_elsets)
		{
		     if (fill_bc_data)
			{

			    char *p = line;
			    while (*p)
				{
				    if (isdigit(*p))
					{
					    long val = strtol(p, &p, 10);
					    element_label[ibc] = val;
					    physical_type[ibc] = boundary_condition;
					    ibc++;
					}
				    else
					{
					    p++;
					}
				}

			    if (ibc > *num_bc)
				{
				    PRINT_ERROR(
						 "Encountered bc element that will not fit in bc array\n");
				    free(line);
				    return 1;
				}
			}

		    ++num_elsets;
		}
	    else if (reading_elsets_volume)
		{
		    if (fill_bc_data)
			{

			    char *p = line;
			    while (*p)
				{
				    if (isdigit(*p))
					{
					    long val = strtol(p, &p, 10);
					    vc_el_label[ivc] = val;
					    vc_el_type[ivc] = boundary_condition;
					    ivc++;
					}
				    else
					{
					    p++;
					}
				}

			    if (ivc > *num_elem)
				{
				    PRINT_ERROR(
						 "Encountered vc element that will not fit in vc array\n");
				    free(line);
				    return 1;
				}
			}
		    
		    ++num_elsets_volume;
		}
	    else if (reading_elements_volume)
		{
		    if (fill_bc_data)
			{
			    long long int v[P4EST_CHILDREN];
			    long long int el_num;
			    int retval;

			    retval = sscanf(line, "%lld, %lld, %lld, %lld, %lld"
#ifdef P4_TO_P8
					    ", %lld, %lld, %lld, %lld"
#endif
					    ,
					    &el_num, &v[0], &v[1], &v[2], &v[3]
#ifdef P4_TO_P8
					    ,
					    &v[4], &v[5], &v[6], &v[7]
#endif
					    );

			    if (retval != P4EST_CHILDREN + 1)
				{
				    PRINT_ERROR("Premature end of file");
				    free(line);
				    return 1;
				}

			    vc_el_label_list[num_elements_volume] = el_num;
			}

		    ++num_elements_volume;
		}

	    ++lines_free;
	    free(line);
	}

    if (fill_bc_data)
	{
	    for (int i = 0; i < *num_bc; i++)
		{
		    for (int j = 0; j < *num_bc; j++)
			{
			    if (element_label[i] == bc_el_label[j])
				{
				    bc_physical_type[j] = physical_type[i];
				    continue;
				}
			}
		}

	    for (int i = 0; i < ivc; i++) // could be shorter loop
		{
		    for (int j = 0; j < *num_elem; j++)
			{
			    if (vc_el_label[i] == vc_el_label_list[j])
				{
				    vc_el_label[i] = j;
				    continue;
				}
			}
		}
		}
    *num_bc = num_elements;
    *num_elem = num_elements_volume;

    if (num_elsets == 0 || num_elements == 0)
	{
	    PRINT_ERROR("No elements or bcs found in mesh file.\n");
	    return -1;
	}
    else
	{
	    return 0;
	}
}
*/

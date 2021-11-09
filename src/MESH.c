#include <stdio.h>
#include <stdlib.h>

#include "MESH.h"
#include "PRINT.h"

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

/*
 * Allocate memory
 */
int MEMORY_ALLOCATE(int flag);
void lgl_allocate(int flag, st_lgl lgl, size_t sz);


/*
 * Free memory
 */
int MEMORY_DEALLOCATE(int flag);
void lgl_deallocate(int flag, st_lgl lgl, size_t sz);


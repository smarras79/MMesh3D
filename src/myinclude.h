#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include<errno.h>
#include<time.h>

#include<mpi.h>

#include "mystructs.h"

#include "visit_writer.h"
#include "nrutil.h"
#include "PRINT.h"
#include "READ_INPUT.h"
#include "BUILD_GRID.h"
#include "BUILD_CONN.h"
#include "GRID_COORD.h"
#include "GRID2CONN.h"
#include "linspace.h"
#include "parabola.h"
#include "gaussj.h"
#include "INTERPOLATE.h"
#include "minmaxval.h"
#include "topo_user_function.h"
#include "SURFACES.h"
#include "TOPOfromTXT.h"
#include "BUILD_LGL.h"
#include "MEMORY.h"
#include "elliptic_solver.h"
#include "mympi_init.h"
#include "DOMAIN_DECOMP.h"





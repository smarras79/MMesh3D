/*
 * Here the parameters are defined
 */

#define nsd                   3

#define HEXA                  8
#define NELFACES              6
#define WEDGE                 6

#define PROB_ENTRIES         14
#define MAX_INPUTS           20
#define NBDY_NODES            8
#define NUMBER_OF_PARAMETERS 11

#define PI                    3.1415926535897932

#define ticks_per_ms          CLOCKS_PER_SEC/1000

#define VISIT_VERTEX         1
#define VISIT_LINE           3
#define VISIT_TRIANGLE       5
#define VISIT_QUAD           9
#define VISIT_TETRA         10
#define VISIT_HEXAHEDRON    12
#define VISIT_WEDGE         13
#define VISIT_PYRAMID       14

/* Physical quantities */
#define earth_radius        6.37122e+6
#define earth_gravity       9.81

#define max(a,b)				\
    ({ __typeof__ (a) _a = (a);			\
	__typeof__ (b) _b = (b);		\
	_a > _b ? _a : _b; })

#define min(a,b)				\
    ({ __typeof__ (a) _a = (a);			\
	__typeof__ (b) _b = (b);		\
	_a < _b ? _a : _b; })


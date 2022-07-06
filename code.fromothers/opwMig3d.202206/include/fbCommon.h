/* fbCommon.h - common include files */

#ifndef FB_COMMON_H
#define FB_COMMON_H

/* INCLUDES */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<limits.h>
#include<float.h>
#include<math.h>

// OpenMP
#include<omp.h>

/* DEFINES */
#ifndef PI
#define PI 3.14159265358979
#endif

#ifndef FLOAT_SIZE_BYTES
#define FLOAT_SIZE_BYTES 4
#endif

#ifndef LONG_SIZE_BYTES
#define LONG_SIZE_BYTES 8
#endif

#ifndef FILE_NAME_MAX_LENGTH
#define FILE_NAME_MAX_LENGTH 4096
#endif

#ifndef FILE_NUMBER_MAX
#define FILE_NUMBER_MAX 256
#endif

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

#ifndef R_table_Max
#define R_table_Max 2
#endif

#ifndef S_table_Max
#define S_table_Max 2
#endif

#include "fballoc.h"
#include "fbsegy_std.h"
#include "fbrw.h"
#include "fbsurw.h"

extern FILE     *cmpGatherFilefp;   // Global variable: file pointer for operating the cmp-gather file.
extern float    Coor_S_R_Scale;     // Global variable: sx = sx/Coor_S_R_Scale; gx = gx/Coor_S_R_Scale

#endif /* FB_COMMON_H */

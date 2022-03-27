#ifndef PIHM_HEADER
#define PIHM_HEADER

#include <stdio.h>
#include <stdlib.h>

#include <math.h>   

#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <errno.h>
#include <ctype.h>
#if defined(_WIN32) || defined(_WIN64)
# include <direct.h>
# include <io.h>
# include <windows.h>
#else
# include <unistd.h>
#endif
#include <stdarg.h>
#if defined(_OPENMP)
# include <omp.h>
#endif

#define VERSION  "0.10.1-alpha"
//#define _SUNDIALS_v3 1
//#define _SUNDIALS_v4 0
/*
 * SUNDIAL Header Files
 */
#if defined(_SUNDIALS_v3)
#include "cvode/cvode.h" 
#elif defined(_SUNDIALS_v4)
#include "cvode/cvode.h" 
#else   /* sundials v2.7*/
#include "cvode.h"
#endif

#if defined(_SUNDIALS_v3)
#include "sunlinsol/sunlinsol_spgmr.h"
#include "cvode/cvode_spils.h"  
#elif defined(_SUNDIALS_v4)
#include "sunlinsol/sunlinsol_spgmr.h"
#include "cvode/cvode_spils.h"
#else
/* CVSPGMR linear header file */
#include "cvode_spgmr.h" 
#endif

/* Definition of type N_Vector */
#if defined(_CVODE_OMP)
# include "nvector/nvector_openmp.h"
#else
# include "nvector/nvector_serial.h"
#endif

#if defined(_SUNDIALS_v3)
#include "sundials/sundials_types.h"  
#include "sundials/sundials_math.h"  
#elif defined(_SUNDIALS_v4)
#include "sundials/sundials_types.h"  
#include "sundials/sundials_math.h"  
#else
/* UnitRoundoff, RSqrt, SQR functions */
#include "sundials/sundials_math.h"
#endif

/* CVDENSE header file */
#if defined(_SUNDIALS_v3)
#include "sunlinsol/sunlinsol_dense.h"
#elif defined(_SUNDIALS_v4)
#include "sunlinsol/sunlinsol_dense.h"
#else
#include "cvode_dense.h"   
#endif

#if defined(_NOAH_)
# include "spa.h"
#endif

#include "custom_io.h"

#include "pihm_const.h"
#include "pihm_input_struct.h"
#include "elem_struct.h"
#include "river_struct.h"
#include "pihm_struct.h"
#include "pihm_func.h"

#endif
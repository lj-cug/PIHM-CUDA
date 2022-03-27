#ifndef PIHM_STRUCT_HEADER_D
#define PIHM_STRUCT_HEADER_D

/* Global calibration coefficients in GPU, only used in is_sm_et.cu */
typedef struct calib_struct_device
{
// 仅在IntcpSnowEt_kernel中用到了下面3个calibration参数
    realtype  ec;
    realtype  ett;
    realtype  edir;

} calib_struct_d;


/* Model control parameters in GPU, only used in hydrol_elem.cu */
typedef struct ctrl_struct_device
{
// 仅在 ChanFlowRiverToRiver核函数中使用以下3个ctrl参数.									
    int surf_mode;  /* surface overland flow formulation
                     * 1 = kinematic, 2 = diffusion */
    int riv_mode;   /* river routing formulation:
                     * 1 = kinematic, 2 = diffusion */
    int stepsize;   /* model step size (s) */
	
} ctrl_struct_d;                                           

// pihm_d作为用户数据的结构体
struct pihm_struct_device
{ 
    int    nelem;        /* Number of Element*/
    int    nriver;       /* Number of River Segments */	
	int    NEQ;          /* Number of equations (state variables) */	
	
//    calib_struct_d  cal_d;     /* Calibration Parameters struct */
//    ctrl_struct_d   ctrl_d;	 /* Contrl Parameters struct */
// 率定和控制参数也直接包含在pihm_struct_device结构体当中
    realtype  ec;
    realtype  ett;
    realtype  edir;
	
    int surf_mode;  /* surface overland flow formulation
                     * 1 = kinematic, 2 = diffusion */
    int riv_mode;   /* river routing formulation:
                     * 1 = kinematic, 2 = diffusion */
    int stepsize;   /* model step size (s) */

// 水文过程计算中, elem_d和river_d的数组, 仅包含1D数组或指针数组(such as Noah LSM etc.)
#include "pihm_device_array.c"

};
typedef pihm_struct_device *pihm_struct_d;    // UserData struct in GPU

#endif
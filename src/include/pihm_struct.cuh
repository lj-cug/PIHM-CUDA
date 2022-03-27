#ifndef PIHM_STRUCT_HEADER_D
#define PIHM_STRUCT_HEADER_D

/* Global calibration coefficients in GPU, only used in is_sm_et.cu */
typedef struct calib_struct_device
{
// ����IntcpSnowEt_kernel���õ�������3��calibration����
    realtype  ec;
    realtype  ett;
    realtype  edir;

} calib_struct_d;


/* Model control parameters in GPU, only used in hydrol_elem.cu */
typedef struct ctrl_struct_device
{
// ���� ChanFlowRiverToRiver�˺�����ʹ������3��ctrl����.									
    int surf_mode;  /* surface overland flow formulation
                     * 1 = kinematic, 2 = diffusion */
    int riv_mode;   /* river routing formulation:
                     * 1 = kinematic, 2 = diffusion */
    int stepsize;   /* model step size (s) */
	
} ctrl_struct_d;                                           

// pihm_d��Ϊ�û����ݵĽṹ��
struct pihm_struct_device
{ 
    int    nelem;        /* Number of Element*/
    int    nriver;       /* Number of River Segments */	
	int    NEQ;          /* Number of equations (state variables) */	
	
//    calib_struct_d  cal_d;     /* Calibration Parameters struct */
//    ctrl_struct_d   ctrl_d;	 /* Contrl Parameters struct */
// �ʶ��Ϳ��Ʋ���Ҳֱ�Ӱ�����pihm_struct_device�ṹ�嵱��
    realtype  ec;
    realtype  ett;
    realtype  edir;
	
    int surf_mode;  /* surface overland flow formulation
                     * 1 = kinematic, 2 = diffusion */
    int riv_mode;   /* river routing formulation:
                     * 1 = kinematic, 2 = diffusion */
    int stepsize;   /* model step size (s) */

// ˮ�Ĺ��̼�����, elem_d��river_d������, ������1D�����ָ������(such as Noah LSM etc.)
#include "pihm_device_array.c"

};
typedef pihm_struct_device *pihm_struct_d;    // UserData struct in GPU

#endif
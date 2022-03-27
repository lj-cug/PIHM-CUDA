#include "pihm.h"
#include "pihm.cuh"

#define max(a,b)(a>b?a:b)
#define min(a,b)(a<b?a:b)
/*
   2021.04.14: We have the problem of race condition in global memory.
   2021.04.15: 核函数的形参大小以及个数的限制
   __global__ function parameters are passed to the device via constant memory and are
   limited to 4 KB.
   2021.04.16: Check NaN:  isnan()   MonthlyLai_cuda等需要定义为__constant__变量，不能是__shared__
   2021.04.18: 与CPU代码对比，确定了void IntcpSnowEt()和void IntcpSnowEt_kernel计算的正确性.
*/


// 简单的截留-融雪-蒸发模块
#if !defined(_NOAH_)
void IntcpSnowEt(int t, int stepsize0, const pihm_struct_d pihm_d)
{
	int nelem = pihm_d->nelem;
	realtype stepsize = (realtype)stepsize0;

	// 初始化3个率定参数
	realtype cal_ec = pihm_d->ec;
	realtype cal_ett = pihm_d->ett;
	realtype cal_edir = pihm_d->edir;

	// pihm_d中的数组首地址赋值给elem_d和river_d的数组
#include "IntcpSnowEt_device_injections.c" 

	int numThreads_elem = min(32, nelem);
	int numBlocks_elem = (nelem + numThreads_elem - 1) / numThreads_elem;

	pihm_t_struct   pihm_time;

	pihm_time = PIHMTime(t);      // 时间t(second)转换为PIHMtime结构体格式的时间
	int month = pihm_time.month;  // get the month used in searching for LAI etc.

	// The following kernel calling has problem of an illegal access was encountered. Check it!
	IntcpSnowEt_kernel << <numThreads_elem, numBlocks_elem >> >
	(month, stepsize, nelem,
	cal_ec, cal_ett, cal_edir,                 // 3个率定参数
	#include "IntcpSnowEt_kernel_arguments.c"  // 类似 elem_d_wf_pcpdrp  
	);

	cudaDeviceSynchronize();

#if 0
	// check for error
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess)
	{
		// print the CUDA error message and exit
		printf("CUDA error in ODE: %s\n", cudaGetErrorString(error));
		exit(-1);
	}
#endif

}

#endif

__constant__ realtype    TSNOW = -3.0;
__constant__ realtype    TRAIN = 1.0;
__constant__ realtype    T0 = 0.0;

__constant__ realtype MonthlyLai_cuda[40][12] = {
	/* Evergreen Needleleaf Forest */
	{ 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0 },
	/* Evergreen Broadleaf Forest */
	{ 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5 },
	/* Deciduous Needleleaf Forest */
	{ 0.0, 0.0, 0.0, 0.6, 1.2, 2.0, 2.6, 1.7, 1.0, 0.5, 0.2, 0.0 },
	/* Deciduous Broadleaf Forest */
	{ 0.0, 0.0, 0.3, 1.2, 3.0, 4.7, 4.5, 3.4, 1.2, 0.3, 0.0, 0.0 },
	/* Mixed Forest */
	{ 2.0, 2.0, 2.2, 2.6, 3.5, 4.3, 4.3, 3.7, 2.6, 2.2, 2.0, 2.0 },
	/* Closed Shrubland */
	{ 0.0, 0.0, 0.3, 0.9, 2.2, 3.5, 3.5, 2.5, 0.9, 0.3, 0.0, 0.0 },
	/* Open Shrubland */
	{ 0.0, 0.0, 0.2, 0.6, 1.5, 2.3, 2.3, 1.7, 0.6, 0.2, 0.0, 0.0 },
	/* Woody Savanna */
	{ 0.2, 0.2, 0.4, 1.0, 2.4, 4.1, 4.1, 2.7, 1.0, 0.4, 0.2, 0.2 },
	/* Savanna */
	{ 0.3, 0.3, 0.5, 0.8, 1.8, 3.6, 3.8, 2.1, 0.9, 0.5, 0.3, 0.3 },
	/* Grassland */
	{ 0.4, 0.5, 0.6, 0.7, 1.2, 3.0, 3.5, 1.5, 0.7, 0.6, 0.5, 0.4 },
	/* Permanent Wetland */
	{ 0.2, 0.3, 0.3, 0.5, 1.5, 2.9, 3.5, 2.7, 1.2, 0.3, 0.3, 0.2 },
	/* Cropland */
	{ 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 1.5, 0.0, 0.0, 0.0 },
	/* Urban and Built-Up */
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	/* Cropland/Natural Veg. Mosaic */
	{ 0.2, 0.3, 0.3, 0.4, 1.1, 2.5, 3.2, 2.2, 1.1, 0.3, 0.3, 0.2 },
	/* Permanent Snow */
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	/* Barren/Sparsely Vegetated */
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	/* IGBP Water */
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	/* Unclassified */
	{ 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999 },
	/* Fill Value */
	{ 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999 },
	/* Unclassified */
	{ 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999 },
	/* Open Water */
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	/* Perennial Ice/Snow */
	{ 1.06, 1.11, 1.20, 1.29, 1.45, 1.67, 1.71, 1.63, 1.45, 1.27, 1.13,
	1.06 },
	/* Developed Open Space */
	{ 0.83, 0.94, 1.06, 1.18, 1.86, 3.72, 4.81, 4.26, 2.09, 1.29, 1.06,
	0.94 },
	/* Developed Low Intensity */
	{ 0.96, 1.07, 1.20, 1.35, 2.04, 3.83, 4.87, 4.35, 2.30, 1.46, 1.18,
	1.06 },
	/* Developed Medium Intensity */
	{ 1.11, 1.22, 1.36, 1.54, 2.26, 3.97, 4.94, 4.46, 2.55, 1.66, 1.34,
	1.20 },
	/* Developed High Intensity */
	{ 1.24, 1.34, 1.50, 1.71, 2.45, 4.09, 5.00, 4.54, 2.76, 1.82, 1.47,
	1.32 },
	/* Barren Land */
	{ 0.03, 0.03, 0.03, 0.02, 0.02, 0.03, 0.04, 0.06, 0.09, 0.06, 0.04,
	0.03 },
	/* Deciduous Forest */
	{ 0.62, 0.67, 0.92, 1.71, 3.42, 5.53, 6.22, 5.60, 3.83, 1.79, 0.92,
	0.67 },
	/* Evergreen Forest */
	{ 3.38, 3.43, 3.47, 3.52, 3.78, 4.54, 4.98, 4.76, 3.87, 3.56, 3.47,
	3.43 },
	/* Mixed Forest */
	{ 3.10, 3.26, 3.61, 4.11, 5.17, 6.73, 7.21, 6.71, 5.34, 4.09, 3.41,
	3.14 },
	/* Dwarf Scrub */
	{ 0.24, 0.24, 0.19, 0.13, 0.15, 0.20, 0.26, 0.48, 0.70, 0.48, 0.30,
	0.24 },
	/* Shrub/Scrub */
	{ 0.35, 0.38, 0.38, 0.38, 0.55, 1.06, 1.53, 1.53, 1.04, 0.58, 0.44,
	0.38 },
	/* Grassland/Herbaceous */
	{ 0.70, 0.80, 0.90, 1.00, 1.60, 3.30, 4.30, 3.80, 1.80, 1.10, 0.90,
	0.80 },
	/* Sedge/Herbaceous */
	{ 0.70, 0.80, 0.90, 1.00, 1.60, 3.30, 4.30, 3.80, 1.80, 1.10, 0.90,
	0.80 },
	/* Lichens */
	{ 0.70, 0.80, 0.90, 1.00, 1.60, 3.30, 4.30, 3.80, 1.80, 1.10, 0.90,
	0.80 },
	/* Moss */
	{ 0.70, 0.80, 0.90, 1.00, 1.60, 3.30, 4.30, 3.80, 1.80, 1.10, 0.90,
	0.80 },
	/* Pasture/Hay */
	{ 0.47, 0.54, 0.60, 0.67, 1.07, 2.20, 2.87, 2.54, 1.20, 0.74, 0.60,
	0.54 },
	/* Cultivated Crops */
	{ 0.47, 0.54, 0.60, 0.67, 1.07, 2.20, 2.87, 2.54, 1.20, 0.74, 0.60,
	0.54 },
	/* Woody Wetland */
	{ 0.35, 0.38, 0.38, 0.38, 0.55, 1.06, 1.53, 1.53, 1.04, 0.58, 0.44,
	0.38 },
	/* Emergent Herbaceous Wetland */
	{ 0.24, 0.24, 0.19, 0.13, 0.15, 0.20, 0.26, 0.48, 0.70, 0.48, 0.30,
	0.24 }
};

__constant__ realtype MonthlyRl_cuda[40][12] = {
	/* Evergreen Needleleaf Forest */
	{ 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
	0.500, 0.500 },
	/* Evergreen Broadleaf Forest */
	{ 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
	0.500, 0.500 },
	/* Deciduous Needleleaf Forest */
	{ 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
	0.500, 0.500 },
	/* Deciduous Broadleaf Forest */
	{ 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
	0.500, 0.500 },
	/* Mixed Forest */
	{ 0.200, 0.200, 0.200, 0.200, 0.278, 0.367, 0.367, 0.300, 0.200, 0.200,
	0.200, 0.200 },
	/* Closed Shrubland */
	{ 0.010, 0.010, 0.010, 0.015, 0.032, 0.048, 0.048, 0.035, 0.015, 0.010,
	0.010, 0.010 },
	/* Open Shrubland */
	{ 0.010, 0.010, 0.010, 0.010, 0.033, 0.053, 0.053, 0.038, 0.010, 0.010,
	0.010, 0.010 },
	/* Woody Savanna */
	{ 0.010, 0.010, 0.010, 0.016, 0.034, 0.050, 0.050, 0.038, 0.016, 0.010,
	0.010, 0.010 },
	/* Savanna */
	{ 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150,
	0.150, 0.150 },
	/* Grassland */
	{ 0.100, 0.100, 0.101, 0.102, 0.106, 0.120, 0.120, 0.108, 0.102, 0.101,
	0.100, 0.100 },
	/* Permanent Wetland */
	{ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
	0.000, 0.000 },
	/* Cropland */
	{ 0.050, 0.050, 0.050, 0.050, 0.050, 0.061, 0.085, 0.085, 0.050, 0.050,
	0.050, 0.050 },
	/* Urban and Built-Up */
	{ 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
	0.500, 0.500 },
	/* Cropland/Natural Veg. Mosaic */
	{ 0.050, 0.050, 0.050, 0.050, 0.050, 0.059, 0.091, 0.050, 0.050, 0.050,
	0.050, 0.050 },
	/* Permanent Snow */
	{ 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
	0.001, 0.001 },
	/* Barren/Sparsely Vegetated */
	{ 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
	0.010, 0.010 },
	/* IGBP Water */
	{ 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
	0.0001, 0.0001, 0.0001, 0.0001 },
	/* Unclassified */
	{ 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300,
	0.300, 0.300 },
	/* Fill Value */
	{ 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150,
	0.150, 0.150 },
	/* Unclassified */
	{ 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050,
	0.050, 0.050 },
	/* Open Water */
	{ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
	0.000, 0.000 },
	/* Perennial Ice/Snow */
	{ 0.152, 0.151, 0.155, 0.165, 0.169, 0.169, 0.169, 0.170, 0.170, 0.166,
	0.157, 0.152 },
	/* Developed Open Space */
	{ 0.089, 0.089, 0.091, 0.093, 0.095, 0.094, 0.093, 0.094, 0.095, 0.094,
	0.091, 0.089 },
	/* Developed Low Intensity */
	{ 0.119, 0.119, 0.123, 0.132, 0.137, 0.136, 0.135, 0.136, 0.137, 0.133,
	0.124, 0.119 },
	/* Developed Medium Intensity */
	{ 0.154, 0.153, 0.163, 0.179, 0.187, 0.187, 0.186, 0.187, 0.188, 0.180,
	0.163, 0.154 },
	/* Developed High Intensity */
	{ 0.183, 0.183, 0.195, 0.218, 0.229, 0.229, 0.228, 0.229, 0.230, 0.220,
	0.196, 0.183 },
	/* Barren Land */
	{ 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.014, 0.014,
	0.013, 0.013 },
	/* Deciduous Forest */
	{ 0.343, 0.343, 0.431, 0.577, 0.650, 0.657, 0.656, 0.653, 0.653, 0.581,
	0.431, 0.343 },
	/* Evergreen Forest */
	{ 1.623, 1.623, 1.623, 1.623, 1.623, 1.623, 1.622, 1.622, 1.623, 1.623,
	1.623, 1.623 },
	/* Mixed Forest */
	{ 0.521, 0.518, 0.557, 0.629, 0.663, 0.664, 0.665, 0.665, 0.667, 0.633,
	0.562, 0.521 },
	/* Dwarf Scrub */
	{ 0.022, 0.022, 0.021, 0.020, 0.020, 0.020, 0.020, 0.023, 0.025, 0.024,
	0.023, 0.022 },
	/* Shrub/Scrub */
	{ 0.034, 0.034, 0.033, 0.033, 0.033, 0.032, 0.032, 0.034, 0.035, 0.035,
	0.034, 0.034 },
	/* Grassland/Herbaceous */
	{ 0.070, 0.070, 0.070, 0.070, 0.070, 0.069, 0.068, 0.069, 0.070, 0.070,
	0.070, 0.070 },
	/* Sedge/Herbaceous */
	{ 0.070, 0.070, 0.070, 0.070, 0.070, 0.069, 0.068, 0.069, 0.070, 0.070,
	0.070, 0.070 },
	/* Lichens */
	{ 0.070, 0.070, 0.070, 0.070, 0.070, 0.069, 0.068, 0.069, 0.070, 0.070,
	0.070, 0.070 },
	/* Moss */
	{ 0.070, 0.070, 0.070, 0.070, 0.070, 0.069, 0.068, 0.069, 0.070, 0.070,
	0.070, 0.070 },
	/* Pasture/Hay */
	{ 0.047, 0.047, 0.047, 0.047, 0.047, 0.046, 0.046, 0.046, 0.047, 0.047,
	0.047, 0.047 },
	/* Cultivated Crops */
	{ 0.047, 0.047, 0.047, 0.047, 0.047, 0.046, 0.046, 0.046, 0.047, 0.047,
	0.047, 0.047 },
	/* Woody Wetland */
	{ 0.038, 0.038, 0.038, 0.037, 0.037, 0.037, 0.037, 0.039, 0.040, 0.039,
	0.039, 0.038 },
	/* Emergent Herbaceous Wetland */
	{ 0.027, 0.027, 0.026, 0.024, 0.024, 0.024, 0.024, 0.028, 0.029, 0.029,
	0.027, 0.027 }
};

__constant__ realtype MonthlyMf_cuda[12] = {
	0.001308019, 0.001633298, 0.002131198, 0.002632776, 0.003031171,
	0.003197325, 0.003095839, 0.002745240, 0.002260213, 0.001759481,
	0.001373646, 0.001202083
};

__global__ void IntcpSnowEt_kernel(int month, realtype stepsize,
	          int nelem,
			  realtype cal_ec, realtype cal_ett, realtype cal_edir,
#include "IntcpSnowEt_kernel_declarations.c"  // 类似 int *elem_d_attrib_soil_type 的指针形参
)
{
   int tid = blockDim.x * blockIdx.x + threadIdx.x;
   
   // 临时局部设备变量
   realtype satn;
   realtype betas;
   realtype fr;
   realtype alphar;
   realtype etas;
   realtype gammas;
   realtype rs;
   realtype pc;
   realtype delta;
   realtype gamma;
   realtype radnet;
   realtype sfctmp;
   realtype wind;
   realtype rh;
   realtype vp;
   realtype pres;
   realtype lai;
   realtype z0;
   realtype ra;
   realtype qvsat;
   realtype qv;
   realtype etp;
   realtype isval;
   realtype frac_snow;
   realtype snow_rate;
   realtype melt_rate;
   realtype intcp_max;
   realtype meltf;

   if (tid < nelem){

		/* Note the dependence on physical units */
		elem_d_ps_albedo[tid] = 0.5 * (elem_d_lc_albedomin[tid] + elem_d_lc_albedomax[tid]);
		radnet = elem_d_ef_soldn[tid] * (1.0 - elem_d_ps_albedo[tid]);
		sfctmp = elem_d_es_sfctmp[tid] - 273.15;
		wind = elem_d_ps_sfcspd[tid];
		rh = elem_d_ps_rh[tid] / 100.0;

		vp = 611.2 * exp(17.67 * sfctmp / (sfctmp + 243.5)) * rh;
		pres = 101.325 * 1.0e3 *
			pow((293.0 - 0.0065 * elem_d_topo_zmax[tid]) / 293.0, 5.26);
		qv = 0.622 * vp / pres;
		qvsat = 0.622 * (vp / rh) / pres;

		if (elem_d_attrib_lai_type[tid] > 0)
		{
			lai = elem_d_ps_proj_lai[tid];
		}
		else
		{   
			lai = MonthlyLai_cuda[elem_d_attrib_lc_type[tid]-1][month-1];  // from table in forcing.c
		}

		meltf = MonthlyMf_cuda[month-1];                                  // from table in forcing.c

		/* Snow accumulation and snow melt calculation */
		frac_snow = (sfctmp < TSNOW) ?
			         1.0 : ((sfctmp > TRAIN) ? 0.0 : (TRAIN - sfctmp) / (TRAIN - TSNOW));
		snow_rate = frac_snow * elem_d_wf_prcp[tid];
		elem_d_ws_sneqv[tid] += snow_rate * stepsize;  
		melt_rate = (sfctmp > T0) ? (sfctmp - T0) * meltf : 0.0;

		if (elem_d_ws_sneqv[tid] > melt_rate * stepsize)
		{
			elem_d_ws_sneqv[tid] -= melt_rate * stepsize;   
		}
		else
		{
			melt_rate = elem_d_ws_sneqv[tid] / stepsize;
			elem_d_ws_sneqv[tid] = 0.0;
		}

		/* ThroughFall and Evaporation from canopy */
		intcp_max = elem_d_lc_cmcfactr[tid] * lai * elem_d_lc_shdfac[tid];

		z0 = MonthlyRl_cuda[elem_d_attrib_lc_type[tid]-1][month-1];          // from table in forcing.c

		ra = log(elem_d_ps_zlvl_wind[tid] / z0) *
			log(10.0 * elem_d_ps_zlvl_wind[tid] / z0) / (wind * 0.16);

		gamma = 4.0 * 0.7 * SIGMA * RD / CP * pow(sfctmp + 273.15, 4) /
			    (pres / ra) + 1.0;
		delta =
			LVH2O * LVH2O * 0.622 / RV / CP / pow(sfctmp + 273.15, 2) * qvsat;

		etp = (radnet * delta + gamma * (1.2 * LVH2O * (qvsat - qv) / ra)) /
			(1000.0 * LVH2O * (delta + gamma));

		if (elem_d_soil_depth[tid] - elem_d_ws_gw[tid] < elem_d_ps_rzd[tid])
		{
			satn = 1.0;
		}
		else
		{
			satn =
				((elem_d_ws_unsat[tid] /
				(elem_d_soil_depth[tid] - elem_d_ws_gw[tid])) > 1.0) ?
				1.0 :
				((elem_d_ws_unsat[tid] /
				(elem_d_soil_depth[tid] - elem_d_ws_gw[tid])) < 0.0) ?
				0.0 :
				0.5 * (1.0 - cos(3.14 *
				(elem_d_ws_unsat[tid] / (elem_d_soil_depth[tid] - elem_d_ws_gw[tid]))));
		}

		betas = (satn * elem_d_soil_porosity[tid] +
			elem_d_soil_smcmin[tid] - elem_d_soil_smcwlt[tid]) /
			(elem_d_soil_smcref[tid] - elem_d_soil_smcwlt[tid]);
		betas = (betas < 0.0001) ? 0.0001 : ((betas > 1.0) ? 1.0 : betas);
		elem_d_wf_edir[tid] = (1.0 - elem_d_lc_shdfac[tid]) * pow(betas, 2.0) * etp;
		elem_d_wf_edir[tid] *= cal_edir;         
		elem_d_wf_edir[tid] = (elem_d_wf_edir[tid] < 0.0) ? 0.0 : elem_d_wf_edir[tid];  


		/* Note the dependence on physical units */
		if (lai > 0.0)
		{
			elem_d_wf_ec[tid] = elem_d_lc_shdfac[tid] *
				pow(((elem_d_ws_cmc[tid] < 0.0) ? 0.0 :
				((elem_d_ws_cmc[tid] > intcp_max) ? intcp_max : elem_d_ws_cmc[tid])) /
				intcp_max, elem_d_lc_cfactr[tid]) * etp;
			elem_d_wf_ec[tid] *= cal_ec;         
			elem_d_wf_ec[tid] = (elem_d_wf_ec[tid] < 0.0) ? 0.0 : elem_d_wf_ec[tid];  

			fr = 1.1 * radnet / (elem_d_epc_rgl[tid] * lai);
			fr = (fr < 0.0) ? 0.0 : fr;   
			alphar =
				(1.0 + fr) / (fr + (elem_d_epc_rsmin[tid] / elem_d_epc_rsmax[tid]));
			alphar = (alphar > 10000.0) ? 10000.0 : alphar;   
			etas =
				1.0 - 0.0016 * (pow((elem_d_epc_topt[tid] - 273.15 - sfctmp), 2.0));
			etas = (etas < 0.0001) ? 0.0001 : etas; 
			gammas = 1.0 / (1.0 + 0.00025 * (vp / rh - vp));
			gammas = (gammas < 0.01) ? 0.01 : gammas;  
			rs = elem_d_epc_rsmin[tid] * alphar / (betas * lai * etas * gammas);
			rs = (rs > elem_d_epc_rsmax[tid]) ? elem_d_epc_rsmax[tid] : rs; 

			pc = (1.0 + delta / gamma) / (1.0 + rs / ra + delta / gamma);

			elem_d_wf_ett[tid] = elem_d_lc_shdfac[tid] * pc *
				(1.0 - pow((elem_d_ws_cmc[tid] < 0.0) ?
				0.0 :
				((elem_d_ws_cmc[tid] > intcp_max) ?
				intcp_max : elem_d_ws_cmc[tid]) / intcp_max, elem_d_lc_cfactr[tid])) *
				etp;
			elem_d_wf_ett[tid] *= cal_ett;          
			elem_d_wf_ett[tid] = (elem_d_wf_ett[tid] < 0.0) ? 0.0 : elem_d_wf_ett[tid];  
			elem_d_wf_ett[tid] = 
				((elem_d_ws_gw[tid] < (elem_d_soil_depth[tid] - elem_d_ps_rzd[tid]))
				&& elem_d_ws_unsat[tid] <= 0.0) ? 0.0 : elem_d_wf_ett[tid];  

			/* Drip function from Rutter and Morton, 1977, Journal of Applied Eco__logfy
			 * D0 = 3.91E-5 m/min = 6.52E-7 m/s */
			elem_d_wf_drip[tid] = (elem_d_ws_cmc[tid] <= 0.0) ?
				0.0 : 6.52E-7 * intcp_max * exp(3.89 * elem_d_ws_cmc[tid] / intcp_max);
		}
		else
		{
			elem_d_wf_ett[tid] = 0.0;
			elem_d_wf_ec[tid]  = 0.0;
			elem_d_wf_drip[tid] = 0.0;
		}

		if (elem_d_wf_drip[tid] < 0.0)
		{
			elem_d_wf_drip[tid] = 0.0;
		}
		if (elem_d_wf_drip[tid] * stepsize > elem_d_ws_cmc[tid])
		{
			elem_d_wf_drip[tid] = elem_d_ws_cmc[tid] / stepsize;
		}

		isval = elem_d_ws_cmc[tid] +
			(1.0 - frac_snow) * elem_d_wf_prcp[tid] * elem_d_lc_shdfac[tid] * stepsize -
			elem_d_wf_ec[tid] * stepsize - elem_d_wf_drip[tid] * stepsize;

		if (isval > intcp_max)
		{
			elem_d_ws_cmc[tid] = intcp_max;
			elem_d_wf_drip[tid] += (isval - intcp_max) / stepsize;  
		}
		else if (isval < 0.0)
		{
			elem_d_ws_cmc[tid] = 0.0;
			if (elem_d_wf_ec[tid] + elem_d_wf_drip[tid] > 0.0)
			{
				elem_d_wf_ec[tid] =
					elem_d_wf_ec[tid] / (elem_d_wf_ec[tid] + elem_d_wf_drip[tid]) *
					(elem_d_ws_cmc[tid] +
					(1.0 - frac_snow) * elem_d_wf_prcp[tid] * elem_d_lc_shdfac[tid] *
					stepsize);   
				elem_d_wf_drip[tid] =
					elem_d_wf_drip[tid] / (elem_d_wf_ec[tid] + elem_d_wf_drip[tid]) *
					(elem_d_ws_cmc[tid] +
					(1.0 - frac_snow) * elem_d_wf_prcp[tid] * elem_d_lc_shdfac[tid] *
					stepsize);   
			}
		}
		else
		{
			elem_d_ws_cmc[tid] = isval;
		}

		elem_d_wf_pcpdrp[tid] = 
			(1.0 - elem_d_lc_shdfac[tid]) * (1.0 - frac_snow) * elem_d_wf_prcp[tid] +
			elem_d_wf_drip[tid] + melt_rate;

    }  //  if (tid < nelem)

}







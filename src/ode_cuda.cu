#include "pihm.h"
#include "pihm.cuh"

#define max(a,b)(a>b?a:b)
#define min(a,b)(a<b?a:b)

/*
    2021.04.19:  与CPU代码计算结果对比，检验并确定了FrictSlope()计算结果的正确性. 
	             确定所有的核函数计算结果与CPU版本的一致.
    2021.04.23: 函数定义为static，对CV_Y和CV_Ydot的内容有影响.

*/

/* The RHS of ODE used in CVode Solver */
static int ODE(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *pihm_data)
{
	pihm_struct_d pihm_d =(pihm_struct_d) pihm_data;

	/* Extract needed constants from data */
	const sunindextype nelem = pihm_d->nelem;
	const sunindextype nriver = pihm_d->nriver;
	const sunindextype NEQ = pihm_d->NEQ;

// pihm_d中的数组首地址赋值给elem_d和river_d的数组
#include "pihm_device_injection.c"   // -Wno-unused-parameter

// pihm_d中的ctrl_d和cal_d的参数
	const sunindextype surf_mode = pihm_d->surf_mode;
	const sunindextype stepsize = pihm_d->stepsize;
	const sunindextype riv_mode = pihm_d->riv_mode;

	/* Extract pointers to vector data */
	realtype  *y_d = N_VGetDeviceArrayPointer_Cuda(CV_Y);
    realtype  *dy_d = N_VGetDeviceArrayPointer_Cuda(CV_Ydot);

   /* Set the threads and blocks for executing kernels */
    unsigned numThreads_elem = min(32, nelem);;
    unsigned numBlocks_elem = (nelem + numThreads_elem - 1) / numThreads_elem;
    unsigned numThreads_river = min(32, nriver);
    unsigned numBlocks_river = (nriver + numThreads_river - 1) / numThreads_river;

	// 0 Initialization of RHS of ODEs   
	unsigned numThreads_NEQ = min(32, NEQ); 
	unsigned numBlocks_NEQ = (NEQ + numThreads_NEQ - 1) / numThreads_NEQ;
	init_dy <<<numThreads_NEQ, numBlocks_NEQ >>>(dy_d, NEQ);  // dy_d = 0.0

	// 1 Use y_d to update some state variables
	Y_update_elem << <numThreads_elem, numBlocks_elem >> >
		                        (y_d, 
								nelem,
								elem_d_ws_surf,elem_d_ws_unsat, elem_d_ws_gw);
	Y_update_river << <numThreads_river, numBlocks_river >> >
		                        (y_d, 
								nelem, nriver,
								river_d_ws_stage, river_d_ws_gw,
								river_d_wf_rivflow0, river_d_wf_rivflow10);

	// 检测y_d是否更新了相关状态变量	下面的输出导致计算中断!  2021.05.07
#if 0
	cudaDeviceSynchronize();

	FILE *outfile0 = NULL;
	outfile0 = fopen("GPU_WStorage.txt", "w+");
	double temp1, temp2, temp3;
	for (int tid = 0; tid < nelem; tid++){
		temp1 = pihm_d->elem_d_ws_surf[tid];
		temp2 = pihm_d->elem_d_ws_unsat[tid];
		temp3 = pihm_d->elem_d_ws_gw[tid];
		fprintf(outfile0, "%15.12f %15.12f %15.12f\n", temp1, temp2, temp3);
	}
	fclose(outfile0);
	exit(0);
#endif

	// 2 PIHM Hydrology fluxes
	// 2.1 计算流域的摩阻坡度, FrictSlope(), moved from lat_flow.c
	realtype *dhbydx, *dhbydy;
	cudaMalloc((void **)&dhbydx, nelem * sizeof(realtype));
	cudaMalloc((void **)&dhbydy, nelem * sizeof(realtype));

	FrictSlope << <numThreads_elem, numBlocks_elem >> >
		(surf_mode, dhbydx, dhbydy,
#include "FrictSlope_device_arguments.c"	
		 );

// 检查FrictSlope核函数计算结果的正确性.
#if 0
	double *dhbydx_h, *dhbydy_h;   // 用户主机上输出检查的数组
	dhbydx_h = (realtype *)malloc(nelem * sizeof(realtype));
	dhbydy_h = (realtype *)malloc(nelem * sizeof(realtype));

	cudaMemcpy(dhbydx_h, dhbydx, nelem * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(dhbydy_h, dhbydy, nelem * sizeof(realtype), cudaMemcpyDeviceToHost);

	FILE *outfile0 = NULL;
	outfile0 = fopen("GPU_FrictSlope.txt", "w+");
	for (int tid = 0; tid < nelem; tid++){
		fprintf(outfile0, "%15.10f %15.10f\n", dhbydx_h[tid], dhbydy_h[tid]);
	}
	fclose(outfile0);

	free(dhbydx_h);
	free(dhbydx_h);

	exit(0);
#endif

	// 2.2 Calculate Hydrologic Process in elements
	Hydrol_elem << <numThreads_elem, numBlocks_elem >> >
		(nelem,
		dhbydx, dhbydy,
		surf_mode, stepsize,
#include "Hydrol_elem_device_arguments.c"
		);

	// cudaFree the variables in FrictSlope
	cudaFree(dhbydx);
	cudaFree(dhbydy);

	// 2.3 Calculate Hydrologic Process about River Flow
	RiverFlow << <numThreads_river, numBlocks_river >> >
		(nriver, riv_mode,
#include "RiverFlow_device_arguments.c"  
		);

	// 3 build RHS of ODE
	build_rhs_elem << <numThreads_elem, numBlocks_elem >> >
		(nelem,
		 dy_d,         
#include "build_rhs_elem_arguments.c"		
		 );

	// ODEs for river segments
	build_rhs_river << <numThreads_river, numBlocks_river >> >
		(dy_d,
		 nelem, nriver,
#include "build_rhs_river_arguments.c"		
		);

	// Wait for GPU to finish before accessing on host
	cudaDeviceSynchronize();

	// 检查河流计算相关变量的正确性
#if 0
	double *river_h_wf_rivflow0, *river_h_wf_rivflow1, *river_h_wf_rivflow2;
	double *river_h_wf_rivflow3, *river_h_wf_rivflow4, *river_h_wf_rivflow5;
	double *river_h_wf_rivflow6, *river_h_wf_rivflow7, *river_h_wf_rivflow8;
	double *river_h_wf_rivflow9, *river_h_wf_rivflow10;

	river_h_wf_rivflow0 = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_wf_rivflow1 = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_wf_rivflow2 = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_wf_rivflow3 = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_wf_rivflow4 = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_wf_rivflow5 = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_wf_rivflow6 = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_wf_rivflow7 = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_wf_rivflow8 = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_wf_rivflow9 = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_wf_rivflow10 = (realtype *)malloc(nriver * sizeof(realtype));

	cudaMemcpy(river_h_wf_rivflow0, pihm_d->river_d_wf_rivflow[0], nriver * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_wf_rivflow1, pihm_d->river_d_wf_rivflow[1], nriver * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_wf_rivflow2, pihm_d->river_d_wf_rivflow[2], nriver * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_wf_rivflow3, pihm_d->river_d_wf_rivflow[3], nriver * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_wf_rivflow4, pihm_d->river_d_wf_rivflow[4], nriver * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_wf_rivflow5, pihm_d->river_d_wf_rivflow[5], nriver * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_wf_rivflow6, pihm_d->river_d_wf_rivflow[6], nriver * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_wf_rivflow7, pihm_d->river_d_wf_rivflow[7], nriver * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_wf_rivflow8, pihm_d->river_d_wf_rivflow[8], nriver * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_wf_rivflow9, pihm_d->river_d_wf_rivflow[9], nriver * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_wf_rivflow10, pihm_d->river_d_wf_rivflow[10], nriver * sizeof(realtype), cudaMemcpyDeviceToHost);
	

	FILE *outfile0 = NULL;
	outfile0 = fopen("GPU_RiverFlow.txt", "w+");
	double temp0, temp1, temp2;
	for (int tid = 0; tid < nriver; tid++){
		temp0 = river_h_wf_rivflow0[tid];
		temp1 = river_h_wf_rivflow1[tid];
		temp2 = river_h_wf_rivflow2[tid];

		fprintf(outfile0, "%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f\n",
			temp0, temp1, temp2,
			river_h_wf_rivflow3[tid], river_h_wf_rivflow4[tid], river_h_wf_rivflow5[tid],
			river_h_wf_rivflow6[tid], river_h_wf_rivflow7[tid], river_h_wf_rivflow8[tid],
			river_h_wf_rivflow9[tid], river_h_wf_rivflow10[tid]
			);
			
	}
	fclose(outfile0);
	exit(0);
#endif

#if 0
	// 检查EtExtract核函数计算的正确性.
	/*
	double *elem_h_wf_edir_surf, *elem_h_wf_edir_unsat, *elem_h_wf_edir_gw;   // 用户主机上输出检查的数组
	elem_h_wf_edir_surf = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_wf_edir_unsat = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_wf_edir_gw = (realtype *)malloc(nelem * sizeof(realtype));

	cudaMemcpy(elem_h_wf_edir_surf, pihm_d->elem_d_wf_edir_surf, nelem * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_wf_edir_unsat, pihm_d->elem_d_wf_edir_unsat, nelem * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_wf_edir_gw, pihm_d->elem_d_wf_edir_gw, nelem * sizeof(realtype), cudaMemcpyDeviceToHost);
	*/

	// 检查LateralFlow核函数计算的正确性
	//double  *elem_h_wf_rechg, *elem_h_wf_edir_gw, *elem_h_wf_ett_gw;   // 用户主机上输出检查的数组
	double  *elem_h_wf_subsurf0, *elem_h_wf_subsurf1, *elem_h_wf_subsurf2;

	elem_h_wf_subsurf0 = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_wf_subsurf1 = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_wf_subsurf2 = (realtype *)malloc(nelem * sizeof(realtype));

	cudaMemcpy(elem_h_wf_subsurf0, pihm_d->elem_d_wf_subsurf[0], nelem * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_wf_subsurf1, pihm_d->elem_d_wf_subsurf[1], nelem * sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_wf_subsurf2, pihm_d->elem_d_wf_subsurf[2], nelem * sizeof(realtype), cudaMemcpyDeviceToHost);

	// 输出中间计算结果
	FILE *outfile0 = NULL;
	//outfile0 = fopen("GPU_SubFlowElemToElem.txt", "w+");
	//outfile0 = fopen("GPU_OvlFlowElemToElem.txt", "w+");
	//outfile0 = fopen("GPU_VerticalFlow.txt", "w+");
	//outfile0 = fopen("GPU_check_dy[offset3].txt", "w+");
	outfile0 = fopen("GPU_check_subsurface.txt", "w+");

	double temp0, temp1, temp2;
	for (int tid = 0; tid < nelem; tid++){
		/*
		temp0 = elem_h_wf_rechg[tid];
		temp1 = elem_h_wf_edir_gw[tid];
		temp2 = elem_h_wf_ett_gw[tid];
		*/
		temp0 = elem_h_wf_subsurf0[tid];
		temp1 = elem_h_wf_subsurf1[tid];
		temp2 = elem_h_wf_subsurf2[tid];

		fprintf(outfile0, "%15.12f %15.12f %15.12f\n", temp0, temp1, temp2);
	}
	fclose(outfile0);
	exit(0);
#endif

// 检查RiverFlow是否修改了河道两侧单元的相关数组
#if 0
	realtype *elem_h_gw;
	elem_h_gw = (realtype *)malloc(nelem * sizeof(realtype));
	cudaMemcpy(elem_h_gw, pihm_d->elem_d_ws_gw, nelem * sizeof(realtype), cudaMemcpyDeviceToHost);

	FILE *outfile1 = NULL;
	outfile1 = fopen("GPU_elem_gw_after_riverflow.txt", "w+");  // 与CPU的计算结果符合!
	for (int i = 0; i < nelem; i++){
		fprintf(outfile1, "%15.12f\n", elem_h_gw[i]);
	}
	exit(0);
#endif


// 检查dy的计算结果      2021.04.23 通过检查!
#if 0
	realtype *dy_h, *y_h;
	dy_h = (realtype *)malloc(NEQ * sizeof(realtype));
	y_h = (realtype *)malloc(NEQ * sizeof(realtype));

    /* 执行cudaMemcpy数据拷贝，无法更新CV_Ydot的内容；
	   执行get_nvector_cuda数据拷贝，才能更新CV_Ydot的内容
	   问题是：这样的数据拷贝降低计算效率   -- 2021.05.06
	*/
	//cudaMemcpy(dy_h, dy_d, NEQ * sizeof(realtype), cudaMemcpyDeviceToHost);  
	get_nvector_cuda(CV_Ydot, dy_h, NEQ);

	//cudaMemcpy(y_h, y_d, NEQ * sizeof(realtype), cudaMemcpyDeviceToHost);
	get_nvector_cuda(CV_Y, y_h, NEQ);

	FILE *outfile0 = NULL;
	outfile0 = fopen("GPU_y_d_in_ODE.txt", "w+");
	for (int i = 0; i < NEQ; i++){
		fprintf(outfile0, "%15.12f\n",y_h[i]);   // y_d 内容正确!
	}
	fclose(outfile0);

	FILE *outfile1 = NULL;
	outfile1 = fopen("GPU_dy_d_in_ODE.txt", "w+");
	for (int i = 0; i < NEQ; i++){
		fprintf(outfile1, "%15.12f\n", dy_h[i]);  // dy_d 内容正确!
	}
	fclose(outfile1);
#endif

// 检查CV_Ydot的内容: 确定经过CVode求解后，dy_d和CV_Ydot的内容更新并是一致的  ~2020.05.07
#if 0
	FILE *outfile2 = NULL;
	outfile2 = fopen("GPU_CV_Y_in_CVode.txt", "w+");
	N_VPrintFile_Cuda(CV_Y, outfile2);     // 内容不对!?
	fclose(outfile2);

	// GPU和CPU计算的CV_Y和CV_Ydot都是一样的，说明涉及的所有核函数都是正确的.   2021.05.07
	FILE *outfile3 = NULL;
	outfile3 = fopen("GPU_CV_Ydot_in_CVode.txt", "w+");  // 需要执行 get_nvector_cuda
	N_VPrintFile_Cuda(CV_Ydot, outfile3);  
	fclose(outfile3);

//	exit(0);
#endif


// 核函数执行的错误检测
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

	return(0);
}


// 初始化dy   
__global__ static void init_dy(realtype  *dy, const int NEQ)
{
	int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid < NEQ){
		dy[tid] = 0.0;
	}	   
}

/*  Updating the elements' variables using y (CV_Y) */
__global__ static void Y_update_elem(const  realtype *y,
	                           const int nelem,
	realtype *elem_d_ws_surf, realtype *elem_d_ws_unsat, realtype *elem_d_ws_gw)
{
	int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid < nelem){

		// 使用elem_d变量，更新计算CV_Y (y_d)
		int offset1 = tid;
		int offset2 = tid + nelem;
		int offset3 = tid + 2 * nelem;

		// 三个部分的水源：地表水(surface), 土壤非饱和水(unsaturated)和地下水(groundwater)
		elem_d_ws_surf[tid] = (y[offset1] >= 0.0) ? y[offset1] : 0.0;
		elem_d_ws_unsat[tid] = (y[offset2] >= 0.0) ? y[offset2] : 0.0;
		elem_d_ws_gw[tid] = (y[offset3] >= 0.0) ? y[offset3] : 0.0;
	}
}

/*  Updating the rivers' variables using y (CV_Y) */
__global__  static void Y_update_river(const  realtype *y,
	                            const int nelem, const int nriver,
	realtype *river_d_ws_stage, realtype *river_d_ws_gw,
	realtype *river_d_wf_rivflow0, realtype *river_d_wf_rivflow10)
{
	int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid < nriver){

		int offset1 = tid + 3 * nelem;
		int offset2 = tid + 3 * nelem + nriver;

		river_d_ws_stage[tid] = (y[offset1] >= 0.0) ? y[offset1] : 0.0;
		river_d_ws_gw[tid] = (y[offset2] >= 0.0) ? y[offset2] : 0.0;

		river_d_wf_rivflow0[tid] = 0.0;   
		river_d_wf_rivflow10[tid] = 0.0; 
	}
}

/*  Calculate hydrologic process in elements */
__global__  static void Hydrol_elem(const int nelem, const realtype *dhbydx, const realtype *dhbydy,
	                         const int surf_mode, const int stepsize,
#include "Hydrol_elem_device_declarations.c"
)
{
	int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid < nelem){

		/* Calculate actual surface water depth */
		elem_d_ws_surfh[tid] = SurfH(elem_d_ws_surf[tid]);

		/* Determine which layers does ET extract water from */
		EtExtract(tid,
#include "EtExtract_device_arguments.c"			
			);

		/* Water flow */
		// 侧向上，单元到单元的流动，坡面流动力学
		LateralFlow(nelem, surf_mode, tid,
			dhbydx, dhbydy,
#include "LateralFlow_river_device_arguments.c"
			);
		// we need stream synchnization after LaterFlow(), because of 3 cuda streams

		// 垂向上，土壤层之间的流动，渗流
		VerticalFlow(stepsize, tid,
#include "VerticalFlow_device_arguments.c"
			);
	}
}

__device__ realtype SurfH(realtype surfeqv)
{
	realtype surfh;

	if (DEPRSTG == 0.0)
	{
		return surfeqv;
	}
	else
	{
		if (surfeqv < 0.0)
		{
			surfh = 0.0;
		}
		else if (surfeqv <= 0.5 * DEPRSTG)
		{
			surfh = sqrt(2.0 * DEPRSTG * surfeqv);
		}
		else
		{
			surfh = DEPRSTG + (surfeqv - 0.5 * DEPRSTG);
		}

		return surfh;
	}

}

/* Build rhs of ODE about Basin Elements */
__global__  static void build_rhs_elem(const int nelem, realtype  *dy,
#include "build_rhs_elem_declarations.c"
)
{
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	
    if (tid < nelem){

		int offset1 = tid;
		int offset2 = tid + nelem;
		int offset3 = tid + 2 * nelem;

/* 使用新的elem_d结构体变量，计算RHS的CV_Ydot (dy_d) */
		dy[offset1] += elem_d_wf_pcpdrp[tid] - elem_d_wf_infil[tid] - elem_d_wf_edir_surf[tid];
		dy[offset2] += elem_d_wf_infil[tid] - elem_d_wf_rechg[tid] - elem_d_wf_edir_unsat[tid] -
	    	           elem_d_wf_ett_unsat[tid];
		dy[offset3] += elem_d_wf_rechg[tid] - elem_d_wf_edir_gw[tid] - elem_d_wf_ett_gw[tid];
	    
	    /* Horizontal water fluxes, 需要RiverFlow的计算结果 */
		dy[offset1] -= (elem_d_wf_ovlflow0[tid] + elem_d_wf_ovlflow1[tid] +
			            elem_d_wf_ovlflow2[tid]) / elem_d_topo_area[tid];	

		dy[offset3] -= (elem_d_wf_subsurf0[tid] + elem_d_wf_subsurf1[tid] +
			            elem_d_wf_subsurf2[tid] )/elem_d_topo_area[tid];   
	    
		dy[offset2] /= elem_d_soil_porosity[tid];
		dy[offset3] /= elem_d_soil_porosity[tid];

    } // if (tid < nelem)

	__syncthreads();
}


/* Build rhs of ODE about River Flow */
__global__  static void build_rhs_river(realtype  *dy,
	                             const int nelem, const int nriver,
#include "build_rhs_river_declarations.c"
)
{
	int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid < nriver){

		int offset1 = tid + 3 * nelem;
		int offset2 = tid + 3 * nelem + nriver;

//-----------使用新的 river_d结构体变量，计算RHS的CV_Ydot (dy_d)---------------------------
			/* Note the limitation due to
			*  d(v) / dt = a * dy / dt + y * da / dt
			*  for cs other than rectangle */
		dy[offset1] -= (river_d_wf_rivflow0[tid] + river_d_wf_rivflow1[tid] +
			            river_d_wf_rivflow2[tid] + river_d_wf_rivflow3[tid] +
						river_d_wf_rivflow4[tid] + river_d_wf_rivflow5[tid] +
						river_d_wf_rivflow6[tid]) / river_d_topo_area[tid];
		
		dy[offset2] += -river_d_wf_rivflow7[tid] -
			            river_d_wf_rivflow8[tid] -
			            river_d_wf_rivflow9[tid] -
			            river_d_wf_rivflow10[tid] +
			            river_d_wf_rivflow6[tid];

		dy[offset2] /= river_d_matl_porosity[tid] * river_d_topo_area[tid];

    }  // if (tid < nriver)

	__syncthreads();
}

// Set CVode parameters in CUDA.
void SetCVodeParam(pihm_struct pihm, pihm_struct_d pihm_d, void *cvode_mem, N_Vector CV_Y)	               
{
	int             cv_flag;
	static int      reset;

	N_Vector abstol = NULL;   // for CUDA

	pihm->ctrl.maxstep = pihm->ctrl.stepsize;

	/* Set the pointer to user-defined data */
	cv_flag = CVodeSetUserData(cvode_mem, pihm_d);

	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMprintf(VL_VERBOSE, "\n\nFailed to do CVodeSetUserData ... \n\n");
		PIHMexit(EXIT_FAILURE);
	}

	if (reset)
	{
		/* When model spins-up and recycles forcing, use CVodeReInit to reset
		* solver time, which does not allocates memory */
		cv_flag = CVodeReInit(cvode_mem, 0.0, CV_Y);
		if (!CheckCVodeFlag(cv_flag))
		{
			PIHMexit(EXIT_FAILURE);
		}
	}
	else
	{
		cv_flag = CVodeInit(cvode_mem, ODE, 0.0, CV_Y);  // 初始化ODE, 从 t=0 开始
		if (!CheckCVodeFlag(cv_flag))
		{
			PIHMprintf(VL_VERBOSE, "\n\nFailed to do CVodeInit ... \n\n");
			PIHMexit(EXIT_FAILURE);
		}
		reset = 1;
	}

	/* Setting abstol for CUDA */
	int size = NumStateVar();
	realtype *abstol_CPU;

	abstol = N_VNew_Cuda(size);   // abstol on GPU
	abstol_CPU = (realtype *)malloc(size * sizeof(realtype));

	for (int i = 0; i < size; i++)
		abstol_CPU[i] = (realtype)pihm->ctrl.abstol;

	set_nvector_cuda(abstol, abstol_CPU, size);
	free(abstol_CPU);

	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	* and scalar absolute tolerance */
	// In CUDA, abstol should also be in GPU.  2021.04.23
	cv_flag = CVodeSVtolerances(cvode_mem, (realtype)pihm->ctrl.reltol,
		                                    abstol);  
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMprintf(VL_VERBOSE, "\n\nFailed to do CVodeSVtolerances ... \n\n");
		PIHMexit(EXIT_FAILURE);
	}

	/* Set Initial time step for CVode */
	cv_flag = CVodeSetInitStep(cvode_mem, (realtype )pihm->ctrl.initstep);
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMprintf(VL_VERBOSE, "\n\nFailed to do CVodeSetInitStep ... \n\n");
		PIHMexit(EXIT_FAILURE);
	}

	cv_flag = CVodeSetStabLimDet(cvode_mem, TRUE);
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMprintf(VL_VERBOSE, "\n\nFailed to do CVodeSetStabLimDet ... \n\n");
		PIHMexit(EXIT_FAILURE);
	}

	cv_flag = CVodeSetMaxStep(cvode_mem, (realtype )pihm->ctrl.maxstep);
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMprintf(VL_VERBOSE, "\n\nFailed to do CVodeSetMaxStep ... \n\n");
		PIHMexit(EXIT_FAILURE);
	}

	cv_flag = CVodeSetMaxNumSteps(cvode_mem, pihm->ctrl.stepsize * 10);
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMprintf(VL_VERBOSE, "\n\nFailed to do CVodeSetMaxNumSteps ... \n\n");
		PIHMexit(EXIT_FAILURE);
	}

#if defined(_SUNDIALS_v4)
	/*
SUNLinearSolver Linsol = NULL;
Linsol = SUNLinSol_SPGMR(CV_Y, PREC_NONE, 0);     

cv_flag = CVodeSetLinearSolver(cvode_mem, Linsol, NULL);     
if (!CheckCVodeFlag(cv_flag))
{
	PIHMprintf(VL_VERBOSE, "\n\nFailed to do CVSpilsSetLinearSolver ... \n\n");
	PIHMexit(EXIT_FAILURE);
}
*/

#elif defined(_SUNDIALS_v3)
	SUNLinearSolver Linsol = NULL;

	// Create SPGMR solver structure without preconditioning
	// and the maximum Krylov dimension maxl */
	Linsol = SUNSPGMR(CV_Y, PREC_NONE, 0); 

	//Set CVSpils linear solver to LS 
	cv_flag = CVSpilsSetLinearSolver(cvode_mem, Linsol);  
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMprintf(VL_VERBOSE, "\n\nFailed to do CVSpilsSetLinearSolver ... \n\n");
		PIHMexit(EXIT_FAILURE);
	}
#else
	cv_flag = CVSpgmr(cvode_mem, PREC_NONE, 0);
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMexit(EXIT_FAILURE);
	}
#endif
}

/* Solve ODE in CVode library whithin CUDA */
// 必须将更新的CV_Y定义为另一个NVector CV_Yout   ~2015.05.06
void SolveCVode(int starttime, int *t, int nextptr, 
	realtype cputime, void *cvode_mem, N_Vector CV_Yout)  
{
	realtype         solvert;
	realtype         tout;
	pihm_t_struct    pihm_time;
	int              cv_flag;

	tout = (realtype )(nextptr - starttime);

	cv_flag = CVodeSetStopTime(cvode_mem, tout);
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMprintf(VL_VERBOSE, "\n\nFailed to do CVodeSetStopTime ... \n\n");
		PIHMexit(EXIT_FAILURE);
	}

	// Call CVode
//	printf("tout=%f\n\n",tout);
	cv_flag = CVode(cvode_mem, tout, CV_Yout, &solvert, CV_NORMAL);   
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMprintf(VL_VERBOSE, "\n\nFailed to do CVode ... \n\n");
		PIHMexit(EXIT_FAILURE);
	}

	*t = (int)round(solvert) + starttime;
//	printf("t=%d\n", *t);

	pihm_time = PIHMTime(*t);

	if (debug_mode)
	{
		PIHMprintf(VL_NORMAL, " Step = %s (%d)\n", pihm_time.str, *t);
	}
	else if (spinup_mode)
	{
		if (pihm_time.t % DAYINSEC == 0)
		{
			PIHMprintf(VL_NORMAL, " Step = %s\n", pihm_time.str);
		}
	}
	else if (pihm_time.t % 3600 == 0)
	{
		PIHMprintf(VL_NORMAL,
			" Step = %s (cputime %f)\n", pihm_time.str, cputime);
	}

}

void AdjCVodeMaxStep(void *cvode_mem, ctrl_struct *ctrl)
{
	/* Variable CVODE max step (to reduce oscillations) */
	long int        nst;
	long int        ncfn;
	long int        nni;
	static long int nst0;
	static long int ncfn0;
	static long int nni0;
	int             cv_flag;
	realtype          nsteps;
	realtype          nfails;
	realtype          niters;

	cv_flag = CVodeGetNumSteps(cvode_mem, &nst);
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMexit(EXIT_FAILURE);
	}

	cv_flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMexit(EXIT_FAILURE);
	}

	cv_flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMexit(EXIT_FAILURE);
	}

	nsteps = (double)(nst - nst0);
	nfails = (double)(ncfn - ncfn0) / nsteps;
	niters = (double)(nni - nni0) / nsteps;

	if (nfails > ctrl->nncfn || niters >= ctrl->nnimax)
	{
		ctrl->maxstep /= ctrl->decr;
	}

	if (nfails == 0.0 && niters <= ctrl->nnimin)
	{
		ctrl->maxstep *= ctrl->incr;
	}

	ctrl->maxstep = (ctrl->maxstep < ctrl->stepsize) ?
		ctrl->maxstep : ctrl->stepsize;
	ctrl->maxstep = (ctrl->maxstep > ctrl->stmin) ?
		ctrl->maxstep : ctrl->stmin;

	cv_flag = CVodeSetMaxStep(cvode_mem, (realtype )ctrl->maxstep);
	if (!CheckCVodeFlag(cv_flag))
	{
		PIHMexit(EXIT_FAILURE);
	}

	nst0 = nst;
	ncfn0 = ncfn;
	nni0 = nni;
}



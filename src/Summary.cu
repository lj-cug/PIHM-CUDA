#include "pihm.h"
#include "pihm.cuh"

/* 采用CVode求解得到的CV_Yout，更新相关状态变量值 */
void Summary(N_Vector CV_Y, void *pihm_data, int stepsize0)
{
	 pihm_struct_d pihm_d = static_cast<pihm_struct_d>(pihm_data);

	 realtype stepsize = (realtype)stepsize0;
	 int nelem = pihm_d->nelem;
	 int nriver = pihm_d->nriver;
	 int NEQ = 3 * nelem + 2 * nriver;

	 // pihm_d中的数组首地址赋值给elem_d和river_d的数组
#include "update_device_injections.c"

	realtype *dy_d = N_VGetDeviceArrayPointer_Cuda(CV_Y);

// 检查dy_d的内容
#if 0
	realtype *dy_h;
	dy_h = (realtype *)malloc(NEQ * sizeof(realtype));
	cudaMemcpy(dy_h, dy_d, NEQ * sizeof(realtype), cudaMemcpyDeviceToHost);

	FILE *outfile1 = NULL;
	outfile1 = fopen("dy_d_in_Summary.txt", "w+");
	for (int i = 0; i < NEQ; i++){
		fprintf(outfile1, "%15.12f\n", dy_h[i]);  
	}
	fclose(outfile1);
	exit(0);
#endif


	// (1) execute the kernel function updating elements
	unsigned int  flags = cudaMemAttachGlobal;
	realtype *subrunoff;
	cudaMallocManaged(&subrunoff, nelem * sizeof(double), flags);

	int numThreads_elem = min(32, nelem);
	int numBlocks_elem = static_cast<int>(ceil((double)nelem / (double)numThreads_elem));
	update_elem << <numThreads_elem, numBlocks_elem >> >	
		(dy_d, stepsize, nelem, subrunoff,
#include "update_elem_device_arguments.c"
		);
	cudaFree(subrunoff);


	// (2) execute the kernel function updating river segments
	int numThreads_river = min(32, nriver);
	int numBlocks_river = static_cast<int>(ceil((double)nriver / (double)numThreads_river));
	update_river <<<numThreads_river, numBlocks_river >>>
		(dy_d, stepsize, nelem, nriver, river_d_ws_stage, river_d_ws_gw);

	cudaDeviceSynchronize();

#if 0 
	// check for error
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess)
	{
		// print the CUDA error message and exit
		printf("CUDA error in Update: %s\n", cudaGetErrorString(error));
		exit(-1);
	}
#endif

	return;
}


__global__ void update_elem(realtype *dy, realtype stepsize, int nelem, realtype *subrunoff,
#include "update_elem_device_declarations.c" 
)
{
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int offset1 = tid;
	int offset2 = tid + nelem;
	int offset3 = tid + 2 * nelem;

	if (tid < nelem){

	// Basin Elements
	elem_d_ws_surf[tid] = dy[offset1];
	elem_d_ws_unsat[tid] = dy[offset2];
	elem_d_ws_gw[tid] = dy[offset3];
	
	MassBalance(subrunoff, stepsize, tid,
#include "mass_balance_arguments.c"
	);
		
	// update the initial ws0
	elem_d_ws0_surf[tid] = elem_d_ws_surf[tid];
	elem_d_ws0_unsat[tid] = elem_d_ws_unsat[tid];
	elem_d_ws0_gw[tid] = elem_d_ws_gw[tid];
	elem_d_ws0_sneqv[tid] = elem_d_ws_sneqv[tid];
	elem_d_ws0_cmcmax[tid] = elem_d_ws_cmcmax[tid];
	elem_d_ws0_cmc[tid] = elem_d_ws_cmc[tid];
	elem_d_ws0_surfh[tid] = elem_d_ws_surfh[tid];
	}
}

__device__ void MassBalance(realtype *subrunoff, realtype stepsize, int tid,
#include "mass_balance_declarations.c"
)
{
	realtype soilw0, soilw1;

	/*
	* Calculate infiltration based on mass conservation
	*/
	soilw0 = elem_d_ws0_gw[tid] + elem_d_ws0_unsat[tid];
	soilw0 = (soilw0 > elem_d_soil_depth[tid]) ? elem_d_soil_depth[tid] : soilw0;
	soilw0 = (soilw0 < 0.0) ? 0.0 : soilw0;

	soilw1 = elem_d_ws_gw[tid] + elem_d_ws_unsat[tid];
	soilw1 = (soilw1 > elem_d_soil_depth[tid]) ? elem_d_soil_depth[tid] : soilw1;
	soilw1 = (soilw1 < 0.0) ? 0.0 : soilw1;

	/* Subsurface runoff rate */
	subrunoff[tid] = 0.0;
	subrunoff[tid] += (elem_d_wf_subsurf0[tid] + elem_d_wf_subsurf1[tid] + elem_d_wf_subsurf2[tid])
			           / elem_d_topo_area[tid];


	elem_d_wf_infil[tid] = (soilw1 - soilw0) * elem_d_soil_porosity[tid] / stepsize + subrunoff[tid] +
		elem_d_wf_edir_unsat[tid] + elem_d_wf_edir_gw[tid] +
		elem_d_wf_ett_unsat[tid] + elem_d_wf_ett_gw[tid];

	if (elem_d_wf_infil[tid] < 0.0)
	{
		subrunoff[tid] -= elem_d_wf_infil[tid];
		elem_d_wf_infil[tid] = 0.0;
	}
}


__global__ void update_river(realtype *dy, realtype stepsize, int nelem, int nriver, 
	realtype *river_d_ws_stage, realtype *river_d_ws_gw)
{
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int offset1 = tid + 3 * nelem;
	int offset2 = offset1 + nriver;
	if (tid < nriver){
	river_d_ws_stage[tid] = dy[offset1];
	river_d_ws_gw[tid] = dy[offset2];

	// river_d_ws0 = river_d_ws;  // river->ws0 was not used.
	}
}



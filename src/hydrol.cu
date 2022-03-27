#include "pihm.h"
#include "pihm.cuh"

/* 由doule *y (NVector CV_Y) 更新 elem_d 和river_d 结构体中的数组 */

void hydrol(double t, double *y_d, void *pihm_data,
	elem_struct_d elem_d, river_struct_d river_d, ctrl_struct_d ctrl_d)
{
	pihm_struct_d pihm_d = static_cast<pihm_struct_d>(pihm_data);
	
	// 从CPU端传入hydrol_elemt和hydrol_river
	//y_d = (double *)malloc(NEQ*sizeof(double));
	//unsigned int  flags = cudaMemAttachGlobal;
	//cudaMallocManaged((void **)&y_d, NEQ*sizeof(double), flags);


	// 从设备端执行y_d代入
   // cudaMalloc((void **)&y_d, NEQ*sizeof(double));
	//y_h = N_VGetHostArrayPointer_Cuda(CV_Y);
	//cudaMemcpy(y_d, y_h, sizeof(double) * NEQ, cudaMemcpyHostToDevice);



	//  (1) 执行计算流域单元的CV_Ydot的CUDA核函数
	int numThreads_elem = 32;  // min(32, nelem);
	int numBlocks_elem = static_cast<int>(ceil((double)nelem / (double)numThreads_elem));
	hydrol_elemt << <numThreads_elem, numBlocks_elem >> > (y_d, pihm_d, elem_d, river_d, ctrl_d);
		           

	//  (2) 执行计算河流单元的CV_Ydot的CUDA核函数
	int numThreads_river = 32;  // min(32, nriver);
	int numBlocks_river = static_cast<int>(ceil((double)nriver / (double)numThreads_river));
	hydrol_river << <numThreads_river, numBlocks_river >> >(y_d, pihm_d, elem_d, river_d, ctrl_d);
		            
	cudaDeviceSynchronize();
	
}

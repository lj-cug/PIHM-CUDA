#include "pihm.h"
#include "pihm.cuh"

int Summary(elem_struct_d *elem, river_struct_d *river, N_Vector CV_Y,
	        int stepsize, int nelem, int nriver)
{
	double stepsize_d;
	stepsize_d = (double)stepsize;

	//printf("Ready to update state variables in cuda. \n\n");
	/*
	FILE *outfile0  = NULL;
	outfile0 = fopen("CV_Y.txt", "w+");
	N_VPrintFile_Cuda(CV_Y, outfile0);
	fclose(outfile0);
	*/

	N_VCopyFromDevice_Cuda(CV_Y);
	double *y_d = N_VGetDeviceArrayPointer_Cuda(CV_Y);

#ifdef _DEBUG_
	printf("Try to output the y_d to verify. \n\n");
	FILE *outfile = NULL;
	outfile = fopen("y_d.txt", "w+");
	double *y_h;
	y_h = (double*)malloc(size*sizeof(double));
	cudaMemcpy(y_h, y_d, size*sizeof(double), cudaMemcpyDeviceToHost);
	printf("Copy y data from device to host for verification. \n\n");

	for (int i = 0; i < size; i++) {
		fprintf(outfile,"%10.6f\n",y_h[i]);
	}
	fclose(outfile);
	exit(0);
#endif

	// (1) execute the kernel function updating elements
	int numThreads_elem = min(32, nelem);
	int numBlocks_elem = static_cast<int>(ceil((double)nelem / (double)numThreads_elem));
	update_elem << <numThreads_elem, numBlocks_elem >> >(elem, y_d, stepsize_d, nelem);

	// (2) execute the kernel function updating river segments
	int numThreads_river = min(32, nriver);
	int numBlocks_river = static_cast<int>(ceil((double)nriver / (double)numThreads_river));
	update_river <<<numThreads_river, numBlocks_river >>>(elem, river, y_d, stepsize_d, nelem, nriver);

	cudaDeviceSynchronize();

	return(0);
}


__global__ void update_elem(elem_struct_d *elem, double *y, double stepsize, int nelem)
{
	double subrunoff;

	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int offset1 = tid;
	int offset2 = tid + nelem;
	int offset3 = tid + 2 * nelem;

	// Basin Elements
	elem->ws.surf[tid] = y[offset1];
	elem->ws.unsat[tid] = y[offset2];
	elem->ws.gw[tid] = y[offset3];
	
	MassBalance(&elem->ws, &elem->ws0, &elem->wf, &subrunoff,
		            &elem->soil, elem->topo.area[tid], stepsize, tid);

	elem->ws0.surf[tid] = elem->ws.surf[tid];
	elem->ws0.unsat[tid] = elem->ws.unsat[tid];
	elem->ws0.gw[tid] = elem->ws.gw[tid];
	elem->ws0.sneqv[tid] = elem->ws.sneqv[tid];
	elem->ws0.cmcmax[tid] = elem->ws.cmcmax[tid];
	elem->ws0.cmc[tid] = elem->ws.cmc[tid];
	elem->ws0.surfh[tid] = elem->ws.surfh[tid];

}

__device__ void MassBalance(const wstate_struct_d *ws, const wstate_struct_d *ws0,
	wflux_struct_d *wf, double *subrunoff, soil_struct_d *soil,
	double area, double stepsize, int tid)
{

	int j;
	double soilw0, soilw1;

	/*
	* Calculate infiltration based on mass conservation
	*/
	soilw0 = ws0->gw[tid] + ws0->unsat[tid];
	soilw0 = (soilw0 > soil->depth[tid]) ? soil->depth[tid] : soilw0;
	soilw0 = (soilw0 < 0.0) ? 0.0 : soilw0;

	soilw1 = ws->gw[tid] + ws->unsat[tid];
	soilw1 = (soilw1 > soil->depth[tid]) ? soil->depth[tid] : soilw1;
	soilw1 = (soilw1 < 0.0) ? 0.0 : soilw1;

	/* Subsurface runoff rate */
	*subrunoff = 0.0;
	for (j = 0; j < NUM_EDGE; j++)
	{
		*subrunoff += wf->subsurf[j][tid] / area;
	}

	wf->infil[tid] = (soilw1 - soilw0) * soil->porosity[tid] / stepsize + *subrunoff +
		wf->edir_unsat[tid] + wf->edir_gw[tid] + wf->ett_unsat[tid] + wf->ett_gw[tid];

	if (wf->infil[tid] < 0.0)
	{
		*subrunoff -= wf->infil[tid];
		wf->infil[tid] = 0.0;
	}
}


__global__ void update_river(elem_struct_d *elem, river_struct_d *river, 
	double *y, double stepsize, int nelem, int nriver)
{

	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int offset1 = tid + 3 * nelem;
	int offset2 = offset1 + nriver;

	// River Segements
	river->ws.stage[tid] = y[offset1];
	river->ws.gw[tid] = y[offset2];

	// river->ws0 = river->ws;  // river->ws0 was not used.
}



#include "pihm.h"
#include "pihm.cuh"

/*
*y��*dy�Ѿ���ȡ��CV_Y��CV_Ydot��ָ�룬������NV_DATA������
*/
__global__ void hydrol_elemt(double *y, void *pihm_data, 
	                         elem_struct_d *elem_d, river_struct_d *river_d,
							 ctrl_struct_d *ctrl_d)
{
	pihm_struct_d pihm_d = static_cast<pihm_struct_d>(pihm_data);

	/* Loop over elements */
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int offset1 = tid;
	int offset2 = tid + pihm_d->nelem;
	int offset3 = tid + 2 * pihm_d->nelem;

	// �������ֵ�ˮԴ���ر�ˮ(surface), �����Ǳ���ˮ(unsaturated)�͵���ˮ(groundwater)
	elem_d->ws.surf[tid] = (y[offset1] >= 0.0) ? y[offset1] : 0.0;
	elem_d->ws.unsat[tid] = (y[offset2] >= 0.0) ? y[offset2] : 0.0;
	elem_d->ws.gw[tid] = (y[offset3] >= 0.0) ? y[offset3] : 0.0;

	/* Calculate actual surface water depth */
	elem_d->ws.surfh[tid] = SurfH(elem_d->ws.surf[tid]);

	/* ����,��hydrol.c�����дΪCUDA�˺���������4���Ӻ�����
	  EtExtract();
	  LateralFlow();   
	  VerticalFlow(); 
	  RiverFlow();      
	*/

	int nelem, surf_mode;
	double stepsize;
	nelem = pihm_d->nelem;
	surf_mode = ctrl_d->surf_mode;
	stepsize = (double)ctrl_d->stepsize;

	/* Determine which layers does ET extract water from */
	EtExtract(elem_d, tid);

	/* Water flow */
	LateralFlow(elem_d, river_d, nelem, surf_mode, tid);    // �����ϣ���Ԫ����Ԫ������������������ѧ
	VerticalFlow(elem_d, stepsize, tid);                  // �����ϣ�������֮�������������
  //RiverFlow(elem, river, ctrl->riv_mode);               // Move to hydrol_river.cu
}


// Determine which layers does ET extract water from 
__device__ void EtExtract(elem_struct_d *elem, int tid)
{
	if (elem->ws.surfh[tid] >= DEPRSTG)
	{
		elem->wf.edir_surf[tid] = elem->wf.edir[tid];
		elem->wf.edir_unsat[tid] = 0.0;
		elem->wf.edir_gw[tid] = 0.0;
	}
	else if (elem->ws.gw[tid] >  elem->soil.depth[tid] - elem->soil.dinf[tid])
	{
		elem->wf.edir_surf[tid] = 0.0;
		elem->wf.edir_unsat[tid] = 0.0;
		elem->wf.edir_gw[tid] = elem->wf.edir[tid];
	}
	else
	{
		elem->wf.edir_surf[tid] = 0.0;
		elem->wf.edir_unsat[tid] = elem->wf.edir[tid];
		elem->wf.edir_gw[tid] = 0.0;
	}

	/* Source of transpiration */
	if (elem->ws.gw[tid] > elem->soil.depth[tid] - elem->ps.rzd[tid])
	{
		elem->wf.ett_unsat[tid] = 0.0;
		elem->wf.ett_gw[tid] = elem->wf.ett[tid];
	}
	else
	{
		elem->wf.ett_unsat[tid] = elem->wf.ett[tid];
		elem->wf.ett_gw[tid] = 0.0;
	}
}

__device__ double SurfH(double surfeqv)
{
	double          surfh;

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
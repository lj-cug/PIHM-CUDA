#include "pihm.h"
#include "pihm.cuh"

// Determine which layers does ET extract water from 
__device__ void EtExtract (int tid, 
#include "EtExtract_device_declarations.c" 
)
{
	if (elem_d_ws_surfh[tid] >= DEPRSTG)
	{
		elem_d_wf_edir_surf[tid] = elem_d_wf_edir[tid];
		elem_d_wf_edir_unsat[tid] = 0.0;
		elem_d_wf_edir_gw[tid] = 0.0;
	}
	else if (elem_d_ws_gw[tid] >  elem_d_soil_depth[tid] - elem_d_soil_dinf[tid])
	{
		elem_d_wf_edir_surf[tid] = 0.0;
		elem_d_wf_edir_unsat[tid] = 0.0;
		elem_d_wf_edir_gw[tid] = elem_d_wf_edir[tid];
	}
	else
	{
		elem_d_wf_edir_surf[tid] = 0.0;
		elem_d_wf_edir_unsat[tid] = elem_d_wf_edir[tid];
		elem_d_wf_edir_gw[tid] = 0.0;
	}

	/* Source of transpiration */
	if (elem_d_ws_gw[tid] > elem_d_soil_depth[tid] - elem_d_ps_rzd[tid])
	{
		elem_d_wf_ett_unsat[tid] = 0.0;
		elem_d_wf_ett_gw[tid] = elem_d_wf_ett[tid];
	}
	else
	{
		elem_d_wf_ett_unsat[tid] = elem_d_wf_ett[tid];
		elem_d_wf_ett_gw[tid] = 0.0;
	}
}

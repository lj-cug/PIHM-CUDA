#include "pihm.h"
#include "pihm.cuh"

/* 
2021.04.19  3个相邻单元上的计算，相互之间没有依赖性，可以使用不同的CUDA stream    
*/
__device__ void LateralFlow(int nelem, int surf_mode, int tid,
	const realtype *dhbydx, const realtype *dhbydy,
#include "LateralFlow_river_device_declarations.c"
)
{
	realtype        avg_sf;
	int             nabr;

// j =0  第1个相邻单元
	if (elem_d_nabr0[tid] > 0)
	{
		nabr = elem_d_nabr0[tid] - 1;   // 第tid号单元的相邻单元编号

		// Subsurface flow between triangular elements
		elem_d_wf_subsurf0[tid] = 
			SubFlowElemToElem0(tid, nabr,
#include "SubFlowElemToElem0_device_arguments.c" 
			                  );

		// Surface flux between triangular elements 
		avg_sf = 0.5 *
			(sqrt(dhbydx[tid] * dhbydx[tid] + dhbydy[tid] * dhbydy[tid]) +
			 sqrt(dhbydx[nabr - 1] * dhbydx[nabr - 1] + dhbydy[nabr - 1] * dhbydy[nabr - 1]));
		elem_d_wf_ovlflow0[tid] = 
			OvlFlowElemToElem0(tid, nabr, avg_sf, surf_mode,
#include "OvlFlowElemToElem0_device_arguments.c" 
			                   );
	}
	else if (elem_d_nabr0[tid] < 0)
	{
		// Do nothing. River-element interactions will be calculated
		// in river_flow.cu 
	}
	else    // Boundary condition flux 
	{
		BoundFluxElem0(elem_d_attrib_bc_type0[tid], tid,
#include "BoundFluxElem0_device_arguments.c" 
			          );
	}

// j =1  第2个相邻单元
	if (elem_d_nabr1[tid] > 0)
	{
		nabr = elem_d_nabr1[tid] - 1;   // 第tid号单元的相邻单元编号

		// Subsurface flow between triangular elements
		elem_d_wf_subsurf1[tid] =
			SubFlowElemToElem1(tid, nabr,
#include "SubFlowElemToElem1_device_arguments.c" 
			                   );

		// Surface flux between triangular elements 
		avg_sf = 0.5 *
			(sqrt(dhbydx[tid] * dhbydx[tid] + dhbydy[tid] * dhbydy[tid]) +
			 sqrt(dhbydx[nabr - 1] * dhbydx[nabr - 1] + dhbydy[nabr - 1] * dhbydy[nabr - 1]));
		elem_d_wf_ovlflow1[tid] = 
			OvlFlowElemToElem1(tid, nabr, avg_sf, surf_mode,
#include "OvlFlowElemToElem1_device_arguments.c" 
			                   );
	}
	else if (elem_d_nabr1[tid] < 0)
	{
		// Do nothing. River-element interactions will be calculated
		// in river_flow.cu 
	}
	else    // Boundary condition flux 
	{
		BoundFluxElem1(elem_d_attrib_bc_type1[tid], tid,
#include "BoundFluxElem1_device_arguments.c" 
			          );
	}


// j =2   第3个相邻单元
	if (elem_d_nabr2[tid] > 0)
	{
		nabr = elem_d_nabr2[tid] - 1;   // 第tid号单元的相邻单元编号

		// Subsurface flow between triangular elements
		elem_d_wf_subsurf2[tid] =
			SubFlowElemToElem2(tid, nabr,
#include "SubFlowElemToElem2_device_arguments.c" 
			                   );

		// Surface flux between triangular elements 
		avg_sf = 0.5 *
			(sqrt(dhbydx[tid] * dhbydx[tid] + dhbydy[tid] * dhbydy[tid]) +
			sqrt(dhbydx[nabr - 1] * dhbydx[nabr - 1] + dhbydy[nabr - 1] * dhbydy[nabr - 1]));
		elem_d_wf_ovlflow2[tid] = 
			OvlFlowElemToElem2(tid, nabr, avg_sf, surf_mode,
#include "OvlFlowElemToElem2_device_arguments.c" 
			                   );
	}
	else if (elem_d_nabr2[tid] < 0)
	{
		// Do nothing. River-element interactions will be calculated
		// in river_flow.cu 
	}
	else    // Boundary condition flux 
	{
		BoundFluxElem2(elem_d_attrib_bc_type2[tid], tid,
#include "BoundFluxElem2_device_arguments.c" 
			          );
	}

}


__device__ realtype AvgHsurf(realtype diff, realtype hsurf, realtype hnabr)
{
	realtype          avg_h;

	if (diff > 0.0)
	{
		if (hsurf > DEPRSTG)
		{
			avg_h = 1.0 * (hsurf - DEPRSTG);
		}
		else
		{
			avg_h = 0.0;
		}
	}
	else
	{
		if (hnabr > DEPRSTG)
		{
			avg_h = 1.0 * (hnabr - DEPRSTG);
		}
		else
		{
			avg_h = 0.0;
		}
	}
	return avg_h;
}


__device__ realtype AvgH(realtype diff, realtype hsub, realtype hnabr)
{
	realtype          avg_h = 0.0;

	if (diff > 0.0)
	{
		if (hsub > 0.0)
		{
			avg_h = hsub;
		}
	}
	else
	{
		if (hnabr > 0.0)
		{
			avg_h = hnabr;
		}
	}
	return avg_h;
}

__device__ realtype EffKh(int ind, realtype gw, 
#include "EffKh_device_declarations.c" 
)
{
	realtype          k1, k2;
	realtype          d1, d2;

	gw = (gw > 0.0) ? gw : 0.0;
	if (gw > elem_d_soil_depth[ind] - elem_d_soil_dmac[ind])
	{
		k1 = elem_d_soil_kmach[ind] * elem_d_soil_areafv[ind] +
			    elem_d_soil_ksath[ind] * (1.0 - elem_d_soil_areafv[ind]);
		k2 = elem_d_soil_ksath[ind];

		if (gw > elem_d_soil_depth[ind])
		{
			d1 = elem_d_soil_dmac[ind];
			d2 = elem_d_soil_depth[ind] - elem_d_soil_dmac[ind];
		}
		else
		{
			d1 = gw - (elem_d_soil_depth[ind] - elem_d_soil_dmac[ind]);
			d2 = elem_d_soil_depth[ind] - elem_d_soil_dmac[ind];
		}

		return (k1 * d1 + k2 * d2) / (d1 + d2);
	}
	else
	{
		return elem_d_soil_ksath[ind];
	}
}


__device__ realtype OverLandFlow(realtype avg_h, realtype grad_h, 
	                realtype avg_sf, realtype crossa, realtype avg_rough)
{
	return crossa * pow(avg_h, 0.6666667) * grad_h / (sqrt(avg_sf) * avg_rough);
}


__device__ realtype SubFlowElemToElem0(int tid, int nabr,
#include "SubFlowElemToElem0_device_declarations.c" 
)
{
	realtype          diff_h;
	realtype          avg_h;
	realtype          grad_h;
	realtype          effk, effk_nabr;
	realtype          avg_ksat;

	/*
	* Subsurface lateral flux calculation between triangular
	* elements
	*/
	diff_h = (elem_d_ws_gw[tid] + elem_d_topo_zmin[tid]) - (elem_d_ws_gw[nabr] + elem_d_topo_zmin[nabr]);
	avg_h = AvgH(diff_h, elem_d_ws_gw[tid], elem_d_ws_gw[nabr]);
	grad_h = diff_h / elem_d_topo_nabrdist0[tid];

	/* Take into account macropore effect */
	effk = EffKh(tid, elem_d_ws_gw[tid], 
#include "EffKh_device_arguments.c" 
		);
	effk_nabr = EffKh(nabr, elem_d_ws_gw[nabr], 
#include "EffKh_device_arguments.c" 
		);
	avg_ksat = 0.5 * (effk + effk_nabr);

	/* Groundwater flow modeled by Darcy's Law */
	return avg_ksat * grad_h * avg_h * elem_d_topo_edge0[tid];
}

__device__ realtype SubFlowElemToElem1(int tid, int nabr,
#include "SubFlowElemToElem1_device_declarations.c" 
)
{
	realtype          diff_h;
	realtype          avg_h;
	realtype          grad_h;
	realtype          effk, effk_nabr;
	realtype          avg_ksat;

	/*
	* Subsurface lateral flux calculation between triangular
	* elements
	*/
	diff_h = (elem_d_ws_gw[tid] + elem_d_topo_zmin[tid]) - (elem_d_ws_gw[nabr] + elem_d_topo_zmin[nabr]);
	avg_h = AvgH(diff_h, elem_d_ws_gw[tid], elem_d_ws_gw[nabr]);
	grad_h = diff_h / elem_d_topo_nabrdist1[tid];

	/* Take into account macropore effect */
	effk = EffKh(tid, elem_d_ws_gw[tid],
#include "EffKh_device_arguments.c" 
		);
	effk_nabr = EffKh(nabr, elem_d_ws_gw[nabr],
#include "EffKh_device_arguments.c" 
		);
	avg_ksat = 0.5 * (effk + effk_nabr);

	/* Groundwater flow modeled by Darcy's Law */
	return avg_ksat * grad_h * avg_h * elem_d_topo_edge1[tid];
}


__device__ realtype SubFlowElemToElem2(int tid, int nabr,
#include "SubFlowElemToElem2_device_declarations.c" 
)
{
	realtype          diff_h;
	realtype          avg_h;
	realtype          grad_h;
	realtype          effk, effk_nabr;
	realtype          avg_ksat;

	/*
	* Subsurface lateral flux calculation between triangular
	* elements
	*/
	diff_h = (elem_d_ws_gw[tid] + elem_d_topo_zmin[tid]) - (elem_d_ws_gw[nabr] + elem_d_topo_zmin[nabr]);
	avg_h = AvgH(diff_h, elem_d_ws_gw[tid], elem_d_ws_gw[nabr]);
	grad_h = diff_h / elem_d_topo_nabrdist2[tid];

	/* Take into account macropore effect */
	effk = EffKh(tid, elem_d_ws_gw[tid],
#include "EffKh_device_arguments.c" 
		         );
	effk_nabr = EffKh(nabr, elem_d_ws_gw[nabr],
#include "EffKh_device_arguments.c" 
		             );
	avg_ksat = 0.5 * (effk + effk_nabr);

	/* Groundwater flow modeled by Darcy's Law */
	return avg_ksat * grad_h * avg_h * elem_d_topo_edge2[tid];
}


__device__ realtype OvlFlowElemToElem0(int tid, int nabr,
	                                   realtype avg_sf, int surf_mode,
#include "OvlFlowElemToElem0_device_declarations.c" 
)
{
	realtype          diff_h;
	realtype          avg_h;
	realtype          grad_h;
	realtype          avg_rough;
	realtype          crossa;

	diff_h = (surf_mode == KINEMATIC) ?
		elem_d_topo_zmax[tid] - elem_d_topo_zmax[nabr] :
		(elem_d_ws_surfh[tid] + elem_d_topo_zmax[tid]) - (elem_d_ws_surfh[nabr] + elem_d_topo_zmax[nabr]);

	avg_h = AvgHsurf(diff_h, elem_d_ws_surfh[tid], elem_d_ws_surfh[nabr]);
	grad_h = diff_h / elem_d_topo_nabrdist0[tid];
	if (surf_mode == KINEMATIC)
	{
		avg_sf = (grad_h > 0.0) ? grad_h : GRADMIN;
	}
	else
	{
		avg_sf = (avg_sf > GRADMIN) ? avg_sf : GRADMIN;
	}
	/* Weighting needed */
	avg_rough = 0.5 * (elem_d_lc_rough[tid] + elem_d_lc_rough[nabr]);
	crossa = avg_h * elem_d_topo_edge0[tid];

	return  OverLandFlow(avg_h, grad_h, avg_sf, crossa, avg_rough);
}


__device__ realtype OvlFlowElemToElem1(int tid, int nabr,
	realtype avg_sf, int surf_mode,
#include "OvlFlowElemToElem1_device_declarations.c" 
)
{
	realtype          diff_h;
	realtype          avg_h;
	realtype          grad_h;
	realtype          avg_rough;
	realtype          crossa;

	diff_h = (surf_mode == KINEMATIC) ?
		elem_d_topo_zmax[tid] - elem_d_topo_zmax[nabr] :
		(elem_d_ws_surfh[tid] + elem_d_topo_zmax[tid]) - (elem_d_ws_surfh[nabr] + elem_d_topo_zmax[nabr]);

	avg_h = AvgHsurf(diff_h, elem_d_ws_surfh[tid], elem_d_ws_surfh[nabr]);
	grad_h = diff_h / elem_d_topo_nabrdist1[tid];
	if (surf_mode == KINEMATIC)
	{
		avg_sf = (grad_h > 0.0) ? grad_h : GRADMIN;
	}
	else
	{
		avg_sf = (avg_sf > GRADMIN) ? avg_sf : GRADMIN;
	}
	/* Weighting needed */
	avg_rough = 0.5 * (elem_d_lc_rough[tid] + elem_d_lc_rough[nabr]);
	crossa = avg_h * elem_d_topo_edge1[tid];

	return  OverLandFlow(avg_h, grad_h, avg_sf, crossa, avg_rough);
}


__device__ realtype OvlFlowElemToElem2(int tid, int nabr,
	realtype avg_sf, int surf_mode,
#include "OvlFlowElemToElem2_device_declarations.c" 
)
{
	realtype          diff_h;
	realtype          avg_h;
	realtype          grad_h;
	realtype          avg_rough;
	realtype          crossa;

	diff_h = (surf_mode == KINEMATIC) ?
		elem_d_topo_zmax[tid] - elem_d_topo_zmax[nabr] :
		(elem_d_ws_surfh[tid] + elem_d_topo_zmax[tid]) - (elem_d_ws_surfh[nabr] + elem_d_topo_zmax[nabr]);

	avg_h = AvgHsurf(diff_h, elem_d_ws_surfh[tid], elem_d_ws_surfh[nabr]);
	grad_h = diff_h / elem_d_topo_nabrdist2[tid];
	if (surf_mode == KINEMATIC)
	{
		avg_sf = (grad_h > 0.0) ? grad_h : GRADMIN;
	}
	else
	{
		avg_sf = (avg_sf > GRADMIN) ? avg_sf : GRADMIN;
	}
	/* Weighting needed */
	avg_rough = 0.5 * (elem_d_lc_rough[tid] + elem_d_lc_rough[nabr]);
	crossa = avg_h * elem_d_topo_edge2[tid];

	return  OverLandFlow(avg_h, grad_h, avg_sf, crossa, avg_rough);
}


__device__ void BoundFluxElem0(int bc_type, int tid,
#include "BoundFluxElem0_device_declarations.c" 
)
{
	realtype          diff_h;
	realtype          avg_h;
	realtype          effk;
	realtype          avg_ksat;
	realtype          grad_h;

	/* No flow (natural) boundary condition is default */
	if (bc_type == NO_FLOW)
	{
		elem_d_wf_ovlflow0[tid] = 0.0;
		elem_d_wf_subsurf0[tid] = 0.0;
	}
	/* Note: ideally different boundary conditions need to be
	* incorporated for surf and subsurf respectively */
	else if (bc_type > 0)
	{
		/* Note: the formulation assumes only Dirichlet TS right now */
		/* note the assumption here is no flow for surface */
		elem_d_wf_ovlflow0[tid] = 0.0;

		diff_h = elem_d_ws_gw[tid] + elem_d_topo_zmin[tid] - elem_d_bc_head0[tid];
		avg_h = AvgH(diff_h, elem_d_ws_gw[tid], elem_d_bc_head0[tid] - elem_d_topo_zmin[tid]);

		/* Minimum distance from circumcenter to the edge of the triangle
		* on which boundary condition is defined */
		effk = EffKh(tid, elem_d_ws_gw[tid], 
#include "EffKh_device_arguments.c" 
			         );
		avg_ksat = effk;
		grad_h = diff_h / elem_d_topo_nabrdist0[tid];
		elem_d_wf_subsurf0[tid] = avg_ksat * grad_h * avg_h * elem_d_topo_edge0[tid];
	}
	else
	{
		/* Neumann bc (note: md->ele[i].bc[j] value has to be
		* = 2+(index of Neumann boundary ts) */
		elem_d_wf_ovlflow0[tid] = 0.0;
		/* Negative sign is added so the positive numbers in forcing time series
		* represents source */
		elem_d_wf_subsurf0[tid] = - elem_d_bc_flux0[tid];
	}
}


__device__ void BoundFluxElem1(int bc_type, int tid,
#include "BoundFluxElem1_device_declarations.c" 
)
{
	realtype          diff_h;
	realtype          avg_h;
	realtype          effk;
	realtype          avg_ksat;
	realtype          grad_h;

	/* No flow (natural) boundary condition is default */
	if (bc_type == NO_FLOW)
	{
		elem_d_wf_ovlflow1[tid] = 0.0;
		elem_d_wf_subsurf1[tid] = 0.0;
	}
	/* Note: ideally different boundary conditions need to be
	* incorporated for surf and subsurf respectively */
	else if (bc_type > 0)
	{
		/* Note: the formulation assumes only Dirichlet TS right now */
		/* note the assumption here is no flow for surface */
		elem_d_wf_ovlflow1[tid] = 0.0;

		diff_h = elem_d_ws_gw[tid] + elem_d_topo_zmin[tid] - elem_d_bc_head1[tid];
		avg_h = AvgH(diff_h, elem_d_ws_gw[tid], elem_d_bc_head1[tid] - elem_d_topo_zmin[tid]);

		/* Minimum distance from circumcenter to the edge of the triangle
		* on which boundary condition is defined */
		effk = EffKh(tid, elem_d_ws_gw[tid],
#include "EffKh_device_arguments.c" 
			);
		avg_ksat = effk;
		grad_h = diff_h / elem_d_topo_nabrdist1[tid];
		elem_d_wf_subsurf1[tid] = avg_ksat * grad_h * avg_h * elem_d_topo_edge1[tid];
	}
	else
	{
		/* Neumann bc (note: md->ele[i].bc[j] value has to be
		* = 2+(index of Neumann boundary ts) */
		elem_d_wf_ovlflow1[tid] = 0.0;
		/* Negative sign is added so the positive numbers in forcing time series
		* represents source */
		elem_d_wf_subsurf1[tid] = -elem_d_bc_flux1[tid];
	}
}

__device__ void BoundFluxElem2(int bc_type, int tid,
#include "BoundFluxElem2_device_declarations.c" 
)
{
	realtype          diff_h;
	realtype          avg_h;
	realtype          effk;
	realtype          avg_ksat;
	realtype          grad_h;

	/* No flow (natural) boundary condition is default */
	if (bc_type == NO_FLOW)
	{
		elem_d_wf_ovlflow2[tid] = 0.0;
		elem_d_wf_subsurf2[tid] = 0.0;
	}
	/* Note: ideally different boundary conditions need to be
	* incorporated for surf and subsurf respectively */
	else if (bc_type > 0)
	{
		/* Note: the formulation assumes only Dirichlet TS right now */
		/* note the assumption here is no flow for surface */
		elem_d_wf_ovlflow2[tid] = 0.0;

		diff_h = elem_d_ws_gw[tid] + elem_d_topo_zmin[tid] - elem_d_bc_head2[tid];
		avg_h = AvgH(diff_h, elem_d_ws_gw[tid], elem_d_bc_head2[tid] - elem_d_topo_zmin[tid]);

		/* Minimum distance from circumcenter to the edge of the triangle
		* on which boundary condition is defined */
		effk = EffKh(tid, elem_d_ws_gw[tid],
#include "EffKh_device_arguments.c" 
			);
		avg_ksat = effk;
		grad_h = diff_h / elem_d_topo_nabrdist2[tid];
		elem_d_wf_subsurf2[tid] = avg_ksat * grad_h * avg_h * elem_d_topo_edge2[tid];
	}
	else
	{
		/* Neumann bc (note: md->ele[i].bc[j] value has to be
		* = 2+(index of Neumann boundary ts) */
		elem_d_wf_ovlflow2[tid] = 0.0;
		/* Negative sign is added so the positive numbers in forcing time series
		* represents source */
		elem_d_wf_subsurf2[tid] = -elem_d_bc_flux2[tid];
	}
}



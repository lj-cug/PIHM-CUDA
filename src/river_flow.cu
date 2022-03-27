#include "pihm.h"
#include "pihm.cuh"
/*
  2021.04.16:  仔细检查每个设备函数的形参，尽可能降低个数, 进行优化.
  2021.04.22:  检查一下是否修改了河道两侧单元的相关数组数值？
*/

/*  Hydrologic processes related with river flow 
    Moved from hydrol.c
*/
__global__ void RiverFlow(const int nriver, const int riv_mode,
#include "RiverFlow_device_declarations.c"  
)
{
	int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid < nriver){

		// 第1个  for (i = 0; i < nriver; i++)
		int             left, right;
		realtype        effk_nabr;
		realtype        effk;
		int             down;      // 下段河道的分段编号    for SoA

		if (river_d_down[tid] > 0)
		{
			//
			// Boundary conditions
			//
			// When a downstream segment is present, boundary conditions are
			// always applied to the upstream node
			//
			if (river_d_attrib_riverbc_type[tid] != 0)
			{
				river_d_wf_rivflow0[tid] +=
					BoundFluxRiver(river_d_attrib_riverbc_type[tid], tid,
#include "BoundFluxRiver_device_arguments.c" 
					);
			}

			down = river_d_down[tid] - 1;       // for SoA

			// Channel flow between river-river segments
			river_d_wf_rivflow1[tid] =
				ChanFlowRiverToRiver(tid, down, riv_mode,
#include "ChanFlowRiverToRiver_device_arguments.c" 
				);

			// Subsurface flow between river-river segments
			left = river_d_leftele[tid] - 1;   //左边的单元编号
			right = river_d_rightele[tid] - 1;  //右边的单元编号
			effk = 0.5 *
				(EffKh(left, elem_d_ws_gw[left],
#include "EffKh_device_arguments.c" 
				) +
				EffKh(right, elem_d_ws_gw[right],
#include "EffKh_device_arguments.c" 
				)
				);

			left = river_d_leftele[down] - 1;   //左边的单元编号
			right = river_d_rightele[down] - 1; //右边的单元编号
			effk_nabr = 0.5 *
				(EffKh(left, elem_d_ws_gw[left],
#include "EffKh_device_arguments.c" 
				) +
				EffKh(right, elem_d_ws_gw[right],
#include "EffKh_device_arguments.c" 
				)
				);

			river_d_wf_rivflow9[tid] =
				SubFlowRiverToRiver(effk, effk_nabr, tid, down,
#include "SubFlowRiverToRiver_device_arguments.c" 
				);
		}

		else
		{
			// Outlet flux
			river_d_wf_rivflow1[tid] =
				OutletFlux(tid, river_d_down[tid],
#include "OutletFlux_device_arguments.c" 
				);

			/* Note: boundary condition for subsurface element can be changed.
			* Assumption: no flow condition */
			river_d_wf_rivflow9[tid] = 0.0;
		}

		// Flux between river segments and triangular elements
		left = river_d_leftele[tid] - 1;    // 第tid段河段左侧单元的编号
		right = river_d_rightele[tid] - 1;  // 第tid段河段右侧单元的编号

		RiverToElem(tid, left, right,
#include "RiverToElem_device_arguments.c" 		
			);  // tid -河段编号; left, right - 单元编号

		// Flux between river channel and subsurface
		river_d_wf_rivflow6[tid] =
			ChanLeak(tid,
#include "ChanLeak_device_arguments.c" 
			);

		// 第2个  for (i = 0; i < nriver; i++)
		/*
		* Accumulate to get in-flow for down segments
		*
		* NOTE: Upstream flux summation must be calculated outside OMP to avoid
		* different threads accessing the same variable at the same time
		*/
		if (river_d_down[tid] > 0)
		{
			down = river_d_down[tid] - 1;

			river_d_wf_rivflow0[down]  -= river_d_wf_rivflow1[tid];
			river_d_wf_rivflow10[down] -= river_d_wf_rivflow9[tid];
		}

	} // if (tid < nriver)
}


__device__ realtype BoundFluxRiver(int riverbc_type, int tid, 
#include "BoundFluxRiver_device_declarations.c" 
)
{
	realtype          total_h;
	realtype          total_h_down;
	realtype          distance;
	realtype          grad_h;
	realtype          avg_h;
	realtype          avg_perim;
	realtype          crossa;
	realtype          flux = 0.0;

	if (riverbc_type > 0)
	{
		// Dirichlet boundary condition 
		total_h = river_d_ws_stage[tid] + river_d_topo_zbed[tid];
		total_h_down = river_d_bc_head[tid];
		distance = 0.5 * river_d_shp_length[tid];
		grad_h = (total_h - total_h_down) / distance;
		avg_h = AvgH(grad_h, river_d_ws_stage[tid],
			((river_d_bc_head[tid] - river_d_topo_zbed[tid] > 0.0) ?
			 (river_d_bc_head[tid] - river_d_topo_zbed[tid]) : 0.0));
		avg_perim = RiverPerim(river_d_shp_intrpl_ord[tid], river_d_ws_stage[tid], river_d_shp_coeff[tid]);
		crossa = RiverCroSectArea(river_d_shp_intrpl_ord[tid], river_d_ws_stage[tid], river_d_shp_coeff[tid]);
		avg_h = (avg_perim == 0.0) ? 0.0 : (crossa / avg_perim);
		flux = OverLandFlow(avg_h, grad_h, grad_h, crossa, river_d_matl_rough[tid]);
	}
	else if (riverbc_type < 0)
	{
		// Neumann boundary condition 
		flux = -river_d_bc_flux[tid];
	}
 	
	return flux;
}


__device__ realtype ChanFlowRiverToRiver(int tid, int down, int riv_mode,
#include "ChanFlowRiverToRiver_device_declarations.c" 	
)
{
	realtype          total_h;
	realtype          perim;
	realtype          total_h_down;
	realtype          perim_down;
	realtype          avg_perim;
	realtype          avg_rough;
	realtype          distance;
	realtype          diff_h;
	realtype          grad_h;
	realtype          avg_sf;
	realtype          crossa;
	realtype          crossa_down;
	realtype          avg_crossa;
	realtype          avg_h;

	total_h = river_d_ws_stage[tid] + river_d_topo_zbed[tid];
	perim =
		RiverPerim(river_d_shp_intrpl_ord[tid], river_d_ws_stage[tid], river_d_shp_coeff[tid]);

	total_h_down = river_d_ws_stage[down] + river_d_topo_zbed[down];
	perim_down =
		RiverPerim(river_d_shp_intrpl_ord[down], river_d_ws_stage[down], river_d_shp_coeff[down]);

	avg_perim = (perim + perim_down) / 2.0;
	avg_rough = (river_d_matl_rough[tid] + river_d_matl_rough[down]) / 2.0;
	distance = 0.5 * (river_d_shp_length[tid] + river_d_shp_length[down]);


	diff_h = (riv_mode == KINEMATIC) ?  
		(river_d_topo_zbed[tid] - river_d_topo_zbed[tid]) : (total_h - total_h_down);

	grad_h = diff_h / distance;
	avg_sf = (grad_h > 0.0) ? grad_h : RIVGRADMIN;
	crossa =
		RiverCroSectArea(river_d_shp_intrpl_ord[tid], river_d_ws_stage[tid], river_d_shp_coeff[tid]);		
	crossa_down =
		RiverCroSectArea(river_d_shp_intrpl_ord[down], river_d_ws_stage[down], river_d_shp_coeff[down]);
	avg_crossa = 0.5 * (crossa + crossa_down);
	avg_h = (avg_perim == 0.0) ? 0.0 : (avg_crossa / avg_perim);

	return OverLandFlow(avg_h, grad_h, avg_sf, crossa, avg_rough);
}


__device__ realtype SubFlowRiverToRiver(realtype effk, realtype effk_nabr, int tid, int down,
#include "SubFlowRiverToRiver_device_declarations.c" 
)
{
	realtype          total_h;
	realtype          total_h_down;
	realtype          avg_wid;
	realtype          diff_h;
	realtype          avg_h;
	realtype          distance;
	realtype          grad_h;
//	realtype          aquifer_depth;
	realtype          avg_ksat;

	/* Lateral flux calculation between element beneath river (ebr) * and ebr */
	total_h = river_d_ws_gw[tid] + river_d_topo_zmin[tid];
	total_h_down = river_d_ws_gw[down] + river_d_topo_zmin[down];
	avg_wid = (river_d_shp_width[tid] + river_d_shp_width[down]) / 2.0;
	diff_h = total_h - total_h_down;
	avg_h = AvgH(diff_h, river_d_ws_gw[tid], river_d_ws_gw[down]);
	distance = 0.5 * (river_d_shp_length[tid] + river_d_shp_length[down]);
	grad_h = diff_h / distance;
//	aquifer_depth = river_d_topo_zbed[tid] - river_d_topo_zmin[tid];  // not used!
#if defined(_ARITH_)
	avg_ksat = 0.5 * (effk + effk_nabr);
#else
	avg_ksat = 2.0 / (1.0 / effk + 1.0 / effk_nabr);
#endif
	/* Groundwater flow modeled by Darcy's law */
	return avg_ksat * grad_h * avg_h * avg_wid;
}


__device__ realtype OutletFlux(int tid, int down,
#include "OutletFlux_device_declarations.c" 
)
{
	realtype          total_h;
	realtype          total_h_down;
	realtype          distance;
	realtype          grad_h;
	realtype          avg_h;
	realtype          avg_perim;
	realtype          crossa;
	realtype          discharge = 0.0;

	switch (down)
	{
	case DIRICHLET:
		/* Dirichlet boundary condition */
		total_h = river_d_ws_stage[tid] + river_d_topo_zbed[tid];
		total_h_down = river_d_bc_head[tid];
		distance = 0.5 * river_d_shp_length[tid];
		grad_h = (total_h - total_h_down) / distance;
		avg_h = AvgH(grad_h, river_d_ws_stage[tid],
			((river_d_bc_head[tid] - (river_d_topo_node_zmax[tid] - river_d_shp_depth[tid]) > 0.0) ?
			river_d_bc_head[tid] - (river_d_topo_node_zmax[tid] - river_d_shp_depth[tid]) : 0.0));
		avg_perim = RiverPerim(river_d_shp_intrpl_ord[tid], river_d_ws_stage[tid], river_d_shp_coeff[tid]);
		crossa = RiverCroSectArea(river_d_shp_intrpl_ord[tid], river_d_ws_stage[tid], river_d_shp_coeff[tid]);
		avg_h = (avg_perim == 0.0) ? 0.0 : (crossa / avg_perim);
		discharge =
			OverLandFlow(avg_h, grad_h, grad_h, crossa, river_d_matl_rough[tid]);
		break;
	case NEUMANN:
		/* Neumann boundary condition */
		discharge = -river_d_bc_flux[tid];
		break;
	case ZERO_DPTH_GRAD:
		/* Zero-depth-gradient boundary conditions */
		distance = 0.5 * river_d_shp_length[tid];
		grad_h = (river_d_topo_zbed[tid] - (river_d_topo_node_zmax[tid] - river_d_shp_depth[tid])) / distance;
		avg_h = river_d_ws_stage[tid];
		avg_perim = RiverPerim(river_d_shp_intrpl_ord[tid], river_d_ws_stage[tid], river_d_shp_coeff[tid]);
		crossa = RiverCroSectArea(river_d_shp_intrpl_ord[tid], river_d_ws_stage[tid], river_d_shp_coeff[tid]);
		discharge = sqrt(grad_h) * crossa * ((avg_perim > 0.0) ?
			pow(crossa / avg_perim, 2.0 / 3.0) : 0.0) / river_d_matl_rough[tid];
		break;
	case CRIT_DPTH:
		/* Critical depth boundary conditions */
		crossa = RiverCroSectArea(river_d_shp_intrpl_ord[tid], river_d_ws_stage[tid], river_d_shp_coeff[tid]);
		discharge = crossa * sqrt(GRAV * river_d_ws_stage[tid]);
		break;
	default:
/*		PIHMprintf(VL_ERROR,
			"Error: River routing boundary condition type (%d) "
			"is not recognized.\n", down);
		PIHMexit(EXIT_FAILURE);   */
		break;
	}

	return discharge;
}


__device__ void RiverToElem(int tid, int left, int right, 
#include "RiverToElem_device_declarations.c" 
)
{
	realtype    effk_left, effk_right;

	/* Lateral surface flux calculation between river-triangular element */
	if (river_d_leftele[tid] > 0)
	{
		river_d_wf_rivflow2[tid] = OvlFlowElemToRiver(tid, left,
#include "OvlFlowElemToRiver_device_arguments.c" 
			                                          );
	}


	if (river_d_rightele[tid] > 0)
	{
		river_d_wf_rivflow3[tid] = OvlFlowElemToRiver(tid, right,
#include "OvlFlowElemToRiver_device_arguments.c" 
			);
	}

	effk_left = EffKh(left, elem_d_ws_gw[left], 
#include "EffKh_device_arguments.c" 
		);
	effk_right = EffKh(right, elem_d_ws_gw[right],
#include "EffKh_device_arguments.c" 
		);

	/* Lateral subsurface flux calculation between river-triangular element */
	if (river_d_leftele[tid] > 0)
	{
		river_d_wf_rivflow4[tid] =
			ChanFlowElemToRiver(left, effk_left, tid, river_d_topo_dist_left[tid], 
#include "ChanFlowElemToRiver_device_arguments.c" 
			);
	}
	if (river_d_rightele[tid] > 0)
	{
		river_d_wf_rivflow5[tid] =
			ChanFlowElemToRiver(right, effk_right, tid,river_d_topo_dist_right[tid],
#include "ChanFlowElemToRiver_device_arguments.c" 
			);			
	}

	/* Lateral flux between rectangular element (beneath river) and triangular
	* element */
	if (river_d_leftele[tid] > 0)
	{
		river_d_wf_rivflow7[tid] =
			SubFlowElemToRiver(left, effk_left, tid,
			         0.5 * (effk_left + effk_right), river_d_topo_dist_left[tid], 
#include "SubFlowElemToRiver_device_arguments.c" 
					 );
	}
	if (river_d_rightele[tid] > 0)
	{
		//river_d_wf_rivflow[RIGHT_AQUIF2AQUIF][tid] =
		river_d_wf_rivflow8[tid] =
			SubFlowElemToRiver(right, effk_right, tid,
			         0.5 * (effk_left + effk_right), river_d_topo_dist_right[tid],
#include "SubFlowElemToRiver_device_arguments.c" 
					  );
	}

	/* Replace flux term */
	/* Left */
	if (elem_d_nabr0[left] == -river_d_ind[tid])
	{
		elem_d_wf_ovlflow0[left] = -river_d_wf_rivflow2[tid];
		elem_d_wf_subsurf0[left] = -(river_d_wf_rivflow4[tid] + river_d_wf_rivflow7[tid]);			
		goto flag0;
	}
	if (elem_d_nabr1[left] == -river_d_ind[tid])
	{
		elem_d_wf_ovlflow1[left] = -river_d_wf_rivflow2[tid];
		elem_d_wf_subsurf1[left] = -(river_d_wf_rivflow4[tid] + river_d_wf_rivflow7[tid]);
		goto flag0;
	}
	if (elem_d_nabr2[left] == -river_d_ind[tid])
	{
		elem_d_wf_ovlflow2[left] = -river_d_wf_rivflow2[tid];
		elem_d_wf_subsurf2[left] = -(river_d_wf_rivflow4[tid] + river_d_wf_rivflow7[tid]);
		goto flag0;
	}
  
	/* Right */
flag0: if (elem_d_nabr0[right] == -river_d_ind[tid])
	{
		elem_d_wf_ovlflow0[right] = -river_d_wf_rivflow3[tid];
		elem_d_wf_subsurf0[right] = -(river_d_wf_rivflow5[tid] + river_d_wf_rivflow8[tid]);			
		goto flag1;
	}
	if (elem_d_nabr1[right] == -river_d_ind[tid])
	{
		elem_d_wf_ovlflow1[right] = -river_d_wf_rivflow3[tid];
		elem_d_wf_subsurf1[right] = -(river_d_wf_rivflow5[tid] + river_d_wf_rivflow8[tid]);
		goto flag1;
	}
	if (elem_d_nabr2[right] == -river_d_ind[tid])
	{
		elem_d_wf_ovlflow2[right] = -river_d_wf_rivflow3[tid];
		elem_d_wf_subsurf2[right] = -(river_d_wf_rivflow5[tid] + river_d_wf_rivflow8[tid]);
		goto flag1;
	}

flag1:  return;
}


__device__ realtype OvlFlowElemToRiver(int tid, int ind, 
#include "OvlFlowElemToRiver_device_declarations.c" 
)
// 这里, tid是河段编号, ind是单元编号, 计算河道与相邻单元的水体交换
{
	realtype          zbank;
	realtype          flux;
	realtype          elem_h;
	realtype          rivseg_h;

	zbank = (river_d_topo_zmax[tid] > elem_d_topo_zmax[ind]) ?
		     river_d_topo_zmax[tid] : elem_d_topo_zmax[ind];

	elem_h = elem_d_topo_zmax[ind] + elem_d_ws_surfh[ind];
	rivseg_h = river_d_topo_zbed[tid] + river_d_ws_stage[tid];

	/*
	* Panday and Hyakorn 2004 AWR Eqs. (23) and (24)
	*/
	if (rivseg_h > elem_h)
	{
		if (elem_h > zbank)
		{
			/* Submerged weir */
			flux = river_d_matl_cwr[tid] * 2.0 * sqrt(2.0 * GRAV) *
				river_d_shp_length[tid] * sqrt(rivseg_h - elem_h) *
				(rivseg_h - zbank) / 3.0;
		}
		else
		{
			if (zbank < rivseg_h)
			{
				/* Free-flowing weir */
				flux = river_d_matl_cwr[tid] * 2.0 * sqrt(2.0 * GRAV) *
					river_d_shp_length[tid] * sqrt(rivseg_h - zbank) *
					(rivseg_h - zbank) / 3.0;
			}
			else
			{
				flux = 0.0;
			}
		}
	}

	else if (elem_d_ws_surfh[ind] > DEPRSTG)
	{
		if (rivseg_h > zbank)
		{
			/* Submerged weir */
			flux = -river_d_matl_cwr[tid] * 2.0 * sqrt(2.0 * GRAV) *
				river_d_shp_length[tid] * sqrt(elem_h - rivseg_h) *
				(elem_h - zbank) / 3.0;
		}
		else
		{
			if (zbank < elem_h)
			{
				/* Free-flowing weir */
				flux = -river_d_matl_cwr[tid] * 2.0 * sqrt(2.0 * GRAV) *
					river_d_shp_length[tid] * sqrt(elem_h - zbank) *
					(elem_h - zbank) / 3.0;
			}
			else
			{
				flux = 0.0;
			}
		}
	}
	else
	{
		flux = 0.0;
	}
	return flux;
}


__device__ realtype ChanFlowElemToRiver(int ind, realtype effk, int tid, realtype distance, 
#include "ChanFlowElemToRiver_device_declarations.c" 
)
// 这里, tid是河段编号, ind是单元编号, 计算河道与相邻单元的水体交换
{
	realtype          diff_h;
	realtype          avg_h;
	realtype          grad_h;
	realtype          avg_ksat;

	diff_h = (river_d_ws_stage[tid] + river_d_topo_zbed[tid]) -
		(elem_d_ws_gw[ind] + elem_d_topo_zmin[ind]);

	/* This is head in neighboring cell representation */
	if (elem_d_topo_zmin[ind] > river_d_topo_zbed[tid])
	{
		avg_h = elem_d_ws_gw[ind];
	}
	else if (elem_d_topo_zmin[ind] + elem_d_ws_gw[ind] > river_d_topo_zbed[tid])
	{
		avg_h = elem_d_topo_zmin[ind] + elem_d_ws_gw[ind] - river_d_topo_zbed[tid];
	}
	else
	{
		avg_h = 0.0;
	}
	avg_h = AvgH(diff_h, river_d_ws_stage[tid], avg_h);

	grad_h = diff_h / distance;

	avg_ksat = 0.5 * (effk + river_d_matl_ksath[tid]);

	return  river_d_shp_length[tid] * avg_ksat * grad_h * avg_h;
}

__device__ realtype SubFlowElemToRiver(int ind, realtype effk, 
	                                   int tid, realtype effk_riv, realtype distance,
#include "SubFlowElemToRiver_device_declarations.c" 
)
// 这里, tid是河段编号, ind是单元编号, 计算河道与相邻单元的水体交换
{
	realtype          diff_h;
	realtype          avg_h;
//	realtype          aquifer_depth;    // not used!
	realtype          avg_ksat;
	realtype          grad_h;

	diff_h = (river_d_ws_gw[tid] + river_d_topo_zmin[tid]) -
		(elem_d_ws_gw[ind] + elem_d_topo_zmin[ind]);

	/* This is head in neighboring cell represention */
	if (elem_d_topo_zmin[ind] > river_d_topo_zbed[tid])
	{
		avg_h = 0.0;
	}
	else if (elem_d_topo_zmin[ind] + elem_d_ws_gw[ind] > river_d_topo_zbed[tid])
	{
		avg_h = river_d_topo_zbed[tid] - elem_d_topo_zmin[ind];
	}
	else
	{
		avg_h = elem_d_ws_gw[ind];
	}
	avg_h = AvgH(diff_h, river_d_ws_gw[tid], avg_h);

#if defined(_ARITH_)
	avg_ksat = 0.5 * (effk + effk_riv);
#else
	avg_ksat = 2.0 / (1.0 / effk + 1.0 / effk_riv);
#endif
	grad_h = diff_h / distance;

	return river_d_shp_length[tid] * avg_ksat * grad_h * avg_h;
}


__device__ realtype ChanLeak(int tid, 
#include "ChanLeak_device_declarations.c" 
)
{
	realtype          diff_h;
	realtype          grad_h;

	if (river_d_topo_zbed[tid] - (river_d_ws_gw[tid] + river_d_topo_zmin[tid]) > 0.0)
	{
		diff_h = river_d_ws_stage[tid];
	}
	else
	{
		diff_h = river_d_ws_stage[tid] + river_d_topo_zbed[tid] - 
			     (river_d_ws_gw[tid] + river_d_topo_zmin[tid]);
	}

	grad_h = diff_h / river_d_matl_bedthick[tid];

	return river_d_matl_ksatv[tid] * river_d_shp_width[tid] * river_d_shp_length[tid] * grad_h;
}


__device__ realtype RiverCroSectArea(int order, realtype depth, realtype coeff)
{
	realtype          cs_area = 0.0;

	depth = (depth > 0.0) ? depth : 0.0;

	switch (order)
	{
	case RECTANGLE:
		cs_area = depth * coeff;
		break;
	case TRIANGLE:
		cs_area = depth * depth / coeff;
		break;
	case QUADRATIC:
		cs_area = 4.0 * depth * sqrt(depth) / (3.0 * sqrt(coeff));
		break;
	case CUBIC:
		cs_area =
			3.0 * pow(depth, 4.0 / 3.0) / (2.0 * pow(coeff, 1.0 / 3.0));
		break;
	default:
	/*	PIHMprintf(VL_ERROR, "Error: River order %d is not defined.\n",
			order);
		PIHMexit(EXIT_FAILURE);  */
		break;
	}

	return cs_area;
}

__device__ realtype RiverPerim(int order, realtype depth, realtype coeff)
{
	realtype          perim = 0.0;

	depth = (depth > 0.0) ? depth : 0.0;

	switch (order)
	{
	case RECTANGLE:
		perim = 2.0 * depth + coeff;
		break;
	case TRIANGLE:
		perim = 2.0 * depth * sqrt(1.0 + coeff * coeff) / coeff;
		break;
	case QUADRATIC:
		perim = sqrt(depth * (1.0 + 4.0 * coeff * depth) / coeff) +
			log(2.0 * sqrt(coeff * depth) +
			sqrt(1.0 + 4.0 * coeff * depth)) / (2.0 * coeff);
		break;
	case CUBIC:
		perim = 2.0 * ((pow(depth * (1.0 + 9.0 * pow(coeff, 2.0 / 3.0) *
			depth), 0.5) / 3.0) +
			(log(3.0 * pow(coeff, 1.0 / 3.0) * sqrt(depth) +
			pow(1.0 + 9.0 * pow(coeff, 2.0 / 3.0) * depth, 0.5)) /
			(9.0 * pow(coeff, 1.0 / 3.0))));
		break;
	default:
	/*	PIHMprintf(VL_ERROR, "Error: River order %d is not defined.\n",
			order);
		PIHMexit(EXIT_FAILURE);  
	*/
		break;
	}

	return perim;
}
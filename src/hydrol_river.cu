#include "pihm.h"
#include "pihm.cuh"

/* 河流计算的核函数 */
__global__ void hydrol_river(double *y, void *pihm_data,
	                         elem_struct_d *elem, river_struct_d *river,
							 ctrl_struct_d *ctrl_d)
{
	pihm_struct_d pihm_d = static_cast<pihm_struct_d>(pihm_data);

	/* Loop over river segments */
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int offset1 = tid + 3 * pihm_d->nelem;
	int offset2 = offset1 + pihm_d->nriver;

	river->ws.stage[tid] = (y[offset1] >= 0.0) ? y[offset1] : 0.0;
	river->ws.gw[tid] = (y[offset2] >= 0.0) ? y[offset2] : 0.0;


	river->wf.rivflow[UP_CHANL2CHANL][tid] = 0.0;
	river->wf.rivflow[UP_AQUIF2AQUIF][tid] = 0.0;

	RiverFlow(elem, river, ctrl_d->riv_mode, tid);		
}

__device__ void RiverFlow(elem_struct_d *elem, river_struct_d *river, int riv_mode, int tid)
{
	// 第1个  for (i = 0; i < nriver; i++)
	int left, right;
	double          effk_nabr;
	double          effk;
	int             down;      // 下段河道的分段编号    for SoA

	if (river->down[tid] > 0)
	{
		/*
		* Boundary conditions
		*
		* When a downstream segment is present, boundary conditions are
		* always applied to the upstream node
		*/
		if (river->attrib.riverbc_type[tid] != 0)
		{
			river->wf.rivflow[UP_CHANL2CHANL][tid] +=
				BoundFluxRiver(river->attrib.riverbc_type[tid], &(river->ws),
				&(river->topo), &(river->shp), &(river->matl),
				&(river->bc),tid);
		}

		 down = river->down[tid] - 1;       // for SoA
		/*
		* Channel flow between river-river segments
		*/
		river->wf.rivflow[DOWN_CHANL2CHANL][tid] = ChanFlowRiverToRiver(river, tid, down, riv_mode);
			
		/*
		* Subsurface flow between river-river segments
		*/
		left = river->leftele[tid] - 1;   //左边的单元编号
		right = river->rightele[tid] - 1; //右边的单元编号
		effk = 0.5 *
			(EffKh(&elem->soil, left, elem->ws.gw[left]) +
			EffKh(&elem->soil, right, elem->ws.gw[right]));

		//left = &elem[down->leftele - 1];
		//right = &elem[down->rightele - 1];
		left = river->leftele[down] - 1;   //左边的单元编号
		right = river->rightele[down] - 1; //右边的单元编号
		effk_nabr = 0.5 *
			(EffKh(&elem->soil, left, elem->ws.gw[left]) +
			 EffKh(&elem->soil, right, elem->ws.gw[right]));

		river->wf.rivflow[DOWN_AQUIF2AQUIF][tid] =
			SubFlowRiverToRiver(river, effk, effk_nabr, tid, down);
	}
	else
	{
		/*
		* Outlet flux
		*/
		river->wf.rivflow[DOWN_CHANL2CHANL][tid] =
			OutletFlux(tid, river->down[tid], &river->ws, &river->topo,
			               &river->shp, &river->matl, &river->bc);
		/* Note: boundary condition for subsurface element can be changed.
		* Assumption: no flow condition */
		river->wf.rivflow[DOWN_AQUIF2AQUIF][tid] = 0.0;
	}

	/*
	* Flux between river segments and triangular elements
	*/
	left = river->leftele[tid] - 1;    // 第tid段河段左侧单元的编号
	right = river->rightele[tid] - 1;  // 第tid段河段右侧单元的编号

	RiverToElem(elem, river, tid, left, right);  // tid -河段编号; left, right - 单元编号

	/*
	* Flux between river channel and subsurface
	*/
	river->wf.rivflow[CHANL_LKG][tid] =
		ChanLeak(&river->ws, &river->topo, &river->shp, &river->matl, tid);
		
	// 第2个  for (i = 0; i < nriver; i++)
	// river_struct   *down;
	if (river->down[tid] > 0)
	{
		// down = &river[river[i].down - 1];
		down = river->down[tid] - 1;

		river->wf.rivflow[UP_CHANL2CHANL][down] -=
			river->wf.rivflow[DOWN_CHANL2CHANL][tid];

		river->wf.rivflow[UP_AQUIF2AQUIF][down] -=
			river->wf.rivflow[DOWN_AQUIF2AQUIF][tid];
	}

}


__device__ double BoundFluxRiver(int riverbc_type, const river_wstate_struct_d *ws,
	const river_topo_struct_d *topo, const shp_struct_d *shp,
	const  matl_struct_d *matl, const river_bc_struct_d *bc, int tid)
{
	double          total_h;
	double          total_h_down;
	double          distance;
	double          grad_h;
	double          avg_h;
	double          avg_perim;
	double          crossa;
	double          flux = 0.0;

	if (riverbc_type > 0)
	{
		/* Dirichlet boundary condition */
		total_h = ws->stage[tid] + topo->zbed[tid];
		total_h_down = bc->head[tid];
		distance = 0.5 * shp->length[tid];
		grad_h = (total_h - total_h_down) / distance;
		avg_h = AvgH(grad_h, ws->stage[tid],
			((bc->head[tid] - topo->zbed[tid] > 0.0) ? (bc->head[tid] - topo->zbed[tid]) : 0.0));
		avg_perim = RiverPerim(shp->intrpl_ord[tid], ws->stage[tid], shp->coeff[tid]);
		crossa = RiverCroSectArea(shp->intrpl_ord[tid], ws->stage[tid], shp->coeff[tid]);
		avg_h = (avg_perim == 0.0) ? 0.0 : (crossa / avg_perim);
		flux = OverLandFlow(avg_h, grad_h, grad_h, crossa, matl->rough[tid]);			
	}
	else if (riverbc_type < 0)
	{
		/* Neumann boundary condition */
		flux = -bc->flux[tid];
	}

	return flux;
}


__device__ double ChanFlowRiverToRiver(const river_struct_d *river, int tid, int down, int riv_mode)
{
	double          total_h;
	double          perim;
	double          total_h_down;
	double          perim_down;
	double          avg_perim;
	double          avg_rough;
	double          distance;
	double          diff_h;
	double          grad_h;
	double          avg_sf;
	double          crossa;
	double          crossa_down;
	double          avg_crossa;
	double          avg_h;

	total_h = river->ws.stage[tid] + river->topo.zbed[tid];
	perim =
		RiverPerim(river->shp.intrpl_ord[tid], river->ws.stage[tid], river->shp.coeff[tid]);

	total_h_down = river->ws.stage[down] + river->topo.zbed[down];
	perim_down =
		RiverPerim(river->shp.intrpl_ord[down], river->ws.stage[down], river->shp.coeff[down]);

	avg_perim = (perim + perim_down) / 2.0;
	avg_rough = (river->matl.rough[tid] + river->matl.rough[down]) / 2.0;
	distance = 0.5 * (river->shp.length[tid] + river->shp.length[down]);

	diff_h = (riv_mode == KINEMATIC) ?
		(river->topo.zbed[tid] - river->topo.zbed[tid]) : (total_h - total_h_down);
	grad_h = diff_h / distance;
	avg_sf = (grad_h > 0.0) ? grad_h : RIVGRADMIN;
	crossa =
		RiverCroSectArea(river->shp.intrpl_ord[tid], river->ws.stage[tid],
		river->shp.coeff[tid]);
	crossa_down =
		RiverCroSectArea(river->shp.intrpl_ord[down], river->ws.stage[down], river->shp.coeff[down]);
	avg_crossa = 0.5 * (crossa + crossa_down);
	avg_h = (avg_perim == 0.0) ? 0.0 : (avg_crossa / avg_perim);

	return OverLandFlow(avg_h, grad_h, avg_sf, crossa, avg_rough);
}

__device__ double SubFlowRiverToRiver(const river_struct_d *river, double effk, double effk_nabr,
	                                  int tid, int down)
{
	double          total_h;
	double          total_h_down;
	double          avg_wid;
	double          diff_h;
	double          avg_h;
	double          distance;
	double          grad_h;
//	double          aquifer_depth;
	double          avg_ksat;

	/* Lateral flux calculation between element beneath river (ebr) * and ebr */
	total_h = river->ws.gw[tid] + river->topo.zmin[tid];
	total_h_down = river->ws.gw[down] + river->topo.zmin[down];
	avg_wid = (river->shp.width[tid] + river->shp.width[down]) / 2.0;
	diff_h = total_h - total_h_down;
	avg_h = AvgH(diff_h, river->ws.gw[tid], river->ws.gw[down]);
	distance = 0.5 * (river->shp.length[tid] + river->shp.length[down]);
	grad_h = diff_h / distance;
//	aquifer_depth = river->topo.zbed[tid] - river->topo.zmin[tid];  // not used!
#if defined(_ARITH_)
	avg_ksat = 0.5 * (effk + effk_nabr);
#else
	avg_ksat = 2.0 / (1.0 / effk + 1.0 / effk_nabr);
#endif
	/* Groundwater flow modeled by Darcy's law */
	return avg_ksat * grad_h * avg_h * avg_wid;
}


__device__ double OutletFlux(int tid, int down, const river_wstate_struct_d *ws,
	const river_topo_struct_d *topo, const shp_struct_d *shp,
	const matl_struct_d *matl, const river_bc_struct_d *bc)
{
	double          total_h;
	double          total_h_down;
	double          distance;
	double          grad_h;
	double          avg_h;
	double          avg_perim;
	double          crossa;
	double          discharge = 0.0;

	switch (down)
	{
	case DIRICHLET:
		/* Dirichlet boundary condition */
		total_h = ws->stage[tid] + topo->zbed[tid];
		total_h_down = bc->head[tid];
		distance = 0.5 * shp->length[tid];
		grad_h = (total_h - total_h_down) / distance;
		avg_h = AvgH(grad_h, ws->stage[tid],
			((bc->head[tid] - (topo->node_zmax[tid] - shp->depth[tid]) > 0.0) ?
			bc->head[tid] - (topo->node_zmax[tid] - shp->depth[tid]) : 0.0));
		avg_perim = RiverPerim(shp->intrpl_ord[tid], ws->stage[tid], shp->coeff[tid]);
		crossa = RiverCroSectArea(shp->intrpl_ord[tid], ws->stage[tid], shp->coeff[tid]);
		avg_h = (avg_perim == 0.0) ? 0.0 : (crossa / avg_perim);
		discharge =
			OverLandFlow(avg_h, grad_h, grad_h, crossa, matl->rough[tid]);
		break;
	case NEUMANN:
		/* Neumann boundary condition */
		discharge = -bc->flux[tid];
		break;
	case ZERO_DPTH_GRAD:
		/* Zero-depth-gradient boundary conditions */
		distance = 0.5 * shp->length[tid];
		grad_h = (topo->zbed[tid] - (topo->node_zmax[tid] - shp->depth[tid])) / distance;
		avg_h = ws->stage[tid];
		avg_perim = RiverPerim(shp->intrpl_ord[tid], ws->stage[tid], shp->coeff[tid]);
		crossa = RiverCroSectArea(shp->intrpl_ord[tid], ws->stage[tid], shp->coeff[tid]);
		discharge = sqrt(grad_h) * crossa * ((avg_perim > 0.0) ?
			pow(crossa / avg_perim, 2.0 / 3.0) : 0.0) / matl->rough[tid];
		break;
	case CRIT_DPTH:
		/* Critical depth boundary conditions */
		crossa = RiverCroSectArea(shp->intrpl_ord[tid], ws->stage[tid], shp->coeff[tid]);
		discharge = crossa * sqrt(GRAV * ws->stage[tid]);
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


__device__ void RiverToElem(elem_struct_d *elem, river_struct_d *river, int tid, int left, int right)
{
	double          effk_left, effk_right;
	int             j;

	/* Lateral surface flux calculation between river-triangular element */
	if (river->leftele[tid] > 0)
	{
		river->wf.rivflow[LEFT_SURF2CHANL][tid] = OvlFlowElemToRiver(elem, river, tid, left);
	}
	if (river->rightele[tid] > 0)
	{
		river->wf.rivflow[RIGHT_SURF2CHANL][tid] = OvlFlowElemToRiver(elem, river, tid, right);
	}

	effk_left = EffKh(&elem->soil, left, elem->ws.gw[left]);
	effk_right = EffKh(&elem->soil, right, elem->ws.gw[right]);

	/* Lateral subsurface flux calculation between river-triangular element */
	if (river->leftele[tid] > 0)
	{
		river->wf.rivflow[LEFT_AQUIF2CHANL][tid] =
			ChanFlowElemToRiver(elem, left, effk_left, river, tid, river->topo.dist_left[tid]);
	}
	if (river->rightele[tid] > 0)
	{
		river->wf.rivflow[RIGHT_AQUIF2CHANL][tid] =
			ChanFlowElemToRiver(elem, right, effk_right, river, tid,river->topo.dist_right[tid]);			
	}

	/* Lateral flux between rectangular element (beneath river) and triangular
	* element */
	if (river->leftele[tid] > 0)
	{
		river->wf.rivflow[LEFT_AQUIF2AQUIF][tid] =
			SubFlowElemToRiver(elem, left, effk_left, river, tid,
			0.5 * (effk_left + effk_right), river->topo.dist_left[tid]);
	}
	if (river->rightele > 0)
	{
		river->wf.rivflow[RIGHT_AQUIF2AQUIF][tid] =
			SubFlowElemToRiver(elem, right, effk_right, river, tid,
			0.5 * (effk_left + effk_right), river->topo.dist_right[tid]);
	}

	/* Replace flux term */
	/* Left */
	for (j = 0; j < NUM_EDGE; j++)
	{
		if (elem->nabr[j][left] == -river->ind[tid])
		{
			elem->wf.ovlflow[j][left] = -river->wf.rivflow[LEFT_SURF2CHANL][tid];
			elem->wf.subsurf[j][left] = -(river->wf.rivflow[LEFT_AQUIF2CHANL][tid] +
				river->wf.rivflow[LEFT_AQUIF2AQUIF][tid]);
			break;
		}
	}

	/* Right */
	for (j = 0; j < NUM_EDGE; j++)
	{
		if (elem->nabr[j][right] == -river->ind[tid])
		{
			elem->wf.ovlflow[j][right] = -river->wf.rivflow[RIGHT_SURF2CHANL][tid];
			elem->wf.subsurf[j][right] = -(river->wf.rivflow[RIGHT_AQUIF2CHANL][tid] +
				river->wf.rivflow[RIGHT_AQUIF2AQUIF][tid]);
			break;
		}
	}
}


__device__ double OvlFlowElemToRiver(const elem_struct_d *elem, const river_struct_d *river, int tid, int ind)
// 这里, tid是河段编号, ind是单元编号, 计算河道与相邻单元的水体交换
{
	double          zbank;
	double          flux;
	double          elem_h;
	double          rivseg_h;

	zbank = (river->topo.zmax[tid] > elem->topo.zmax[ind]) ?
		     river->topo.zmax[tid] : elem->topo.zmax[ind];

	elem_h = elem->topo.zmax[ind] + elem->ws.surfh[ind];
	rivseg_h = river->topo.zbed[tid] + river->ws.stage[tid];

	/*
	* Panday and Hyakorn 2004 AWR Eqs. (23) and (24)
	*/
	if (rivseg_h > elem_h)
	{
		if (elem_h > zbank)
		{
			/* Submerged weir */
			flux = river->matl.cwr[tid] * 2.0 * sqrt(2.0 * GRAV) *
				river->shp.length[tid] * sqrt(rivseg_h - elem_h) *
				(rivseg_h - zbank) / 3.0;
		}
		else
		{
			if (zbank < rivseg_h)
			{
				/* Free-flowing weir */
				flux = river->matl.cwr[tid] * 2.0 * sqrt(2.0 * GRAV) *
					river->shp.length[tid] * sqrt(rivseg_h - zbank) *
					(rivseg_h - zbank) / 3.0;
			}
			else
			{
				flux = 0.0;
			}
		}
	}
	else if (elem->ws.surfh[ind] > DEPRSTG)
	{
		if (rivseg_h > zbank)
		{
			/* Submerged weir */
			flux = -river->matl.cwr[tid] * 2.0 * sqrt(2.0 * GRAV) *
				river->shp.length[tid] * sqrt(elem_h - rivseg_h) *
				(elem_h - zbank) / 3.0;
		}
		else
		{
			if (zbank < elem_h)
			{
				/* Free-flowing weir */
				flux = -river->matl.cwr[tid] * 2.0 * sqrt(2.0 * GRAV) *
					river->shp.length[tid] * sqrt(elem_h - zbank) *
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


__device__ double ChanFlowElemToRiver(const elem_struct_d *elem, int ind, double effk,
	const river_struct_d *river, int tid, double distance)
// 这里, tid是河段编号, ind是单元编号, 计算河道与相邻单元的水体交换
{
	double          diff_h;
	double          avg_h;
	double          grad_h;
	double          avg_ksat;

	diff_h = (river->ws.stage[tid] + river->topo.zbed[tid]) -
		(elem->ws.gw[ind] + elem->topo.zmin[ind]);

	/* This is head in neighboring cell representation */
	if (elem->topo.zmin[ind] > river->topo.zbed[tid])
	{
		avg_h = elem->ws.gw[ind];
	}
	else if (elem->topo.zmin[ind] + elem->ws.gw[ind] > river->topo.zbed[tid])
	{
		avg_h = elem->topo.zmin[ind] + elem->ws.gw[ind] - river->topo.zbed[tid];
	}
	else
	{
		avg_h = 0.0;
	}
	avg_h = AvgH(diff_h, river->ws.stage[tid], avg_h);

	grad_h = diff_h / distance;

	avg_ksat = 0.5 * (effk + river->matl.ksath[tid]);

	return  river->shp.length[tid] * avg_ksat * grad_h * avg_h;
}

__device__ double SubFlowElemToRiver(const elem_struct_d *elem, int ind, double effk,
	const river_struct_d *river, int tid, double effk_riv, double distance)
// 这里, tid是河段编号, ind是单元编号, 计算河道与相邻单元的水体交换
{
	double          diff_h;
	double          avg_h;
//	double          aquifer_depth;    // not used!
	double          avg_ksat;
	double          grad_h;

	diff_h = (river->ws.gw[tid] + river->topo.zmin[tid]) -
		(elem->ws.gw[ind] + elem->topo.zmin[ind]);

	/* This is head in neighboring cell represention */
	if (elem->topo.zmin[ind] > river->topo.zbed[tid])
	{
		avg_h = 0.0;
	}
	else if (elem->topo.zmin[ind] + elem->ws.gw[ind] > river->topo.zbed[tid])
	{
		avg_h = river->topo.zbed[tid] - elem->topo.zmin[ind];
	}
	else
	{
		avg_h = elem->ws.gw[ind];
	}
	avg_h = AvgH(diff_h, river->ws.gw[tid], avg_h);
//	aquifer_depth = river->topo.zbed[tid] - river->topo.zmin[tid];

#if defined(_ARITH_)
	avg_ksat = 0.5 * (effk + effk_riv);
#else
	avg_ksat = 2.0 / (1.0 / effk + 1.0 / effk_riv);
#endif
	grad_h = diff_h / distance;

	return river->shp.length[tid] * avg_ksat * grad_h * avg_h;
}


__device__ double ChanLeak(const river_wstate_struct_d *ws, const river_topo_struct_d *topo,
	const shp_struct_d *shp, const matl_struct_d *matl, int tid)
{
	double          diff_h;
	double          grad_h;

	if (topo->zbed[tid] - (ws->gw[tid] + topo->zmin[tid]) > 0.0)
	{
		diff_h = ws->stage[tid];
	}
	else
	{
		diff_h = ws->stage[tid] + topo->zbed[tid] - (ws->gw[tid] + topo->zmin[tid]);
	}

	grad_h = diff_h / matl->bedthick[tid];

	return matl->ksatv[tid] * shp->width[tid] * shp->length[tid] * grad_h;
}


__device__ double RiverCroSectArea(int order, double depth, double coeff)
{
	double          cs_area = 0.0;

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

__device__ double RiverPerim(int order, double depth, double coeff)
{
	double          perim = 0.0;

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
		PIHMexit(EXIT_FAILURE);  */
		break;
	}

	return perim;
}
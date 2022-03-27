#include "pihm.h"
#include "pihm.cuh"

// 计算流域的摩阻坡度, 将原始的LateralFlow中的FrictSlope()移到这里.
__global__ void FrictSlope(int surf_mode, realtype *dhbydx, realtype *dhbydy,
#include "FrictSlope_device_declarations.c"
)
{
	int tid = blockDim.x * blockIdx.x + threadIdx.x;

	int             nabr, rivnabr;
	realtype        surfh[NUM_EDGE];
	realtype        topo_x[NUM_EDGE], topo_y[NUM_EDGE];

	if (surf_mode == DIFF_WAVE)
	{
		//   j = 0 (单元的第1条边)
		if (elem_d_nabr0[tid] > 0)        // elements side not linked with river segment
		{
			nabr = elem_d_nabr0[tid] - 1;

			surfh[0] = elem_d_topo_zmax[nabr] + elem_d_ws_surfh[nabr];
		}
		else if (elem_d_nabr0[tid] < 0)   // interconnect with river
		{
			rivnabr = -elem_d_nabr0[tid] - 1;   // 与第tid个单元的第j条边衔接的河段编号

			if (river_d_ws_stage[rivnabr] > river_d_shp_depth[rivnabr])
			{
				surfh[0] = river_d_topo_zbed[rivnabr] + river_d_ws_stage[rivnabr];
			}
			else
			{
				surfh[0] = river_d_topo_zmax[rivnabr];
			}
		}
		else     // = 0 rhe element side is on the basin boundary
		{
			if (elem_d_attrib_bc_type0[tid] == NO_FLOW)
			{
				surfh[0] = elem_d_topo_zmax[tid] + elem_d_ws_surfh[tid];
			}
			else
			{
				surfh[0] = elem_d_bc_head0[tid];
			}
		}

		//   j = 1 (单元的第2条边)
		if (elem_d_nabr1[tid] > 0)        // elements side not linked with river segment
		{
			nabr = elem_d_nabr1[tid] - 1;
			surfh[1] = elem_d_topo_zmax[nabr] + elem_d_ws_surfh[nabr];
		}
		else if (elem_d_nabr1[tid] < 0)   // interconnect with river
		{
			rivnabr = -elem_d_nabr1[tid] - 1;   // 与第tid个单元的第j条边衔接的河段编号

			if (river_d_ws_stage[rivnabr] > river_d_shp_depth[rivnabr])
			{
				surfh[1] = river_d_topo_zbed[rivnabr] + river_d_ws_stage[rivnabr];
			}
			else
			{
				surfh[1] = river_d_topo_zmax[rivnabr];
			}
		}
		else     // = 0 rhe element side is on the basin boundary
		{
			if (elem_d_attrib_bc_type1[tid] == NO_FLOW)
			{
				surfh[1] = elem_d_topo_zmax[tid] + elem_d_ws_surfh[tid];
			}
			else
			{
				surfh[1] = elem_d_bc_head1[tid];
			}
		}

		//   j = 2 (单元的第3条边)
		if (elem_d_nabr2[tid] > 0)        // elements side not linked with river segment
		{
			nabr = elem_d_nabr2[tid] - 1;
			surfh[2] = elem_d_topo_zmax[nabr] + elem_d_ws_surfh[nabr];
		}
		else if (elem_d_nabr2[tid] < 0)   // interconnect with river
		{
			rivnabr = -elem_d_nabr2[tid] - 1;   // 与第tid个单元的第j条边衔接的河段编号

			if (river_d_ws_stage[rivnabr] > river_d_shp_depth[rivnabr])
			{
				surfh[2] = river_d_topo_zbed[rivnabr] + river_d_ws_stage[rivnabr];
			}
			else
			{
				surfh[2] = river_d_topo_zmax[rivnabr];
			}
		}
		else     // = 0 rhe element side is on the basin boundary
		{
			if (elem_d_attrib_bc_type2[tid] == NO_FLOW)
			{
				surfh[2] = elem_d_topo_zmax[tid] + elem_d_ws_surfh[tid];
			}
			else
			{
				surfh[2] = elem_d_bc_head2[tid];
			}
		}

		topo_x[0] = elem_d_topo_nabr_x0[tid];
		topo_x[1] = elem_d_topo_nabr_x1[tid];
		topo_x[2] = elem_d_topo_nabr_x2[tid];
		topo_y[0] = elem_d_topo_nabr_y0[tid];
		topo_y[1] = elem_d_topo_nabr_y1[tid];
		topo_y[2] = elem_d_topo_nabr_y2[tid];

		dhbydx[tid] = DhByDl(topo_y, topo_x, surfh);
		dhbydy[tid] = DhByDl(topo_x, topo_y, surfh);

	}  // if (surf_mode == DIFF_WAVE)

	__syncthreads();
}


__device__  realtype DhByDl(const realtype *l1,
	const realtype *l2,
	const realtype *surfh)
{
	return -1.0 *
		(l1[2] * (surfh[1] - surfh[0]) + l1[1] * (surfh[0] - surfh[2]) +
		l1[0] * (surfh[2] - surfh[1])) /
		(l2[2] * (l1[1] - l1[0]) + l2[1] * (l1[0] - l1[2]) +
		l2[0] * (l1[2] - l1[1]));
}
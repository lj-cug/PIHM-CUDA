#include "pihm.h"
#include "pihm.cuh"

void NoahHydrol(pihm_struct_d pihm_d, realtype dt)
{
	int nelem = pihm_d->nelem;

#include "NoahHydrol_device_injections.c" 

	int numThreads_elem = min(32, nelem);
	int numBlocks_elem = (nelem + numThreads_elem - 1) / numThreads_elem;

	NoahHydrol_kernel << < numThreads_elem, numBlocks_elem >> >(nelem, dt,
#include "NoahHydrol_kernel_arguments.c" 
		);

	cudaDeviceSynchronize();

#if 0
	// check for error
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess)
	{
		// print the CUDA error message and exit
		printf("CUDA error in Noah: %s\n", cudaGetErrorString(error));
		exit(-1);
	}
#endif

}


__global__ void NoahHydrol_kernel(int nelem, realtype dt,
#include "NoahHydrol_kernel_declarations.c" 
)
{
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int             k;

	if (tid < nelem){

		/* Find water table position 
		cuda版本中, elem_d_ps_sldpth elem_d_ps_satdpth  是一维数组
		*/
		/*
		elem_d_ps_nwtbl[tid] = FindWaterTable(elem_d_ps_sldpth[tid], elem_d_ps_nsoil[tid],
			                   elem_d_ws_gw[tid], elem_d_ps_satdpth[tid]);  
		*/

		// 现在根据最大土壤分层数(11)，全部计算所有分层的elem_d_ps_sldpth elem_d_ps_satdpth
		// CPU版本的代码，调用FindWaterTable，根据不同单元上的土壤分层计算2个变量
		int             layer = -999;
		double          depth = 0.0;
		double          dsum = 0.0;

		// Initialization 
		elem_d_ps_satdpth0[tid] = 0.0;
		elem_d_ps_satdpth1[tid] = 0.0;
		elem_d_ps_satdpth2[tid] = 0.0;
		elem_d_ps_satdpth3[tid] = 0.0;
		elem_d_ps_satdpth4[tid] = 0.0;
		elem_d_ps_satdpth5[tid] = 0.0;
		elem_d_ps_satdpth6[tid] = 0.0;
		elem_d_ps_satdpth7[tid] = 0.0;
		elem_d_ps_satdpth8[tid] = 0.0;
		elem_d_ps_satdpth9[tid] = 0.0;
		elem_d_ps_satdpth10[tid] = 0.0;

		depth += elem_d_ps_sldpth0[tid] + elem_d_ps_sldpth1[tid] + elem_d_ps_sldpth2[tid] +
			elem_d_ps_sldpth3[tid] + elem_d_ps_sldpth4[tid] + elem_d_ps_sldpth5[tid] +
			elem_d_ps_sldpth6[tid] + elem_d_ps_sldpth7[tid] + elem_d_ps_sldpth8[tid] +
			elem_d_ps_sldpth9[tid] + elem_d_ps_sldpth10[tid];

		    if (elem_d_ws_gw[tid] <= 0.0)
			{
				layer = elem_d_ps_nsoil[tid];

				elem_d_ps_sldpth0[tid] = 1.0E-3;
				elem_d_ps_sldpth1[tid] = 1.0E-3;
				elem_d_ps_sldpth2[tid] = 1.0E-3;
				elem_d_ps_sldpth3[tid] = 1.0E-3;
				elem_d_ps_sldpth4[tid] = 1.0E-3;
				elem_d_ps_sldpth5[tid] = 1.0E-3;
				elem_d_ps_sldpth6[tid] = 1.0E-3;
				elem_d_ps_sldpth7[tid] = 1.0E-3;
				elem_d_ps_sldpth8[tid] = 1.0E-3;
				elem_d_ps_sldpth9[tid] = 1.0E-3;
				elem_d_ps_sldpth10[tid] = 1.0E-3;
			}
			else if (elem_d_ws_gw[tid] > depth)
			{
				layer = 0;

				elem_d_ps_satdpth0[tid] = elem_d_ps_sldpth0[tid];
				elem_d_ps_satdpth1[tid] = elem_d_ps_sldpth1[tid];
				elem_d_ps_satdpth2[tid] = elem_d_ps_sldpth2[tid];
				elem_d_ps_satdpth3[tid] = elem_d_ps_sldpth3[tid];
				elem_d_ps_satdpth4[tid] = elem_d_ps_sldpth4[tid];
				elem_d_ps_satdpth5[tid] = elem_d_ps_sldpth5[tid];
				elem_d_ps_satdpth6[tid] = elem_d_ps_sldpth6[tid];
				elem_d_ps_satdpth7[tid] = elem_d_ps_sldpth7[tid];
				elem_d_ps_satdpth8[tid] = elem_d_ps_sldpth8[tid];
				elem_d_ps_satdpth9[tid] = elem_d_ps_sldpth9[tid];
				elem_d_ps_satdpth10[tid] = elem_d_ps_sldpth10[tid];
			}
			else
			{
// j = 10
				if (dsum + elem_d_ps_sldpth10[tid] > elem_d_ws_gw[tid])
				{
					elem_d_ps_satdpth10[tid] = elem_d_ws_gw[tid] - dsum;
					layer = 10 + 1;
				}
				else
				{
					elem_d_ps_satdpth10[tid] = elem_d_ps_sldpth10[tid];
					dsum += elem_d_ps_sldpth10[tid];
				}
// j = 9
				if (dsum + elem_d_ps_sldpth9[tid] > elem_d_ws_gw[tid])
				{
					elem_d_ps_satdpth9[tid] = elem_d_ws_gw[tid] - dsum;
					layer = 9 + 1;
				}
				else
				{
					elem_d_ps_satdpth9[tid] = elem_d_ps_sldpth9[tid];
					dsum += elem_d_ps_sldpth9[tid];
				}
// j = 8
				if (dsum + elem_d_ps_sldpth8[tid] > elem_d_ws_gw[tid])
				{
					elem_d_ps_satdpth8[tid] = elem_d_ws_gw[tid] - dsum;
					layer = 8 + 1;
				}
				else
				{
					elem_d_ps_satdpth8[tid] = elem_d_ps_sldpth8[tid];
					dsum += elem_d_ps_sldpth8[tid];
				}
// j = 7
				if (dsum + elem_d_ps_sldpth7[tid] > elem_d_ws_gw[tid])
				{
					elem_d_ps_satdpth7[tid] = elem_d_ws_gw[tid] - dsum;
					layer = 7 + 1;
				}
				else
				{
					elem_d_ps_satdpth7[tid] = elem_d_ps_sldpth7[tid];
					dsum += elem_d_ps_sldpth7[tid];
				}
// j = 6
				if (dsum + elem_d_ps_sldpth6[tid] > elem_d_ws_gw[tid])
				{
					elem_d_ps_satdpth6[tid] = elem_d_ws_gw[tid] - dsum;
					layer = 6 + 1;
				}
				else
				{
					elem_d_ps_satdpth6[tid] = elem_d_ps_sldpth6[tid];
					dsum += elem_d_ps_sldpth6[tid];
				}
// j = 5
				if (dsum + elem_d_ps_sldpth5[tid] > elem_d_ws_gw[tid])
				{
					elem_d_ps_satdpth5[tid] = elem_d_ws_gw[tid] - dsum;
					layer = 5 + 1;
				}
				else
				{
					elem_d_ps_satdpth5[tid] = elem_d_ps_sldpth5[tid];
					dsum += elem_d_ps_sldpth5[tid];
				}
// j = 4
				if (dsum + elem_d_ps_sldpth4[tid] > elem_d_ws_gw[tid])
				{
					elem_d_ps_satdpth4[tid] = elem_d_ws_gw[tid] - dsum;
					layer = 4 + 1;
				}
				else
				{
					elem_d_ps_satdpth4[tid] = elem_d_ps_sldpth4[tid];
					dsum += elem_d_ps_sldpth4[tid];
				}
// j =3
				if (dsum + elem_d_ps_sldpth3[tid] > elem_d_ws_gw[tid])
				{
					elem_d_ps_satdpth3[tid] = elem_d_ws_gw[tid] - dsum;
					layer = 3 + 1;
				}
				else
				{
					elem_d_ps_satdpth3[tid] = elem_d_ps_sldpth3[tid];
					dsum += elem_d_ps_sldpth3[tid];
				}
// j =2
				if (dsum + elem_d_ps_sldpth2[tid] > elem_d_ws_gw[tid])
				{
					elem_d_ps_satdpth2[tid] = elem_d_ws_gw[tid] - dsum;
					layer = 2 + 1;
				}
				else
				{
					elem_d_ps_satdpth2[tid] = elem_d_ps_sldpth2[tid];
					dsum += elem_d_ps_sldpth2[tid];
				}
// j =1
				if (dsum + elem_d_ps_sldpth1[tid] > elem_d_ws_gw[tid])
				{
					elem_d_ps_satdpth1[tid] = elem_d_ws_gw[tid] - dsum;
					layer = 1 + 1;
				}
				else
				{
					elem_d_ps_satdpth1[tid] = elem_d_ps_sldpth1[tid];
					dsum += elem_d_ps_sldpth1[tid];
				}
// j = 0
				if (dsum + elem_d_ps_sldpth0[tid] > elem_d_ws_gw[tid])
				{
					elem_d_ps_satdpth0[tid] = elem_d_ws_gw[tid] - dsum;
					layer = 0 + 1;
				}
				else
				{
					elem_d_ps_satdpth0[tid] = elem_d_ps_sldpth0[tid];
					dsum += elem_d_ps_sldpth0[tid];
				}

			}

			elem_d_ps_nwtbl[tid] = layer;

		/*
		for (k = 0; k < elem_d_ps_nsoil[tid]; k++)
		{
			elem_d_ws_smc[k][tid] =
				(elem_d_ws_smc[k][tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
				elem_d_ws_smc[k][tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
			elem_d_ws_smc[k][tid] =
				(elem_d_ws_smc[k][tid] < elem_d_soil_smcmax[tid]) ?
				elem_d_ws_smc[k][tid] : elem_d_soil_smcmax[tid];
			elem_d_ws_sh2o[k][tid] =
				(elem_d_ws_sh2o[k][tid] < elem_d_ws_smc[k][tid]) ?
				elem_d_ws_sh2o[k][tid] : elem_d_ws_smc[k][tid];
		}
		*/
		// k==0
		elem_d_ws_smc0[tid] =
			(elem_d_ws_smc0[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc0[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc0[tid] =
			(elem_d_ws_smc0[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc0[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o0[tid] =
			(elem_d_ws_sh2o0[tid] < elem_d_ws_smc0[tid]) ?
			elem_d_ws_sh2o0[tid] : elem_d_ws_smc0[tid];
		// k==1
		elem_d_ws_smc1[tid] =
			(elem_d_ws_smc1[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc1[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc1[tid] =
			(elem_d_ws_smc1[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc1[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o1[tid] =
			(elem_d_ws_sh2o1[tid] < elem_d_ws_smc1[tid]) ?
			elem_d_ws_sh2o1[tid] : elem_d_ws_smc1[tid];
		// k==2
		elem_d_ws_smc2[tid] =
			(elem_d_ws_smc2[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc2[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc2[tid] =
			(elem_d_ws_smc2[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc2[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o2[tid] =
			(elem_d_ws_sh2o2[tid] < elem_d_ws_smc2[tid]) ?
			elem_d_ws_sh2o2[tid] : elem_d_ws_smc2[tid];
		// k==3
		elem_d_ws_smc3[tid] =
			(elem_d_ws_smc3[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc3[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc3[tid] =
			(elem_d_ws_smc3[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc3[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o3[tid] =
			(elem_d_ws_sh2o3[tid] < elem_d_ws_smc3[tid]) ?
			elem_d_ws_sh2o3[tid] : elem_d_ws_smc3[tid];
		// k==4
		elem_d_ws_smc4[tid] =
			(elem_d_ws_smc4[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc4[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc4[tid] =
			(elem_d_ws_smc4[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc4[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o4[tid] =
			(elem_d_ws_sh2o4[tid] < elem_d_ws_smc4[tid]) ?
			elem_d_ws_sh2o4[tid] : elem_d_ws_smc4[tid];
		// k==5
		elem_d_ws_smc5[tid] =
			(elem_d_ws_smc5[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc5[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc5[tid] =
			(elem_d_ws_smc5[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc5[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o5[tid] =
			(elem_d_ws_sh2o5[tid] < elem_d_ws_smc5[tid]) ?
			elem_d_ws_sh2o5[tid] : elem_d_ws_smc5[tid];
		// k==6
		elem_d_ws_smc6[tid] =
			(elem_d_ws_smc6[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc6[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc6[tid] =
			(elem_d_ws_smc6[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc6[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o6[tid] =
			(elem_d_ws_sh2o6[tid] < elem_d_ws_smc6[tid]) ?
			elem_d_ws_sh2o6[tid] : elem_d_ws_smc6[tid];
		// k==7
		elem_d_ws_smc7[tid] =
			(elem_d_ws_smc7[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc7[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc7[tid] =
			(elem_d_ws_smc7[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc7[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o7[tid] =
			(elem_d_ws_sh2o7[tid] < elem_d_ws_smc7[tid]) ?
			elem_d_ws_sh2o7[tid] : elem_d_ws_smc7[tid];
		// k==8
		elem_d_ws_smc8[tid] =
			(elem_d_ws_smc8[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc8[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc8[tid] =
			(elem_d_ws_smc8[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc8[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o8[tid] =
			(elem_d_ws_sh2o8[tid] < elem_d_ws_smc8[tid]) ?
			elem_d_ws_sh2o8[tid] : elem_d_ws_smc8[tid];
		// k==9
		elem_d_ws_smc9[tid] =
			(elem_d_ws_smc9[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc9[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc9[tid] =
			(elem_d_ws_smc9[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc9[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o9[tid] =
			(elem_d_ws_sh2o9[tid] < elem_d_ws_smc9[tid]) ?
			elem_d_ws_sh2o9[tid] : elem_d_ws_smc9[tid];
		// k==10
		elem_d_ws_smc10[tid] =
			(elem_d_ws_smc10[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc10[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc10[tid] =
			(elem_d_ws_smc10[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc10[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o10[tid] =
			(elem_d_ws_sh2o10[tid] < elem_d_ws_smc10[tid]) ?
			elem_d_ws_sh2o10[tid] : elem_d_ws_smc10[tid];


#if defined(_CYCLES_)
		SmFlx();  // Solute transport
#else
		SmFlx(tid, dt,
#include "SmFlx_kernel_arguments.c"		
			 );
#endif

#if defined(_CYCLES_)
		/*
		 * Calcluate vertical transport of solute
		 */

#endif

	}  // if (tid < nelem)

}

/*   已经不用了
__device__ int FindWaterTable(double *sldpth, int nsoil, double gw, double *satdpth)
{
	int             layer = -999;
	int             j;
	double          dsum = 0.0;
	double          depth;

	for (j = 0; j < MAXLYR; j++)
	{
		satdpth[j] = 0.0;
	}

	depth = 0.0;
	for (j = 0; j < nsoil; j++)
	{
		depth += sldpth[j];
	}

	if (gw <= 0.0)
	{
		layer = nsoil;
		satdpth[nsoil - 1] = 1.0E-3;
	}
	else if (gw > depth)
	{
		layer = 0;
		for (j = 0; j < nsoil; j++)
		{
			satdpth[j] = sldpth[j];
		}
	}
	else
	{
		for (j = nsoil - 1; j >= 0; j--)
		{
			if (dsum + sldpth[j] > gw)
			{
				satdpth[j] = gw - dsum;
				layer = j + 1;
				break;
			}
			else
			{
				satdpth[j] = sldpth[j];
				dsum += sldpth[j];
			}
		}
	}

	return layer;
}
*/





#if defined(_CYCLES_)
void SmFlx(const soil_struct *soil, const cstate_struct *cs, double dt,
	pstate_struct *ps, wstate_struct *ws, wflux_struct *wf)
#else
__device__ void SmFlx(int tid, double dt,
#include "SmFlx_kernel_declarations.c" 
)
#endif
{
	/*
	* Function SmFlx
	*
	* Calculate soil moisture flux. The soil moisture content (smc - a per unit
	* volume measurement) is a dependent variable that is updated with
	* prognostic eqns.
	* Frozen ground version: new states added: sh2o, and frozen ground
	* correction factor, frzfact and parameter slope.
	*/
	int             i;
	double          ai[MAXLYR], bi[MAXLYR], ci[MAXLYR];
	double          rhstt[MAXLYR];
	double          sice[MAXLYR];

	/* Store ice content at each soil layer before calling SRT and SStep */
	/*
	for (i = 0; i < elem_d_ps_nsoil[tid]; i++)
	{
		sice[i] = elem_d_ws_smc[i][tid] - elem_d_ws_sh2o[i][tid];
	}
	*/
	sice[0] = elem_d_ws_smc0[tid] - elem_d_ws_sh2o0[tid];
	sice[1] = elem_d_ws_smc1[tid] - elem_d_ws_sh2o1[tid];
	sice[2] = elem_d_ws_smc2[tid] - elem_d_ws_sh2o2[tid];
	sice[3] = elem_d_ws_smc3[tid] - elem_d_ws_sh2o3[tid];
	sice[4] = elem_d_ws_smc4[tid] - elem_d_ws_sh2o4[tid];
	sice[5] = elem_d_ws_smc5[tid] - elem_d_ws_sh2o5[tid];
	sice[6] = elem_d_ws_smc6[tid] - elem_d_ws_sh2o6[tid];
	sice[7] = elem_d_ws_smc7[tid] - elem_d_ws_sh2o7[tid];
	sice[8] = elem_d_ws_smc8[tid] - elem_d_ws_sh2o8[tid];
	sice[9] = elem_d_ws_smc9[tid] - elem_d_ws_sh2o9[tid];
	sice[10] = elem_d_ws_smc10[tid] - elem_d_ws_sh2o10[tid];

	/* Call subroutines SRT and SStep to solve the soil moisture tendency
	* equations.
	* Call the SRT/SStep function in the manner of time scheme "d" (implicit
	* state, explicit coefficient) of Section 2 of Kalnay and Kanamitsu pcpdrp
	* is units of m/s, zsoil is negative depth in m
	* According to Dr. Ken Mitchell's suggestion, add the second contraint to
	* remove numerical instability of runoff and soil moisture
	*
	* Frozen ground version:
	* smc states replaced by sh2o states in SRT subr. sh2o & sice states
	* included in SStep subr. Frozen ground correction factor, frzfact added.
	* All water balance calculations using unfrozen water */
	if (0 == elem_d_ps_nwtbl[tid])
	{
		/* Special case: all soil layers are saturated */
		/*
		for (i = 0; i < elem_d_ps_nsoil[tid]; i++)
		{
			elem_d_ws_smc[i][tid] = elem_d_soil_smcmax[tid];
			elem_d_ws_sh2o[i][tid] = elem_d_ws_smc[i][tid] - sice[i];
		}
		*/
		elem_d_ws_smc0[tid] = elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o0[tid] = elem_d_ws_smc0[tid] - sice[0];
		elem_d_ws_smc1[tid] = elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o1[tid] = elem_d_ws_smc1[tid] - sice[1];
		elem_d_ws_smc2[tid] = elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o2[tid] = elem_d_ws_smc2[tid] - sice[2];
		elem_d_ws_smc3[tid] = elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o3[tid] = elem_d_ws_smc3[tid] - sice[3];
		elem_d_ws_smc4[tid] = elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o4[tid] = elem_d_ws_smc4[tid] - sice[4];
		elem_d_ws_smc5[tid] = elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o5[tid] = elem_d_ws_smc5[tid] - sice[5];
		elem_d_ws_smc6[tid] = elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o6[tid] = elem_d_ws_smc6[tid] - sice[6];
		elem_d_ws_smc7[tid] = elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o7[tid] = elem_d_ws_smc7[tid] - sice[7];
		elem_d_ws_smc8[tid] = elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o8[tid] = elem_d_ws_smc8[tid] - sice[8];
		elem_d_ws_smc9[tid] = elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o9[tid] = elem_d_ws_smc9[tid] - sice[9];
		elem_d_ws_smc10[tid] = elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o10[tid] = elem_d_ws_smc10[tid] - sice[10];
	}
	else
	{
#if defined(_CYCLES_)
		SRT(soil, cs, dt, ps, ws, wf, rhstt, sice, ai, bi, ci);
#else
		SRT(tid, rhstt, sice, ai, bi, ci,
#include "SRT_kernel_arguments.c"				
		   );
#endif
		SStep(tid, rhstt, sice, ai, bi, ci, dt,
#include "SStep_kernel_arguments.c"	
             );
	}
}



#if defined(_CYCLES_)
void SRT(const soil_struct *soil, const cstate_struct *cs, double dt,
	pstate_struct *ps, wstate_struct *ws, wflux_struct *wf, double *rhstt,
	double *sice, double *ai, double *bi, double *ci)
#else
__device__ void SRT(int tid, double *rhstt, double *sice, double *ai,
	                double *bi, double *ci,
#include "SRT_kernel_declarations.c"
)
#endif
{
	/*
	* Function SRT
	* Calculate the right hand side of the time tendency term of the soil water
	* diffusion equation. Also to compute (prepare) the matrix coefficients for
	* the tri-diagonal matrix of the implicit time scheme.
	*/
	int             iohinf;
	int             j, jj, k, ks;
	double          ddz;
	double          ddz2;
	double          denom;
	double          denom2;
	double          numer;
	double          pddum;
	double          mxsmc, mxsmc2;
	double          sicemax;
	double          wcnd;
	double          wcnd2;
	double          wdf;
	double          wdf2;
	double          dsmdz, dsmdz2;
	double          dice;
	const double    CVFRZ = 3.0;
	double          acrt;
	double          sum;
	double          ialp1;
	/* Frozen ground version:
	* Reference frozen ground parameter, cvfrz, is a shape parameter of areal
	* distribution function of soil ice content which equals 1/cv.
	* Cv is a coefficient of spatial variation of soil ice content. Based on
	* field data cv depends on areal mean of frozen depth, and it close to
	* constant = 0.6 if areal mean frozen depth is above 20 cm. That is why
	* parameter cvfrz = 3 (int{1/0.6*0.6}). Current logic doesn't allow cvfrz
	* be bigger than 3 */

	/* Let sicemax be the greatest, if any, frozen water content within soil
	* layers. */
	iohinf = 1;
	sicemax = 0.0;

	for (ks = 0; ks < elem_d_ps_nsoil[tid]; ks++)
	{
		sicemax = (sice[ks] > sicemax) ? sice[ks] : sicemax;
	}

	/* Calculate infiltration reduction factor due to frozen soil */
	dice = -elem_d_ps_zsoil[0][tid] * sice[0];
	for (ks = 1; ks < elem_d_ps_nsoil[tid]; ks++)
	{
		dice += (elem_d_ps_zsoil[ks - 1][tid] - elem_d_ps_zsoil[ks][tid]) * sice[ks];
	}

	elem_d_ps_fcr[tid] = 1.0;

	if (dice > 1.0e-2)
	{
		acrt = CVFRZ * elem_d_ps_frzx[tid] / dice;
		sum = 1.0;
		ialp1 = (int)round(CVFRZ) - 1;
		for (j = 1; j < ialp1 + 1; j++)
		{
			k = 1;
			for (jj = j + 1; jj < ialp1; jj++)
			{
				k *= jj;
			}
			sum += pow(acrt, CVFRZ - (double)j) / (double)k;
		}

		elem_d_ps_fcr[tid] = 1.0 - exp(-acrt) * sum;
	}

	/* Determine rainfall infiltration rate and runoff */
#if defined(_CYCLES_)
	ResidueWetting(ps, cs, dt, ws, wf);
#endif

	pddum = elem_d_wf_infil[tid];

	mxsmc = elem_d_ws_sh2o[0][tid];

	dsmdz = (elem_d_ws_sh2o[0][tid] - elem_d_ws_sh2o[1][tid]) / (-0.5 * elem_d_ps_zsoil[1][tid]);
	WDfCnd(tid, &wdf, &wcnd, mxsmc, sicemax, 
#include "WDfCnd_kernel_arguments.c"		
		  );

	/* Calc the matrix coefficients ai, bi, and ci for the top layer */
	ddz = 1.0 / (-0.5 * elem_d_ps_zsoil[1][tid]);
	ai[0] = 0.0;
	bi[0] = wdf * ddz / (-elem_d_ps_zsoil[0][tid]);
	ci[0] = -bi[0];
	/* Calc rhstt for the top layer after calc'ng the vertical soil moisture
	* gradient btwn the top and next to top layers. */
	rhstt[0] = (wdf * dsmdz + wcnd - pddum + elem_d_wf_edir[tid] + elem_d_wf_et[0][tid]) /
		        elem_d_ps_zsoil[0][tid];

	rhstt[0] += elem_d_wf_runoff2_lyr[0][tid] / elem_d_ps_zsoil[0][tid];

	/* Loop thru the remaining soil layers, repeating the abv process */
	/* Initialize ddz2 */
	ddz2 = 0.0;
	for (k = 1; k < elem_d_ps_nsoil[tid]; k++)
	{
		denom2 = (elem_d_ps_zsoil[k - 1][tid] - elem_d_ps_zsoil[k][tid]);
		if (k < elem_d_ps_nsoil[tid] - 1)
		{
			mxsmc2 = elem_d_ws_sh2o[k][tid];
			denom = elem_d_ps_zsoil[k - 1][tid] - elem_d_ps_zsoil[k + 1][tid];
			dsmdz2 = (elem_d_ws_sh2o[k][tid] - elem_d_ws_sh2o[k + 1][tid]) / (denom * 0.5);
			WDfCnd(tid, &wdf2, &wcnd2, mxsmc2, sicemax,
#include "WDfCnd_kernel_arguments.c"
				   );

			/* Calc some partial products for later use in calc'ng rhstt
			* Calc the matrix coef, ci, after calc'ng its partial product */
			ddz2 = 2.0 / denom;
			ci[k] = -wdf2 * ddz2 / denom2;
		}
		else
		{
			/* Retrieve the soil water diffusivity and hydraulic conductivity
			* for this layer */
			wdf2 = 0.0;
			wcnd2 = 0.0;

			/* Calc a partial product for later use in calc'ng rhstt */
			dsmdz2 = 0.0;
			/* Set matrix coef ci to zero */
			ci[k] = 0.0;
		}

		/* Calc rhstt for this layer after calc'ng its numerator */
		numer = (wdf2 * dsmdz2) + wcnd2 - (wdf * dsmdz) - wcnd + elem_d_wf_et[k][tid];

		numer += elem_d_wf_runoff2_lyr[k][tid];

		rhstt[k] = numer / (-denom2);

		/* Calc matrix coefs, ai, and bi for this layer */
		ai[k] = -wdf * ddz / denom2;
		bi[k] = -(ai[k] + ci[k]);

		/* Reset values of wdf, wcnd, dsmdz, and ddz for loop to next lyr */
		if (k != elem_d_ps_nsoil[tid] - 1)
		{
			wdf = wdf2;
			wcnd = wcnd2;
			dsmdz = dsmdz2;
			ddz = ddz2;
		}
	}
}

__device__ void WDfCnd(int tid, double *wdf, double *wcnd, double smc, double sicemax,
#include "WDfCnd_kernel_declarations.c"
)
{
	/*
	* Function WDfCnd
	*
	* Calculate soil water diffusivity and soil hydraulic conductivity.
	* Flux-PIHM: using van Genuchten parameters
	*/
	double          expon;
	double          factr1;
	double          factr2;
	double          vkwgt;
	double          satkfunc;
	double          dpsidsm;

	/* Calc the ratio of the actual to the max psbl soil h2o content */
	factr1 = 0.05 / (elem_d_soil_smcmax[tid] - elem_d_soil_smcmin[tid]);
	factr2 = (smc - elem_d_soil_smcmin[tid]) / (elem_d_soil_smcmax[tid] - elem_d_soil_smcmin[tid]);

	/* Factr2 should avoid to be 0 or 1 */
	factr2 = (factr2 > 1.0 - 5.0e-4) ? 1.0 - 5.0e-4 : factr2;
	factr2 = (factr2 < 0.0 + 5.0e-4) ? 5.0e-4 : factr2;

	factr1 = (factr1 < factr2) ? factr1 : factr2;
	expon = 1.0 - 1.0 / elem_d_soil_beta[tid];

	satkfunc = KrFunc_cuda(elem_d_soil_beta[tid], factr2);
	dpsidsm =
		(1.0 - expon) / elem_d_soil_alpha[tid] / expon / (elem_d_soil_smcmax[tid] - elem_d_soil_smcmin[tid]) *
		pow(pow(factr2, -1.0 / expon) - 1.0, -expon) *
		pow(factr2, -(1.0 / expon + 1.0));

	*wcnd = elem_d_soil_ksatv[tid] * satkfunc;

	*wdf = *wcnd * dpsidsm;

	if (sicemax > 0.0)
	{
		vkwgt = 1.0 / (1.0 + pow(500.0 * sicemax, 3.0));
		satkfunc = KrFunc_cuda(elem_d_soil_beta[tid], factr1);
		dpsidsm = (1.0 - expon) / elem_d_soil_alpha[tid] / expon /
			(elem_d_soil_smcmax[tid] - elem_d_soil_smcmin[tid]) *
			pow(pow(factr1, -1.0 / expon) - 1.0, -expon) *
			pow(factr1, -(1.0 / expon + 1.0));
		*wdf = vkwgt * *wdf + (1.0 - vkwgt) * dpsidsm * satkfunc * elem_d_soil_ksatv[tid];
	}
}


__device__ void SStep(int tid, double *rhstt, double *sice, double *ai,
	double *bi, double *ci, double dt,
#include "SStep_kernel_declarations.c"
)
{
	/*
	* Function SStep
	*
	* Calculate/update soil moisture content values and canopy moisture content
	* values.
	*/
	int             k;
	double          rhsttin[MAXLYR];
	double          ciin[MAXLYR];
	double          sh2o0[MAXLYR];

	/* Create 'amount' values of variables to be input to the tri-diagonal
	* matrix routine. */
	for (k = 0; k < elem_d_ps_nsoil[tid]; k++)
	{
		rhstt[k] *= dt;
		ai[k] *= dt;
		bi[k] = 1.0 + bi[k] * dt;
		ci[k] *= dt;
		sh2o0[k] = elem_d_ws_sh2o[k][tid];
		elem_d_wf_smflxv[k][tid] = 0.0;
	}

	/* Copy values for input variables before call to Rosr12 */
	for (k = 0; k < elem_d_ps_nsoil[tid]; k++)
	{
		rhsttin[k] = rhstt[k];
		ciin[k] = ci[k];
	}

	/* Call Rosr12 to solve the tri-diagonal matrix */
	Rosr12(ci, ai, bi, ciin, rhsttin, rhstt, elem_d_ps_nsoil[tid]);  /* 显示错误，编译没有问题 */

	/* Sum the previous smc value and the matrix solution to get a new value. */
	for (k = 0; k < elem_d_ps_nsoil[tid]; k++)
	{
		elem_d_ws_sh2o[k][tid] += ci[k];
	}

	AdjSmProf(tid, sice, dt, 
#include "AdjSmProf_kernel_arguments.c"
		     );

	/* Calculate soil moisture flux within soil layers */
	for (k = elem_d_ps_nsoil[tid] - 1; k > 0; k--)
	{
		/* Positive smflxv[k] is flux out of soil layer k */
		elem_d_wf_smflxv[k - 1][tid] =
			(elem_d_ws_sh2o[k][tid] - sh2o0[k]) * elem_d_ps_sldpth[k][tid] / dt +
			 elem_d_wf_runoff2_lyr[k][tid] + elem_d_wf_et[k][tid] + elem_d_wf_smflxv[k][tid];
	}
}


__device__ void AdjSmProf(int tid, const realtype *sice, realtype dt,
#include "AdjSmProf_kernel_declarations.c"
)
{
	/* Min allowable value of smc will be SH2OMIN. */
	/* In Flux-PIHM, the soil layers are gone thru twice:
	* 1. From bottom to top, to make sure all layers below water table is
	* saturated;
	* 2. From top to bottom, to make sure soil moisture from all layers are
	* within plausible ranges */
	int             k;
	realtype          ddz;
	realtype          sh2omid[MAXLYR];
	realtype          wplus;
	realtype          stot;
	realtype          stotmin;

	/* Runoff3: runoff within soil layers */
	wplus = 0.0;
	elem_d_wf_runoff3[tid] = 0.0;

	for (k = elem_d_ps_nsoil[tid] - 1; k >= 0; k--)
	{
		if (k != 0)
		{
			ddz = elem_d_ps_zsoil[k - 1][tid] - elem_d_ps_zsoil[k][tid];
		}
		else
		{
			ddz = -elem_d_ps_zsoil[0][tid];
		}

		sh2omid[k] = elem_d_ws_sh2o[k][tid] + wplus / ddz;
		stot = sh2omid[k] + sice[k];

		if (stot > elem_d_soil_smcmax[tid])
		{
			elem_d_ws_smc[k][tid] = elem_d_soil_smcmax[tid];
			wplus = (stot - elem_d_soil_smcmax[tid]) * ddz;
		}
		else
		{
			stotmin = (elem_d_ps_satdpth[k][tid] * elem_d_soil_smcmax[tid] +
				(elem_d_ps_sldpth[k][tid] - elem_d_ps_satdpth[k][tid]) *
				(elem_d_soil_smcmin[tid] + SH2OMIN)) /
				elem_d_ps_sldpth[k][tid];
			stotmin = (stotmin > elem_d_soil_smcmax[tid]) ? elem_d_soil_smcmax[tid] : stotmin;
			stotmin = (stotmin < elem_d_soil_smcmin[tid] + SH2OMIN) ?
				(elem_d_soil_smcmin[tid] + SH2OMIN) : stotmin;

			if (stot < stotmin)
			{
				elem_d_ws_smc[k][tid] = stotmin;
				wplus = (stot - stotmin) * ddz;
			}
			else
			{
				elem_d_ws_smc[k][tid] = stot;
				wplus = 0.0;
			}
		}

		sh2omid[k] = elem_d_ws_smc[k][tid] - sice[k];
	}

	ddz = -elem_d_ps_zsoil[0][tid];
	for (k = 0; k < elem_d_ps_nsoil[tid]; k++)
	{
		if (k != 0)
		{
			ddz = elem_d_ps_zsoil[k - 1][tid] - elem_d_ps_zsoil[k][tid];
		}
		elem_d_ws_sh2o[k][tid] = sh2omid[k] + wplus / ddz;
		stot = elem_d_ws_sh2o[k][tid] + sice[k];
		if (stot > elem_d_soil_smcmax[tid])
		{
			wplus = (stot - elem_d_soil_smcmax[tid]) * ddz;
		}
		else
		{
			wplus = 0.0;
		}

		elem_d_ws_smc[k][tid] = (stot < elem_d_soil_smcmax[tid]) ? stot : elem_d_soil_smcmax[tid];
		elem_d_ws_smc[k][tid] = (elem_d_ws_smc[k][tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc[k][tid] : (elem_d_soil_smcmin[tid] + SH2OMIN);
		elem_d_ws_sh2o[k][tid] = elem_d_ws_smc[k][tid] - sice[k];
		elem_d_ws_sh2o[k][tid] = (elem_d_ws_sh2o[k][tid] > 0.0) ? elem_d_ws_sh2o[k][tid] : 0.0;
	}

	elem_d_wf_runoff3[tid] = wplus / dt;
}


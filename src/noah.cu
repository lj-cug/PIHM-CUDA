#include "pihm.h"
#include "pihm.cuh"

/*
   Noah LSM模型中，由于核函数中循环计算展开的需要，土壤垂向分层数 NSOIL，选取固定的5层   2021.05.12

   MM-PIHM-0.11.1的Noah模型考虑了冰川和冰盖的模块  glacier  ice; MM-PIHM-0.10没有冰川模块

*/
void Noah(pihm_struct_d pihm_d, realtype dt)
{
	int nelem = pihm_d->nelem;

#include "Noah_device_injections.c"   // Noah模型的形参变量

	int numThreads_elem = min(32, nelem);
	int numBlocks_elem = (nelem + numThreads_elem - 1) / numThreads_elem;

/* 
   把原本的3个子函数fused成一个核函数 Noah_kernel(),
   包括：CalHum_kernel + SFlx_kernel + some state variables' calculation
*/
	Noah_kernel << < numThreads_elem, numBlocks_elem >> >(nelem, dt,
#include "Noah_kernel_arguments.c" 
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

// Some constans used in CalHum
__constant__  realtype    A2 = 17.67;
__constant__  realtype    A3 = 273.15;
__constant__  realtype    A4 = 29.65;
__constant__  realtype    ELWV = 2.501e6;
__constant__  realtype    E0 = 611.0;
__constant__  realtype    RVV = 461.0;
__constant__  realtype    EPSILON = 0.622;
__constant__  realtype    SVP1 = 611.2;
__constant__  realtype    SVP2 = 17.67;
__constant__  realtype    SVP3 = 29.65;
__constant__  realtype    SVPT0 = 273.15;
__constant__  int         IZ0TLND = 0;     // originally in SFlx_kernel

__global__ void Noah_kernel(int nelem, realtype dt,
#include "Noah_kernel_declarations.c" 
)
{
/*
   __global__ void CalHum_kernel( )
*/
	int tid = blockDim.x * blockIdx.x + threadIdx.x;

	int               j;
	realtype          e;
	realtype          esat;
	realtype          svp;
	realtype          a23m4;
	realtype          t2v;
	realtype          rho;
	realtype          rh;

	if (tid < nelem){

// Fused CalHum_kernel here.   2021.05.13

		rh = elem_d_ps_rh[tid] / 100.0;

		svp = SVP1 * exp(SVP2 * (elem_d_es_sfctmp[tid] - SVPT0) /
			                    (elem_d_es_sfctmp[tid] - SVP3));
		e = rh * svp;

		elem_d_ps_q2[tid] = (0.622 * e) / (elem_d_ps_sfcprs[tid] - (1.0 - 0.622) * e);

		elem_d_es_th2[tid] = elem_d_es_sfctmp[tid] + (0.0098 * elem_d_ps_zlvl[tid]);
		//elem_d_es_t1v[tid] = elem_d_es_t1[tid] * (1.0 + 0.61 * elem_d_ps_q2[tid]);
		//elem_d_es_th2v[tid] = elem_d_es_th2[tid] * (1.0 + 0.61 * elem_d_ps_q2[tid]);

		t2v = elem_d_es_sfctmp[tid] * (1.0 + 0.61 * elem_d_ps_q2[tid]);
		rho = elem_d_ps_sfcprs[tid] / (RD * t2v);

		a23m4 = A2 * (A3 - A4);

		esat = E0 * exp(ELWV / RVV * (1.0 / A3 - 1.0 / elem_d_es_sfctmp[tid]));

		elem_d_ps_q2sat[tid] = EPSILON * esat / (elem_d_ps_sfcprs[tid] - (1.0 - EPSILON) * esat);

		elem_d_ps_dqsdt2[tid] = elem_d_ps_q2sat[tid] * a23m4 /
			                    ((elem_d_es_sfctmp[tid] - A4) * (elem_d_es_sfctmp[tid] - A4));
	
// Fused FrozRain here.
		elem_d_ps_ffrozp[tid] =
			FrozRain(elem_d_wf_prcp[tid], elem_d_es_sfctmp[tid]);
	
// Fused some state variablesin Noah here
		elem_d_ps_alb[tid] = BADVAL;

		if (elem_d_ps_q1[tid] == BADVAL)
		{
			elem_d_ps_q1[tid] = elem_d_ps_q2[tid];
		}

		elem_d_ef_solnet[tid] = elem_d_ef_soldn[tid] * (1.0 - elem_d_ps_albedo[tid]);
		elem_d_ef_lwdn[tid] = elem_d_ef_longwave[tid] * elem_d_ps_emissi[tid];

 //     #pragma unroll   // 展开的次数未知，不能展开
		/*
		for (j = 0; j < elem_d_ps_nsoil[tid]; j++)  
		{
			elem_d_ws_smc[j][tid] =
				(elem_d_ws_smc[j][tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
				elem_d_ws_smc[j][tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
			elem_d_ws_smc[j][tid] = (elem_d_ws_smc[j][tid] < elem_d_soil_smcmax[tid]) ?
				elem_d_ws_smc[j][tid] : elem_d_soil_smcmax[tid];
			elem_d_ws_sh2o[j][tid] = (elem_d_ws_sh2o[j][tid] < elem_d_ws_smc[j][tid]) ?
				elem_d_ws_sh2o[j][tid] : elem_d_ws_smc[j][tid];
		} 
		*/

		// nsoil 是变化的. 但CUDA版本下，统一展开，用array_cuda[MAXLYR]
// j==0
		elem_d_ws_smc0[tid] =
			(elem_d_ws_smc0[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc0[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc0[tid] = (elem_d_ws_smc0[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc0[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o0[tid] = (elem_d_ws_sh2o0[tid] < elem_d_ws_smc0[tid]) ?
			elem_d_ws_sh2o0[tid] : elem_d_ws_smc0[tid];
// j==1
		elem_d_ws_smc1[tid] =
			(elem_d_ws_smc1[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc1[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc1[tid] = (elem_d_ws_smc1[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc1[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o1[tid] = (elem_d_ws_sh2o1[tid] < elem_d_ws_smc1[tid]) ?
			elem_d_ws_sh2o1[tid] : elem_d_ws_smc1[tid];
// j==2
		elem_d_ws_smc2[tid] =
			(elem_d_ws_smc2[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc2[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc2[tid] = (elem_d_ws_smc2[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc2[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o2[tid] = (elem_d_ws_sh2o2[tid] < elem_d_ws_smc2[tid]) ?
			elem_d_ws_sh2o2[tid] : elem_d_ws_smc2[tid];
// j==3
		elem_d_ws_smc3[tid] =
			(elem_d_ws_smc3[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc3[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc3[tid] = (elem_d_ws_smc3[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc3[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o3[tid] = (elem_d_ws_sh2o3[tid] < elem_d_ws_smc3[tid]) ?
			elem_d_ws_sh2o3[tid] : elem_d_ws_smc3[tid];
// j==4
		elem_d_ws_smc4[tid] =
			(elem_d_ws_smc4[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc4[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc4[tid] = (elem_d_ws_smc4[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc4[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o4[tid] = (elem_d_ws_sh2o4[tid] < elem_d_ws_smc4[tid]) ?
			elem_d_ws_sh2o4[tid] : elem_d_ws_smc4[tid];
// j==5
		elem_d_ws_smc5[tid] =
			(elem_d_ws_smc5[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc5[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc5[tid] = (elem_d_ws_smc5[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc5[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o5[tid] = (elem_d_ws_sh2o5[tid] < elem_d_ws_smc5[tid]) ?
			elem_d_ws_sh2o5[tid] : elem_d_ws_smc5[tid];
// j==6
		elem_d_ws_smc6[tid] =
			(elem_d_ws_smc6[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc6[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc6[tid] = (elem_d_ws_smc6[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc6[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o6[tid] = (elem_d_ws_sh2o6[tid] < elem_d_ws_smc6[tid]) ?
			elem_d_ws_sh2o6[tid] : elem_d_ws_smc6[tid];
// j==7
		elem_d_ws_smc7[tid] =
			(elem_d_ws_smc7[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc7[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc7[tid] = (elem_d_ws_smc7[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc7[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o7[tid] = (elem_d_ws_sh2o7[tid] < elem_d_ws_smc7[tid]) ?
			elem_d_ws_sh2o7[tid] : elem_d_ws_smc7[tid];
// j==8
		elem_d_ws_smc8[tid] =
			(elem_d_ws_smc8[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc8[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc8[tid] = (elem_d_ws_smc8[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc8[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o8[tid] = (elem_d_ws_sh2o8[tid] < elem_d_ws_smc8[tid]) ?
			elem_d_ws_sh2o8[tid] : elem_d_ws_smc8[tid];
// j==9
		elem_d_ws_smc9[tid] =
			(elem_d_ws_smc9[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc9[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc9[tid] = (elem_d_ws_smc9[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc9[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o9[tid] = (elem_d_ws_sh2o9[tid] < elem_d_ws_smc9[tid]) ?
			elem_d_ws_sh2o9[tid] : elem_d_ws_smc9[tid];
// j==10
		elem_d_ws_smc10[tid] =
			(elem_d_ws_smc10[tid] > elem_d_soil_smcmin[tid] + SH2OMIN) ?
			elem_d_ws_smc10[tid] : elem_d_soil_smcmin[tid] + SH2OMIN;
		elem_d_ws_smc10[tid] = (elem_d_ws_smc10[tid] < elem_d_soil_smcmax[tid]) ?
			elem_d_ws_smc10[tid] : elem_d_soil_smcmax[tid];
		elem_d_ws_sh2o10[tid] = (elem_d_ws_sh2o10[tid] < elem_d_ws_smc10[tid]) ?
			elem_d_ws_sh2o10[tid] : elem_d_ws_smc10[tid];

// Fused SFlx_kernel here.   2021.05.13

		/*
          Run Noah LSM
		  __global__ void SFlx_kernel( )
		*/
			/*
			* Subroutine SFlx - unified noahlsm version 1.0 july 2007
			*
			* Sub-driver for "Noah LSM" family of physics subroutines for a
			* soil/veg/snowpack land-surface model to update soil moisture, soil ice,
			* soil temperature, skin temperature, snowpack water content, snowdepth,
			* and all terms of the surface energy balance and surface water balance
			* (excluding input atmospheric forcings of downward radiation and precip)
			*/
			int             frzgra, snowng;
			realtype        df1;
			realtype        df1a;
			realtype        dsoil;
			realtype        dtot;
			realtype        frcsno, frcsoi;
			realtype        t1v;
			realtype        th2v;
			realtype        t2v;
			realtype        t24;
			realtype        interp_fraction;
			realtype        sn_new;
			realtype        prcpf;
			realtype        soilwm;
			realtype        soilww;
			realtype        smav[MAXLYR];
			int             k;

			/*
			* Initialization
			*/
			elem_d_wf_snomlt[tid] = 0.0;

			elem_d_wf_pcpdrp[tid] = 0.0;

			/* Urban */
			if (elem_d_lc_isurban[tid])
			{
				elem_d_lc_shdfac[tid] = 0.05;
#if !defined(_CYCLES_)
				elem_d_epc_rsmin[tid] = 400.0;
#endif
				elem_d_soil_smcmax[tid] = 0.45;
				elem_d_soil_smcmin[tid] = 0.0;
				elem_d_soil_smcref[tid] = 0.42;
				elem_d_soil_smcwlt[tid] = 0.40;
				elem_d_soil_smcdry[tid] = 0.40;
			}

#if defined(_CYCLES_)
			elem_d_lc_shdfac[tid] = CommRadIntcp(crop);  //
#endif

			/* Set minimum LAI for non-barren land cover to improve performance */
			if (elem_d_lc_shdfac[tid] > 0.0)
			{
				elem_d_ps_proj_lai[tid] = (elem_d_ps_proj_lai[tid] > 0.5) ? elem_d_ps_proj_lai[tid] : 0.5;
			}

			elem_d_ws_cmcmax[tid] = elem_d_lc_shdfac[tid] * elem_d_lc_cmcfactr[tid] * elem_d_ps_proj_lai[tid];

			/* Flux-PIHM uses LAI as a forcing variable.
			* Vegetation fraction is calculated from LAI following Noah-MP */
			if (elem_d_ps_proj_lai[tid] >= elem_d_lc_laimax[tid])
			{
				elem_d_ps_embrd[tid] = elem_d_lc_emissmax[tid];
				elem_d_ps_alb[tid] = elem_d_lc_albedomin[tid];
				elem_d_ps_z0brd[tid] = elem_d_lc_z0max[tid];
			}
			else if (elem_d_ps_proj_lai[tid] <= elem_d_lc_laimin[tid])
			{
				elem_d_ps_embrd[tid] = elem_d_lc_emissmin[tid];
				elem_d_ps_alb[tid] = elem_d_lc_albedomax[tid];
				elem_d_ps_z0brd[tid] = elem_d_lc_z0min[tid];
			}
			else
			{
				if (elem_d_lc_laimax[tid] > elem_d_lc_laimin[tid])
				{
					interp_fraction =
						(elem_d_ps_proj_lai[tid] - elem_d_lc_laimin[tid]) /
						(elem_d_lc_laimax[tid] - elem_d_lc_laimin[tid]);

					/* Bound interp_fraction between 0 and 1 */
					interp_fraction = (interp_fraction < 1.0) ? interp_fraction : 1.0;
					interp_fraction = (interp_fraction > 0.0) ? interp_fraction : 0.0;

					/* Scale emissivity and LAI between emissmin and emissmax by
					* interp_fraction */
					elem_d_ps_embrd[tid] = ((1.0 - interp_fraction) * elem_d_lc_emissmin[tid]) +
						(interp_fraction * elem_d_lc_emissmax[tid]);
					elem_d_ps_alb[tid] = ((1.0 - interp_fraction) * elem_d_lc_albedomax[tid]) +
						(interp_fraction * elem_d_lc_albedomin[tid]);
					elem_d_ps_z0brd[tid] = ((1.0 - interp_fraction) * elem_d_lc_z0min[tid]) +
						(interp_fraction * elem_d_lc_z0max[tid]);
				}
				else
				{
					elem_d_ps_embrd[tid] = 0.5 * elem_d_lc_emissmin[tid] + 0.5 * elem_d_lc_emissmax[tid];
					elem_d_ps_alb[tid] = 0.5 *  elem_d_lc_albedomin[tid] + 0.5 * elem_d_lc_albedomax[tid];
					elem_d_ps_z0brd[tid] = 0.5 * elem_d_lc_z0min[tid] + 0.5 * elem_d_lc_z0max[tid];
				}
			}

			/* Initialize precipitation logicals. */
			snowng = 0;
			frzgra = 0;

			/* If input snowpack is nonzero, then compute snow density "sndens" and
			* snow thermal conductivity "sncond" subroutine */
			if (elem_d_ws_sneqv[tid] <= 1.0e-7)    /* Safer if KMH (2008/03/25) */
			{
				elem_d_ws_sneqv[tid] = 0.0;
				elem_d_ps_sndens[tid] = 0.0;
				elem_d_ps_snowh[tid] = 0.0;
				elem_d_ps_sncond[tid] = 1.0;
			}
			else
			{
				elem_d_ps_sndens[tid] = elem_d_ws_sneqv[tid] / elem_d_ps_snowh[tid];

#if defined(_USE_CPU_)
				if (elem_d_ps_sndens[tid] > 1.0)
				{
					PIHMprintf(VL_ERROR,
						"Error: Physical snow depth is less than snow water equiv.\n");
					PIHMexit(EXIT_FAILURE);
				}
#endif

				elem_d_ps_sncond[tid] = CSnow(elem_d_ps_sndens[tid]);
			}

			/* Determine if it's precipitating and what kind of precip it is.
			* If it's prcping and the air temp is colder than 0 C, it's snowing!
			* If it's prcping and the air temp is warmer than 0 C, but the grnd temp is
			* colder than 0 C, freezing rain is presumed to be falling. */
			if (elem_d_wf_prcp[tid] > 0.0)
			{
				/* Snow defined when fraction of frozen precip (ffrozp) > 0.5, passed in
				* from model microphysics.  */
				if (elem_d_ps_ffrozp[tid] > 0.5)
				{
					snowng = 1;
				}
				else
				{
					if (elem_d_es_t1[tid] <= TFREEZ)
					{
						frzgra = 1;
					}
				}
			}

			/* If either prcp flag is set, determine new snowfall and add it to the
			* existing snowpack.
			* Note that since all precip is added to snowpack, no precip infiltrates
			* into the soil so that prcpf is set to zero. */
			if (snowng || frzgra)
			{
				sn_new = elem_d_wf_prcp[tid] * dt;
				elem_d_ws_sneqv[tid] += sn_new;
				prcpf = 0.0;

				/* Update snow density based on new snowfall, using old and new snow.
				* Update snow thermal conductivity */
				SnowNew(tid, sn_new,
					elem_d_ps_snowh, elem_d_es_sfctmp, elem_d_ps_sndens);
				elem_d_ps_sncond[tid] = CSnow(elem_d_ps_sndens[tid]);
			}
			else
			{
				/* Precip is liquid (rain), hence save in the precip variable (along
				* with any canopy "drip" added to this later) */
				prcpf = elem_d_wf_prcp[tid];
			}

			/*
			* Determine snowcover and albedo over land.
			*/
			if (elem_d_ws_sneqv[tid] == 0.0)
			{
				/* If snow depth = 0, set snow fraction = 0, albedo = snow free albedo.
				*/
				elem_d_ps_sncovr[tid] = 0.0;
				elem_d_ps_albedo[tid] = elem_d_ps_alb[tid];
				elem_d_ps_emissi[tid] = elem_d_ps_embrd[tid];
			}
			else
			{
				/* Determine snow fractional coverage.
				* Determine surface albedo modification due to snowdepth state. */
				elem_d_ps_sncovr[tid] = SnFrac(elem_d_ws_sneqv[tid], elem_d_lc_snup[tid], elem_d_ps_salp[tid]);

				elem_d_ps_sncovr[tid] = (elem_d_ps_sncovr[tid] < 0.98) ? elem_d_ps_sncovr[tid] : 0.98;

				AlCalc(tid, dt, snowng,
#include "AlCalc_kernel_arguments.c" 	
					);
			}

			/*
			* Next calculate the subsurface heat flux, which first requires
			* calculation of the thermal diffusivity. Treatment of the latter
			* follows that on Pages 148-149 from "Heat transfer in cold climates", by
			* V. J. Lunardini (published in 1981 by van Nostrand Reinhold Co.) i.e.
			* treatment of two contiguous "plane parallel" mediums (namely here the
			* first soil layer and the snowpack layer, if any). This diffusivity
			* treatment behaves well for both zero and nonzero snowpack, including the
			* limit of very thin snowpack.  this treatment also eliminates the need to
			* impose an arbitrary upper bound on subsurface heat flux when the snowpack
			* becomes extremely thin.
			*
			* First calculate thermal diffusivity of top soil layer, using both the
			* frozen and liquid soil moisture, following the soil thermal diffusivity
			* function of Peters-Lidard et al. (1998, JAS, Vol 55, 1209-1224), which
			* requires the specifying the quartz content of the given soil class (see
			* routine RedPrm)
			*
			* Next add subsurface heat flux reduction effect from the overlying green
			* canopy, adapted from Section 2.1.2 of Peters-Lidard et al. (1997, JGR,
			* Vol 102(D4))
			*/
			df1 = TDfCnd(elem_d_ws_smc0[tid], elem_d_soil_quartz[tid], elem_d_soil_smcmax[tid],
				elem_d_soil_smcmin[tid], elem_d_ws_sh2o0[tid]);

			/* Urban */
			if (elem_d_lc_isurban[tid])
			{
				df1 = 3.24;
			}

			df1 *= exp(elem_d_ps_sbeta[tid] * elem_d_lc_shdfac[tid]);

			/* KMH 09/03/2006
			* KMH 03/25/2008  Change sncovr threshold to 0.97 */
			if (elem_d_ps_sncovr[tid] > 0.97)
			{
				df1 = elem_d_ps_sncond[tid];
			}

			/* Finally "plane parallel" snowpack effect following V. J. Linardini
			* reference cited above. Note that dtot is combined depth of snowdepth and
			* thickness of first soil layer */
			dsoil = -(0.5 * elem_d_ps_zsoil0[tid]);
			if (elem_d_ws_sneqv[tid] == 0.0)
			{
				elem_d_ef_ssoil[tid] = df1 * (elem_d_es_t1[tid] - elem_d_es_stc0[tid]) / dsoil;
			}
			else
			{
				dtot = elem_d_ps_snowh[tid] + dsoil;
				frcsno = elem_d_ps_snowh[tid] / dtot;
				frcsoi = dsoil / dtot;

				/* 1. harmonic mean (series flow) */
				//df1h = (elem_d_ps_sncond * df1) / (frcsoi * elem_d_ps_sncond + frcsno * df1);

				/* 2. arithmetic mean (parallel flow) */
				df1a = frcsno * elem_d_ps_sncond[tid] + frcsoi * df1;

				/* 3. geometric mean (intermediate between harmonic and arithmetic
				* mean) */
				//df1 = pow (sncond, frcsno) * pow(df1, frcsoi);
				/* weigh df by snow fraction */
				//df1 = df1h * sncovr + df1a * (1.0-sncovr);
				//df1 = df1h * sncovr + df1 * (1.0-sncovr);

				/* Calculate subsurface heat flux, ssoil, from final thermal
				* diffusivity of surface mediums, df1 above, and skin temperature and
				* top mid-layer soil temperature */
				df1 = df1a * elem_d_ps_sncovr[tid] + df1 * (1.0 - elem_d_ps_sncovr[tid]);

				elem_d_ef_ssoil[tid] = df1 * (elem_d_es_t1[tid] - elem_d_es_stc0[tid]) / dtot;
			}

			/*
			* Determine surface roughness over snowpack using snow condition from the
			* previous timestep.
			*/
			if (elem_d_ps_sncovr[tid] > 0.0)
			{
				elem_d_ps_z0[tid] = Snowz0(elem_d_ps_sncovr[tid], elem_d_ps_z0brd[tid], elem_d_ps_snowh[tid]);
			}
			else
			{
				elem_d_ps_z0[tid] = elem_d_ps_z0brd[tid];
			}

			/*
			* Next call function SfcDif to calculate the sfc exchange coef (ch) for
			* heat and moisture.
			*
			* Note !!!
			* Do not call SfcDif until after above call to RedPrm, in case alternative
			* values of roughness length (z0) and Zilintinkevich coef (czil) are set
			* there via namelist i/o.
			*
			* Note !!!
			* Function SfcDif returns a ch that represents the wind spd times the
			* "original" nondimensional "ch" typical in literature. Hence the ch
			* returned from SfcDif has units of m/s. The important companion
			* coefficient of ch, carried here as "rch", is the ch from sfcdif times air
			* density and parameter "CP". "rch" is computed in "Penman".
			* rch rather than ch is the coeff usually invoked later in eqns.
			*
			* Note !!!
			* SfcDif also returns the surface exchange coefficient for momentum, cm,
			* also known as the surface drag coefficient. Needed as a state variable
			* for iterative/implicit solution of ch in SfcDif
			*/
			t1v = elem_d_es_t1[tid] * (1.0 + 0.61 * elem_d_ps_q2[tid]);
			th2v = elem_d_es_th2[tid] * (1.0 + 0.61 * elem_d_ps_q2[tid]);

			SfcDifOff(tid, t1v, th2v, IZ0TLND,
#include "SfcDifOff_kernel_arguments.c"		
				);

			/*
			* Call Penman function to calculate potential evaporation (ETP), and other
			* partial products and sums save in common/rite for later calculations.
			*/

			/* Calculate total downward radiation (solar plus longwave) needed in
			* Penman ep subroutine that follows */
			elem_d_ef_fdown[tid] = elem_d_ef_solnet[tid] + elem_d_ef_lwdn[tid];

			/* Calc virtual temps and virtual potential temps needed by Penman. */
			t2v = elem_d_es_sfctmp[tid] * (1.0 + 0.61 * elem_d_ps_q2[tid]);

			Penman(tid, &t24, t2v, snowng, frzgra,
#include "Penman_kernel_arguments.c"	
				  );

			/*
			* Call CanRes to calculate the canopy resistance and convert it into pc if
			* nonzero greenness fraction
			*/
			if (elem_d_lc_shdfac[tid] > 0.0)
			{
#if defined(_CYCLES_)
				CanRes(es, ps);   // under coding
#else
				CanRes(tid,
#include "CanRes_kernel_arguments.c"		
					);
#endif
			}
			else
			{
				elem_d_ps_rc[tid] = 0.0;
			}

			/*
			* Now decide major pathway branch to take depending on whether snowpack
			* exists or not
			*/
			elem_d_wf_esnow[tid] = 0.0;

			if (elem_d_ws_sneqv[tid] == 0.0)
			{
#if defined(_CYCLES_)
				NoPac(soil, lc, cs, dt, t24, crop, ps, ws, wf, es, ef);
#else
				NoPac(tid, dt, t24,
#include "NoPac_kernel_arguments.c"		
					 );
#endif
				elem_d_ps_eta_kinematic[tid] = elem_d_wf_eta[tid] * 1000.0;
			}
			else
			{
#if defined(_CYCLES_)
				SnoPac(soil, lc, cs, snowng, dt, t24, prcpf, df1, crop, ps, ws, wf, es,
					ef);
#else
				SnoPac(tid, snowng, dt, t24, prcpf, df1,
#include "SnoPac_kernel_arguments.c"
					   );
#endif
				elem_d_ps_eta_kinematic[tid] = (elem_d_wf_esnow[tid] + elem_d_wf_etns[tid]) * 1000.0;
			}

			/* Calculate effective mixing ratio at grnd level (skin) */
			elem_d_ps_q1[tid] = elem_d_ps_q2[tid] + elem_d_ps_eta_kinematic[tid] * CP / elem_d_ps_rch[tid];

			/* Determine sensible heat (H) in energy units (W m-2) */
			elem_d_ef_sheat[tid] = -(elem_d_ps_ch[tid] * CP * elem_d_ps_sfcprs[tid]) / 
				                (RD * t2v) * (elem_d_es_th2[tid] - elem_d_es_t1[tid]);

			/* Convert evap terms from rate (m s-1) to energy units (w m-2) */
			elem_d_ef_edir[tid] = elem_d_wf_edir[tid] * 1000.0 * LVH2O;
			elem_d_ef_ec[tid] = elem_d_wf_ec[tid] * 1000.0 * LVH2O;

			/*
			for (k = 0; k < elem_d_ps_nsoil[tid]; k++)
			{
				elem_d_ef_et[k][tid] = elem_d_wf_et[k][tid] * 1000.0 * LVH2O;
			}
			*/
			elem_d_ef_et0[tid] = elem_d_wf_et0[tid] * 1000.0 * LVH2O;
			elem_d_ef_et1[tid] = elem_d_wf_et1[tid] * 1000.0 * LVH2O;
			elem_d_ef_et2[tid] = elem_d_wf_et2[tid] * 1000.0 * LVH2O;
			elem_d_ef_et3[tid] = elem_d_wf_et3[tid] * 1000.0 * LVH2O;
			elem_d_ef_et4[tid] = elem_d_wf_et4[tid] * 1000.0 * LVH2O;
			elem_d_ef_et5[tid] = elem_d_wf_et5[tid] * 1000.0 * LVH2O;
			elem_d_ef_et6[tid] = elem_d_wf_et6[tid] * 1000.0 * LVH2O;
			elem_d_ef_et7[tid] = elem_d_wf_et7[tid] * 1000.0 * LVH2O;
			elem_d_ef_et8[tid] = elem_d_wf_et8[tid] * 1000.0 * LVH2O;
			elem_d_ef_et9[tid] = elem_d_wf_et9[tid] * 1000.0 * LVH2O;
			elem_d_ef_et10[tid] = elem_d_wf_et10[tid] * 1000.0 * LVH2O;

			elem_d_ef_ett[tid] = elem_d_wf_ett[tid] * 1000.0 * LVH2O;
			elem_d_ef_esnow[tid] = elem_d_wf_esnow[tid] * 1000.0 * LSUBS;
			elem_d_ef_etp[tid] =
				elem_d_wf_etp[tid] * 1000.0 * ((1.0 - elem_d_ps_sncovr[tid]) * LVH2O + elem_d_ps_sncovr[tid] * LSUBS);
			if (elem_d_ef_etp[tid] > 0.0)
			{
				elem_d_ef_eta[tid] = elem_d_ef_edir[tid] + elem_d_ef_ec[tid] + elem_d_ef_ett[tid] + elem_d_ef_esnow[tid];
#if defined(_CYCLES_)
				elem_d_ef_eta[tid] += elem_d_wf_eres * RHOH2O * LVH2O;
#endif
			}
			else
			{
				elem_d_ef_eta[tid] = elem_d_ef_etp[tid];
			}

			/* Determine beta (ratio of actual to potential evap) */
			elem_d_ps_beta[tid] = (elem_d_ef_etp[tid] == 0.0) ? 0.0 : (elem_d_ef_eta[tid] / elem_d_ef_etp[tid]);

			/* Convert the sign of soil heat flux so that:
			*   ssoil>0: warm the surface  (night time)
			*   ssoil<0: cool the surface  (day time) */
			elem_d_ef_ssoil[tid] *= -1.0;

			elem_d_ws_soilm[tid] = -1.0 * elem_d_ws_smc0[tid] * elem_d_ps_zsoil0[tid];
			/*
			for (k = 1; k < elem_d_ps_nsoil[tid]; k++)
			{
				elem_d_ws_soilm[tid] += elem_d_ws_smc[k][tid] * 
					(elem_d_ps_zsoil[k - 1][tid] - elem_d_ps_zsoil[k][tid]);
			}
			*/
			elem_d_ws_soilm[tid] +=
				elem_d_ws_smc1[tid] *
				(elem_d_ps_zsoil0[tid] - elem_d_ps_zsoil1[tid]) +
				elem_d_ws_smc2[tid] *
				(elem_d_ps_zsoil1[tid] - elem_d_ps_zsoil2[tid]) +
				elem_d_ws_smc3[tid] *
				(elem_d_ps_zsoil2[tid] - elem_d_ps_zsoil3[tid]) +
				elem_d_ws_smc4[tid] *
				(elem_d_ps_zsoil3[tid] - elem_d_ps_zsoil4[tid]) +
				elem_d_ws_smc5[tid] *
				(elem_d_ps_zsoil4[tid] - elem_d_ps_zsoil5[tid]) +
				elem_d_ws_smc6[tid] *
				(elem_d_ps_zsoil5[tid] - elem_d_ps_zsoil6[tid]) +
				elem_d_ws_smc7[tid] *
				(elem_d_ps_zsoil6[tid] - elem_d_ps_zsoil7[tid]) +
				elem_d_ws_smc8[tid] *
				(elem_d_ps_zsoil7[tid] - elem_d_ps_zsoil8[tid]) +
				elem_d_ws_smc9[tid] *
				(elem_d_ps_zsoil8[tid] - elem_d_ps_zsoil9[tid]) +
				elem_d_ws_smc10[tid] *
				(elem_d_ps_zsoil9[tid] - elem_d_ps_zsoil10[tid]);

			soilwm = -1.0 * (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]) * elem_d_ps_zsoil0[tid];
			soilww = -1.0 * (elem_d_ws_smc0[tid] - elem_d_soil_smcwlt[tid]) * elem_d_ps_zsoil0[tid];

			/*
			for (k = 0; k < elem_d_ps_nsoil[tid]; k++)
			{
				smav[k] = (elem_d_ws_smc[k][tid] - elem_d_soil_smcwlt[tid]) / 
					      (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]);
			}
			*/
			realtype	temp;
			temp = elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid];
			smav[0] = (elem_d_ws_smc0[tid] - elem_d_soil_smcwlt[tid]) / temp;
			smav[1] = (elem_d_ws_smc1[tid] - elem_d_soil_smcwlt[tid]) / temp;
			smav[2] = (elem_d_ws_smc2[tid] - elem_d_soil_smcwlt[tid]) / temp;
			smav[3] = (elem_d_ws_smc3[tid] - elem_d_soil_smcwlt[tid]) / temp;
			smav[4] = (elem_d_ws_smc4[tid] - elem_d_soil_smcwlt[tid]) / temp;
			smav[5] = (elem_d_ws_smc5[tid] - elem_d_soil_smcwlt[tid]) / temp;
			smav[6] = (elem_d_ws_smc6[tid] - elem_d_soil_smcwlt[tid]) / temp;
			smav[7] = (elem_d_ws_smc7[tid] - elem_d_soil_smcwlt[tid]) / temp;
			smav[8] = (elem_d_ws_smc8[tid] - elem_d_soil_smcwlt[tid]) / temp;
			smav[9] = (elem_d_ws_smc9[tid] - elem_d_soil_smcwlt[tid]) / temp;
			smav[10] = (elem_d_ws_smc10[tid] - elem_d_soil_smcwlt[tid]) / temp;

			// 根系分层也是变化的，在CVode中又不能用二维数组，咋整???
			if (elem_d_ps_nroot[tid] > 1)
			{
		        /*
				for (k = 1; k < elem_d_ps_nroot[tid]; k++)
				{
					soilwm += (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil[k - 1][tid] - elem_d_ps_zsoil[k][tid]);
					soilww += (elem_d_ws_smc[k][tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil[k - 1][tid] - elem_d_ps_zsoil[k][tid]);
				}
			     */
				if (elem_d_ps_nroot[tid] == 2){
					soilwm += (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil0[tid] - elem_d_ps_zsoil1[tid]);
					soilwm += (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil1[tid] - elem_d_ps_zsoil2[tid]);
					soilww += (elem_d_ws_smc1[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil0[tid] - elem_d_ps_zsoil1[tid]);
					soilww += (elem_d_ws_smc2[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil1[tid] - elem_d_ps_zsoil2[tid]);
				}
				else if (elem_d_ps_nroot[tid] == 3){
					soilwm += (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil0[tid] - elem_d_ps_zsoil1[tid]);
					soilwm += (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil1[tid] - elem_d_ps_zsoil2[tid]);
					soilwm += (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil2[tid] - elem_d_ps_zsoil3[tid]);
					soilww += (elem_d_ws_smc1[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil0[tid] - elem_d_ps_zsoil1[tid]);
					soilww += (elem_d_ws_smc2[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil1[tid] - elem_d_ps_zsoil2[tid]);
					soilww += (elem_d_ws_smc3[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil2[tid] - elem_d_ps_zsoil3[tid]);
				}
				else if (elem_d_ps_nroot[tid] == 4){
					soilwm += (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil0[tid] - elem_d_ps_zsoil1[tid]);
					soilwm += (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil1[tid] - elem_d_ps_zsoil2[tid]);
					soilwm += (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil2[tid] - elem_d_ps_zsoil3[tid]);
					soilwm += (elem_d_soil_smcmax[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil3[tid] - elem_d_ps_zsoil4[tid]);
					soilww += (elem_d_ws_smc1[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil0[tid] - elem_d_ps_zsoil1[tid]);
					soilww += (elem_d_ws_smc2[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil1[tid] - elem_d_ps_zsoil2[tid]);
					soilww += (elem_d_ws_smc3[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil2[tid] - elem_d_ps_zsoil3[tid]);
					soilww += (elem_d_ws_smc4[tid] - elem_d_soil_smcwlt[tid]) *
						(elem_d_ps_zsoil3[tid] - elem_d_ps_zsoil4[tid]);
				}
				else{
					/* 还需要计算其他根系分层的情况 */


				}
			}
			
			if (soilwm < 1.0e-6)
			{
				soilwm = 0.0;
				elem_d_ps_soilw[tid] = 0.0;
				elem_d_ws_soilm[tid] = 0.0;
			}
			else
			{
				elem_d_ps_soilw[tid] = soilww / soilwm;
			}

// Fused some state variables here.
			/* ET: convert from W m-2 to m s-1 */
			elem_d_wf_ec[tid] = elem_d_ef_ec[tid] / LVH2O / 1000.0;
			elem_d_wf_ett[tid] = elem_d_ef_ett[tid] / LVH2O / 1000.0;
			elem_d_wf_edir[tid] = elem_d_ef_edir[tid] / LVH2O / 1000.0;
	} // if (tid < nelem)
}
//-------- The end of Noah_kernel() ----------------------------


__device__ void AlCalc(int tid, realtype dt, int snowng,
#include "AlCalc_kernel_declarations.c" 	
)
{
	/*
	* Calculate albedo including snow effect (0 -> 1)
	*/
	/* snoalb is argument representing maximum albedo over deep snow, as passed
	* into SFlx, and adapted from the satellite-based maximum snow albedo
	* fields provided by D. Robinson and G. Kukla (1985, JCAM, Vol 24,
	* 402-411) */
	realtype          snoalb2;
	realtype          snoalb1;
	const realtype    SNACCA = 0.94;
	const realtype    SNACCB = 0.58;

	/* Turn of vegetation effect */
	//elem_d_ps_albedo = elem_d_ps_alb + (1.0 - (elem_d_lc_shdfac - elem_d_lc_shdmin)) * elem_d_ps_sncovr *
	//  (elem_d_ps_snoalb - elem_d_ps_alb);
	//elem_d_ps_albedo = (1.0 - elem_d_ps_sncovr) * elem_d_ps_alb +
	//  elem_d_ps_sncovr * elem_d_ps_snoalb;    /* this is equivalent to below */

	elem_d_ps_albedo[tid] = elem_d_ps_alb[tid] + elem_d_ps_sncovr[tid] * 
		                    (elem_d_ps_snoalb[tid] - elem_d_ps_alb[tid]);
	elem_d_ps_emissi[tid] = elem_d_ps_embrd[tid] + elem_d_ps_sncovr[tid] * 
		                    (EMISSI_S - elem_d_ps_embrd[tid]);

	/* Formulation by livneh
	* snoalb is considered as the maximum snow albedo for new snow, at a value
	* of 85%. Snow albedo curve defaults are from Bras P.263. should not be
	* changed except for serious problems with snow melt.
	* To implement accumulating parameters, SNACCA and SNACCB, assert that it
	* is indeed accumulation season. i.e. that snow surface temp is below zero
	* and the date falls between October and February */
	snoalb1 = elem_d_ps_snoalb[tid] + elem_d_ps_lvcoef[tid] * (0.85 - elem_d_ps_snoalb[tid]);
	snoalb2 = snoalb1;

	/* Initial lstsnw */
	if (snowng)
	{
		elem_d_ps_snotime1[tid] = 0.0;
	}
	else
	{
		elem_d_ps_snotime1[tid] += dt;
		snoalb2 = snoalb1 * pow(SNACCA, pow(elem_d_ps_snotime1[tid] / 86400.0, SNACCB));
	}

	snoalb2 = (snoalb2 > elem_d_ps_alb[tid]) ? snoalb2 : elem_d_ps_alb[tid];
	elem_d_ps_albedo[tid] = elem_d_ps_alb[tid] + elem_d_ps_sncovr[tid] * (snoalb2 - elem_d_ps_alb[tid]);
	elem_d_ps_albedo[tid] = (elem_d_ps_albedo[tid] > snoalb2) ? snoalb2 : elem_d_ps_albedo[tid];
}


__constant__ realtype    SLV = 2.501000e6;	
__constant__ realtype    RC = 70.0;
#if defined(_CYCLES_)
void CanRes(const estate_struct *es, pstate_struct *ps)
#else
__device__ void CanRes(int tid,
#include "CanRes_kernel_declarations.c"
)
#endif
{
	/*
	* Function CanRes
	*
	* Calculate canopy resistance which depends on incoming solar radiation,
	* air temperature, atmospheric water vapor pressure deficit at the lowest
	* model level, and soil moisture (preferably unfrozen soil moisture rather
	* than total)
	*
	* Source:  Jarvis (1976), Noilhan and Planton (1989, MWR), Jacquemin and
	* Noilhan (1990, BLM)
	* See also:  Chen et al. (1996, JGR, Vol 101(D3), 7251-7268), Eqns 12-14
	* and Table 2 of Sec. 3.1.2
	*/
	realtype          delta;
	realtype          rr;
#if !defined(_CYCLES_)
	realtype          ff;
	realtype          gx;
	int               k;
	realtype          part[MAXLYR];
#endif

#if defined(_CYCLES_)
	elem_d_ps_rc[tid] = RC;
#else
	/* Initialize canopy resistance multiplier terms. */
	elem_d_ps_rcs[tid] = 0.0;
	elem_d_ps_rct[tid] = 0.0;
	elem_d_ps_rcq[tid] = 0.0;
	elem_d_ps_rcsoil[tid] = 0.0;

	elem_d_ps_rc[tid] = 0.0;

	/* Contribution due to incoming solar radiation */
	ff = 0.55 * 2.0 * elem_d_ef_soldn[tid] / (elem_d_epc_rgl[tid] * elem_d_ps_proj_lai[tid]);
	elem_d_ps_rcs[tid] = (ff + elem_d_epc_rsmin[tid] / elem_d_epc_rsmax[tid]) / (1.0 + ff);
	elem_d_ps_rcs[tid] = (elem_d_ps_rcs[tid] > 0.0001) ? elem_d_ps_rcs[tid] : 0.0001;

	/* Contribution due to air temperature at first model level above ground rct
	* expression from Noilhan and Planton (1989, MWR). */
	elem_d_ps_rct[tid] = 1.0 - 0.0016 * pow(elem_d_epc_topt[tid] - elem_d_es_sfctmp[tid], 2.0);
	elem_d_ps_rct[tid] = (elem_d_ps_rct[tid] > 0.0001) ? elem_d_ps_rct[tid] : 0.0001;

	/* Contribution due to vapor pressure deficit at first model level.
	* rcq expression from ssib */
	elem_d_ps_rcq[tid] = 1.0 / (1.0 + elem_d_epc_hs[tid] * (elem_d_ps_q2sat[tid] - elem_d_ps_q2[tid]));
	elem_d_ps_rcq[tid] = (elem_d_ps_rcq[tid] > 0.01) ? elem_d_ps_rcq[tid] : 0.01;

	/* Contribution due to soil moisture availability.
	* Determine contribution from each soil layer, then add them up. */
	gx = (elem_d_ws_sh2o[0][tid] - elem_d_soil_smcwlt[tid]) / 
		(elem_d_soil_smcref[tid] - elem_d_soil_smcwlt[tid]);
	gx = (gx > 1.0) ? 1.0 : gx;
	gx = (gx < 0.0) ? 0.0 : gx;

	/* Use root distribution as weighting factor */
	//part[0] = rtdis[0] * gx;
	/* Use soil depth as weighting factor */
	part[0] = (elem_d_ps_zsoil[0][tid] / elem_d_ps_zsoil[elem_d_ps_nroot[tid] - 1][tid]) * gx;

	for (k = 1; k < elem_d_ps_nroot[tid]; k++)   // 其他分层土壤的part
	{
		// Frozen ground extension: total soil water "smc" was replaced by
		// unfrozen soil water "sh2o"
		gx = (elem_d_ws_sh2o[k][tid] - elem_d_soil_smcwlt[tid]) / 
			(elem_d_soil_smcref[tid] - elem_d_soil_smcwlt[tid]);
		gx = (gx > 1.0) ? 1.0 : gx;
		gx = (gx < 0.0) ? 0.0 : gx;

		// Use root distribution as weighting factor 
		// part[k] = rtdis[k] * gx;
		// Use soil depth as weighting factor 
		part[k] = ((elem_d_ps_zsoil[k][tid] - elem_d_ps_zsoil[k - 1][tid]) /
			elem_d_ps_zsoil[elem_d_ps_nroot[tid] - 1][tid]) * gx;
	}
	
	for (k = 0; k < elem_d_ps_nroot[tid]; k++){
		elem_d_ps_rcsoil[tid] += part[k];
	}
	elem_d_ps_rcsoil[tid] = (elem_d_ps_rcsoil[tid] > 0.0001) ? elem_d_ps_rcsoil[tid] : 0.0001;

	/* Determine canopy resistance due to all factors.  convert canopy
	* resistance (rc) to plant coefficient (pc) to be used with potential evap
	* in determining actual evap. pc is determined by:
	*   pc * linearized Penman potential evap =
	*   Penman-Monteith actual evaporation (containing rc term). */
	elem_d_ps_rc[tid] = elem_d_epc_rsmin[tid] /
		(elem_d_ps_proj_lai[tid] * elem_d_ps_rcs[tid] * elem_d_ps_rct[tid] *
		 elem_d_ps_rcq[tid] * elem_d_ps_rcsoil[tid]);
#endif

	rr = (4.0 * elem_d_ps_emissi[tid] * SIGMA * RD / CP) *
		pow(elem_d_es_sfctmp[tid], 4.0) / (elem_d_ps_sfcprs[tid] * elem_d_ps_ch[tid]) + 1.0;

	delta = (SLV / CP) * elem_d_ps_dqsdt2[tid];

	elem_d_ps_pc[tid] = (rr + delta) / (rr * (1.0 + elem_d_ps_rc[tid] * elem_d_ps_ch[tid]) + delta);
}

__device__ realtype CSnow(realtype dsnow)
{
	/*
	* Calculate snow thermal conductivity
	*/
	realtype          c;
	realtype          sncond;
	const realtype    UNIT = 0.11631;

	/* Sncond in units of cal/(cm*hr*c), returned in W/(m*C)
	* Csnow in units of cal/(cm*hr*c), returned in W/(m*C)
	* Basic version is Dyachkova equation (1960), for range 0.1-0.4 */
	c = 0.328 * pow(10.0, 2.25 * dsnow);

	/* De Vaux equation (1933), in range 0.1-0.6 */
	//sncond = 0.0293 * (1.0 + 100.0 * dsnow * dsnow);
	//csnow = 0.0293 * (1.0 + 100.0 * dsnow * dsnow);

	/* E. Andersen from Flerchinger */
	//sncond = 0.021 + 2.51 * dsnow * dsnow;
	//csnow = 0.021 + 2.51 * dsnow * dsnow;

	//sncond = UNIT * c;

	/* Double snow thermal conductivity */
	sncond = 2.0 * UNIT * c;

	return sncond;
}

#if defined(_CYCLES_)
void Evapo(const soil_struct *soil, const lc_struct *lc,
	const pstate_struct *ps, const estate_struct *es,
	const cstate_struct *cs, realtype dt, crop_struct crop[], wstate_struct *ws,
	wflux_struct *wf)
#else
__device__ void Evapo(int tid, realtype dt,
#include "Evapo_kernel_declarations.c"
)
#endif
{
	/*
	* Function Evapo
	*
	* Calculate soil moisture flux. The soil moisture content (smc - a per unit
	* volume measurement) is a dependent variable that is updated with
	* prognostic eqns. The canopy moisture content (cmc) is also updated.
	* Frozen ground version: new states added: sh2o, and frozen ground
	* correction factor, frzfact and parameter slope.
	*/
	int               k;
	realtype          cmc2ms;

	/* Executable code begins here if the potential evapotranspiration is
	* greater than zero. */
	elem_d_wf_edir[tid] = 0.0;
	elem_d_wf_ec[tid] = 0.0;
	elem_d_wf_ett[tid] = 0.0;
	
	/*
	for (k = 0; k < elem_d_ps_nsoil[tid]; k++){
		elem_d_wf_et[k][tid] = 0.0;
	}
	*/
	elem_d_wf_et0[tid] = 0.0;
	elem_d_wf_et1[tid] = 0.0;
	elem_d_wf_et2[tid] = 0.0;
	elem_d_wf_et3[tid] = 0.0;
	elem_d_wf_et4[tid] = 0.0;
	elem_d_wf_et5[tid] = 0.0;
	elem_d_wf_et6[tid] = 0.0;
	elem_d_wf_et7[tid] = 0.0;
	elem_d_wf_et8[tid] = 0.0;
	elem_d_wf_et9[tid] = 0.0;
	elem_d_wf_et10[tid] = 0.0;

	if (elem_d_wf_etp[tid] > 0.0)
	{
		if (elem_d_lc_shdfac[tid] < 1.0)
		{
			/* Retrieve direct evaporation from soil surface. Call this function
			* only if veg cover not complete.
			* Frozen ground version:  sh2o states replace smc states. */
			DEvap(tid,
#include "DEvap_kernel_arguments.c"	
				 );
		}

#if defined(_CYCLES_)
		/* Evaporation from residue (Cycles function) */
		ResidueEvaporation(elem_d_wf_etp * RHOH2O * dt, dt, elem_d_ps_sncovr, crop, ps, cs,
			ws, wf);
#endif

		if (elem_d_lc_shdfac[tid] > 0.0)
		{
			/* Initialize plant total transpiration, retrieve plant
			* transpiration, and accumulate it for all soil layers. */
#if defined(_CYCLES_)
			WaterUptake(soil, es, ps, dt, crop, ws, wf);
#else
			Transp(tid,
#include "Transp_kernel_arguments.c"				
				  );
#endif
			/*
			for (k = 0; k < elem_d_ps_nsoil[tid]; k++)
				elem_d_wf_ett[tid] += elem_d_wf_et[k][tid];
			*/
			elem_d_wf_ett[tid] += elem_d_wf_et0[tid] + elem_d_wf_et1[tid] +
				elem_d_wf_et2[tid] + elem_d_wf_et3[tid] + elem_d_wf_et4[tid] + elem_d_wf_et5[tid] +
				elem_d_wf_et6[tid] + elem_d_wf_et7[tid] + elem_d_wf_et8[tid] + elem_d_wf_et9[tid] +
				elem_d_wf_et10[tid];

			/* Calculate canopy evaporation.
			* If statements to avoid tangent linear problems near cmc = 0.0. */
			if (elem_d_ws_cmc[tid] > 0.0)
			{
				elem_d_wf_ec[tid] =
					elem_d_lc_shdfac[tid] * pow((elem_d_ws_cmc[tid] / elem_d_ws_cmcmax[tid] > 1.0) ?
					1.0 : elem_d_ws_cmc[tid] / elem_d_ws_cmcmax[tid], elem_d_lc_cfactr[tid]) *
					elem_d_wf_etp[tid];
			}
			else
			{
				elem_d_wf_ec[tid] = 0.0;
			}

			/* Ec should be limited by the total amount of available water on
			* the canopy. F.Chen, 18-oct-1994 */
			cmc2ms = elem_d_ws_cmc[tid] / dt;
			elem_d_wf_ec[tid] = (cmc2ms < elem_d_wf_ec[tid]) ? cmc2ms : elem_d_wf_ec[tid];
		}
	}

	/* Total up evap and transp types to obtain actual evapotransp */
#if defined(_CYCLES_)
	elem_d_wf_etns[tid] = elem_d_wf_edir[tid] + elem_d_wf_ett[tid] + elem_d_wf_ec[tid] + elem_d_wf_eres[tid];
#else
	elem_d_wf_etns[tid] = elem_d_wf_edir[tid] + elem_d_wf_ett[tid] + elem_d_wf_ec[tid];
#endif
}

__device__ void DEvap(int tid,
#include "DEvap_kernel_declarations.c"
)
{
	/*
	* Function DEvap
	*
	* Calculate direct soil evaporation
	*/
	realtype          fx, sratio;

	/* Direct evap a function of relative soil moisture availability, linear
	* when fxexp = 1.
	* fx > 1 represents demand control
	* fx < 1 represents flux control */
	sratio = (elem_d_ws_sh2o[0][tid] - elem_d_soil_smcdry[tid]) /
		     (elem_d_soil_smcmax[tid] - elem_d_soil_smcdry[tid]);
	if (sratio > 0.0)
	{
		fx = pow(sratio, elem_d_ps_fxexp[tid]);
		fx = (fx > 1.0) ? 1.0 : fx;
		fx = (fx < 0.0) ? 0.0 : fx;
	}
	else
	{
		fx = 0.0;
	}

	/* Allow for the direct-evap-reducing effect of shade */
	elem_d_wf_edir[tid] = fx * (1.0 - elem_d_lc_shdfac[tid]) * elem_d_wf_etp[tid];
}


__device__ void Transp(int tid,
#include "Transp_kernel_declarations.c"
)
{
	/*
	* Function Transp
	*
	* Calculate transpiration for the veg class.
	*/
	int               i, k;
	realtype          denom;
	realtype          etpa;
	realtype          gx[MAXLYR];
	realtype          rtx[MAXLYR], sgx;

	/* Initialize plant transp to zero for all soil layers. */
	/*
	for (k = 0; k < elem_d_ps_nsoil[tid]; k++){
		elem_d_wf_et[k][tid] = 0.0;
	}
	*/
	elem_d_wf_et0[tid] = 0.0;
	elem_d_wf_et1[tid] = 0.0;
	elem_d_wf_et2[tid] = 0.0;
	elem_d_wf_et3[tid] = 0.0;
	elem_d_wf_et4[tid] = 0.0;
	elem_d_wf_et5[tid] = 0.0;
	elem_d_wf_et6[tid] = 0.0;
	elem_d_wf_et7[tid] = 0.0;
	elem_d_wf_et8[tid] = 0.0;
	elem_d_wf_et9[tid] = 0.0;
	elem_d_wf_et10[tid] = 0.0;

	/* Calculate an 'adjusted' potential transpiration
	* If statement below to avoid tangent linear problems near zero
	* Note: gx and other terms below redistribute transpiration by layer,
	* et(k), as a function of soil moisture availability, while preserving
	* total etpa. */
	if (elem_d_ws_cmc[tid] != 0.0)
	{
		etpa = elem_d_lc_shdfac[tid] * elem_d_ps_pc[tid] * elem_d_wf_etp[tid] *
			(1.0 - pow(elem_d_ws_cmc[tid] / elem_d_ws_cmcmax[tid], elem_d_lc_cfactr[tid]));
	}
	else
	{
		etpa = elem_d_lc_shdfac[tid] * elem_d_ps_pc[tid] * elem_d_wf_etp[tid];
	}

	sgx = 0.0;

	// CUDA version.
	realtype     temp;
	temp = elem_d_soil_smcref[tid] - elem_d_soil_smcwlt[tid];
	gx[0] = (elem_d_ws_smc0[tid] - elem_d_soil_smcwlt[tid]) / temp;
	gx[1] = (elem_d_ws_smc1[tid] - elem_d_soil_smcwlt[tid]) / temp;
	gx[2] = (elem_d_ws_smc2[tid] - elem_d_soil_smcwlt[tid]) / temp;
	gx[3] = (elem_d_ws_smc3[tid] - elem_d_soil_smcwlt[tid]) / temp;
	gx[4] = (elem_d_ws_smc4[tid] - elem_d_soil_smcwlt[tid]) / temp;
	gx[5] = (elem_d_ws_smc5[tid] - elem_d_soil_smcwlt[tid]) / temp;
	gx[6] = (elem_d_ws_smc6[tid] - elem_d_soil_smcwlt[tid]) / temp;
	gx[7] = (elem_d_ws_smc7[tid] - elem_d_soil_smcwlt[tid]) / temp;
	gx[8] = (elem_d_ws_smc8[tid] - elem_d_soil_smcwlt[tid]) / temp;
	gx[9] = (elem_d_ws_smc9[tid] - elem_d_soil_smcwlt[tid]) / temp;
	gx[10] = (elem_d_ws_smc10[tid] - elem_d_soil_smcwlt[tid]) / temp;

	for (i = 0; i < elem_d_ps_nroot[tid]; i++)    
	{
		//gx[i] = (elem_d_ws_smc[i][tid] - elem_d_soil_smcwlt[tid]) /
		//	     (elem_d_soil_smcref[tid] - elem_d_soil_smcwlt[tid]);
		gx[i] = (gx[i] < 0.0) ? 0.0 : gx[i];
		gx[i] = (gx[i] > 1.0) ? 1.0 : gx[i];
		sgx += gx[i];
	}
	sgx /= (realtype)elem_d_ps_nroot[tid];

	denom = 0.0;
	// CUDA version.
	rtx[0] = elem_d_ps_rtdis0[tid] + gx[i] - sgx;
	rtx[1] = elem_d_ps_rtdis1[tid] + gx[i] - sgx;
	rtx[2] = elem_d_ps_rtdis2[tid] + gx[i] - sgx;
	rtx[3] = elem_d_ps_rtdis3[tid] + gx[i] - sgx;
	rtx[4] = elem_d_ps_rtdis4[tid] + gx[i] - sgx;
	rtx[5] = elem_d_ps_rtdis5[tid] + gx[i] - sgx;
	rtx[6] = elem_d_ps_rtdis6[tid] + gx[i] - sgx;
	rtx[7] = elem_d_ps_rtdis7[tid] + gx[i] - sgx;
	rtx[8] = elem_d_ps_rtdis8[tid] + gx[i] - sgx;
	rtx[9] = elem_d_ps_rtdis9[tid] + gx[i] - sgx;
	rtx[10] = elem_d_ps_rtdis10[tid] + gx[i] - sgx;

	for (i = 0; i < elem_d_ps_nroot[tid]; i++)
	{
	//	rtx[i] = elem_d_ps_rtdis[i][tid] + gx[i] - sgx;  // above code
		gx[i] *= ((rtx[i] > 0.0) ? rtx[i] : 0.0);
		denom += gx[i];
	}

	denom = (denom <= 0.0) ? 1.0 : denom;

	/* 根系在垂向上均匀分布的情况 */
	/*
	for (i = 0; i < elem_d_ps_nroot[tid]; i++){
		elem_d_wf_et[i][tid] = etpa * gx[i] / denom;
	}
	*/
	elem_d_wf_et0[tid] = etpa * gx[0] / denom;
	elem_d_wf_et1[tid] = etpa * gx[1] / denom;
	elem_d_wf_et2[tid] = etpa * gx[2] / denom;
	elem_d_wf_et3[tid] = etpa * gx[3] / denom;
	elem_d_wf_et4[tid] = etpa * gx[4] / denom;
	elem_d_wf_et5[tid] = etpa * gx[5] / denom;
	elem_d_wf_et6[tid] = etpa * gx[6] / denom;
	elem_d_wf_et7[tid] = etpa * gx[7] / denom;
	elem_d_wf_et8[tid] = etpa * gx[8] / denom;
	elem_d_wf_et9[tid] = etpa * gx[9] / denom;
	elem_d_wf_et10[tid] = etpa * gx[10] / denom;


	/* 考虑根系在垂向上不均匀分布的情况 */
#if NOT_YET_IMPLEMENTED
	/* Above code assumes a vertically uniform root distribution
	* Code below tests a variable root distribution */
	elem_d_wf_et[0] = (zsoil[0] / zsoil[elem_d_ps_nroot - 1]) * gx * etpa;
	elem_d_wf_et[0] = (zsoil[0] / zsoil[elem_d_ps_nroot - 1]) * etpa;
	/* Using root distribution as weighting factor */
	elem_d_wf_et[0] = (elem_d_lc_rtdis[0] * etpa);
	elem_d_wf_et[0] = etpa * part[0];
	/* Loop down thru the soil layers repeating the operation above, but using
	* the thickness of the soil layer (rather than the absolute depth of each
	* layer) in the final calculation. */
	for (k = 0; k < elem_d_ps_nroot; k++)
	{
		gx = (elem_d_ws_smc[k] - elem_d_soil_smcwlt) / (elem_d_soil_smcref - elem_d_soil_smcwlt);
		gx = (gx < 0.0) ? 0.0 : gx;
		gx = (gx > 1.0) ? 1.0 : gx;
		/* Test canopy resistance */
		gx = 1.0;
		elem_d_wf_et[k] = ((zsoil[k] - zsoil[k - 1]) /
			zsoil[elem_d_ps_nroot - 1]) * gx * etpa;
		elem_d_wf_et[k] = ((zsoil[k] - zsoil[k - 1]) /
			zsoil[elem_d_ps_nroot - 1]) * etpa;
		/* Using root distribution as weighting factor */
		elem_d_wf_et[k] = elem_d_lc_rtdis[k] * etpa;
		elem_d_wf_et[k] = etpa * part[k];
	}
#endif

}



__device__ realtype FrozRain(realtype prcp, realtype sfctmp)
{
	realtype          ffrozp;

	if (prcp > 0.0 && sfctmp < TFREEZ)
	{
		ffrozp = 1.0;
	}
	else
	{
		ffrozp = 0.0;
	}

	return ffrozp;
}


__constant__ realtype    CK = 8.0;
__constant__ realtype    ERROR0 = 0.005;
__device__ realtype FrH2O(int tid, realtype tkelv, realtype smc, realtype sh2o, 
	realtype *elem_d_soil_alpha,
	realtype *elem_d_soil_beta,
	realtype *elem_d_soil_smcmin,
	realtype *elem_d_soil_smcmax
)
{
	/*
	* Function FrH2O
	*
	* Calculate amount of supercooled liquid soil water content if temperature
	* is below 273.15K (t0). Requires Newton-type iteration to solve the
	* nonlinear implicit equation given in Eqn 17 of Koren et al
	* (1999, JGR, Vol 104(D16), 19569-19585).
	*
	* New version (June 2001): much faster and more accurate Newton iteration
	* achieved by first taking log of eqn cited above -- less than 4
	* (typically 1 or 2) iterations achieves convergence. Also, explicit
	* 1-step solution option for special case of parameter ck = 0, which
	* reduces the original implicit equation to a simpler explicit form, known
	* as the "Flerchinger eqn". Improved handling of solution in the limit of
	* freezing point temperature t0.
	*
	* In Flux-PIHM, van Genuchten parameters are used. See Technical
	* Documentation for details
	*/
	realtype          denom;
	realtype          df;
	realtype          dswl;
	realtype          fk;
	realtype          swl;
	realtype          swlk;
	realtype          satn;
	int               nlog, kcount;
	realtype          mx;
	realtype          freew;

	nlog = 0;
	kcount = 0;

	/* If temperature not significantly below freezing (t0), sh2o = smc */
	if (tkelv > (TFREEZ - 1.0e-3))
	{
		freew = smc;
	}
	else
	{
		/* Option 1: iterated solution for nonzero CK in Koren et al, JGR, 1999,
		* Eqn 17 */
		if (CK != 0.0)
		{
			/* Initial guess for swl (frozen content) */
			swl = smc - sh2o;

			/* Keep within bounds. */
			swl = (swl > smc - SH2OMIN) ? (smc - SH2OMIN) : swl;
			swl = (swl < 0.0) ? 0.0 : swl;

			/* Start of iterations */
			while (nlog < 10 && kcount == 0)
			{
				nlog++;

				satn =
					(smc - swl - elem_d_soil_smcmin[tid]) / 
					(elem_d_soil_smcmax[tid] - elem_d_soil_smcmin[tid]);
				mx = elem_d_soil_beta[tid] / (1.0 - elem_d_soil_beta[tid]);

				df = log(GRAV / elem_d_soil_alpha[tid] / LSUBF) +
					1.0 / elem_d_soil_beta[tid] * log(pow(satn, mx) - 1.0) +
					2.0 * log(1.0 + CK * swl) - log(-(tkelv - TFREEZ) / tkelv);

				denom = 1.0 / (elem_d_soil_beta[tid] - 1.0) /
					(elem_d_soil_smcmax[tid] - elem_d_soil_smcmin[tid]) * pow(satn, mx - 1.0) /
					(pow(satn, mx) - 1.0) + 2.0 * CK / (1.0 + CK * swl);

				swlk = swl - df / denom;

				/* Bounds useful for mathematical solution. */
				swlk = (swlk > smc - SH2OMIN) ? (smc - SH2OMIN) : swlk;
				swlk = (swlk < 0.0) ? 0.0 : swlk;

				/* Mathematical solution bounds applied. */
				dswl = fabs(swlk - swl);

				/* If more than 10 iterations, use explicit method (ck = 0
				* approx.) when dswl less or eq. error, no more iterations
				* required. */
				swl = swlk;
				if (dswl <= ERROR0)
				{
					kcount++;
				}
				/* End of iterations
				* Bounds applied within do-block are valid for physical
				* solution */
			}

			freew = smc - swl;
		}

		/* Option 2: explicit solution for flerchinger eq. i.e. ck = 0 in Koren
		* et al., JGR, 1999, Eqn 17
		* Apply physical bounds to flerchinger solution */
		if (kcount == 0)
		{
			fk = pow(pow(-(tkelv - TFREEZ) / tkelv * elem_d_soil_alpha[tid] *
				LSUBF / GRAV, elem_d_soil_beta[tid]), 1.0 / mx) *
				(elem_d_soil_smcmax[tid] - elem_d_soil_smcmin[tid]) - elem_d_soil_smcmin[tid];
			fk = (fk < SH2OMIN) ? SH2OMIN : fk;

			freew = (fk < smc) ? fk : smc;
		}
	}

	return freew;
}

__constant__ realtype    CH2O = 4.2e6;
// 计算土壤中的热扩散方程的时间项(RHS)
__device__ void HRT(int tid,
	realtype *rhsts, realtype yy, realtype zz1, realtype dt, realtype df1, 
	realtype *ai, realtype *bi, realtype *ci,
#include "HRT_kernel_declarations.c"
)
{
	/*
	* Function HRT
	*
	* Calculate the right hand side of the time tendency term of the soil
	* thermal diffusion equation. Also to compute (prepare) the matrix
	* coefficients for the tri-diagonal matrix of the implicit time scheme.
	*/
	int               itavg;
	int               k;
	realtype          ddz, ddz2;
	realtype          denom;
	realtype          df1k;
	realtype          dtsdz;
	realtype          dtsdz2;
	realtype          hcpct;
	realtype          sice;
	realtype          csoil_loc;
	realtype          df1n;
	realtype          qtot;
	realtype          tavg;
	realtype          tbk;
	realtype          tbk1;
	realtype          tsnsr;
	realtype          tsurf;

	/* Urban */
	csoil_loc = (elem_d_lc_isurban[tid]) ? 3.0e6 : elem_d_soil_csoil[tid];

	/* Initialize logical for soil layer temperature averaging. */
	itavg = 1;

	/* Begin section for top soil layer
	// k==0			// 表层土壤温度
	* Calc the heat capacity of the top soil layer */
	hcpct = elem_d_ws_sh2o0[tid] * CH2O + (1.0 - elem_d_soil_smcmax[tid]) * csoil_loc +
		(elem_d_soil_smcmax[tid] - elem_d_ws_smc0[tid])
		* CP + (elem_d_ws_smc0[tid] - elem_d_ws_sh2o0[tid]) * CPICE;

	/* Calc the matrix coefficients ai, bi, and ci for the top layer */
	ddz = 1.0 / (-0.5 * elem_d_ps_zsoil1[tid]);
	ai[0] = 0.0;
	ci[0] = (df1 * ddz) / (elem_d_ps_zsoil0[tid] * hcpct);

	/* Calculate the vertical soil temp gradient btwn the 1st and 2nd soil
	* layers. Then calculate the subsurface heat flux. use the temp gradient
	* and subsfc heat flux to calc "right-hand side tendency terms", or
	* "rhsts", for top soil layer. */
	bi[0] = -ci[0] + df1 / (0.5 * elem_d_ps_zsoil0[tid] * elem_d_ps_zsoil0[tid] * hcpct * zz1);
	dtsdz = (elem_d_es_stc0[tid] - elem_d_es_stc1[tid]) / (-0.5 * elem_d_ps_zsoil1[tid]);
	elem_d_ef_ssoil[tid] = df1 * (elem_d_es_stc0[tid] - yy) / (0.5 * elem_d_ps_zsoil0[tid] * zz1);
	denom = elem_d_ps_zsoil0[tid] * hcpct;

	rhsts[0] = (df1 * dtsdz - elem_d_ef_ssoil[tid]) / denom;

	/* Next capture the vertical difference of the heat flux at top and bottom
	* of first soil layer for use in heat flux constraint applied to potential
	* soil freezing/thawing in routine SnkSrc. */
	qtot = -1.0 * rhsts[0] * denom;

	/* Calculate frozen water content in 1st soil layer. */
	sice = elem_d_ws_smc0[tid] - elem_d_ws_sh2o0[tid];

	/* If temperature averaging invoked (itavg = true; else skip):
	* Set temp "tsurf" at top of soil column (for use in freezing soil physics
	* later in function subroutine SnkSrc). If snowpack content is zero, then
	* tsurf expression below gives tsurf = skin temp. If snowpack is nonzero
	* (hence argument zz1 = 1), then tsurf expression below yields soil column
	* top temperature under snowpack. Then calculate temperature at bottom
	* interface of 1st soil layer for use later in function subroutine SnkSrc
	*/
	if (itavg)
	{
		tsurf = (yy + (zz1 - 1.0) * elem_d_es_stc0[tid]) / zz1;  

		/* If frozen water present or any of layer-1 mid-point or bounding
		* interface temperatures below freezing, then call SnkSrc to compute
		* heat source/sink (and change in frozen water content) due to possible
		* soil water phase change */
		/*
		tbk = TBnd(elem_d_es_stc0[tid], elem_d_es_stc1[tid], elem_d_ps_zsoil[tid],  
			       elem_d_ps_zbot[tid], 0, elem_d_ps_nsoil[tid]);
		*/
		tbk = TBnd(elem_d_es_stc0[tid], elem_d_es_stc1[tid], elem_d_ps_zbot[tid], 0,// k ==0
			elem_d_ps_zsoil0[tid], elem_d_ps_zsoil1[tid], elem_d_ps_zsoil2[tid],
			elem_d_ps_zsoil3[tid], elem_d_ps_zsoil4[tid], elem_d_ps_zsoil5[tid],
			elem_d_ps_zsoil6[tid], elem_d_ps_zsoil7[tid], elem_d_ps_zsoil8[tid],
			elem_d_ps_zsoil9[tid], elem_d_ps_zsoil10[tid]);

		// 顶层土壤
		if ((sice > 0.0) || (elem_d_es_stc0[tid] < TFREEZ) || (tsurf < TFREEZ) ||
			(tbk < TFREEZ))
		{
		//	tavg = TmpAvg(tsurf, elem_d_es_stc0[tid], tbk, elem_d_ps_zsoil[tid], 0);
			tavg = TmpAvg(tsurf, elem_d_es_stc0[tid], tbk, 0,   // k==0
				elem_d_ps_zsoil0[tid], elem_d_ps_zsoil1[tid], elem_d_ps_zsoil2[tid],
				elem_d_ps_zsoil3[tid], elem_d_ps_zsoil4[tid], elem_d_ps_zsoil5[tid],
				elem_d_ps_zsoil6[tid], elem_d_ps_zsoil7[tid], elem_d_ps_zsoil8[tid],
				elem_d_ps_zsoil9[tid], elem_d_ps_zsoil10[tid]);

			SnkSrc(tid, &tsnsr, tavg, elem_d_ws_smc0[tid], &elem_d_ws_sh2o0[tid], 
				elem_d_ps_zsoil0[tid], elem_d_ps_zsoil1[tid], elem_d_ps_zsoil2[tid],
				elem_d_ps_zsoil3[tid], elem_d_ps_zsoil4[tid], elem_d_ps_zsoil5[tid],
				elem_d_ps_zsoil6[tid], elem_d_ps_zsoil7[tid], elem_d_ps_zsoil8[tid],
				elem_d_ps_zsoil9[tid], elem_d_ps_zsoil10[tid],
					dt, 0, qtot,  // k==0
					elem_d_soil_alpha,
					elem_d_soil_beta,
					elem_d_soil_smcmin,
					elem_d_soil_smcmax
					);
				
			rhsts[0] -= tsnsr / denom;
		}
	}
	else
	{
		if ((sice > 0.0) || (elem_d_es_stc0[tid] < TFREEZ))
		{
			SnkSrc(tid, &tsnsr, elem_d_es_stc0[tid], elem_d_ws_smc0[tid], &elem_d_ws_sh2o0[tid],
				elem_d_ps_zsoil0[tid], elem_d_ps_zsoil1[tid], elem_d_ps_zsoil2[tid],
				elem_d_ps_zsoil3[tid], elem_d_ps_zsoil4[tid], elem_d_ps_zsoil5[tid],
				elem_d_ps_zsoil6[tid], elem_d_ps_zsoil7[tid], elem_d_ps_zsoil8[tid],
				elem_d_ps_zsoil9[tid], elem_d_ps_zsoil10[tid],	
				dt, 0, qtot,
				elem_d_soil_alpha,
				elem_d_soil_beta,
				elem_d_soil_smcmin,
				elem_d_soil_smcmax
				);

			rhsts[0] -= tsnsr / denom;
		}    /* This ends section for top soil layer. */
	}  // 顶层土壤温度的平均

	/* Initialize ddz2 */
	ddz2 = 0.0;
	df1k = df1;

	/* Loop thru the remaining soil layers, repeating the above process (except
	* subsfc or "ground" heat flux not repeated in lower layers)
	* Calculate heat capacity for this soil layer. */
	// CPU version   for()
//#include "HRT_soil_unroll.c"

	// CUDA version     unroll for
// k == 1
	hcpct = elem_d_ws_sh2o1[tid] * CH2O + (1.0 - elem_d_soil_smcmax[tid]) * csoil_loc +
		(elem_d_soil_smcmax[tid] - elem_d_ws_smc1[tid]) * CP +
		(elem_d_ws_smc1[tid] - elem_d_ws_sh2o1[tid]) * CPICE;

	if (1 != elem_d_ps_nsoil[tid] - 1)
	{
		/* Calc the vertical soil temp gradient thru this layer */
		df1n = TDfCnd(elem_d_ws_smc1[tid], elem_d_soil_quartz[tid], elem_d_soil_smcmax[tid],
			elem_d_soil_smcmin[tid], elem_d_ws_sh2o1[tid]);

		/* Urban */
		df1n = (elem_d_lc_isurban[tid]) ? 3.24 : df1n;

		denom = 0.5 * (elem_d_ps_zsoil0[tid] - elem_d_ps_zsoil2[tid]);

		/* Calc the matrix coef, ci, after calc'ng its partial product */
		dtsdz2 = (elem_d_es_stc1[tid] - elem_d_es_stc2[tid]) / denom;
		ddz2 = 2.0 / (elem_d_ps_zsoil0[tid] - elem_d_ps_zsoil2[tid]);

		/* If temperature averaging invoked (itavg = true; else skip)
		* Calculate temp at bottom of layer. */
		ci[k] = -df1n * ddz2 / ((elem_d_ps_zsoil0[tid] - elem_d_ps_zsoil1[tid]) * hcpct);
		if (itavg)
		{
			tbk1 = TBnd(elem_d_es_stc1[tid], elem_d_es_stc2[tid], elem_d_ps_zbot[tid], 1,  // k==1
				elem_d_ps_zsoil0[tid], elem_d_ps_zsoil1[tid], elem_d_ps_zsoil2[tid],
				elem_d_ps_zsoil3[tid], elem_d_ps_zsoil4[tid], elem_d_ps_zsoil5[tid],
				elem_d_ps_zsoil6[tid], elem_d_ps_zsoil7[tid], elem_d_ps_zsoil8[tid],
				elem_d_ps_zsoil9[tid], elem_d_ps_zsoil10[tid]);
		}
	}  
	else
	{
		/* Special case of bottom soil layer
		* Calculate thermal diffusivity for bottom layer. */
		df1n = TDfCnd(elem_d_ws_smc1[tid], elem_d_soil_quartz[tid], elem_d_soil_smcmax[tid],
			elem_d_soil_smcmin[tid],
			elem_d_ws_sh2o1[tid]);

		/* Urban */
		df1n = (elem_d_lc_isurban[tid]) ? 3.24 : df1n;

		/* Calc the vertical soil temp gradient thru bottom layer. */
		denom = 0.5 * (elem_d_ps_zsoil0[tid] + elem_d_ps_zsoil1[tid]) - elem_d_ps_zbot[tid];
		dtsdz2 = (elem_d_es_stc1[tid] - elem_d_ps_tbot[tid]) / denom;

		/* Set matrix coef, ci to zero if bottom layer. */
		ci[k] = 0.0;

		/* If temperature averaging invoked (itavg = true; else skip)
		* Calculate temp at bottom of last layer. */
		if (itavg)
		{
			tbk1 = TBnd(elem_d_es_stc1[tid], elem_d_ps_tbot[tid], elem_d_ps_zbot[tid], 1,  // k==1
				elem_d_ps_zsoil0[tid], elem_d_ps_zsoil1[tid], elem_d_ps_zsoil2[tid],
				elem_d_ps_zsoil3[tid], elem_d_ps_zsoil4[tid], elem_d_ps_zsoil5[tid],
				elem_d_ps_zsoil6[tid], elem_d_ps_zsoil7[tid], elem_d_ps_zsoil8[tid],
				elem_d_ps_zsoil9[tid], elem_d_ps_zsoil10[tid]);

		}    /* This ends special loop for bottom layer. */
	}  // if( k != elem_d_ps_nsoil[tid] - 1)

	denom = (elem_d_ps_zsoil1[tid] - elem_d_ps_zsoil0[tid]) * hcpct;
	rhsts[1] = (df1n * dtsdz2 - df1k * dtsdz) / denom;
	qtot = -1.0 * denom * rhsts[1];

	sice = elem_d_ws_smc1[tid] - elem_d_ws_sh2o1[tid];

	if (itavg)
	{
		tavg = TmpAvg(tbk, elem_d_es_stc1[tid], tbk1, 1,
			elem_d_ps_zsoil0[tid], elem_d_ps_zsoil1[tid], elem_d_ps_zsoil2[tid],
			elem_d_ps_zsoil3[tid], elem_d_ps_zsoil4[tid], elem_d_ps_zsoil5[tid],
			elem_d_ps_zsoil6[tid], elem_d_ps_zsoil7[tid], elem_d_ps_zsoil8[tid],
			elem_d_ps_zsoil9[tid], elem_d_ps_zsoil10[tid]);

		if ((sice > 0.0) || (elem_d_es_stc1[tid] < TFREEZ) || (tbk < TFREEZ) || (tbk1 < TFREEZ))
		{
			SnkSrc(tid, &tsnsr, tavg, elem_d_ws_smc1[tid], &elem_d_ws_sh2o1[tid],
				elem_d_ps_zsoil0[tid], elem_d_ps_zsoil1[tid], elem_d_ps_zsoil2[tid],
				elem_d_ps_zsoil3[tid], elem_d_ps_zsoil4[tid], elem_d_ps_zsoil5[tid],
				elem_d_ps_zsoil6[tid], elem_d_ps_zsoil7[tid], elem_d_ps_zsoil8[tid],
				elem_d_ps_zsoil9[tid], elem_d_ps_zsoil10[tid],
				dt, 1, qtot,
				elem_d_soil_alpha,
				elem_d_soil_beta,
				elem_d_soil_smcmin,
				elem_d_soil_smcmax
				);
			rhsts[1] = rhsts[1] - tsnsr / denom;
		}
	}
	else
	{
		if ((sice > 0.0) || (elem_d_es_stc1[tid] < TFREEZ))
		{
			SnkSrc(tid, &tsnsr, elem_d_es_stc1[tid], elem_d_ws_smc1[tid], &elem_d_ws_sh2o1[tid],
				elem_d_ps_zsoil0[tid], elem_d_ps_zsoil1[tid], elem_d_ps_zsoil2[tid],
				elem_d_ps_zsoil3[tid], elem_d_ps_zsoil4[tid], elem_d_ps_zsoil5[tid],
				elem_d_ps_zsoil6[tid], elem_d_ps_zsoil7[tid], elem_d_ps_zsoil8[tid],
				elem_d_ps_zsoil9[tid], elem_d_ps_zsoil10[tid],
				dt, 1, qtot,
				elem_d_soil_alpha,
				elem_d_soil_beta,
				elem_d_soil_smcmin,
				elem_d_soil_smcmax
				);
			rhsts[1] = rhsts[1] - tsnsr / denom;
		}
	}

	/* Calc matrix coefs, ai, and bi for this layer. */
	ai[1] = -df1k * ddz / ((elem_d_ps_zsoil0[tid] - elem_d_ps_zsoil1[tid]) * hcpct);
	bi[1] = -(ai[1] + ci[1]);

	/* Reset values of df1, dtsdz, ddz, and tbk for loop to next soil layer.
	*/
	tbk = tbk1;
	df1k = df1n;
	dtsdz = dtsdz2;
	ddz = ddz2;
//------------ end of (k==1) --------------------------

// k==2









}








// 计算/更新 土壤层的温度分布
__device__ void HStep(int tid, realtype *rhsts, realtype dt, realtype *ai,
	realtype *bi, realtype *ci, 
	realtype *elem_d_es_stc0, realtype *elem_d_es_stc1, realtype *elem_d_es_stc2,
	realtype *elem_d_es_stc3, realtype *elem_d_es_stc4, realtype *elem_d_es_stc5,
	realtype *elem_d_es_stc6, realtype *elem_d_es_stc7, realtype *elem_d_es_stc8,
	realtype *elem_d_es_stc9, realtype *elem_d_es_stc10)	
{
	/*
	* Subroutine HStep
	*
	* Calculate/update the soil temperature field.
	*/
	int             k;
	realtype        rhstsin[MAXLYR];
	realtype        ciin[MAXLYR];

	/* Create finite difference values for use in Rosr12 routine */
#pragma unroll
	for (k = 0; k < MAXLYR; k++)
	{
		rhsts[k] *= dt;
		ai[k] *= dt;
		bi[k] = 1.0 + bi[k] * dt;
		ci[k] *= dt;
	}

	/* Copy values for input variables before call to Rosr12 */
#pragma unroll
	for (k = 0; k < MAXLYR; k++)
	{
		rhstsin[k] = rhsts[k];
	}

#pragma unroll
	for (k = 0; k <MAXLYR; k++)
	{
		ciin[k] = ci[k];
	}

	/* Solve the tri-diagonal matrix equation */
	Rosr12(ci, ai, bi, ciin, rhstsin, rhsts);

	/* Calc/update the soil temps using matrix solution */
	/*
	for (k = 0; k < MAXLYR; k++)
	{
		elem_d_es_stc[k][tid] += ci[k];
	}
	*/
	elem_d_es_stc0[tid] += ci[0];
	elem_d_es_stc1[tid] += ci[1];
	elem_d_es_stc2[tid] += ci[2];
	elem_d_es_stc3[tid] += ci[3];
	elem_d_es_stc4[tid] += ci[4];
	elem_d_es_stc5[tid] += ci[5];
	elem_d_es_stc6[tid] += ci[6];
	elem_d_es_stc7[tid] += ci[7];
	elem_d_es_stc8[tid] += ci[8];
	elem_d_es_stc9[tid] += ci[9];
	elem_d_es_stc10[tid] += ci[10];
}


__device__ void Rosr12(realtype *p, const realtype *a, const realtype *b, realtype *c,
	const realtype *d, realtype *delta)
{
	/*
	* Function Rosr12
	*
	* Invert (solve) the tri-diagonal matrix problem shown below:
	* ###                                            ### ###  ###   ###  ###
	* #b[0], c[0],  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
	* #a[1], b[1], c[1],  0  ,  0  ,   . . .  ,    0   # #      #   #      #
	* # 0  , a[2], b[2], c[2],  0  ,   . . .  ,    0   # #      #   # d[2] #
	* # 0  ,  0  , a[3], b[3], c[3],   . . .  ,    0   # # p[3] #   # d[3] #
	* # 0  ,  0  ,  0  , a[4], b[4],   . . .  ,    0   # # p[4] #   # d[4] #
	* # .                                          .   # #  .   # = #   .  #
	* # .                                          .   # #  .   #   #   .  #
	* # .                                          .   # #  .   #   #   .  #
	* # 0  , . . . , 0 , a[m-2], b[m-2], c[m-2],   0   # #p(m-2)#   #d[m-2]#
	* # 0  , . . . , 0 ,   0   , a[m-1], b[m-1], c[m-1]# #p(m-1)#   #d[m-1]#
	* # 0  , . . . , 0 ,   0   ,   0   ,  a(m) ,  b[m] # # p(m) #   # d[m] #
	* ###                                            ### ###  ###   ###  ###
	*/
	int             k, kk;

	/* Initialize eqn coef c for the lowest soil layer */
	c[MAXLYR - 1] = 0.0;
	p[0] = -c[0] / b[0];

	/* Solve the coefs for the 1st soil layer */
	delta[0] = d[0] / b[0];

	/* Solve the coefs for soil layers 2 thru nsoil */
#pragma unroll
	for (k = 1; k < MAXLYR; k++)
	{
		p[k] = -c[k] * (1.0 / (b[k] + a[k] * p[k - 1]));
		delta[k] =
			(d[k] - a[k] * delta[k - 1]) * (1.0 / (b[k] + a[k] * p[k - 1]));
	}

	/* Set p to delta for lowest soil layer */
	p[MAXLYR - 1] = delta[MAXLYR - 1];

	/* Adjust p for soil layers 2 thru nsoil */
#pragma unroll
	for (k = 1; k < MAXLYR; k++)
	{
		kk = MAXLYR - k - 1;
		p[kk] = p[kk] * p[kk + 1] + delta[kk];
	}
}




#if defined(_CYCLES_)
void NoPac(const soil_struct *soil, const lc_struct *lc,
	const cstate_struct *cs, realtype dt, realtype t24, crop_struct crop[],
	pstate_struct *ps, wstate_struct *ws, wflux_struct *wf,
	estate_struct *es, eflux_struct *ef)
#else
__device__ void NoPac(int tid, realtype dt, realtype t24,
#include "NoPac_kernel_declarations.c"	
)
#endif
{
	/*
	* Function NoPac
	*
	* Calculate soil moisture and heat flux values and update soil moisture
	* content and soil heat content values for the case when no snow pack is
	* present.
	*/
	int               k;
	realtype          df1;
	realtype          yy;
	realtype          zz1;
	realtype          yynum;
	realtype          prcpf;

	prcpf = elem_d_wf_prcp[tid];

	/* Initialize dew */
	elem_d_wf_dew[tid] = 0.0;

	/* Initialize evap terms */
	elem_d_wf_edir[tid] = 0.0;
	elem_d_wf_ec[tid] = 0.0;

	/*
	for (k = 0; k < elem_d_ps_nsoil[tid]; k++)
		elem_d_wf_et[k][tid] = 0.0;
	*/
	elem_d_wf_et0[tid] = 0.0;
	elem_d_wf_et1[tid] = 0.0;
	elem_d_wf_et2[tid] = 0.0;
	elem_d_wf_et3[tid] = 0.0;
	elem_d_wf_et4[tid] = 0.0;
	elem_d_wf_et5[tid] = 0.0;
	elem_d_wf_et6[tid] = 0.0;
	elem_d_wf_et7[tid] = 0.0;
	elem_d_wf_et8[tid] = 0.0;
	elem_d_wf_et9[tid] = 0.0;
	elem_d_wf_et10[tid] = 0.0;

	elem_d_wf_ett[tid] = 0.0;

	if (elem_d_wf_etp[tid] > 0.0)
	{
#if defined(_CYCLES_)
		Evapo(soil, lc, ps, es, cs, dt, crop, ws, wf);
#else
		Evapo(tid, dt,
#include "Evapo_kernel_arguments.c"			
			 );
#endif

		elem_d_wf_eta[tid] = elem_d_wf_etns[tid];
	}
	else
	{
		/* If etp < 0, assume dew forms (transform etp into dew and reinitialize
		*  etp to zero). */
		elem_d_wf_dew[tid] = -elem_d_wf_etp[tid];

		/* Add dew amount to prcp */
		prcpf += elem_d_wf_dew[tid];
	}

	PcpDrp(tid, prcpf, dt,
#include "PcpDrp_kernel_arguments.c"			
		  );

	/* Based on etp and e values, determine beta */
	if (elem_d_wf_etp[tid] <= 0.0)
	{
		elem_d_ps_beta[tid] = 0.0;
		elem_d_wf_eta[tid] = elem_d_wf_etp[tid];
		if (elem_d_wf_etp[tid] < 0.0)
		{
			elem_d_ps_beta[tid] = 1.0;
		}
	}
	else
	{
		elem_d_ps_beta[tid] = elem_d_wf_eta[tid] / elem_d_wf_etp[tid];
	}

	/* Get soil thermal diffusivity/conductivity for top soil lyr, Calc.
	* adjusted top lyr soil temp and adjusted soil flux, then call ShFlx to
	* compute/update soil heat flux and soil temps. */
	df1 = TDfCnd(elem_d_ws_smc0[tid], elem_d_soil_quartz[tid], 
		         elem_d_soil_smcmax[tid], elem_d_soil_smcmin[tid],
		         elem_d_ws_sh2o0[tid]);

	/* Urban */
	df1 = (elem_d_lc_isurban[tid]) ? 3.24 : df1;

	/* Vegetation greenness fraction reduction in subsurface heat flux via
	* reduction factor, which is convenient to apply here to thermal
	* diffusivity that is later used in HRT to compute sub sfc heat flux (see
	* additional comments on veg effect sub-sfc heat flx in function SFlx) */
	df1 *= exp(elem_d_ps_sbeta[tid] * elem_d_lc_shdfac[tid]);

	/* Compute intermediate terms passed to routine HRT (via routine ShFlx
	* below) for use in computing subsurface heat flux in HRT */
	yynum = elem_d_ef_fdown[tid] - elem_d_ps_emissi[tid] * SIGMA * t24;
	yy = elem_d_es_sfctmp[tid] +
		(yynum / elem_d_ps_rch[tid] + elem_d_es_th2[tid] - elem_d_es_sfctmp[tid] -
		 elem_d_ps_beta[tid] * elem_d_ps_epsca[tid]) /elem_d_ps_rr[tid];

	zz1 = df1 / (-0.5 * elem_d_ps_zsoil0[tid] * elem_d_ps_rch[tid] * elem_d_ps_rr[tid]) + 1.0;

	ShFlx(tid, dt, yy, zz1, df1,
#include "ShFlx_kernel_arguments.c"		
		 );

	/* Set flx1 and flx3 (snopack phase change heat fluxes) to zero since they
	* are not used here in SnoPac. flx2 (freezing rain heat flux) was similarly
	* initialized in the Penman routine. */
	elem_d_ef_flx1[tid] = CPH2O * elem_d_wf_prcp[tid] * 1000.0 * 
		               (elem_d_es_t1[tid] - elem_d_es_sfctmp[tid]);
	elem_d_ef_flx3[tid] = 0.0;
}

__constant__ realtype    KD = 6.54e-7;
__constant__ realtype    BFACTR = 3.89;
__device__ void PcpDrp(int tid, realtype prcp, realtype dt,
#include "PcpDrp_kernel_declarations.c"	
)
{
	/*
	* Function PcpDrp
	*
	* Separated from Noah SmFlx function
	* The canopy moisture content (cmc) is updated.
	*/
	realtype          excess;
	realtype          rhsct;
	realtype          trhsct;

	/* Compute the right hand side of the canopy eqn term (rhsct)
	* Convert rhsct (a rate) to trhsct (an amount) and add it to existing cmc.
	* If resulting amt exceeds max capacity, it becomes drip and will fall to
	* the grnd. */
	rhsct = elem_d_lc_shdfac[tid] * prcp - elem_d_wf_ec[tid];
	elem_d_wf_drip[tid] = 0.0;
	trhsct = dt * rhsct;
	excess = elem_d_ws_cmc[tid] + trhsct;

	/* Pcpdrp is the combined prcp and drip (from cmc) that goes into the soil
	*/
	if (excess > 0.0)
	{
		/* PIHM drip calculation following Rutter and Mortan (1977 JAE) */
		if (excess >= elem_d_ws_cmcmax[tid])
		{
			elem_d_wf_drip[tid] =
				(KD * elem_d_ws_cmcmax[tid] * exp(BFACTR)) + (excess - elem_d_ws_cmcmax[tid]) / dt;
		}
		else
		{
			elem_d_wf_drip[tid] = (KD * elem_d_ws_cmcmax[tid] * exp(BFACTR * excess / elem_d_ws_cmcmax[tid]));
		}

		elem_d_wf_drip[tid] = (elem_d_wf_drip[tid] > excess / dt) ? excess / dt : elem_d_wf_drip[tid];
	}

	rhsct -= elem_d_wf_drip[tid];

	elem_d_wf_pcpdrp[tid] = (1.0 - elem_d_lc_shdfac[tid]) * prcp + elem_d_wf_drip[tid];

	/* Update canopy water content/interception (cmc). Convert rhsct to an
	* "amount" value and add to previous cmc value to get new cmc. */
	elem_d_ws_cmc[tid] += dt * rhsct;
	if (elem_d_ws_cmc[tid] < 1.0e-20)
	{
		elem_d_ws_cmc[tid] = 0.0;
	}

	elem_d_ws_cmc[tid] = (elem_d_ws_cmc[tid] < elem_d_ws_cmcmax[tid]) ? elem_d_ws_cmc[tid] : elem_d_ws_cmcmax[tid];
}

__constant__ realtype    ELCP = 2.4888e+3;
__constant__ realtype    LSUBC = 2.501000e+6;
__device__ void Penman(int tid,  realtype *t24, realtype t2v, int snowng, int frzgra,
#include "Penman_kernel_declarations.c"	
)
{
	/*
	* Function Penman
	*
	* Calculate potential evaporation for the current point. Various partial
	* sums/products are also calculated and passed back to the calling routine
	* for later use.
	*/
	realtype          a;
	realtype          delta;
	realtype          fnet;
	realtype          rad;
	realtype          rho;
	realtype          emissi;
	realtype          elcp1;
	realtype          lvs;

	/* Prepare partial quantities for Penman equation. */
	emissi = elem_d_ps_emissi[tid];
	elcp1 = (1.0 - elem_d_ps_sncovr[tid]) * ELCP + elem_d_ps_sncovr[tid] * ELCP * LSUBS / LSUBC;
	lvs = (1.0 - elem_d_ps_sncovr[tid]) * LSUBC + elem_d_ps_sncovr[tid] * LSUBS;

	elem_d_ef_flx2[tid] = 0.0;
	delta = elcp1 * elem_d_ps_dqsdt2[tid];
	*t24 = elem_d_es_sfctmp[tid] * elem_d_es_sfctmp[tid] * elem_d_es_sfctmp[tid] * elem_d_es_sfctmp[tid];
	elem_d_ps_rr[tid] = emissi * *t24 * 6.48e-8 / (elem_d_ps_sfcprs[tid] * elem_d_ps_ch[tid]) + 1.0;
	rho = elem_d_ps_sfcprs[tid] / (RD * t2v);

	elem_d_ps_rch[tid] = rho * CP * elem_d_ps_ch[tid];

	/* Adjust the partial sums/products with the latent heat effects caused by
	* falling precipitation. */
	if (!snowng)
	{
		if (elem_d_wf_prcp[tid] > 0.0)
		{
			elem_d_ps_rr[tid] += CPH2O * elem_d_wf_prcp[tid] * 1000.0 / elem_d_ps_rch[tid];
		}
	}
	else
	{
		elem_d_ps_rr[tid] += CPICE * elem_d_wf_prcp[tid] * 1000.0 / elem_d_ps_rch[tid];
	}

	/* Include the latent heat effects of frzng rain converting to ice on impact
	* in the calculation of flx2 and fnet. */
	fnet = elem_d_ef_fdown[tid] - emissi * SIGMA * *t24 - elem_d_ef_ssoil[tid];
	if (frzgra)
	{
		elem_d_ef_flx2[tid] = -LSUBF * elem_d_wf_prcp[tid] * 1000.0;
		fnet -= elem_d_ef_flx2[tid];
	}

	/* Finish Penman equation calculations */
	rad = fnet / elem_d_ps_rch[tid] + elem_d_es_th2[tid] - elem_d_es_sfctmp[tid];
	a = elcp1 * (elem_d_ps_q2sat[tid] - elem_d_ps_q2[tid]);
	elem_d_ps_epsca[tid] = (a * elem_d_ps_rr[tid] + rad * delta) / (delta + elem_d_ps_rr[tid]);
	elem_d_wf_etp[tid] = elem_d_ps_epsca[tid] * elem_d_ps_rch[tid] / lvs / 1000.0;
}

__constant__ realtype    WWST = 1.2;
__constant__ realtype    VKRM = 0.40;
__constant__ realtype    EXCM = 0.001;
__constant__ realtype    BETA = 1.0 / 270.0;
__constant__ realtype    WOLD = 0.15;
__constant__ int       ITRMX = 5;
__constant__ realtype    EPSU2 = 1.0e-4;
__constant__ realtype    EPSUST = 0.07;
__constant__ realtype    ZTMIN = -5.0;
__constant__ realtype    ZTMAX = 1.0;
__constant__ realtype    HPBL = 1000.0;
__constant__ realtype    SQVISC = 258.2;
__device__ void SfcDifOff(int tid, realtype t1v, realtype th2v, int iz0tlnd,
#include "SfcDifOff_kernel_declarations.c"	
)
{
	/*
	* Calculate surface layer exchange coefficients via iterative process.
	* See Chen et al. (1997, BLM)
	*
	* This routine sfcdif can handle both over open water (sea, ocean) and
	* over solid surface (land, sea-ice).
	*/
	realtype          zilfc;
	realtype          zu;
	realtype          zt;
	realtype          rdz;
	realtype          cxch;
	realtype          dthv;
	realtype          du2;
	realtype          btgh;
	realtype          wstar2;
	realtype          ustar;
	realtype          zslu;
	realtype          zslt;
	realtype          rlogu;
	realtype          rlogt;
	realtype          rlmo;
	realtype          zetalt;
	realtype          zetalu;
	realtype          zetau;
	realtype          zetat;
	realtype          xlu4;
	realtype          xlt4;
	realtype          xu4;
	realtype          xt4;
	realtype          xlu;
	realtype          xlt;
	realtype          xu;
	realtype          xt;
	realtype          psmz;
	realtype          simm;
	realtype          pshz;
	realtype          simh;
	realtype          ustark;
	realtype          rlmn;
	realtype          rlma;
	int             ilech;
	int             itr;
	realtype          wwst2;
	realtype          btg;
	realtype          elfc;
	realtype          wnew;

	wwst2 = WWST * WWST;
	btg = BETA * GRAV;
	elfc = VKRM * btg;
	wnew = 1.0 - WOLD;

	/* czil: constant C in Zilitinkevich, S. S.1995 */

	ilech = 0;

	if (iz0tlnd == 0 || elem_d_lc_isurban[tid])
	{
		/* Just use the original Czil value. */
		zilfc = -elem_d_ps_czil[tid] * VKRM * SQVISC;
	}
	else
	{
		/* Modify czil according to Chen & Zhang, 2009
		* czil = 10 ** -0.40 h, ( where h = 10*zo ) */
		elem_d_ps_czil[tid] = pow(10.0, -0.4 * (elem_d_ps_z0[tid] / 0.07));
		zilfc = -elem_d_ps_czil[tid] * VKRM * SQVISC;
	}

	zu = elem_d_ps_z0[tid];
	rdz = 1.0 / elem_d_ps_zlvl_wind[tid];
	cxch = EXCM * rdz;
	dthv = th2v - t1v;

	/* Beljars correction of ustar */
	du2 = elem_d_ps_sfcspd[tid] * elem_d_ps_sfcspd[tid];
	du2 = (du2 > EPSU2) ? du2 : EPSU2;

	btgh = btg * HPBL;
	/* If statements to avoid tangent linear problems near zero */
	if (btgh * elem_d_ps_ch[tid] * dthv != 0.0)
	{
		wstar2 = wwst2 * pow(fabs(btgh * elem_d_ps_ch[tid] * dthv), 2.0 / 3.0);
	}
	else
	{
		wstar2 = 0.0;
	}

	ustar = sqrt(elem_d_ps_cm[tid] * sqrt(du2 + wstar2));
	ustar = (ustar > EPSUST) ? ustar : EPSUST;

	/* Zilitinkevitch approach for zt */
	zt = exp(zilfc * sqrt(ustar * elem_d_ps_z0[tid])) * elem_d_ps_z0[tid];
	zslu = elem_d_ps_zlvl_wind[tid] + zu;

	zslt = elem_d_ps_zlvl[tid] + zt;

	rlogu = log(zslu / zu);

	rlogt = log(zslt / zt);

	rlmo = elfc * elem_d_ps_ch[tid] * dthv / pow(ustar, 3);

	for (itr = 0; itr < ITRMX; itr++)
	{
		zetalt = zslt * rlmo;
		zetalt = zetalt > ZTMIN ? zetalt : ZTMIN;
		rlmo = zetalt / zslt;
		zetalu = zslu * rlmo;
		zetau = zu * rlmo;

		zetat = zt * rlmo;

		/* 1. Monin-Obukkhov length-scale */
		if (ilech == 0)
		{
			if (rlmo < 0.0)
			{
				xlu4 = 1.0 - 16.0 * zetalu;
				xlt4 = 1.0 - 16.0 * zetalt;
				xu4 = 1.0 - 16.0 * zetau;
				xt4 = 1.0 - 16.0 * zetat;

				xlu = sqrt(sqrt(xlu4));
				xlt = sqrt(sqrt(xlt4));
				xu = sqrt(sqrt(xu4));
				xt = sqrt(sqrt(xt4));

				psmz = Pspmu(xu);
				simm = Pspmu(xlu) - psmz + rlogu;
				pshz = Psphu(xt);
				simh = Psphu(xlt) - pshz + rlogt;
			}
			else
			{
				zetalu = (zetalu < ZTMAX) ? zetalu : ZTMAX;
				zetalt = (zetalt < ZTMAX) ? zetalt : ZTMAX;

				psmz = Pspms(zetau);
				simm = Pspms(zetalu) - psmz + rlogu;
				pshz = Psphs(zetat);
				simh = Psphs(zetalt) - pshz + rlogt;
			}
		}
		else
		{
			/* Lech's functions */
			if (rlmo < 0.0)
			{
				psmz = Pslmu(zetau);
				simm = Pslmu(zetalu) - psmz + rlogu;
				pshz = Pslhu(zetat);
				simh = Pslhu(zetalt) - pshz + rlogt;
			}
			else
			{
				zetalu = zetalu < ZTMAX ? zetalu : ZTMAX;
				zetalt = zetalt < ZTMAX ? zetalt : ZTMAX;

				psmz = Pslms(zetau);
				simm = Pslms(zetalu) - psmz + rlogu;
				pshz = Pslhs(zetat);
				simh = Pslhs(zetalt) - pshz + rlogt;
			}
		}

		/* Beljaars correction for ustar */
		ustar = sqrt(elem_d_ps_cm[tid] * sqrt(du2 + wstar2));
		ustar = (ustar > EPSUST) ? ustar : EPSUST;

		/* Zilitinkevitch fix for zt */
		zt = exp(zilfc * sqrt(ustar * elem_d_ps_z0[tid])) * elem_d_ps_z0[tid];
		zslt = elem_d_ps_zlvl[tid] + zt;

		rlogt = log(zslt / zt);
		ustark = ustar * VKRM;

		elem_d_ps_cm[tid] = (ustark / simm > cxch) ? ustark / simm : cxch;
		elem_d_ps_ch[tid] = (ustark / simh > cxch) ? ustark / simh : cxch;

		/* If statements to avoid tangent linear problems near zero */
		if (btgh * elem_d_ps_ch[tid] * dthv != 0.0)
		{
			wstar2 = wwst2 * pow(fabs(btgh * elem_d_ps_ch[tid] * dthv), 2.0 / 3.0);
		}
		else
		{
			wstar2 = 0.0;
		}

		rlmn = elfc * elem_d_ps_ch[tid] * dthv / pow(ustar, 3.0);

		rlma = rlmo * WOLD + rlmn * wnew;

		rlmo = rlma;
	}
}


__device__ realtype SnFrac(realtype sneqv, realtype snup, realtype salp)
{
	/*
	* Function SnFrac
	*
	* Calculate snow fraction (0 -> 1)
	*/
	realtype          rsnow;
	realtype          sncovr;

	/* Snup is veg-class dependent snowdepth threshold above which snocvr = 1.
	*/
	if (sneqv < snup)
	{
		rsnow = sneqv / snup;
		sncovr = 1.0 - (exp(-salp * rsnow) - rsnow * exp(-salp));
	}
	else
	{
		sncovr = 1.0;
	}

	/* Formulation of Dickinson et al. 1986 */
	//z0n = 0.035;

	//sncovr = snowh / (snowh + 5.0 * z0n);

	/* Formulation of Marshall et al. 1994 */
	//sncovr = sneqv / (sneqv + 2.0 * z0n);

	return sncovr;
}

__device__ void ShFlx(int tid,
	            realtype dt, realtype yy, realtype zz1, realtype df1,
#include "ShFlx_kernel_declarations.c"
)
{
	/*
	* Function ShFlx
	*
	* Update the temperature state of the soil column based on the thermal
	* diffusion equation and update the frozen soil moisture content based on
	* the temperature.
	*/
	realtype          ai[MAXLYR], bi[MAXLYR], ci[MAXLYR];
	realtype          rhsts[MAXLYR];

	/* HRT routine calcs the right hand side of the soil temp dif eqn */
	/* Land case */
	HRT(tid, rhsts, yy, zz1, dt, df1, ai, bi, ci,
#include "HRT_kernel_arguments.c"    // 1D数组  CUDA
		);

	// elem_d_ps_nsoil[tid]不确定, elem_d_es_stc是二级指针形参  2021.06.08
	HStep(tid, rhsts, dt, ai, bi, ci, 
		elem_d_es_stc0, elem_d_es_stc1, elem_d_es_stc2, elem_d_es_stc3,
		elem_d_es_stc4, elem_d_es_stc5, elem_d_es_stc6, elem_d_es_stc7,
		elem_d_es_stc8, elem_d_es_stc9, elem_d_es_stc10	);

	/* In the no snowpack case (via routine NoPac branch,) update the grnd
	* (skin) temperature here in response to the updated soil temperature
	* profile above. (Note: inspection of routine SnoPac shows that t1 below is
	* a dummy variable only, as skin temperature is updated differently in
	* routine SnoPac) */
	elem_d_es_t1[tid] = (yy + (zz1 - 1.0) * elem_d_es_stc[0][tid]) / zz1;

	/* Calculate surface soil heat flux */
	elem_d_ef_ssoil[tid] = df1 * (elem_d_es_stc0[tid] - elem_d_es_t1[tid]) / 
		                  (0.5 * elem_d_ps_zsoil0[tid]);
}

__device__ void SnkSrc(int tid, realtype *tsnsr, realtype tavg, realtype smc, realtype *sh2o,
	realtype zsoil0, realtype zsoil1, realtype zsoil2, realtype zsoil3, realtype zsoil4,
	realtype zsoil5, realtype zsoil6, realtype zsoil7, realtype zsoil8, realtype zsoil9,realtype zsoil10,
	realtype dt, int k, realtype qtot,
	realtype *elem_d_soil_alpha,
	realtype *elem_d_soil_beta,
	realtype *elem_d_soil_smcmin,
	realtype *elem_d_soil_smcmax
)
{
	/*
	* Subroutine SnkSrc
	*
	* Calculate sink/source term of the thermal diffusion equation. (sh2o) is
	* available liquid water.
	*/
	realtype          dz;
	realtype          freew;
	realtype          xh2o;

	realtype          DH2O = 1.0000e3;

	if (0 == k)
	{
		dz = -zsoil0;
	}
	else if (1 == k){

		dz = zsoil0 - zsoil1;
	}
	else if (2 == k){

		dz = zsoil1 - zsoil2;
	}
	else if (3 == k){

		dz = zsoil2 - zsoil3;
	}
	else if (4 == k){

		dz = zsoil3 - zsoil4;
	}
	else if (5 == k){

		dz = zsoil4 - zsoil5;
	}
	else if (6 == k){

		dz = zsoil5 - zsoil6;
	}
	else if (7 == k){

		dz = zsoil6 - zsoil7;
	}
	else if (8 == k){

		dz = zsoil7 - zsoil8;
	}
	else if (9 == k){

		dz = zsoil8 - zsoil9;
	}
	else if (10 == k){

		dz = zsoil9 - zsoil10;
	}
	else
	{
	//	dz = zsoil[k - 1] - zsoil[k];
	}

	/* Via function FrH2O, compute potential or 'equilibrium' unfrozen
	* supercooled free water for given soil type and soil layer temperature.
	* Function FrH2O invokes Eqn (17) from V. Koren et al (1999, JGR, Vol. 104,
	* Pg 19573). (Aside: latter eqn in journal in centigrade units.
	* routine FrH2O use form of eqn in kelvin units.) */
	freew = FrH2O(tid, tavg, smc, *sh2o, 
		elem_d_soil_alpha,
		elem_d_soil_beta,
		elem_d_soil_smcmin,
		elem_d_soil_smcmax
		);

	/* In next block of code, invoke Eqn 18 of V. Koren et al (1999, JGR,
	* Vol. 104, Pg 19573.) that is, first estimate the new amount of liquid
	* water, 'xh2o', implied by the sum of (1) the liquid water at the begin of
	* current time step, and (2) the freeze of thaw change in liquid water
	* implied by the heat flux 'qtot' passed in from routine HRT. Second,
	* determine if xh2o needs to be bounded by 'freew' (equil amt) or if
	* 'freew' needs to be bounded by xh2o. */
	xh2o = *sh2o + qtot * dt / (DH2O * LSUBF * dz);

	/* First, if freezing and remaining liquid less than lower bound, then
	* reduce extent of freezing, thereby letting some or all of heat flux qtot
	* cool the soil temp later in routine HRT. */
	if (xh2o < *sh2o && xh2o < freew)
	{
		xh2o = (freew > *sh2o) ? *sh2o : freew;
	}

	/* Second, if thawing and the increase in liquid water greater than upper
	* bound, then reduce extent of thaw, thereby letting some or all of heat
	* flux qtot warm the soil temp later in routine HRT. */
	if (xh2o > *sh2o && xh2o > freew)
	{
		xh2o = (freew < *sh2o) ? *sh2o : freew;
	}

	xh2o = (xh2o < 0.0) ? 0.0 : xh2o;
	xh2o = (xh2o > smc) ? smc : xh2o;

	/* Calculate phase-change heat source/sink term for use in routine HRT and
	* update liquid water to reflect final freeze/thaw increment. */
	*tsnsr = -DH2O * LSUBF * dz * (xh2o - *sh2o) / dt;

	*sh2o = xh2o;
}


#if defined(_CYCLES_)
void SnoPac(const soil_struct *soil, const lc_struct *lc,
	const cstate_struct *cs, int snowng, double dt, double t24, double prcpf,
	double df1, crop_struct crop[], pstate_struct *ps, wstate_struct *ws,
	wflux_struct *wf, estate_struct *es, eflux_struct *ef)
#else
__device__ void SnoPac(int tid, int snowng, double dt, 
	                   double t24, double prcpf, double df1,
#include "SnoPac_kernel_declarations.c"
)
#endif
{
	/*
	* Function SnoPac
	*
	* Calculate soil moisture and heat flux values & update soil moisture
	* content and soil heat content values for the case when a snow pack is
	* present.
	*/
	int             k;

	double          denom;
	double          dsoil;
	double          dtot;
	double          esnow1;
	double          esnow2;
	double          etp3;
	double          etanrg;
	double          ex;
	double          seh;
	double          sncond;
	double          t12;
	double          t12a;
	double          t12b;
	double          t14;
	double          ssoil1;
	double          t11;
	double          yy;
	double          zz1;
	const double    ESDMIN = 1.0e-6;
	const double    SNOEXP = 2.0;

	/* Initialize evap terms. */
	elem_d_wf_dew[tid] = 0.0;
	elem_d_wf_edir[tid] = 0.0;
	elem_d_wf_ec[tid] = 0.0;
	/*
	for (k = 0; k < elem_d_ps_nsoil[tid]; k++)
	{
		elem_d_wf_et[k][tid] = 0.0;
	}
	*/
	elem_d_wf_et0[tid] = 0.0;
	elem_d_wf_et1[tid] = 0.0;
	elem_d_wf_et2[tid] = 0.0;
	elem_d_wf_et3[tid] = 0.0;
	elem_d_wf_et4[tid] = 0.0;
	elem_d_wf_et5[tid] = 0.0;
	elem_d_wf_et6[tid] = 0.0;
	elem_d_wf_et7[tid] = 0.0;
	elem_d_wf_et8[tid] = 0.0;
	elem_d_wf_et9[tid] = 0.0;
	elem_d_wf_et10[tid] = 0.0;

	elem_d_wf_ett[tid] = 0.0;
	elem_d_wf_esnow[tid] = 0.0;
	esnow1 = 0.0;
	esnow2 = 0.0;

	elem_d_ps_beta[tid] = 1.0;

	/* If etp < 0 (downward) then dewfall (= frostfall in this case). */
	if (elem_d_wf_etp[tid] <= 0.0)
	{
		if ((elem_d_ps_ribb[tid] >= 0.1) && (elem_d_ef_fdown[tid] > 150.0))
		{
			elem_d_wf_etp[tid] =
				(((elem_d_wf_etp[tid] * (1.0 - elem_d_ps_ribb[tid]) < 0.0) ?
				elem_d_wf_etp[tid] * (1.0 - elem_d_ps_ribb[tid]) : 0.0)
				* elem_d_ps_sncovr[tid] / 0.980 + elem_d_wf_etp[tid] * (0.980 - elem_d_ps_sncovr[tid])) / 0.980;
		}

		if (elem_d_wf_etp[tid] == 0.0)
		{
			elem_d_ps_beta[tid] = 0.0;
		}

		elem_d_wf_dew[tid] = -elem_d_wf_etp[tid];
		esnow2 = elem_d_wf_etp[tid] * dt;
		etanrg = elem_d_wf_etp[tid] * 1000.0 *
			((1.0 - elem_d_ps_sncovr[tid]) * LVH2O + elem_d_ps_sncovr[tid] * LSUBS);
	}
	else
	{
		/* Land case */
		if (elem_d_ps_sncovr[tid] < 1.0)
		{
#if defined(_CYCLES_)
			Evapo(soil, lc, ps, es, cs, dt, crop, ws, wf);
#else
			Evapo(tid, dt,
#include "Evapo_kernel_arguments.c"				
			   	  );
#endif
			elem_d_wf_edir[tid] *= (1.0 - elem_d_ps_sncovr[tid]);
			elem_d_wf_ec[tid] *= (1.0 - elem_d_ps_sncovr[tid]);
			/*
			for (k = 0; k < elem_d_ps_nsoil[tid]; k++)
			{
				elem_d_wf_et[k][tid] = elem_d_wf_et[k][tid] * (1.0 - elem_d_ps_sncovr[tid]);
			}
			*/
			realtype		temp;
			temp = 1.0 - elem_d_ps_sncovr[tid];
			elem_d_wf_et0[tid] = elem_d_wf_et0[tid] * temp;
			elem_d_wf_et1[tid] = elem_d_wf_et1[tid] * temp;
			elem_d_wf_et2[tid] = elem_d_wf_et2[tid] * temp;
			elem_d_wf_et3[tid] = elem_d_wf_et3[tid] * temp;
			elem_d_wf_et4[tid] = elem_d_wf_et4[tid] * temp;
			elem_d_wf_et5[tid] = elem_d_wf_et5[tid] * temp;
			elem_d_wf_et6[tid] = elem_d_wf_et6[tid] * temp;
			elem_d_wf_et7[tid] = elem_d_wf_et7[tid] * temp;
			elem_d_wf_et8[tid] = elem_d_wf_et8[tid] * temp;
			elem_d_wf_et9[tid] = elem_d_wf_et9[tid] * temp;
			elem_d_wf_et10[tid] = elem_d_wf_et10[tid] * temp;

			elem_d_wf_ett[tid] *= (1.0 - elem_d_ps_sncovr[tid]);
			elem_d_wf_etns[tid] *= (1.0 - elem_d_ps_sncovr[tid]);
		}
		elem_d_wf_esnow[tid] = elem_d_wf_etp[tid] * elem_d_ps_sncovr[tid];
		esnow2 = elem_d_wf_esnow[tid] * dt;
		etanrg = elem_d_wf_esnow[tid] * 1000.0 * LSUBS + elem_d_wf_etns[tid] * 1000.0 * LVH2O;
	}

	/* If precip is falling, calculate heat flux from snow sfc to newly
	* accumulating precip.  note that this reflects the flux appropriate for
	* the not-yet-updated skin temperature (t1).  assumes temperature of the
	* snowfall striking the ground is = sfctmp (lowest model level air temp) */
	elem_d_ef_flx1[tid] = 0.0;
	if (snowng)
	{
		elem_d_ef_flx1[tid] = CPICE * elem_d_wf_prcp[tid] * 1000.0 * (elem_d_es_t1[tid] - elem_d_es_sfctmp[tid]);
	}
	else
	{
		if (elem_d_wf_prcp[tid] > 0.0)
		{
			elem_d_ef_flx1[tid] = CPH2O * elem_d_wf_prcp[tid] * 1000.0 * 
				                 (elem_d_es_t1[tid] - elem_d_es_sfctmp[tid]);
		}
	}

	dsoil = -0.5 * elem_d_ps_zsoil0[tid];
	dtot = elem_d_ps_snowh[tid] + dsoil;
	denom = 1.0 + df1 / (dtot * elem_d_ps_rr[tid] * elem_d_ps_rch[tid]);

	/* Calculate an 'effective snow-grnd sfc temp' (t12) based on heat fluxes
	* between the snow pack and the soil and on net radiation.
	* Include flx1 (precip-snow sfc) and flx2 (freezing rain latent heat)
	* fluxes. flx1 from above, flx2 brought in via common block rite.
	* flx2 reflects freezing rain latent heat flux using t1 calculated in
	* Penman. */
	/* Surface emissivity weighted by snow cover fraction */
	t12a = ((elem_d_ef_fdown[tid] - elem_d_ef_flx1[tid] - elem_d_ef_flx2[tid] -
		     elem_d_ps_emissi[tid] * SIGMA * t24) /
			 elem_d_ps_rch[tid] + elem_d_es_th2[tid] - elem_d_es_sfctmp[tid] -
			  etanrg / elem_d_ps_rch[tid]) / elem_d_ps_rr[tid];
	t12b = df1 * elem_d_es_stc0[tid] / (dtot * elem_d_ps_rr[tid] * elem_d_ps_rch[tid]);

	t12 = (elem_d_es_sfctmp[tid] + t12a + t12b) / denom;

	/* If the 'effective snow-grnd sfc temp' is at or below freezing, no snow
	* melt will occur. Set the skin temp to this effective temp. Reduce (by
	* sublimation) or increase (by frost) the depth of the snowpack, depending
	* on sign of etp.
	* Update soil heat flux (ssoil) using new skin temperature (t1) since no
	* snowmelt, set accumulated snowmelt to zero, set 'effective' precip from
	* snowmelt to zero, set phase-change heat flux from snowmelt to zero. */
	if (t12 <= TFREEZ)
	{
		/* Sub-freezing block */
		elem_d_es_t1[tid] = t12;
		elem_d_ef_ssoil[tid] = df1 * (elem_d_es_t1[tid] - elem_d_es_stc0[tid]) / dtot;

		if (elem_d_ws_sneqv[tid] - esnow2 > 0.0)
		{
			elem_d_ws_sneqv[tid] -= esnow2;
		}
		else
		{
			elem_d_ws_sneqv[tid] = 0.0;
			esnow2 = elem_d_ws_sneqv[tid];
			elem_d_wf_esnow[tid] = esnow2 / dt;
		}

		elem_d_ef_flx3[tid] = 0.0;
		ex = 0.0;

		elem_d_wf_snomlt[tid] = 0.0;
	}
	else
	{
		/* If the 'effective snow-grnd sfc temp' is above freezing, snow melt
		* will occur. Call the snow melt rate, ex and amt, snomlt. Revise the
		* effective snow depth. Revise the skin temp because it would have chgd
		* due to the latent heat released by the melting. Calc the latent heat
		* released, flx3.
		* Set the effective precip, prcp1 to the snow melt rate, ex for use in
		* SmFlx. Adjustment to t1 to account for snow patches.
		* Calculate qsat valid at freezing point. Note that esat (saturation
		* vapor pressure) value of 6.11e+2 used here is that valid at freezing
		* point.
		* Note that etp from call Penman in sflx is ignored here in favor of
		* bulk etp over 'open water' at freezing temp. Update soil heat flux
		* (s) using new skin temperature (t1) */
		/* Above freezing block */
		elem_d_es_t1[tid] = TFREEZ * pow(elem_d_ps_sncovr[tid], SNOEXP) +
			t12 * (1.0 - pow(elem_d_ps_sncovr[tid], SNOEXP));
		elem_d_ps_beta[tid] = 1.0;

		elem_d_ef_ssoil[tid] = df1 * (elem_d_es_t1[tid] - elem_d_es_stc0[tid]) / dtot;

		if (elem_d_ws_sneqv[tid] - esnow2 <= ESDMIN)
		{
			/* If potential evap (sublimation) greater than depth of snowpack.
			* beta < 1
			* snowpack has sublimated away, set depth to zero. */
			elem_d_ws_sneqv[tid] = 0.0;
			ex = 0.0;
			elem_d_wf_snomlt[tid] = 0.0;
			elem_d_ef_flx3[tid] = 0.0;
		}
		else
		{
			/* Sublimation less than depth of snowpack
			* Snowpack (esd) reduced by esnow2 (depth of sublimated snow) */
			elem_d_ws_sneqv[tid] -= esnow2;
			etp3 = elem_d_wf_etp[tid] * 1000.0 * LVH2O;
			seh = elem_d_ps_rch[tid] * (elem_d_es_t1[tid] - elem_d_es_th2[tid]);
			t14 = elem_d_es_t1[tid] * elem_d_es_t1[tid];
			t14 = t14 * t14;
			elem_d_ef_flx3[tid] =
				elem_d_ef_fdown[tid] - elem_d_ef_flx1[tid] - elem_d_ef_flx2[tid] - 
				elem_d_ps_emissi[tid] * SIGMA * t14 -
				elem_d_ef_ssoil[tid] - seh - etanrg;
			elem_d_ef_flx3[tid] = (elem_d_ef_flx3[tid] <= 0.0) ? 0.0 : elem_d_ef_flx3[tid];

			ex = elem_d_ef_flx3[tid] * 0.001 / LSUBF;

			/* Snowmelt reduction depending on snow cover */
			elem_d_wf_snomlt[tid] = ex;

			/* ESDMIN represents a snowpack depth threshold value below which
			* we choose not to retain any snowpack, and instead include it in
			* snowmelt. */
			if (elem_d_ws_sneqv[tid] - elem_d_wf_snomlt[tid] * dt >= ESDMIN)
			{
				elem_d_ws_sneqv[tid] -= elem_d_wf_snomlt[tid] * dt;
			}
			else
			{
				/* Snowmelt exceeds snow depth */
				ex = elem_d_ws_sneqv[tid] / dt;
				elem_d_ef_flx3[tid] = ex * 1000.0 * LSUBF;
				elem_d_wf_snomlt[tid] = elem_d_ws_sneqv[tid] / dt;

				elem_d_ws_sneqv[tid] = 0.0;
			}
		}

		/* If non-glacial land, add snowmelt rate (ex) to precip rate to be used
		* in subroutine SmFlx (soil moisture evolution) via infiltration.
		* Runoff/baseflow later near the end of sflx (after return from call to
		* subroutine SnoPac) */
		prcpf += ex;

		/* Set the effective potnl evapotransp (etp1) to zero since this is snow
		* case, so surface evap not calculated from edir, ec, or ett in SmFlx
		* (below).
		* SmFlx returns updated soil moisture values for non-glacial land. */
	}

	PcpDrp(tid, prcpf, dt,
#include "PcpDrp_kernel_arguments.c"			
		   );

	/* Before call ShFlx in this snowpack case, set zz1 and yy arguments to
	* special values that ensure that ground heat flux calculated in ShFlx
	* matches that already computer for below the snowpack, thus the sfc heat
	* flux to be computed in ShFlx will effectively be the flux at the snow top
	* surface.  t11 is a dummy argument so we will not use the skin temp value
	* as revised by ShFlx. */
	zz1 = 1.0;
	yy = elem_d_es_stc0[tid] - 0.5 * elem_d_ef_ssoil[tid] * elem_d_ps_zsoil0[tid] * zz1 / df1;

	t11 = elem_d_es_t1[tid];
	ssoil1 = elem_d_ef_ssoil[tid];

	/* ShFlx will calc/update the soil temps. Note: the sub-sfc heat flux
	* (ssoil1) and the skin temp (t11) output from this ShFlx call are not used
	* in any subsequent calculations. Rather, they are dummy variables here in
	* the SnoPac case, since the skin temp and sub-sfc heat flux are updated
	* instead near the beginning of the call to SnoPac. */
	ShFlx(tid, dt, yy, zz1, df1,
#include "ShFlx_kernel_arguments.c"	
		  );

	elem_d_es_t1[tid] = t11;
	elem_d_ef_ssoil[tid] = ssoil1;

	/* Snow depth and density adjustment based on snow compaction. yy is assumed
	* to be the soil temperature at the top of the soil column. */
	if (elem_d_ws_sneqv[tid] > 0.0)
	{
		SnowPack(elem_d_ws_sneqv[tid], dt, &elem_d_ps_snowh[tid],
			    &elem_d_ps_sndens[tid], elem_d_es_t1[tid], yy);
	}
	else
	{
		elem_d_ws_sneqv[tid] = 0.0;
		elem_d_ps_snowh[tid] = 0.0;
		elem_d_ps_sndens[tid] = 0.0;
		sncond = 1.0;
		elem_d_ps_sncovr[tid] = 0.0;
	}
}

__device__ void SnowPack(double esd, double dtsec, double *snowh, double *sndens,
	double tsnow, double tsoil)
{
	/*
	* Function SnowPack
	*
	* Calculate compaction of snowpack under conditions of increasing snow
	* density, as obtained from an approximate solution of E. Anderson's
	* differential equation (3.29), NOAA technical report NWS 19, by Victor
	* Koren, 03/25/95.
	*
	* esd     water equivalent of snow (m)
	* dtsec   time step (sec)
	* snowh   snow depth (m)
	* sndens  snow density (g/cm3=dimensionless fraction of h2o density)
	* tsnow   snow surface temperature (k)
	* tsoil   soil surface temperature (k)
	* Function will return new values of snowh and sndens
	*/
	int             ipol;
	int             j;
	double          bfac;
	double          dsx;
	double          dthr;
	double          dw;
	double          snowhc;
	double          pexp;
	double          tavgc;
	double          tsnowc;
	double          tsoilc;
	double          esdc;
	double          esdcx;
	const double    C1 = 0.01;
	const double    C2 = 21.0;

	/* Conversion into simulation units */
	snowhc = *snowh * 100.0;
	esdc = esd * 100.0;
	dthr = dtsec / 3600.0;
	tsnowc = tsnow - 273.15;
	tsoilc = tsoil - 273.15;

	/* Calculating of average temperature of snow pack */
	tavgc = 0.5 * (tsnowc + tsoilc);

	esdcx = (esdc > 1.0e-2) ? esdc : 1.0e-2;

	/* Calculating of snow depth and density as a result of compaction
	*  sndens = ds0 * (exp (bfac * esd) - 1.) / (bfac * esd)
	*  bfac = dthr * c1 * exp (0.08 * tavgc - c2 * ds0)
	* Note: bfac * esd in sndens eqn above has to be carefully treated
	* numerically below:
	*   c1 is the fractional increase in density (1/(cm*hr))
	*   c2 is a constant (cm3/g) kojima estimated as 21 cms/g */
	bfac = dthr * C1 * exp(0.08 * tavgc - C2 * *sndens);

	/* The function of the form (e**x - 1) / x embedded in above expression for
	* dsx was causing numerical difficulties when the denominator "x"
	* (i.e. bfac*esdc) became zero or approached zero (despite the fact that
	* the analytical function (e**x-1)/x has a well defined limit as "x"
	* approaches zero), hence below we replace the (e**x-1)/x expression with
	* an equivalent, numerically well-behaved polynomial expansion.
	* Number of terms of polynomial expansion, and hence its accuracy, is
	* governed by iteration limit "ipol".
	*      ipol greater than 9 only makes a difference on double precision
	*      (relative errors given in percent %).
	*       ipol=9, for rel.error <~ 1.6 e-6 % (8 significant digits)
	*       ipol=8, for rel.error <~ 1.8 e-5 % (7 significant digits)
	*       ipol=7, for rel.error <~ 1.8 e-4 % ... */
	ipol = 4;
	pexp = 0.0;
	for (j = ipol; j > 0; j--)
	{
		pexp = (1.0 + pexp) * bfac * esdcx / (double)(j + 1);
	}

	pexp += 1.0;

	/* Set upper/lower limit on snow density */
	dsx = *sndens * pexp;
	dsx = (dsx > 0.40) ? 0.40 : dsx;
	dsx = (dsx < 0.05) ? 0.05 : dsx;

	/* Update of snow depth and density depending on liquid water during
	* snowmelt. Assumed that 13% of liquid water can be stored in snow per day
	* during snowmelt till snow density 0.40. */
	*sndens = dsx;
	if (tsnowc >= 0.0)
	{
		dw = 0.13 * dthr / 24.0;
		*sndens = *sndens * (1.0 - dw) + dw;
		if (*sndens >= 0.40)
		{
			*sndens = 0.40;
		}
	}
	/* Calculate snow depth (cm) from snow water equivalent and snow density.
	* Change snow depth units to meters */
	snowhc = esdc / *sndens;
	*snowh = snowhc * 0.01;
}

__constant__ double    Z0S = 0.001;    /* Snow roughness length (m) */
__device__ double Snowz0(double sncovr, double z0brd, double snowh)
{
	/*
	* Function Snowz0
	*
	* Calculate total roughness length over snow
	*/
	double          burial;
	double          z0eff;
	double          z0;

	burial = 7.0 * z0brd - snowh;
	z0eff = (burial <= 0.0007) ? Z0S : (burial / 7.0);
	z0 = (1.0 - sncovr) * z0brd + sncovr * z0eff;

	return z0;
}

__device__ void SnowNew(int tid, double newsn, 
	realtype *elem_d_ps_snowh,
	realtype *elem_d_es_sfctmp,
	realtype *elem_d_ps_sndens
)
{
	/*
	* Function SnowNew
	*
	* Calculate snow depth and density to account for the new snowfall.
	* New values of snow depth & density returned.
	*/
	double          dsnew;
	double          hnewc;
	double          snowhc;
	double          newsnc;
	double          tempc;

	/* Conversion into simulation units */
	snowhc = elem_d_ps_snowh[tid] * 100.0;
	newsnc = newsn * 100.0;

	/* Calculating new snowfall density depending on temperature equation from
	* Gottlib L. 'A general runoff model for snowcovered and glacierized
	* basin', 6th Nordic Hydrological Conference, Vemadolen, Sweden, 1980,
	* 172-177pp. */
	tempc = elem_d_es_sfctmp[tid] - 273.15;

	if (tempc <= -15.0)
	{
		dsnew = 0.05;
	}
	else
	{
		dsnew = 0.05 + 0.0017 * pow(tempc + 15.0, 1.5);
	}

	/* Adjustment of snow density depending on new snowfall */
	hnewc = newsnc / dsnew;
	if (snowhc + hnewc < 0.001)
	{
		elem_d_ps_sndens[tid] = (dsnew > elem_d_ps_sndens[tid]) ? dsnew : elem_d_ps_sndens[tid];
	}
	else
	{
		elem_d_ps_sndens[tid] = (snowhc * elem_d_ps_sndens[tid] + hnewc * dsnew) / (snowhc + hnewc);
	}
	snowhc = snowhc + hnewc;
	elem_d_ps_snowh[tid] = snowhc * 0.01;
}


// 计算边界上的土壤温度
__device__  realtype TBnd(realtype tu, realtype tb, realtype zbot, int k, 
	realtype zsoil0, realtype zsoil1, realtype zsoil2, realtype zsoil3,
	realtype zsoil4, realtype zsoil5, realtype zsoil6, realtype zsoil7,
	realtype zsoil8, realtype zsoil9, realtype zsoil10)
{
	/*
	* Subroutine TBnd
	*
	* Calculate temperature on the boundary of the layer by interpolation of
	* the middle layer temperatures */
	realtype          zb, zup;
	realtype          tbnd1;

	/* Use surface temperature on the top of the first layer */
	if (k == 0)
	{
		zup = 0.0;
	}
	else if (k == 1){
		zup = zsoil0;
	}
	else if (k == 2){
		zup = zsoil1;
	}
	else if (k == 3){
		zup = zsoil2;
	}
	else if (k == 4){
		zup = zsoil3;
	}
	else if (k == 5){
		zup = zsoil4;
	}
	else if (k == 6){
		zup = zsoil5;
	}
	else if (k == 7){
		zup = zsoil6;
	}
	else if (k == 8){
		zup = zsoil7;
	}
	else if (k == 9){
		zup = zsoil8;
	}
	else if (k == 10){
		zup = zsoil9;
	}
	else
	{
	//	zup = zsoil[k - 1];
	}

	/* Use depth of the constant bottom temperature when interpolate temperature
	* into the last layer boundary */
	// Now nsoil is MAXLYR, that is 11
	if (k == 10)
	{
		zb = 2.0 * zbot - zsoil10;
	}
	else if (k == 9)  // k != nsoil -1
	{
		zb = zsoil10;
	}
	else if (k == 8) 
	{
		zb = zsoil9;
	}
	else if (k == 7)
	{
		zb = zsoil8;
	}
	else if (k == 6)
	{
		zb = zsoil7;
	}
	else if (k == 5)
	{
		zb = zsoil6;
	}
	else if (k == 4)
	{
		zb = zsoil5;
	}
	else if (k == 3)
	{
		zb = zsoil4;
	}
	else if (k == 2)
	{
		zb = zsoil3;
	}
	else if (k == 1)
	{
		zb = zsoil2;
	}
	else if (k == 0)
	{
		zb = zsoil1;
	}
	else {  // k != nsoil -1
		//zb = zsoil[k + 1];
	}

	/* Linear interpolation between the average layer temperatures */
	if (k == 0){
		tbnd1 = tu + (tb - tu) * (zup - zsoil0) / (zup - zb);
	}
	else if (k == 1){
		tbnd1 = tu + (tb - tu) * (zup - zsoil1) / (zup - zb);
	}
	else if (k == 2){
		tbnd1 = tu + (tb - tu) * (zup - zsoil2) / (zup - zb);
	}
	else if (k == 3){
		tbnd1 = tu + (tb - tu) * (zup - zsoil3) / (zup - zb);
	}
	else if (k == 4){
		tbnd1 = tu + (tb - tu) * (zup - zsoil4) / (zup - zb);
	}
	else if (k == 5){
		tbnd1 = tu + (tb - tu) * (zup - zsoil5) / (zup - zb);
	}
	else if (k == 6){
		tbnd1 = tu + (tb - tu) * (zup - zsoil6) / (zup - zb);
	}
	else if (k == 7){
		tbnd1 = tu + (tb - tu) * (zup - zsoil7) / (zup - zb);
	}
	else if (k == 8){
		tbnd1 = tu + (tb - tu) * (zup - zsoil8) / (zup - zb);
	}
	else if (k == 9){
		tbnd1 = tu + (tb - tu) * (zup - zsoil9) / (zup - zb);
	}
	else if (k == 10){
		tbnd1 = tu + (tb - tu) * (zup - zsoil10) / (zup - zb);
	}
	else{
		// No case.
	}

	return tbnd1;
}










__constant__ realtype    THKICE = 2.2;
__constant__ realtype    THKQTZ = 7.7;
__constant__ realtype    THKW = 0.57;
__device__ realtype TDfCnd(realtype smc, realtype qz, realtype smcmax, realtype smcmin, realtype sh2o)
{
	/*
	* Function TDfCnd
	*
	* Calculate thermal diffusivity and conductivity of the soil for a given
	* point and time.
	* Peters-Lidard approach (Peters-Lidard et al., 1998)
	* June 2001 changes: frozen soil condition.
	*/
	realtype          df;
	realtype          ake;
	realtype          gammd;
	realtype          thkdry;
	realtype          thko;
	realtype          thksat;
	realtype          thks;
	realtype          satratio;
	realtype          xu;
	realtype          xunfroz;

	/* We now get quartz as an input argument
	*
	* If the soil has any moisture content compute a partial sum/product
	* otherwise use a constant value which works well with most soils
	*
	*  thkw ......water thermal conductivity
	*  thkqtz ....thermal conductivity for quartz
	*  thko ......thermal conductivity for other soil components
	*  thks ......thermal conductivity for the solids combined(quartz+other)
	*  THKICE ....ice thermal conductivity
	*  smcmax ....porosity (= smcmax)
	*  qz .........quartz content (soil type dependent)
	*
	* Use as in Peters-lidard, 1998 (Modif. from Johansen, 1975).
	* Pablo Grunmann, 08/17/98
	* Refs.:
	*  Farouki, O. T.,1986: Thermal properties of soils. Series on rock and
	*      soil mechanics, Vol. 11, trans tech, 136 pp.
	*  Johansen, O., 1975: Thermal conductivity of soils. Ph.D. thesis,
	*      University of Trondheim
	*  Peters-Lidard, C. D., et al., 1998: The effect of soil thermal
	*      conductivity parameterization on surface energy fluxes and
	*      temperatures. Journal of the Atmospheric Sciences, Vol. 55,
	*      pp. 1209-1224.
	*/
	satratio = (smc - smcmin) / (smcmax - smcmin);

	/* Thermal conductivity of "other" soil components */
	thko = 2.0;

	/* Solids' conductivity */
	thks = pow(THKQTZ, qz) * pow(thko, 1.0 - qz);

	/* Unfrozen fraction (from 1, i.e., 100% liquid, to 0 (100% frozen)) */
	xunfroz = sh2o / smc;

	/* Unfrozen volume for saturation (porosity * xunfroz) */
	xu = xunfroz * smcmax;

	/* Saturated thermal conductivity */
	thksat = pow(thks, 1.0 - smcmax) * pow(THKICE, smcmax - xu) * pow(THKW, xu);

	/* Dry density in kg/m3 */
	gammd = (1.0 - smcmax) * 2700.0;

	/* Dry thermal conductivity in W m-1 K-1 */
	thkdry = (0.135 * gammd + 64.7) / (2700.0 - 0.947 * gammd);

	if ((sh2o + 0.0005) < smc)
	{
		/* Frozen */
		ake = satratio;
	}
	else
	{
		/* Unfrozen
		* range of validity for the kersten number (ake) Kersten number (using
		* "fine" formula, valid for soils containing at least 5% of particles
		* with diameter less than 2.e-6 meters.)
		* (for "coarse" formula, see Peters-Lidard et al., 1998). */
		if (satratio > 0.1)
		{
			ake = log10(satratio) + 1.0;
		}
		else
		{
			/* Use k = kdry */
			ake = 0.0;
		}
	}

	/* Thermal conductivity */
	df = ake * (thksat - thkdry) + thkdry;

	return df;
}

__device__ realtype TmpAvg(realtype tup, realtype tm, realtype tdn, int k,
	realtype zsoil0, realtype zsoil1, realtype zsoil2, realtype zsoil3,
	realtype zsoil4, realtype zsoil5, realtype zsoil6, realtype zsoil7,
	realtype zsoil8, realtype zsoil9, realtype zsoil10)
{
	/*
	* Function TmpAvg
	*
	* Calculate soil layer average temperature (tavg) in freezing/thawing layer
	* using up, down, and middle layer temperatures (tup, tdn, tm), where tup
	* is at top boundary of layer, tdn is at bottom boundary of layer.
	* Tm is layer prognostic state temperature.
	*/
	realtype          dz;
	realtype          dzh;
	realtype          x0;
	realtype          xdn;
	realtype          xup;
	realtype          tavg;

	if (k == 0)
	{
		dz = -zsoil0;
	}
	else if (k==1)
	{
		dz = zsoil0 - zsoil1;
	}
	else if (k == 2)
	{
		dz = zsoil1 - zsoil2;
	}
	else if (k == 3)
	{
		dz = zsoil2 - zsoil3;
	}
	else if (k == 4)
	{
		dz = zsoil3 - zsoil4;
	}
	else if (k == 5)
	{
		dz = zsoil4 - zsoil5;
	}
	else if (k == 6)
	{
		dz = zsoil5 - zsoil6;
	}
	else if (k == 7)
	{
		dz = zsoil6 - zsoil7;
	}
	else if (k == 8)
	{
		dz = zsoil7 - zsoil8;
	}
	else if (k == 9)
	{
		dz = zsoil8 - zsoil9;
	}
	else if (k == 10)
	{
		dz = zsoil9 - zsoil10;
	}
	else{
		dz = 0.0001;   // 超过了最大土壤分层数 MAXLYR
	}

	dzh = dz * 0.5;

	if (tup < TFREEZ)
	{
		if (tm < TFREEZ)
		{
			if (tdn < TFREEZ)
			{
				/* tup, tm, tdn < TFREEZ */
				tavg = (tup + 2.0 * tm + tdn) / 4.0;
			}
			else
			{
				/* tup & tm < TFREEZ, tdn > TFREEZ */
				x0 = (TFREEZ - tm) * dzh / (tdn - tm);
				tavg = 0.5 *
					(tup * dzh + tm * (dzh + x0) + TFREEZ * (2.0 * dzh - x0)) /
					dz;
			}
		}
		else
		{
			if (tdn < TFREEZ)
			{
				/* tup < TFREEZ, tm > TFREEZ, tdn < TFREEZ */
				xup = (TFREEZ - tup) * dzh / (tm - tup);
				xdn = dzh - (TFREEZ - tm) * dzh / (tdn - tm);
				tavg = 0.5 *
					(tup * xup + TFREEZ * (2.0 * dz - xup - xdn) + tdn * xdn) /
					dz;
			}
			else
			{
				/* tup < TFREEZ, tm > TFREEZ, tdn > TFREEZ */
				xup = (TFREEZ - tup) * dzh / (tm - tup);
				tavg = 0.5 * (tup * xup + TFREEZ * (2.0 * dz - xup)) / dz;
			}
		}
	}
	else
	{
		if (tm < TFREEZ)
		{
			if (tdn < TFREEZ)
			{
				/* tup > TFREEZ, tm < TFREEZ, tdn < TFREEZ */
				xup = dzh - (TFREEZ - tup) * dzh / (tm - tup);
				tavg = 0.5 *
					(TFREEZ * (dz - xup) + tm * (dzh + xup) + tdn * dzh) / dz;
			}
			else
			{
				/* tup > TFREEZ, tm < TFREEZ, tdn > TFREEZ */
				xup = dzh - (TFREEZ - tup) * dzh / (tm - tup);
				xdn = (TFREEZ - tm) * dzh / (tdn - tm);
				tavg = 0.5 *
					(TFREEZ * (2.0 * dz - xup - xdn) + tm * (xup + xdn)) / dz;
			}
		}
		else
		{
			if (tdn < TFREEZ)
			{
				/* tup > TFREEZ, tm > TFREEZ, tdn < TFREEZ */
				xdn = dzh - (TFREEZ - tm) * dzh / (tdn - tm);
				tavg = (TFREEZ * (dz - xdn) + 0.5 * (TFREEZ + tdn) * xdn) / dz;
			}
			else
			{
				/* tup > TFREEZ, tm > TFREEZ, tdn > TFREEZ */
				tavg = (tup + 2.0 * tm + tdn) / 4.0;
			}
		}
	}

	return tavg;
}




/* Lech's surface functions */
__device__ realtype Pslmu(realtype zz)
{
	return -0.96 * log(1.0 - 4.5 * zz);
}

__device__ realtype Pslms(realtype zz)
{
	const realtype    RIC = 0.183;
	realtype          rric;

	rric = 1.0 / RIC;
	return zz * rric - 2.076 * (1.0 - 1.0 / (zz + 1.0));
}

__device__ realtype Pslhu(realtype zz)
{
	return -0.96 * log(1.0 - 4.5 * zz);
}

__device__ realtype Pslhs(realtype zz)
{
	const realtype    RIC = 0.183;
	const realtype    FHNEU = 0.8;
	const realtype    RFC = 0.191;
	realtype          rfac;

	rfac = RIC / (FHNEU * RFC * RFC);
	return zz * rfac - 2.076 * (1.0 - exp(-1.2 * zz));
}

/* Paulson's surface functions */
__device__ realtype Pspmu(realtype xx)
{
	return -2.0 * log((xx + 1.0) * 0.5) - log((xx * xx + 1.0) * 0.5) +
		2.0 * atan(xx) - PI * 0.5;
}

__device__ realtype Pspms(realtype yy)
{
	return 5.0 * yy;
}

__device__ realtype Psphu(realtype xx)
{
	return -2.0 * log((xx * xx + 1.0) * 0.5);
}

__device__ realtype Psphs(realtype yy)
{
	return 5.0 * yy;
}

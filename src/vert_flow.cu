#include "pihm.h"
#include "pihm.cuh"

__device__ void VerticalFlow(realtype  dt, int tid, 
#include "VerticalFlow_device_declarations.c" 
)
{
	/* Calculate infiltration rate */
	elem_d_wf_infil[tid] = 
		Infil(dt, tid, 
#include "Infil_device_arguments.c" 
	     );

	/* Calculate recharge rate */
	elem_d_wf_rechg[tid] = 
		Recharge(tid, 
#include "Recharge_device_arguments.c" 
	     );

}

__device__ realtype Infil(realtype  dt, int tid, 
#include "Infil_device_declarations.c" 
)
{
	realtype          applrate;
	realtype          wetfrac;
	realtype          dh_by_dz;
	realtype          satn;
	realtype          satkfunc;
	realtype          infil;
	realtype          infil_max;
	realtype          kinf;
	realtype          deficit;
	realtype          psi_u;
	realtype          h_u;
	
	if (elem_d_ws_unsat[tid] + elem_d_ws_gw[tid] > elem_d_soil_depth[tid])
	{
		infil = 0.0;
	}
	else
	{
		/* Water appliation rate is the sum of net gain of overland flow and net
		* precipitation */
		applrate = 0.0;
		applrate += -(elem_d_wf_ovlflow0[tid] + elem_d_wf_ovlflow1[tid] +
			            elem_d_wf_ovlflow2[tid]) / elem_d_topo_area[tid];
		applrate = (applrate > 0.0) ? applrate : 0.0;
		applrate += elem_d_wf_pcpdrp[tid];

		if (DEPRSTG > 0.0)
		{
			wetfrac = elem_d_ws_surfh[tid] / DEPRSTG;
			wetfrac = (wetfrac > 0.0) ? wetfrac : 0.0;
			wetfrac = (wetfrac < 1.0) ? wetfrac : 1.0;
		}
		else
		{
			wetfrac = (elem_d_ws_surfh[tid] > 0.0) ? 1.0 : 0.0;
		}

		if (elem_d_ws_gw[tid] > elem_d_soil_depth[tid] - elem_d_soil_dinf[tid])
		{
			/* Assumption: Dinf < Dmac */
			dh_by_dz = (elem_d_ws_surfh[tid] + elem_d_topo_zmax[tid] - 
				       (elem_d_ws_gw[tid] + elem_d_topo_zmin[tid])) /
				(0.5 * (elem_d_ws_surfh[tid] + elem_d_soil_dinf[tid]));
			dh_by_dz = (elem_d_ws_surfh[tid] <= 0.0 && dh_by_dz > 0.0) ? 0.0 : dh_by_dz;

			satn = 1.0;
			satkfunc = KrFunc_cuda(elem_d_soil_beta[tid], satn);

			kinf = EffKinf(dh_by_dz, satkfunc, satn, applrate, elem_d_ws_surfh[tid], tid, 
#include "EffKinf_device_arguments.c" 
				           );

			infil = kinf * dh_by_dz;
		}
		else
		{
			deficit = elem_d_soil_depth - elem_d_ws_gw;

			satn = elem_d_ws_unsat[tid] / deficit;
			satn = (satn > 1.0) ? 1.0 : satn;
			satn = (satn < SATMIN) ? SATMIN : satn;

			psi_u = Psi(satn, elem_d_soil_alpha[tid], elem_d_soil_beta[tid]);
			/* Note: for psi calculation using van Genuchten relation, cutting
			* the psi-sat tail at small saturation can be performed for
			* computational advantage. If you do not want to perform this,
			* comment the statement that follows */
			psi_u = (psi_u > PSIMIN) ? psi_u : PSIMIN;

			h_u = psi_u + elem_d_topo_zmax[tid] - 0.5 * elem_d_soil_dinf[tid];
			dh_by_dz = (elem_d_ws_surfh[tid] + elem_d_topo_zmax[tid] - h_u) /
				(0.5 * (elem_d_ws_surfh[tid] + elem_d_soil_dinf[tid]));
			dh_by_dz = (elem_d_ws_surfh[tid] <= 0.0 && dh_by_dz > 0.0) ? 0.0 : dh_by_dz;

			satkfunc = KrFunc_cuda(elem_d_soil_beta[tid], satn);

			kinf = EffKinf(dh_by_dz, satkfunc, satn, applrate, elem_d_ws_surfh[tid], tid, 
#include "EffKinf_device_arguments.c" 
				           );

			infil = kinf * dh_by_dz;
			infil = (infil > 0.0) ? infil : 0.0;
		}

		infil_max = applrate + ((elem_d_ws0_surf[tid] > 0.0) ? elem_d_ws0_surf[tid] / dt : 0.0);

		infil = (infil > infil_max) ? infil_max : infil;

		infil *= wetfrac;
	}

	return infil;
}


__device__ realtype  Recharge(int tid, 
#include "Recharge_device_declarations.c" 
)
{
	realtype          satn;
	realtype          satkfunc;
	realtype          dh_by_dz;
	realtype          psi_u;
	realtype          kavg;
	realtype          deficit;
	realtype          rechg;

	if (elem_d_ws_gw[tid] > elem_d_soil_depth[tid] - elem_d_soil_dinf[tid])
	{
		rechg = elem_d_wf_infil[tid];
	}
	else
	{
		deficit = elem_d_soil_depth[tid] - elem_d_ws_gw[tid];
		satn = elem_d_ws_unsat[tid] / deficit;
		satn = (satn > 1.0) ? 1.0 : satn;
		satn = (satn < SATMIN) ? SATMIN : satn;

		satkfunc = KrFunc_cuda(elem_d_soil_beta[tid], satn);

		psi_u = Psi(satn, elem_d_soil_alpha[tid], elem_d_soil_beta[tid]);

		dh_by_dz =
			(0.5 * deficit + psi_u) / (0.5 * (deficit + elem_d_ws_gw[tid]));

		kavg = AvgKv(deficit, elem_d_ws_gw[tid], satkfunc, tid, 
#include "AvgKv_device_arguments.c" 
			         );

		rechg = kavg * dh_by_dz;

		rechg = (rechg > 0.0 && elem_d_ws_unsat[tid] <= 0.0) ? 0.0 : rechg;
		rechg = (rechg < 0.0 && elem_d_ws_gw[tid] <= 0.0) ? 0.0 : rechg;
	}

	return rechg;
}


__device__ realtype  AvgKv(realtype  deficit, realtype  gw,
	                    realtype  satkfunc, int tid,
#include "AvgKv_device_declarations.c" 
)
{
	realtype          k1, k2, k3;
	realtype          d1, d2, d3;

	if (deficit > elem_d_soil_dmac[tid])
	{
		k1 = satkfunc * elem_d_soil_ksatv[tid];
		d1 = elem_d_soil_dmac[tid];

		k2 = satkfunc * elem_d_soil_ksatv[tid];
		d2 = deficit - elem_d_soil_dmac[tid];

		k3 = elem_d_soil_ksatv[tid];
		d3 = gw;
	}
	else
	{
		k1 = satkfunc * elem_d_soil_ksatv[tid];
		d1 = deficit;

		k2 = (elem_d_soil_areafh[tid] > 0.0) ?
			elem_d_soil_kmacv[tid] * elem_d_soil_areafh[tid] 
			+ elem_d_soil_ksatv[tid] * (1.0 - elem_d_soil_areafh[tid]) :
			elem_d_soil_ksatv[tid];
		d2 = elem_d_soil_dmac[tid] - deficit;

		k3 = elem_d_soil_ksatv[tid];
		d3 = gw - (elem_d_soil_dmac[tid] - deficit);
	}

#if defined(_ARITH_)
	/* Arithmetic mean formulation */
	return (k1 * d1 + k2 * d2 + k3 * d3) / (d1 + d2 + d3);
#else
	return (d1 + d2 + d3) / (d1 / k1 + d2 / k2 + d3 / k3);
#endif
}

// 与soilc中的 KrFunc() 相同
__device__ realtype  KrFunc_cuda(realtype  beta, realtype  satn)
{
    return sqrt(satn) *
		(1.0 - pow(1.0 - pow(satn, beta / (beta - 1.0)), (beta - 1.0) / beta)) *
		(1.0 - pow(1.0 - pow(satn, beta / (beta - 1.0)), (beta - 1.0) / beta));
}


__device__ realtype  EffKinf(realtype dh_by_dz, realtype ksatfunc,
	                         realtype elemsatn, realtype applrate, 
							 realtype surfh, int tid,
#include "EffKinf_device_declarations.c" 
)
{
	/*
	* For infiltration, macropores act as cracks, and are hydraulically
	* effective in rapidly conducting water flow (Chen and Wagenet, 1992,
	* Journal of Hydro__logfy, 130, 105-126).
	* For macropore hydraulic conductivities, use van Genuchten parameters for
	* fractures (Gerke and van Genuchten, 1993, Water Resources Research, 29,
	* 1225-1238).
	*/
	realtype          keff = 0.0;
	realtype          kmax;
	const realtype    BETA_CRACK = 2.0;

	if (elem_d_soil_areafh[tid] == 0.0)
	{
		/* Matrix */
		keff = elem_d_soil_kinfv[tid] * ksatfunc;
	}
	else if (surfh > DEPRSTG)
	{
		/* When surface wet fraction is larger than 1 (surface is totally
		* ponded), i.e., surfh > DEPRSTG, flow situation is macropore control,
		* regardless of the application rate */
		keff = elem_d_soil_kinfv[tid] * (1.0 - elem_d_soil_areafh[tid]) * ksatfunc +
			elem_d_soil_kmacv[tid] * elem_d_soil_areafh[tid];
	}
	else
	{
		if (applrate <= dh_by_dz * elem_d_soil_kinfv[tid] * ksatfunc)
		{
			/* Matrix control */
			keff = elem_d_soil_kinfv[tid] * ksatfunc;
		}
		else
		{
			kmax = dh_by_dz * (elem_d_soil_kmacv[tid] * elem_d_soil_areafh[tid] +
				elem_d_soil_kinfv[tid] * (1.0 - elem_d_soil_areafh[tid]) * ksatfunc);
			if (applrate < kmax)
			{
				/* Application control */
				keff = elem_d_soil_kinfv[tid] * (1.0 - elem_d_soil_areafh[tid]) * ksatfunc +
					   elem_d_soil_kmacv[tid] * elem_d_soil_areafh[tid] *
				KrFunc_cuda(BETA_CRACK, elemsatn);
			}
			else
			{
				/* Macropore control */
				keff = elem_d_soil_kinfv[tid] * (1.0 - elem_d_soil_areafh[tid]) * ksatfunc +
					   elem_d_soil_kmacv[tid] * elem_d_soil_areafh[tid];
			}
		}
	}

	return keff;
}

__device__ realtype Psi(realtype satn, realtype alpha, realtype beta)
{
	satn = (satn < SATMIN) ? SATMIN : satn;

	/* van Genuchten 1980 SSSAJ */
	return -pow(pow(1.0 / satn, beta / (beta - 1.0)) - 1.0, 1.0 / beta) / alpha;
}


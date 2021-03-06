#include "pihm.h"
#if defined(_CVODE_CUDA)
#include <math.h>    // math function on CPU.
#include "pihm.cuh"
#endif

#if defined(_CVODE_CUDA)
void Spinup(pihm_struct pihm, pihm_struct_d pihm_d, 
	        N_Vector CV_Y, N_Vector CV_Yout, void *cvode_mem)
#else
void Spinup(pihm_struct pihm, N_Vector CV_Y, void *cvode_mem)
#endif
{
    int             spinyears = 0;
    int             first_spin_cycle = 1;
    int             steady;
    int             metyears;
    ctrl_struct    *ctrl;

    ctrl = &pihm->ctrl;

    metyears = (ctrl->endtime - ctrl->starttime) / DAYINSEC / 365;

    do
    {
        PIHMprintf(VL_NORMAL, "Spinup year: %6d\n", spinyears + 1);

#if defined(_BGC_)
        ResetSpinupStat(pihm->elem);
#endif

        for (ctrl->cstep = 0; ctrl->cstep < ctrl->nstep; ctrl->cstep++)
        {
#if defined(_CVODE_CUDA)
			PIHM(pihm, pihm_d, cvode_mem, CV_Y, CV_Yout, 0.0);	
#else
            PIHM(pihm, cvode_mem, CV_Y, 0.0);
#endif
        }
#if defined(_CVODE_CUDA)
		SetCVodeParam(pihm, pihm_d, cvode_mem, CV_Y);
#else
        /* Reset solver parameters */
        SetCVodeParam(pihm, cvode_mem, CV_Y);
#endif
        spinyears += metyears;

#if defined(_BGC_)
        steady = CheckSteadyState(pihm->elem, pihm->siteinfo.area,
            first_spin_cycle, ctrl->endtime - ctrl->starttime,
            spinyears);
#else
        steady = CheckSteadyState(pihm->elem, pihm->siteinfo.area,
            first_spin_cycle, spinyears);
#endif

        first_spin_cycle = 0;
    } while (spinyears < ctrl->maxspinyears && (!steady));
}

#if defined(_BGC_)
void ResetSpinupStat(elem_struct *elem)
{
    int             i;

#if defined(_LUMPED_)
    i = LUMPED;
#else
# if defined(_OPENMP)
#  pragma omp parallel for
# endif
    for (i = 0; i < nelem; i++)
#endif
    {
        elem[i].spinup.soilc = 0.0;
        elem[i].spinup.totalc = 0.0;
    }
}
#endif

#if defined(_BGC_)
int CheckSteadyState(const elem_struct *elem, realtype total_area,
    int first_cycle, int totalt, int spinyears)
#else
int CheckSteadyState(const elem_struct *elem, realtype total_area,
    int first_cycle, int spinyears)
#endif
{
    int             i;
    realtype          totalw = 0.0;
    static realtype   totalw_prev = 0.0;
    int             steady;
#if defined(_BGC_)
    realtype          t1 = 0.0;
    realtype          soilc = 0.0;
    realtype          totalc = 0.0;
    static realtype   soilc_prev = 0.0;
#endif
#if defined(_FBR_)
    realtype          fbrgw = 0.0;
    static realtype   fbrgw_prev = 0.0;
#endif

#if defined(_LUMPED_)
    i = LUMPED;
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        totalw += (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity *
                   elem[i].topo.area / total_area;

#if defined(_FBR_)
        totalw += (elem[i].ws.fbr_unsat + elem[i].ws.fbr_gw) *
            elem[i].geol.porosity * elem[i].topo.area / total_area;

        fbrgw += elem[i].ws.fbr_gw * elem[i].topo.area / total_area;
#endif

#if defined(_BGC_)
        /* Convert soilc and totalc to average daily soilc */
        soilc += elem[i].spinup.soilc * elem[i].topo.area /
            (realtype)(totalt / DAYINSEC) / total_area;
        totalc += elem[i].spinup.totalc * elem[i].topo.area /
            (realtype)(totalt / DAYINSEC) / total_area;
        t1 = (soilc - soilc_prev) / (realtype)(totalt / DAYINSEC / 365);
#endif
    }

    if (!first_cycle)
    {
        /* Check if domain reaches steady state */
        steady = (fabs(totalw - totalw_prev) < SPINUP_W_TOLERANCE);
#if defined(_FBR_)
        steady = (steady && (fabs(fbrgw - fbrgw_prev) < SPINUP_W_TOLERANCE));
#endif
#if defined(_BGC_)
        steady = (steady && (fabs(t1) < SPINUP_C_TOLERANCE));
#endif

        PIHMprintf(VL_NORMAL, "spinyears = %d ", spinyears);
        PIHMprintf(VL_NORMAL, "totalw_prev = %lg totalw = %lg wdif = %lg\n",
            totalw_prev, totalw, totalw - totalw_prev);
#if defined(_FBR_)
        PIHMprintf(VL_NORMAL, "fbrgw_prev = %lg fbrgw = %lg wdif = %lg\n",
            fbrgw_prev, fbrgw, fbrgw - fbrgw_prev);
#endif
#if defined(_BGC_)
        PIHMprintf(VL_NORMAL, "soilc_prev = %lg soilc = %lg pdif = %lg\n",
            soilc_prev, soilc, t1);
#endif
    }
    else
    {
        steady = 0;
    }

    totalw_prev = totalw;
#if defined(_FBR_)
    fbrgw_prev = fbrgw;
#endif
#if defined(_BGC_)
    soilc_prev = soilc;
#endif

    if (steady)
    {
        PIHMprintf(VL_NORMAL, "Reaches steady state after %d year.\n",spinyears);            
    }

    return steady;
}

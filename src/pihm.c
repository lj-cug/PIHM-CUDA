#include "pihm.h"
#if defined(_CVODE_CUDA)
#include "pihm.cuh"
#endif

#if defined(_CVODE_CUDA)
void PIHM(pihm_struct pihm, pihm_struct_d pihm_d,
	void *cvode_mem, N_Vector CV_Y, N_Vector CV_Yout, realtype cputime)
#else
void PIHM(pihm_struct pihm, void *cvode_mem, N_Vector CV_Y, realtype cputime)
#endif
{
    int             t;

    t = pihm->ctrl.tout[pihm->ctrl.cstep];

    /* Apply boundary conditions 
	   ÿһ����ʩ�ӱ߽�����
	*/
#if defined(_CVODE_CUDA)
	ApplyBc(&pihm->forc, pihm->elem, pihm->river, pihm_d, t);
#else
    ApplyBc(&pihm->forc, pihm->elem, pihm->river, t);   
#endif

    /*
     * Apply forcing and simulate land surface processes 
	   ��һ����ʱ��Ƶ��ʩ��½�ر���������
     */
    if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)   
    {
        /* Apply forcing */
#if defined(_NOAH_)
#if defined(_CVODE_CUDA)
		ApplyForc(&pihm->forc, pihm->elem, pihm_d, t, pihm->ctrl.rad_mode,
			      &pihm->siteinfo);
#else
        ApplyForc(&pihm->forc, pihm->elem, t , pihm->ctrl.rad_mode,
                  &pihm->siteinfo);
#endif
#else
#if defined(_CVODE_CUDA)
		ApplyForc(&pihm->forc, pihm->elem, pihm_d, t);
#else
        ApplyForc(&pihm->forc, pihm->elem, t); 
#endif
#endif

#if defined(_NOAH_)
#if defined(_CVODE_CUDA)
		//printf("Start to calculate Noah LSM.\n");
		// Noah2d�еĶ���ָ����Ȼ��CVode��Ӱ��!?    ~2021.06.07
		Noah(pihm_d, (realtype)pihm->ctrl.etstep); 
		//printf("Finished to calculate Noah LSM.\n");
#else  /* CPU */
        /* Calculate surface energy balance */
        Noah(pihm->elem, (realtype)pihm->ctrl.etstep);
#endif
#else  /* is_sm_et Model */
#if defined(_CVODE_CUDA)
		// cuda�˺�������elem.wf��elem.ws�ı���
		IntcpSnowEt(t, pihm->ctrl.stepsize, pihm_d);
#else  /* CPU */
        /* Calculate Interception storage and ET */
        IntcpSnowEt(t, (realtype)pihm->ctrl.etstep, pihm->elem, &pihm->cal); 
#endif
#endif

        /* Update print variables for land surface step variables 
		   ��һ����Ƶ��, ������½�ر�����㲽����������������´�����ı���ֵ.
		   �˴���varctrl��tp_varctrl����map_output�л�ȡ��CPU�ô����������ָ��.
		   ���ǣ�ֱ�ӻ�ȡUM�е�����ָ�룬��������.
		*/
#ifndef _CVODE_CUDA
       UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, LS_STEP);
       UpdPrintVar(pihm->print.tp_varctrl, pihm->print.ntpprint, LS_STEP);
	  // printf("Finished to Update Land Surface Variables.\n\n");
#endif
    }  // ÿ��һ��ʱ�䣬ʩ�����������߽�������½�ر������.

    // ʹ��CVode���ODE�����¼�����CV_Y
#if defined(_CVODE_CUDA)
/*
    CV_Yout��CVode���õ����µ�NVector, ʹ��CUDA�汾��CVodeʱ��Ҫ���������  2021.05.07
*/
    SolveCVode(pihm->ctrl.starttime, &t, pihm->ctrl.tout[pihm->ctrl.cstep + 1], 
		       cputime, cvode_mem, CV_Yout);   
#else
	SolveCVode(pihm->ctrl.starttime, &t, pihm->ctrl.tout[pihm->ctrl.cstep + 1],
			cputime, cvode_mem, CV_Y);
#endif
	 // printf("Finished to Solve CVode in CUDA.\n\n");

    /* Use mass balance to calculate model fluxes or variables */
#if defined(_CVODE_CUDA)
// ���ڣ���CV_Yout����elem��river�ṹ���еı���.  2021.05.07
	Summary(CV_Yout, pihm_d, pihm->ctrl.stepsize);
#else
    Summary(pihm->elem, pihm->river, CV_Y, (realtype)pihm->ctrl.stepsize);  
#endif

/* ����Noah LSMģ���е�״̬����ֵ */
#if defined(_NOAH_)
#if defined(_CVODE_CUDA)
	//printf("Start to update Noah LSM.\n");
	NoahHydrol(pihm_d, (realtype)pihm->ctrl.stepsize);
	//printf("Finished to update Noah LSM.\n");
#else
    NoahHydrol(pihm->elem, (realtype)pihm->ctrl.stepsize);
#endif
#endif

#if defined(_CYCLES_)
//    SoluteTransport(pihm->elem, pihm->river, (realtype)pihm->ctrl.stepsize);
#endif

    /* Update print variables for hydrology step variables 
	   ��һ����Ƶ��, ������ˮ�Ĺ��̼��㲽����������������´�����ı���ֵ.
	   �˴���varctrl��tp_varctrl�Ѿ���ȡ��CPU��GPU�ô����������ָ��.
	*/
#ifndef _CVODE_CUDA
    UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, HYDROL_STEP);
    UpdPrintVar(pihm->print.tp_varctrl, pihm->print.ntpprint, HYDROL_STEP);
#endif

#if defined(_DAILY_)
    DailyVar(t, pihm->ctrl.starttime, pihm->elem);

    /*
     * Daily timestep modules
     */
    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
# if defined(_BGC_)
        /* Daily BGC processes */
        DailyBgc(pihm, t - DAYINSEC);
# endif

# if defined(_CYCLES_)
        DailyCycles(t - DAYINSEC, pihm);
# endif

        /* Update print variables for CN (daily) step variables */
        UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, CN_STEP);
        UpdPrintVar(pihm->print.tp_varctrl, pihm->print.ntpprint, CN_STEP);

        /* Initialize daily structures */
        InitDailyStruct(pihm->elem);
    }
#endif

    /*
     * Print outputs on CPU  ���ǽ���������ʱ,ִ�����ݴ�GPU��CPU�Ŀ���
     */
    /* Print water balance (ÿ��dt�����ˮ��ƽ�����)*/
    if (pihm->ctrl.waterbal)  // �Ƿ����ˮ��ƽ�����?
    {
#if defined(_CVODE_CUDA)
	// ��Ҫ����pihm_d�е�ʱ��
		//printf("Start to Printout Water Balance state.\n");
		//PrintWaterBal(pihm->print.watbal_file, t, pihm->ctrl.starttime,
		//	          pihm->ctrl.stepsize, pihm_d);
		//printf("Finished to Printout Water Balance state.\n");
#else
        PrintWaterBal(pihm->print.watbal_file, t, pihm->ctrl.starttime,
                      pihm->ctrl.stepsize, pihm->elem, pihm->river);
#endif
    }

	/*����м������*/
#if defined(_CVODE_CUDA)



#else
    /* Print binary and txt output files ��һ����Ƶ����������� */
    PrintData(pihm->print.varctrl, pihm->print.nprint, t,
              t - pihm->ctrl.starttime, pihm->ctrl.ascii);

    /* Print tecplot output files ��һ����Ƶ��)���Tecplot��ʽ�ļ����� */
    if (tecplot)
    {
        PrintDataTecplot(pihm->print.tp_varctrl, pihm->print.ntpprint, t,
                         t - pihm->ctrl.starttime);
    }
#endif
}

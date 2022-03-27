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
	   每一步都施加边界条件
	*/
#if defined(_CVODE_CUDA)
	ApplyBc(&pihm->forc, pihm->elem, pihm->river, pihm_d, t);
#else
    ApplyBc(&pihm->forc, pihm->elem, pihm->river, t);   
#endif

    /*
     * Apply forcing and simulate land surface processes 
	   按一定的时间频率施加陆地表面驱动力
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
		// Noah2d中的二级指针仍然对CVode有影响!?    ~2021.06.07
		Noah(pihm_d, (realtype)pihm->ctrl.etstep); 
		//printf("Finished to calculate Noah LSM.\n");
#else  /* CPU */
        /* Calculate surface energy balance */
        Noah(pihm->elem, (realtype)pihm->ctrl.etstep);
#endif
#else  /* is_sm_et Model */
#if defined(_CVODE_CUDA)
		// cuda核函数计算elem.wf和elem.ws的变量
		IntcpSnowEt(t, pihm->ctrl.stepsize, pihm_d);
#else  /* CPU */
        /* Calculate Interception storage and ET */
        IntcpSnowEt(t, (realtype)pihm->ctrl.etstep, pihm->elem, &pihm->cal); 
#endif
#endif

        /* Update print variables for land surface step variables 
		   按一定的频率, 更新在陆地表面计算步的输出变量，仅更新待输出的变量值.
		   此处，varctrl和tp_varctrl已在map_output中获取了CPU用待输出变量的指针.
		   但是，直接获取UM中的数组指针，会有问题.
		*/
#ifndef _CVODE_CUDA
       UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, LS_STEP);
       UpdPrintVar(pihm->print.tp_varctrl, pihm->print.ntpprint, LS_STEP);
	  // printf("Finished to Update Land Surface Variables.\n\n");
#endif
    }  // 每隔一段时间，施加驱动力、边界条件和陆地表面过程.

    // 使用CVode求解ODE，更新计算了CV_Y
#if defined(_CVODE_CUDA)
/*
    CV_Yout是CVode求解得到的新的NVector, 使用CUDA版本的CVode时需要定义此向量  2021.05.07
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
// 现在，用CV_Yout更新elem和river结构体中的变量.  2021.05.07
	Summary(CV_Yout, pihm_d, pihm->ctrl.stepsize);
#else
    Summary(pihm->elem, pihm->river, CV_Y, (realtype)pihm->ctrl.stepsize);  
#endif

/* 更新Noah LSM模型中的状态变量值 */
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
	   按一定的频率, 更新在水文过程计算步的输出变量，仅更新待输出的变量值.
	   此处，varctrl和tp_varctrl已经获取了CPU或GPU用待输出变量的指针.
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
     * Print outputs on CPU  考虑仅在输出结果时,执行数据从GPU到CPU的拷贝
     */
    /* Print water balance (每个dt都输出水量平衡情况)*/
    if (pihm->ctrl.waterbal)  // 是否输出水量平衡情况?
    {
#if defined(_CVODE_CUDA)
	// 需要代入pihm_d中的时间
		//printf("Start to Printout Water Balance state.\n");
		//PrintWaterBal(pihm->print.watbal_file, t, pihm->ctrl.starttime,
		//	          pihm->ctrl.stepsize, pihm_d);
		//printf("Finished to Printout Water Balance state.\n");
#else
        PrintWaterBal(pihm->print.watbal_file, t, pihm->ctrl.starttime,
                      pihm->ctrl.stepsize, pihm->elem, pihm->river);
#endif
    }

	/*输出中间计算结果*/
#if defined(_CVODE_CUDA)



#else
    /* Print binary and txt output files 按一定的频率输出计算结果 */
    PrintData(pihm->print.varctrl, pihm->print.nprint, t,
              t - pihm->ctrl.starttime, pihm->ctrl.ascii);

    /* Print tecplot output files 按一定的频率)输出Tecplot格式的计算结果 */
    if (tecplot)
    {
        PrintDataTecplot(pihm->print.tp_varctrl, pihm->print.ntpprint, t,
                         t - pihm->ctrl.starttime);
    }
#endif
}

#include "pihm.h"
int NumStateVar(void);

#if defined(_CVODE_CUDA)
#include "pihm.cuh"
static pihm_struct_d SetUserData(void);
#endif

/* Global variables */
int             verbose_mode;
int             debug_mode;
int             append_mode;
int             corr_mode;
int             spinup_mode;
int             tecplot;
char            project[MAXSTRING];
int             nelem;
int             nriver;
#if defined(_OPENMP)
int             nthreads = 1;    /* Default value */
#endif
#if defined(_BGC_)
int             first_balance;
#endif

int main(int argc, char *argv[])
{
    char            outputdir[MAXSTRING];
    pihm_struct     pihm;
#if defined(_CVODE_CUDA)
	pihm_struct_d   pihm_d;
	N_Vector        CV_Yout = NULL;
#endif
    ctrl_struct    *ctrl;
    N_Vector        CV_Y = NULL;
    void           *cvode_mem = NULL;

#if defined(_OPENMP)
    realtype        start_omp;
#else
    clock_t         start;
#endif
    realtype        cputime, cputime_dt;    /* Time cpu duration */

#if defined(_OPENMP)
    /* Set the number of threads to use */
	// export OMP_NUM_THREADS = 4   // Linux  或者在Windows CMD: set OMP_NUM_THREADS = 4 
    nthreads = omp_get_max_threads();
#endif

    memset(outputdir, 0, MAXSTRING);

    /* Read command line arguments */
    ParseCmdLineParam(argc, argv, outputdir);

    /* Print ASCII art */
    StartupScreen();

	pihm = (pihm_struct)malloc(sizeof(*pihm)); // PIHM struct on CPU, AOS-style

    /* Read PIHM input files */
    ReadAlloc(pihm);

#if defined(_CVODE_CUDA)
	/* Set model parameters */
	pihm_d = SetUserData();

	pihm_d->nelem = nelem;   // 流域单元个数
	pihm_d->nriver = nriver; // 河道分段个数
	pihm_d->NEQ = NumStateVar();  // 变量个数
	printf("Number of basin elements=%d; Number of river segments=%d\n\n", nelem, nriver);

#if defined(_SUNDIALS_v3)
	CV_Y = N_VNew_Cuda(NumStateVar());
	CV_Yout = N_VNew_Cuda(NumStateVar());
#elif defined(_SUNDIALS_v4)
	/* We can also create managed N_Vector, 
	   but we should have higher at least Ver. 4.0.0 SUNDIALS */
	CV_Y = N_VNewManaged_Cuda(NumStateVar()); 
#endif

#else
    CV_Y = N_VNew(NumStateVar());  
#endif
	if (CV_Y == NULL)
    {
        PIHMprintf(VL_ERROR, "Error creating CVODE state variable vector.\n");
        PIHMexit(EXIT_FAILURE);
    }

    /* Initialize PIHM structure on Host 
	   In CUDA case, besides PIHM structure was initialized,  
	   we also initialized the CV_Y on GPU.
	*/
#if defined(_NOAH_)
//	printf("Try to initialize the variables on CPU.\n\n");
	Initialize(pihm, pihm_d, CV_Y);
//	printf("Finished to initialize the variables on CPU.\n\n");
#else
    Initialize(pihm, CV_Y); 
#endif

#if 0
	FILE *outfile0 = NULL;
	outfile0 = fopen("CV_Y_after_initialize.txt", "w+");
	N_VPrintFile_Cuda(CV_Y, outfile0);
	fclose(outfile0);
	exit(0);
#endif

	/* Allocate memory for solver */
#if defined(_SUNDIALS_v3)
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON); 
#elif defined(_SUNDIALS_v4)
	cvode_mem = CVodeCreate(CV_BDF);  
#endif
	if (cvode_mem == NULL)
	{
		PIHMprintf(VL_ERROR, "Error in allocating memory for solver.\n");
		PIHMexit(EXIT_FAILURE);
	}

#if defined(_CVODE_CUDA)
	/* Initialize PIHM structure on Device */
	// Allocate CUDA managed space of the SoA-stype PIHM_d structure on GPU
	Initialize_cuda(pihm_d);

	// Array on CPU are transferred into PIHM_d on GPU
	AoS_to_SoA(pihm, pihm_d);   

#if 0
	// 检查pihm_d中数组内容的正确性
	FILE *outfile0 = NULL;
	outfile0 = fopen("array_in_pihm_d.txt", "w+");
	realtype temp1, temp2;
	for (int tid = 0; tid < nelem; tid++){		
		temp1 = pihm_d->elem_d_wf_pcpdrp[tid];
		temp2 = pihm->elem[tid].wf.pcpdrp;
		fprintf(outfile0,"%10.5f %10.5f\n",temp1, temp2);
	}

	fclose(outfile0);
#endif
#endif

    /* Create output directory */
    CreateOutputDir(outputdir);

    /* Create output structures */
#if defined(_CYCLES_)
    MapOutput(pihm->ctrl.prtvrbl, pihm->ctrl.tpprtvrbl, pihm->epctbl,
              pihm->elem, pihm->river, &pihm->meshtbl, outputdir, &pihm->print);
#else
#ifndef _CVODE_CUDA
    MapOutput(pihm->ctrl.prtvrbl, pihm->ctrl.tpprtvrbl, pihm->elem, pihm->river,
             &pihm->meshtbl, outputdir, &pihm->print);
#endif
#endif

    /* Backup input files */
#if !defined(_MSC_VER)
    if (!append_mode)
    {
        BackupInput(outputdir, &pihm->filename);
    }
#endif

    InitOutputFile(&pihm->print, outputdir, pihm->ctrl.waterbal,
                    pihm->ctrl.ascii);

    /* Set solver parameters */
#if defined(_CVODE_CUDA)
	SetCVodeParam(pihm, pihm_d, cvode_mem, CV_Y);
	//printf("Set solver parameters.\n\n");
#else
    SetCVodeParam(pihm, cvode_mem, CV_Y);
	//printf("Finished to Set CVode Parameters.\n\n");
#endif
#if defined(_BGC_)
    first_balance = 1;
#endif

    /*
     * Run PIHM
     */
#if defined(_OPENMP)
    start_omp = omp_get_wtime();
#else
    start = clock();
#endif

    ctrl = &pihm->ctrl;

    if (spinup_mode)
    {
#if defined(_CVODE_CUDA)
		Spinup(pihm, pihm_d, CV_Y, CV_Yout, cvode_mem);

		PrintInit(pihm_d, outputdir,
			ctrl->endtime, ctrl->starttime,
			ctrl->endtime, ctrl->prtvrbl[IC_CTRL]);
#else
		Spinup(pihm, CV_Y, cvode_mem);  

        /* In spin-up mode, initial conditions are always printed */
        PrintInit(pihm->elem, pihm->river, outputdir,
                  ctrl->endtime, ctrl->starttime,
                  ctrl->endtime, ctrl->prtvrbl[IC_CTRL]);
#endif
#if defined(_BGC_)
        WriteBgcIc(outputdir, pihm->elem, pihm->river);
#endif
    }

    else  // normal running mode
    {
        for (ctrl->cstep = 0; ctrl->cstep < ctrl->nstep; ctrl->cstep++)   
        {
#if defined(_OPENMP)
            RunTime(start_omp, &cputime, &cputime_dt);
#else
            RunTime(start, &cputime, &cputime_dt);
#endif

#if defined(_CVODE_CUDA)
			PIHM(pihm, pihm_d, cvode_mem, CV_Y, CV_Yout, cputime);
#else
            PIHM(pihm, cvode_mem, CV_Y, cputime);
#endif
            /* Adjust CVODE max step to reduce oscillation */
            AdjCVodeMaxStep(cvode_mem, &pihm->ctrl);

            /* Print CVODE performance and statistics */
            if (debug_mode)
            {
                PrintPerf(cvode_mem, ctrl->tout[ctrl->cstep + 1],
                    ctrl->starttime, cputime_dt, cputime,
                    ctrl->maxstep, pihm->print.cvodeperf_file);
            }

            /* Write init files */
            if (ctrl->write_ic)
            {
#if defined(_CVODE_CUDA)
				PrintInit(pihm_d, outputdir,
					ctrl->tout[ctrl->cstep + 1], ctrl->starttime,
					ctrl->endtime, ctrl->prtvrbl[IC_CTRL]);
#else
                PrintInit(pihm->elem, pihm->river, outputdir,
                    ctrl->tout[ctrl->cstep + 1], ctrl->starttime,
                    ctrl->endtime, ctrl->prtvrbl[IC_CTRL]);
#endif
            }
        }  // 时间推进计算

#if defined(_BGC_)
        if (ctrl->write_bgc_restart)
        {
            WriteBgcIc(outputdir, pihm->elem, pihm->river);
        }
#endif

# if TEMP_DISABLED
#if defined(_CYCLES_)
        if (ctrl->write_cycles_restart)
        {
            WriteCyclesIC(pihm->filename.cyclesic, pihm->elem, pihm->river);
        }
# endif
#endif
    }

    if (debug_mode)
    {
        PrintCVodeFinalStats(cvode_mem);
    }

    /* Free memory */
    N_VDestroy(CV_Y);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

	/* free CPU memory */
    FreeMem(pihm);
    free(pihm);

#if defined(_CVODE_CUDA)
	N_VDestroy(CV_Yout);

	/* free GPU memory */
	FreeMem_cuda(pihm_d);
	cudaFree(pihm_d);
#endif

    PIHMprintf(VL_BRIEF, "\nSimulation completed.\n");
    return EXIT_SUCCESS;
}


#if defined(_CVODE_CUDA)
pihm_struct_d SetUserData(void)
{
	/* Allocate user data structure */
	pihm_struct_d Userdata = (pihm_struct_d)malloc(sizeof *Userdata);  // Allocation must be on CPU

	return Userdata;
}
#endif


// 计算状态变量个数，也就是方程个数
int NumStateVar(void)
{
	/*
	* Return number of state variables
	*/
	int             nsv;

	nsv = 3 * nelem + 2 * nriver;

#if defined(_BGC_)
# if defined(_LUMPED_)
	nsv += 1;
# else
	nsv += 2 * nelem + 2 * nriver;
# endif
#endif

#if defined(_CYCLES_)
	nsv += 2 * nelem + 4 * nriver;
#endif

#if defined(_FBR_)
	nsv += 2 * nelem;
#endif

	return nsv;
}



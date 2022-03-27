#include "pihm.h"
#if defined(_CVODE_CUDA)
#include "pihm.cuh"
#endif

#define TEC_HEADER          "VARIABLES = \"X\" \"Y\" \"Zmin\" \"Zmax\" \"h\""
#define WB_HEADER           "VARIABLES = \"TIME (s)\" \"SRC (m)\" \"SNK (m)\" \"STRG (m)\""
#define RIVER_TEC_HEADER2   "ZONE T = \"Water Depth River\""
#define RIVER_TEC_HEADER3   "StrandID=1, SolutionTime="
#define ELEM_TEC_HEADER3    "VARSHARELIST = ([1, 2, 3, 4]=1), CONNECTIVITYSHAREZONE = 1"

void StartupScreen(void)
{
#if defined(_CVODE_CUDA)
	PIHMprintf(VL_NORMAL, "\n");
	PIHMprintf(VL_NORMAL, "\t########    ####   ##     ##   ##     ##\n");
	PIHMprintf(VL_NORMAL, "\t##     ##    ##    ##     ##   ###   ###\n");
	PIHMprintf(VL_NORMAL, "\t##     ##    ##    ##     ##   #### ####\n");
	PIHMprintf(VL_NORMAL, "\t########     ##    #########   ## ### ##    ####  ##   ##  ## #     ##\n");
	PIHMprintf(VL_NORMAL, "\t##           ##    ##     ##   ##     ##   ##     ##   ##  #   #   #  #\n");
	PIHMprintf(VL_NORMAL, "\t##           ##    ##     ##   ##     ##   ##     ##   ##  #   #   ####\n");
	PIHMprintf(VL_NORMAL, "\t##          ####   ##     ##   ##     ##     ###    ####   ## #   #    #\n");

	PIHMprintf(VL_BRIEF, "\n\tThe Penn State Integrated Hydrologic Model in CUDA\n\n");
	PIHMprintf(VL_BRIEF, "\n\tDeveloped by Dr. LiJian at China University of Geosciences (Wuhan)\n\n");
#else
    PIHMprintf(VL_NORMAL, "\n");
    PIHMprintf(VL_NORMAL, "\t########    ####   ##     ##   ##     ##\n");
    PIHMprintf(VL_NORMAL, "\t##     ##    ##    ##     ##   ###   ###\n");
    PIHMprintf(VL_NORMAL, "\t##     ##    ##    ##     ##   #### ####\n");
    PIHMprintf(VL_NORMAL, "\t########     ##    #########   ## ### ##\n");
    PIHMprintf(VL_NORMAL, "\t##           ##    ##     ##   ##     ##\n");
    PIHMprintf(VL_NORMAL, "\t##           ##    ##     ##   ##     ##\n");
    PIHMprintf(VL_NORMAL, "\t##          ####   ##     ##   ##     ##\n");

    PIHMprintf(VL_BRIEF, "\n\tThe Penn State Integrated Hydrologic Model\n\n");
#endif

#if defined(_NOAH_)
    PIHMprintf(VL_BRIEF, "\t* Noah Land surface module turned on.\n");
#endif
#if defined(_RT_)
    PIHMprintf(VL_BRIEF, "\t* Reactive transport module turned on.\n");
#endif
#if defined(_BGC_)
    PIHMprintf(VL_BRIEF, "\t* Biogeochemistry module turned on.\n");
#endif
#if defined(_CYCLES_)
    PIHMprintf(VL_BRIEF, "\t* Crop module turned on.\n");
#endif
#if defined(_FBR_)
    PIHMprintf(VL_BRIEF, "\t* Fractured bedrock module turned on.\n");
#endif
#if defined(_OPENMP)
    PIHMprintf(VL_BRIEF, "\t* OpenMP (# of threads = %d).\n", nthreads);
#endif
    PIHMprintf(VL_BRIEF, "\n");

    if (1 == corr_mode)
    {
        PIHMprintf(VL_NORMAL,
            "\tSurface elevation correction mode turned on.\n");
    }
    if (1 == debug_mode)
    {
        PIHMprintf(VL_NORMAL,
            "\tDebug mode turned on.\n");
    }
    if (1 == tecplot)
    {
        PIHMprintf(VL_NORMAL,
            "\tTecplot output turned on.\n");
    }
    if (VL_BRIEF == verbose_mode)
    {
        PIHMprintf(VL_NORMAL,
            "\tBrief mode turned on.\n");
    }
    if (VL_VERBOSE == verbose_mode)
    {
        PIHMprintf(VL_NORMAL,
            "\tVerbose mode turned on.\n");
    }
    if (1 == append_mode)
    {
        PIHMprintf(VL_NORMAL,
            "\tAppend mode turned on.\n");
    }
}

void InitOutputFile(print_struct *print, const char *outputdir, int watbal,
    int ascii)
{
    char            ascii_fn[MAXSTRING];
    char            dat_fn[MAXSTRING];
    char            watbal_fn[MAXSTRING];
    char            perf_fn[MAXSTRING];
    int             i;
    char            mode[2];

    if (append_mode)
    {
        strcpy(mode, "a");
    }
    else
    {
        strcpy(mode, "w");
    }

    /* Initialize water balance file*/
    if (watbal)
    {
        sprintf(watbal_fn, "%s%s.watbal.txt", outputdir, project);
        print->watbal_file = fopen(watbal_fn, mode);
        CheckFile(print->watbal_file, watbal_fn);
    }

    /* Initialize cvode output files */
    if (debug_mode)
    {
        sprintf(perf_fn, "%s%s.cvode.log", outputdir, project);
        print->cvodeperf_file = fopen(perf_fn, mode);
        CheckFile(print->cvodeperf_file, perf_fn);
        /* Print header lines */
        fprintf(print->cvodeperf_file,
            "%-8s%-8s%-16s%-8s%-8s%-8s%-8s%-8s%-8s\n",
            "step", "cpu_dt", "cputime", "maxstep",
            "nsteps", "nevals", "niters", "ncfails", "nefails");
    }

    /*
     * Initialize model variable output files
     */
    for (i = 0; i < print->nprint; i++)
    {
        sprintf(dat_fn, "%s.dat", print->varctrl[i].name);
        print->varctrl[i].datfile = fopen(dat_fn, mode);
        CheckFile(print->varctrl[i].datfile, dat_fn);

        if (ascii)
        {
            sprintf(ascii_fn, "%s.txt", print->varctrl[i].name);
            print->varctrl[i].txtfile = fopen(ascii_fn, mode);
            CheckFile(print->varctrl[i].txtfile, ascii_fn);
        }
    }

    /* Tecplot files */
    if (tecplot)
    {
        for (i = 0; i < print->ntpprint; i++)
        {
            int             j;

            sprintf(dat_fn, "%s.plt", print->tp_varctrl[i].name);
            print->tp_varctrl[i].datfile = fopen(dat_fn, mode);
            CheckFile(print->tp_varctrl[i].datfile, dat_fn);

            if (print->tp_varctrl[i].intr == 0)
            {
                fprintf(print->tp_varctrl[i].datfile, "%s \n", TEC_HEADER);
                fprintf(print->tp_varctrl[i].datfile,
                    "ZONE T=\"%s\", N=%d, E=%d, DATAPACKING=%s, "
                    "SOLUTIONTIME=%lf, ZONETYPE=%s\n",
                    print->tp_varctrl[i].name, print->tp_varctrl[i].nnodes,
                    print->tp_varctrl[i].nvar, "POINT", 0.0000, "FETRIANGLE");

                for (j = 0; j < print->tp_varctrl[i].nnodes; j++)
                {
                    fprintf(print->tp_varctrl[i].datfile,
                        "%lf %lf %lf %lf %lf\n",
                        print->tp_varctrl[i].x[j], print->tp_varctrl[i].y[j],
                        print->tp_varctrl[i].zmin[j],
                        print->tp_varctrl[i].zmax[j],
                        0.000001);
                }
                for (j = 0; j < print->tp_varctrl[i].nvar; j++)
                {
                    fprintf(print->tp_varctrl[i].datfile, "%d %d %d\n",
                        print->tp_varctrl[i].node0[j],
                        print->tp_varctrl[i].node1[j],
                        print->tp_varctrl[i].node2[j]);
                }
            }
        }
    }
}

// 只更新待输出的变量，已经在map_output.c中获取了数组地址
void UpdPrintVar(varctrl_struct *varctrl, int nprint, int module_step)
{
	int             i;
#if defined(_OPENMP)
# pragma omp parallel for
#endif
	for (i = 0; i < nprint; i++)
	{
		int             j;

		if (varctrl[i].upd_intvl == module_step)
		{
			for (j = 0; j < varctrl[i].nvar; j++)
			{
				varctrl[i].buffer[j] += *varctrl[i].var[j];
			}

			varctrl[i].counter++;
		}
	}
}

#if defined(_CVODE_CUDA)
void PrintInit(const pihm_struct_d pihm_d, 
	           const char *outputdir, int t, int starttime, int endtime, int intvl)	            
#else
void PrintInit(const elem_struct *elem, const river_struct *river,
	           const char *outputdir, int t, int starttime, int endtime, int intvl)
#endif
{
	pihm_t_struct   pihm_time;

	pihm_time = PIHMTime(t);

#if defined(_CVODE_CUDA)
	int nelem = pihm_d->nelem;
	realtype *elem_h_ws_cmc;
	realtype *elem_h_ws_sneqv;
	realtype *elem_h_ws_surf;
	realtype *elem_h_ws_unsat;
	realtype *elem_h_ws_gw;

	elem_h_ws_cmc = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_ws_sneqv = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_ws_surf = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_ws_unsat = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_ws_gw = (realtype *)malloc(nelem * sizeof(realtype));

	cudaMemcpy(elem_h_ws_cmc, pihm_d->elem_d_ws_cmc, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_ws_sneqv, pihm_d->elem_d_ws_sneqv, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_ws_surf, pihm_d->elem_d_ws_surf, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_ws_unsat, pihm_d->elem_d_ws_unsat, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_ws_gw, pihm_d->elem_d_ws_gw, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);

	int nriver = pihm_d->nriver;
	realtype *river_h_ws_stage;
	realtype *river_h_ws_gw;

	river_h_ws_stage = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_ws_gw = (realtype *)malloc(nriver * sizeof(realtype));
	cudaMemcpy(river_h_ws_stage, pihm_d->river_d_ws_stage, nriver*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_ws_gw, pihm_d->river_d_ws_gw, nriver*sizeof(realtype), cudaMemcpyDeviceToHost);
#endif

	if (PrintNow(intvl, t - starttime, &pihm_time) || t == endtime)
	{
		FILE           *init_file;
		char            fn[MAXSTRING];
		int             i;

		sprintf(fn, "%s/restart/%s.%s.ic", outputdir, project,
			pihm_time.strshort);

		init_file = fopen(fn, "wb");
		CheckFile(init_file, fn);

		for (i = 0; i < nelem; i++)
		{
#if defined(_CVODE_CUDA)
			fwrite(&elem_h_ws_cmc[i], sizeof(realtype), 1, init_file);
			fwrite(&elem_h_ws_sneqv[i], sizeof(realtype), 1, init_file);
			fwrite(&elem_h_ws_surf[i], sizeof(realtype), 1, init_file);
			fwrite(&elem_h_ws_unsat[i], sizeof(realtype), 1, init_file);
			fwrite(&elem_h_ws_gw[i], sizeof(realtype), 1, init_file);
#else							
			fwrite(&elem[i].ws.cmc, sizeof(realtype), 1, init_file);
			fwrite(&elem[i].ws.sneqv, sizeof(realtype), 1, init_file);
			fwrite(&elem[i].ws.surf, sizeof(realtype), 1, init_file);
			fwrite(&elem[i].ws.unsat, sizeof(realtype), 1, init_file);
			fwrite(&elem[i].ws.gw, sizeof(realtype), 1, init_file);

#if defined(_FBR_)
			fwrite(&elem[i].ws.fbr_unsat, sizeof(realtype), 1, init_file);
			fwrite(&elem[i].ws.fbr_gw, sizeof(realtype), 1, init_file);
#endif
#if defined(_NOAH_)
			fwrite(&elem[i].es.t1, sizeof(realtype), 1, init_file);
			fwrite(&elem[i].ps.snowh, sizeof(realtype), 1, init_file);

			int             j;

			for (j = 0; j < MAXLYR; j++)
			{
				fwrite(&elem[i].es.stc[j], sizeof(realtype), 1, init_file);
			}
			for (j = 0; j < MAXLYR; j++)
			{
				fwrite(&elem[i].ws.smc[j], sizeof(realtype), 1, init_file);
			}
			for (j = 0; j < MAXLYR; j++)
			{
				fwrite(&elem[i].ws.sh2o[j], sizeof(realtype), 1, init_file);
			}
#endif

#endif
		}  // for (i = 0; i < nelem; i++)

		for (i = 0; i < nriver; i++)
		{
#if defined(_CVODE_CUDA)
			fwrite(&river_h_ws_stage[i], sizeof(realtype), 1, init_file);
			fwrite(&river_h_ws_gw[i], sizeof(realtype), 1, init_file);
#else
			fwrite(&river[i].ws.stage, sizeof(realtype), 1, init_file);
			fwrite(&river[i].ws.gw, sizeof(realtype), 1, init_file);
#endif
		}  // for (i = 0; i < nriver; i++)

		fflush(init_file);
		fclose(init_file);
	}

#if defined(_CVODE_CUDA)
	free(elem_h_ws_cmc);
	free(elem_h_ws_sneqv);
	free(elem_h_ws_surf);
	free(elem_h_ws_unsat);
	free(elem_h_ws_gw);

	free(river_h_ws_stage);
	free(river_h_ws_gw);
#endif
}

/* 水量平衡计算检测 */
#if defined(_CVODE_CUDA)
void PrintWaterBal(FILE *watbal_file, int t, int tstart, int dt,
	               const pihm_struct_d pihm_d)
#else
void PrintWaterBal(FILE *watbal_file, int t, int tstart, int dt,
	               const elem_struct *elem, const river_struct *river)
#endif
{
	int               i;
	realtype          tot_src = 0.0, tot_snk = 0.0, tot_strg = 0.0;  // 水源  汇流  滞留
	static realtype   tot_strg_prev = 0.0;
	static realtype   error = 0.0;

	if (t == tstart + dt)
	{
		fprintf(watbal_file, "%s\n", WB_HEADER);
	}

#if defined(_CVODE_CUDA)
#if defined(_NOT_UM_)
// 拷贝数组，从设备到主机, 没有使用统一内存时，必须显式地拷贝数据   
	/* 流域单元的水量平衡计算 */
	int nelem = pihm_d->nelem;
	realtype *elem_h_wf_prcp;
	realtype *elem_h_wf_edir;
	realtype *elem_h_wf_ett;
	realtype *elem_h_wf_ec;
	realtype *elem_h_topo_area;   // elem结构体中已有  
	realtype *elem_h_ws_cmc;
	realtype *elem_h_ws_sneqv;
	realtype *elem_h_ws_surf;
	realtype *elem_h_ws_unsat;
	realtype *elem_h_ws_gw;
	realtype *elem_h_soil_porosity; // elem结构体中已有

	elem_h_wf_prcp = (realtype *)malloc(nelem * sizeof(realtype));	
	elem_h_wf_edir = (realtype *)malloc(nelem * sizeof(realtype));	 
	elem_h_wf_ett = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_wf_ec = (realtype *)malloc(nelem * sizeof(realtype));
    elem_h_topo_area = (realtype *)malloc(nelem * sizeof(realtype));

	elem_h_ws_cmc = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_ws_sneqv = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_ws_surf = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_ws_unsat = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_ws_gw = (realtype *)malloc(nelem * sizeof(realtype));
	elem_h_soil_porosity = (realtype *)malloc(nelem * sizeof(realtype)); 

	cudaMemcpy(elem_h_wf_prcp, pihm_d->elem_d_wf_prcp, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_wf_edir, pihm_d->elem_d_wf_edir, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_wf_ett, pihm_d->elem_d_wf_ett, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_wf_ec, pihm_d->elem_d_wf_ec, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_topo_area, pihm_d->elem_d_topo_area, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);

	cudaMemcpy(elem_h_ws_cmc, pihm_d->elem_d_ws_cmc, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_ws_sneqv, pihm_d->elem_d_ws_sneqv, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_ws_surf, pihm_d->elem_d_ws_surf, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_ws_unsat, pihm_d->elem_d_ws_unsat, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_ws_gw, pihm_d->elem_d_ws_gw, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(elem_h_soil_porosity, pihm_d->elem_d_soil_porosity, nelem*sizeof(realtype), cudaMemcpyDeviceToHost);

	/* 河流单元的水量平衡计算 */
	int nriver = pihm_d->nriver;
	realtype *river_h_ws_stage;
	realtype *river_h_ws_gw;
	realtype *river_h_matl_porosity;
	realtype *river_h_topo_area;    // river 结构体中已有
	realtype *river_h_wf_rivflow1;

	river_h_ws_stage = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_ws_gw = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_matl_porosity = (realtype *)malloc(nriver * sizeof(realtype));
	river_h_topo_area = (realtype *)malloc(nriver * sizeof(realtype));

	// river_d_wf_rivflow[DOWN_CHANL2CHANL] 
	river_h_wf_rivflow1 = (realtype*)malloc(nriver * sizeof(realtype));  

	cudaMemcpy(river_h_ws_stage, pihm_d->river_d_ws_stage, nriver*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_ws_gw, pihm_d->river_d_ws_gw, nriver*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_matl_porosity, pihm_d->river_d_matl_porosity, nriver*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_topo_area, pihm_d->river_d_topo_area, nriver*sizeof(realtype), cudaMemcpyDeviceToHost);
	cudaMemcpy(river_h_wf_rivflow1, pihm_d->river_d_wf_rivflow[DOWN_CHANL2CHANL], nriver*sizeof(realtype), cudaMemcpyDeviceToHost);
#endif

#if 0
	// 检查数组从设备到主机拷贝中的问题
	cudaError_t error_cuda = cudaGetLastError();
	printf("CUDA error in PrintWaterBal: %s\n", cudaGetErrorString(error_cuda)); // no error
#endif
#endif  /* _CVODE_CUDA */

#if defined(_OPENMP)
# pragma omp parallel for
#endif
	for (i = 0; i < nelem; i++)
	{
		/* 检查 SRC 的正确性 */
#if defined(_CVODE_CUDA)
		tot_src += pihm_d->elem_d_wf_prcp[i] * pihm_d->elem_d_topo_area[i] * dt;
#if defined(_NOAH_)
		tot_src += (pihm_d->elem_d_wf_dew[i] + pihm_d->elem_d_wf_snomlt[i]) *
			        pihm_d->elem_d_topo_area[i] * dt;
#endif
#else
		tot_src += elem[i].wf.prcp * elem[i].topo.area * dt;
#if defined(_NOAH_)
		tot_src += (elem[i].wf.dew + elem[i].wf.snomlt) *
			        elem[i].topo.area * dt;
#endif
#endif

		/* 检查 SNK 的正确性 */
#if defined(_CVODE_CUDA)
		tot_snk += (pihm_d->elem_d_wf_edir[i] + pihm_d->elem_d_wf_ett[i] + pihm_d->elem_d_wf_ec[i]) *
			        pihm_d->elem_d_topo_area[i] * dt;
#if defined(_NOAH_)
		tot_snk += pihm_d->elem_d_wf_esnow[i] * pihm_d->elem_d_topo_area[i] * dt;
#endif
#else
		tot_snk += (elem[i].wf.edir + elem[i].wf.ett + elem[i].wf.ec) *
			        elem[i].topo.area * dt;
#if defined(_NOAH_)
		tot_snk += elem[i].wf.esnow * elem[i].topo.area * dt;
#endif
#endif

		/* 检查 STRG 的正确性 */
#if defined(_CVODE_CUDA)   
		tot_strg += (pihm_d->elem_d_ws_cmc[i] + pihm_d->elem_d_ws_sneqv[i] +
			         pihm_d->elem_d_ws_surf[i] +
					(pihm_d->elem_d_ws_unsat[i] + pihm_d->elem_d_ws_gw[i]) * 
					 pihm_d->elem_d_soil_porosity[i]) *
					 pihm_d->elem_d_topo_area[i];
#else
		tot_strg += (elem[i].ws.cmc + elem[i].ws.sneqv + elem[i].ws.surf +
			        (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity) *
			         elem[i].topo.area;
#endif
	}   // 流域单元

#if defined(_OPENMP)
# pragma omp parallel for
#endif
	for (i = 0; i < nriver; i++)
	{
		/* 检查 STRG 的正确性 */
#if defined(_CVODE_CUDA)
		tot_strg +=
			(pihm_d->river_d_ws_stage[i] + pihm_d->river_d_ws_gw[i] *
			 pihm_d->river_d_matl_porosity[i]) * pihm_d->river_d_topo_area[i];
#else
		tot_strg +=
			(river[i].ws.stage + river[i].ws.gw * river[i].matl.porosity) *
			 river[i].topo.area;
#endif
		/* 检查 SNK 的正确性 */
#if defined(_CVODE_CUDA)
		if (pihm_d->river_d_down[i] < 0)
		{
			tot_snk += pihm_d->river_d_wf_rivflow[DOWN_CHANL2CHANL][i] * dt;
		}
#else
		if (river[i].down < 0)
		{
			tot_snk += river[i].wf.rivflow[DOWN_CHANL2CHANL] * dt;
		}
#endif
	}

    /* 水量平衡误差检验 
	   CUDA计算的tot_snk的误差较大，其他各项误差较小.
	*/
	if (tot_strg_prev != 0.0)
	{
		error += tot_src - tot_snk - (tot_strg - tot_strg_prev);   // 累积的水量计算误差

		fprintf(watbal_file, "%d %lg %lg %lg %lg %lg\n", 
			t - tstart,tot_src, tot_snk, tot_strg - tot_strg_prev,   //前4项的CUDA与CPU计算结果几乎一致
			tot_src - tot_snk - (tot_strg - tot_strg_prev), error);  //但后2项的偏离逐渐增大?!
		fflush(watbal_file);
	}

	tot_strg_prev = tot_strg;  // 上一时刻的滞留水量

	/* 释放临时的主机数组 */
#if defined(_CVODE_CUDA)
#if defined(_NOT_UM_)
	free(elem_h_wf_prcp);
	free(elem_h_wf_edir);
	free(elem_h_wf_ett);
	free(elem_h_wf_ec);
	free(elem_h_topo_area);
	free(elem_h_ws_cmc);
	free(elem_h_ws_sneqv);
	free(elem_h_ws_surf);
	free(elem_h_ws_unsat);
	free(elem_h_ws_gw);
	free(elem_h_soil_porosity);

	free(river_h_ws_stage);
	free(river_h_ws_gw);
	free(river_h_matl_porosity);
	free(river_h_topo_area);
    free(river_h_wf_rivflow1);
#endif
#endif
}

// CPU上输出计算变量
void PrintData(varctrl_struct *varctrl, int nprint, int t, int lapse, int ascii)
{
	int             i;
	pihm_t_struct   pihm_time;

	pihm_time = PIHMTime(t);

#if defined(_OPENMP)
# pragma omp parallel for
#endif
	for (i = 0; i < nprint; i++)
	{
		int             j;
		realtype          outval;
		realtype          outtime;

		if (PrintNow(varctrl[i].intvl, lapse, &pihm_time))
		{
			if (ascii)
			{
				fprintf(varctrl[i].txtfile, "\"%s\"", pihm_time.str);
				for (j = 0; j < varctrl[i].nvar; j++)
				{
					if (varctrl[i].counter > 0)
					{
						fprintf(varctrl[i].txtfile, "\t%lf",
							varctrl[i].buffer[j] / (realtype)varctrl[i].counter);
					}
					else
					{
						fprintf(varctrl[i].txtfile, "\t%lf",
							varctrl[i].buffer[j]);
					}
				}
				fprintf(varctrl[i].txtfile, "\n");
				fflush(varctrl[i].txtfile);
			}

			outtime = (realtype)t;
			fwrite(&outtime, sizeof(realtype), 1, varctrl[i].datfile);
			for (j = 0; j < varctrl[i].nvar; j++)
			{
				if (varctrl[i].counter > 0)
				{
					outval = varctrl[i].buffer[j] / (realtype)varctrl[i].counter;
				}
				else
				{
					outval = varctrl[i].buffer[j];
				}
				fwrite(&outval, sizeof(realtype), 1, varctrl[i].datfile);

				varctrl[i].buffer[j] = 0.0;
			}
			varctrl[i].counter = 0;
			fflush(varctrl[i].datfile);
		}
	}
}

void PrintDataTecplot(varctrl_struct *varctrl, int nprint, int t, int lapse)
{
	int             i;
	pihm_t_struct   pihm_time;
	realtype         *hnodes;    /* h at nodes */
	int            *inodes;

	pihm_time = PIHMTime(t);

	for (i = 0; i < nprint; i++)
	{
		int             j;
		realtype          outval;
		realtype          outtime;

		if (PrintNow(varctrl[i].intvl, lapse, &pihm_time))
		{
			outtime = (realtype)t;

			if (varctrl[i].intr == RIVERVAR)
			{
				/*Print river files */
				fprintf(varctrl[i].datfile, "%s\n", TEC_HEADER);
				fprintf(varctrl[i].datfile, "%s\n", RIVER_TEC_HEADER2);
				fprintf(varctrl[i].datfile, "%s %d\n", RIVER_TEC_HEADER3, t);
				for (j = 0; j < varctrl[i].nvar; j++)
				{
					if (varctrl[i].counter > 0)
					{
						outval = varctrl[i].buffer[j] /
							(realtype)varctrl[i].counter;
					}
					else
					{
						outval = varctrl[i].buffer[j];
					}

					fprintf(varctrl[i].datfile, "%lf %lf %lf %lf %lf\n",
						varctrl[i].x[j], varctrl[i].y[j],
						varctrl[i].zmin[j], varctrl[i].zmax[j], outval);
					varctrl[i].buffer[j] = 0.0;
				}
			}
			else
			{
				/*Print element files */
				hnodes = (realtype *)calloc(varctrl[i].nnodes, sizeof(realtype));
				inodes = (int *)calloc(varctrl[i].nnodes, sizeof(int));
				for (j = 0; j < varctrl[i].nnodes; j++)
				{
					hnodes[j] = 0.0;
					inodes[j] = 0;
				}
				fprintf(varctrl[i].datfile,
					"ZONE T=\"%s\", N=%d, E=%d, DATAPACKING=%s, "
					"SOLUTIONTIME=%lf, ZONETYPE=%s\n",
					varctrl[i].name, varctrl[i].nnodes, varctrl[i].nvar,
					"POINT", outtime, "FETRIANGLE");
				fprintf(varctrl[i].datfile, "%s\n", ELEM_TEC_HEADER3);
				for (j = 0; j < varctrl[i].nvar; j++)
				{
					if (varctrl[i].counter > 0)
					{
						outval = varctrl[i].buffer[j] /
							(realtype)varctrl[i].counter;
					}
					else
					{
						outval = varctrl[i].buffer[j];
					}

					hnodes[varctrl[i].node0[j] - 1] += outval;
					hnodes[varctrl[i].node1[j] - 1] += outval;
					hnodes[varctrl[i].node2[j] - 1] += outval;
					inodes[varctrl[i].node0[j] - 1] += 1;
					inodes[varctrl[i].node1[j] - 1] += 1;
					inodes[varctrl[i].node2[j] - 1] += 1;
					varctrl[i].buffer[j] = 0.0;
				}
				for (j = 0; j < varctrl[i].nnodes; j++)
				{
					if (inodes[j] == 0)
					{
						fprintf(varctrl[i].datfile, "%8.6f\n", 0.0);
					}
					else
					{
						fprintf(varctrl[i].datfile, "%8.6f\n",
							hnodes[j] / inodes[j]);
					}
				}

				free(hnodes);
				free(inodes);
			}

			varctrl[i].counter = 0;
			fflush(varctrl[i].datfile);
		}
	}
}



void PrintPerf(void *cvode_mem, int t, int starttime, realtype cputime_dt,
    realtype cputime, realtype maxstep, FILE *perf_file)
{
    static realtype   dt;
    static long int nst0, nfe0, nni0, ncfn0, netf0;
    long int        nst, nfe, nni, ncfn, netf;
    int             cv_flag;

    cv_flag = CVodeGetNumSteps(cvode_mem, &nst);
    if (!CheckCVodeFlag(cv_flag))
    {
        PIHMexit(EXIT_FAILURE);
    }

    cv_flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    if (!CheckCVodeFlag(cv_flag))
    {
        PIHMexit(EXIT_FAILURE);
    }

    cv_flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    if (!CheckCVodeFlag(cv_flag))
    {
        PIHMexit(EXIT_FAILURE);
    }

    cv_flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    if (!CheckCVodeFlag(cv_flag))
    {
        PIHMexit(EXIT_FAILURE);
    }

    cv_flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    if (!CheckCVodeFlag(cv_flag))
    {
        PIHMexit(EXIT_FAILURE);
    }

    fprintf(perf_file, "%-8d%-8.3f%-16.3f%-8.2f",
        t - starttime, cputime_dt, cputime, maxstep);
    fprintf(perf_file, "%-8ld%-8ld%-8ld%-8ld%-8ld\n",
        nst - nst0, nni - nni0, nfe - nfe0, netf - netf0, ncfn - ncfn0);
    fflush(perf_file);

    dt = 0.0;

    nst0 = nst;
    nni0 = nni;
    nfe0 = nfe;
    netf0 = netf;
    ncfn0 = ncfn;

    dt += cputime_dt;
}


void PrintCVodeFinalStats(void *cvode_mem)
{
    int             cv_flag;
    long int        nst;
    long int        nfe;
    long int        netf;
    long int        nni;
    long int        ncfn;

    cv_flag = CVodeGetNumSteps(cvode_mem, &nst);
    if (!CheckCVodeFlag(cv_flag))
    {
        PIHMexit(EXIT_FAILURE);
    }

    cv_flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    if (!CheckCVodeFlag(cv_flag))
    {
        PIHMexit(EXIT_FAILURE);
    }

    cv_flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    if (!CheckCVodeFlag(cv_flag))
    {
        PIHMexit(EXIT_FAILURE);
    }

    cv_flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    if (!CheckCVodeFlag(cv_flag))
    {
        PIHMexit(EXIT_FAILURE);
    }

    cv_flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    if (!CheckCVodeFlag(cv_flag))
    {
        PIHMexit(EXIT_FAILURE);
    }

    PIHMprintf(VL_NORMAL, "\n");
    PIHMprintf(VL_NORMAL,
        "num of steps = %-6ld num of rhs evals = %-6ld\n", nst, nfe);
    PIHMprintf(VL_NORMAL,
        "num of nonlin solv iters = %-6ld "
        "num of nonlin solv conv fails = %-6ld "
        "num of err test fails = %-6ld\n",
        nni, ncfn, netf);
}

int PrintNow(int intvl, int lapse, const pihm_t_struct *pihm_time)
{
    int             print = 0;

    if (intvl != 0)
    {
        switch (intvl)
        {
            case YEARLY_OUTPUT:
                if (pihm_time->month == 1 && pihm_time->day == 1 &&
                    pihm_time->hour == 0 && pihm_time->minute == 0)
                {
                    print = 1;
                }
                break;
            case MONTHLY_OUTPUT:
                if (pihm_time->day == 1 && pihm_time->hour == 0 &&
                    pihm_time->minute == 0)
                {
                    print = 1;
                }
                break;
            case DAILY_OUTPUT:
                if (pihm_time->hour == 0 && pihm_time->minute == 0)
                {
                    print = 1;
                }
                break;
            case HOURLY_OUTPUT:
                if (pihm_time->minute == 0)
                {
                    print = 1;
                }
                break;
            default:
                if (lapse % intvl == 0)
                {
                    print = 1;
                }
        }
    }

    return print;
}

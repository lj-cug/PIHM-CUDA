#include "pihm.h"
#include "pihm.cuh"

/* 分配状态变量的CUDA空间 */
cudaError_t Initialize_cuda(pihm_struct_d pihm_d)
{
	int  i;
	cudaError_t cudaStatus;
	unsigned int  flags = cudaMemAttachGlobal;   // or cudaMemAttachHost

	// 全部分配为统一内存, "仅分配CUDA核函数(DOE)计算中需要用到的一些数组或变量" !!!
	/*
	 * Initialize PIHM structure on GPU, now we use SoA-style struct
	 */
	// --------------- Basin Element ---------------------
	/* Element geometrics */
	for (i = 0; i < NUM_EDGE; i++){
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_node[i]), nelem * sizeof(realtype), flags);
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_nabr[i]), nelem * sizeof(realtype), flags);
	}
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ind), nelem * sizeof(realtype), flags);

	/* Element attribute */
	for (i = 0; i < NUM_EDGE; i++)
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_attrib_bc_type[i]), nelem * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_attrib_soil_type), nelem * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_attrib_lc_type), nelem * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_attrib_meteo_type), nelem * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_attrib_lai_type), nelem * sizeof(int), flags);

	/* Topography parameters */
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_area), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_x), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_y), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_zmin), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_zmax), nelem * sizeof(realtype), flags);
#if defined(_FBR_)
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_zbed), nelem * sizeof(realtype), flags);
#endif	
	for (i = 0; i < NUM_EDGE; i++){
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_edge[i]), nelem * sizeof(realtype), flags);
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_nabrdist[i]), nelem * sizeof(realtype), flags);
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_nabr_x[i]), nelem * sizeof(realtype), flags);
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_nabr_y[i]), nelem * sizeof(realtype), flags);
	}
#if defined(_NOAH_)
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_slope), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_aspect), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_svf), nelem * sizeof(realtype), flags);
	/*
	for (i = 0; i < 36; i++){
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_h_phi[i]), nelem * sizeof(realtype), flags);
    }
	*/
#endif	
#if defined(_RT_)
	for (i = 0; i < NUM_EDGE; i++){
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_topo_areasub[i]), nelem * sizeof(realtype), flags);
	}
#endif

	/* Soil parameters  soil */
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_depth), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_ksath), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_ksatv), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_kinfv), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_dinf), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_alpha), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_beta), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_porosity), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_smcmax), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_smcmin), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_smcwlt), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_smcref), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_dmac), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_kmach), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_kmacv), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_areafv), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_areafh), nelem * sizeof(realtype), flags);
#if defined(_NOAH_)	
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_csoil), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_quartz), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_soil_smcdry), nelem * sizeof(realtype), flags);
#endif	
#if defined(_CYCLES_)

#endif
#if defined(_FBR_)

#endif
	/* Land cover parameters */
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_shdfac), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_shdmin), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_shdmax), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_laimin), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_laimax), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_snup), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_cfactr), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_emissmax), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_emissmin), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_albedomax), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_albedomin), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_z0max), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_z0min), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_rough), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_cmcfactr), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_bare), nelem * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_lc_isurban), nelem * sizeof(int), flags);

	/* Ecophysiological parameters */
#if !defined(_CYCLES_)	
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_epc_rsmin), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_epc_rgl), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_epc_hs), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_epc_topt), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_epc_rsmax), nelem * sizeof(realtype), flags);
#if defined(_BGC_)


#endif	
#endif	/* Ecophysiological parameters */

	/* Physical states  ps */
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_rzd), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_rc), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_pc), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_proj_lai), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_rcs), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_rct), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_rcq), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_rcsoil), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_albedo), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_zlvl), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_zlvl_wind), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_sfcspd), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_rh), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_sfcprs), nelem * sizeof(realtype), flags);
#if defined(_NOAH_)
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_nroot), nelem * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_nsoil), nelem * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_nmacd), nelem * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_nwtbl), nelem * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_alb), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_snoalb), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_soilw), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_frzk), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_frzx), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_czil), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_emissi), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_ch), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_cm), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_rch), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_z0), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_fcr), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_salp), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_sbeta), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_lvcoef), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_snotime1), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_ribb), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_beta), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_sncovr), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_q1), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_q2), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_ffrozp), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_z0brd), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_embrd), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_q2sat), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_q2d), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_dqsdt2), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_sndens), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_snowh), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_sncond), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_rr), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_epsca), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_eta_kinematic), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_zbot), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_tbot), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ps_gwet), nelem * sizeof(realtype), flags);

//  分配临时的 "一级指针数组" 变量的空间
	for (i = 0; i < MAXLYR; i++){
		cudaMallocManaged(&(pihm_d->elem_d_ps_rtdis[i]), nelem * sizeof(realtype));
		cudaMallocManaged(&(pihm_d->elem_d_ps_sldpth[i]), nelem * sizeof(realtype));
		cudaMallocManaged(&(pihm_d->elem_d_ps_zsoil[i]), nelem * sizeof(realtype));
		cudaMallocManaged(&(pihm_d->elem_d_ps_satdpth[i]), nelem * sizeof(realtype));
	}
#endif
#if defined(_CYCLES_)

#endif
#if defined(_BGC_)

#endif

	/* Water states  ws */
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws_surf), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws_unsat), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws_gw), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws_sneqv), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws_cmcmax), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws_cmc), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws_surfh), nelem * sizeof(realtype), flags);
#if defined(_FBR_)

#endif
#if defined(_NOAH_)
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws_soilm), nelem * sizeof(realtype), flags);

	for (i = 0; i < MAXLYR; i++){
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws_smc[i]),
			nelem * sizeof(realtype), flags);
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws_sh2o[i]),
			nelem * sizeof(realtype), flags);
	}
#endif
#if defined(_CYCLES_)

#endif	
	/* Initial Water states   ws0 仅IS_ET_SM LSM模型需要 */
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws0_surf), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws0_unsat), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws0_gw), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws0_sneqv), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws0_cmcmax), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws0_cmc), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ws0_surfh), nelem * sizeof(realtype), flags);

	/* Water fluxes  wf */
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_prcp), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_pcpdrp), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_infil), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_rechg), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_drip), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_edir), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_ett), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_ec), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_etp), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_eta), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_edir_surf), nelem * sizeof(realtype), flags);
	for (i = 0; i < NUM_EDGE; i++){
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_ovlflow[i]), nelem * sizeof(realtype), flags);
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_subsurf[i]), nelem * sizeof(realtype), flags);
	}
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_edir_unsat), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_edir_gw), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_ett_unsat), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_ett_gw), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_esnow), nelem * sizeof(realtype), flags);
#if defined(_FBR_)	


#endif
#if defined(_NOAH_)
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_runoff2), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_runoff3), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_dew), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_snomlt), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_etns), nelem * sizeof(realtype), flags);

	for (i = 0; i < MAXLYR; i++){
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_et[i]),
			nelem * sizeof(realtype), flags);
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_runoff2_lyr[i]),
			nelem * sizeof(realtype), flags);
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_smflxv[i]),
			nelem * sizeof(realtype), flags);
	}

	/*  三级指针 ?
	for (i = 0; i < NUM_EDGE; i++)
		for(int j = 0; j < MAXLYR; j++)
			// *elem_d_wf_smflxh[NUM_EDGE][MAXLYR];
			cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_wf_smflxh[i][j]), 
			                                 nelem * sizeof(realtype), flags);
	*/
#endif
#if defined(_CYCLES_)


#endif      /* Water fluxes */   

	/* Energy states */
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_es_sfctmp), nelem * sizeof(realtype), flags);
#if defined(_NOAH_)
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_es_t1), nelem * sizeof(realtype), flags);	
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_es_th2), nelem * sizeof(realtype), flags);
	for (i = 0; i < MAXLYR; i++){
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_es_stc[i]), nelem * sizeof(realtype), flags);
	}
#endif     /* Energy states */

	/* Energy fluxes */
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_soldn), nelem * sizeof(realtype), flags);
#if defined(_NOAH_)
	for (i = 0; i < MAXLYR; i++){
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_et[i]), nelem * sizeof(realtype), flags);
	}
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_solnet), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_etp), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_ssoil), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_eta), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_sheat), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_fdown), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_lwdn), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_ec), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_edir), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_ett), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_esnow), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_soldir), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_soldif), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_longwave), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_flx1), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_flx2), nelem * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_ef_flx3), nelem * sizeof(realtype), flags);
#endif
#if defined(_BGC_)


#endif

	/* Boundary conditions */
	for (i = 0; i < NUM_EDGE; i++){
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_bc_head[i]), nelem * sizeof(realtype), flags);
		cudaStatus = cudaMallocManaged(&(pihm_d->elem_d_bc_flux[i]), nelem * sizeof(realtype), flags);
	}

	// Check cudaMallocManaged Status
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Element: cudaMallocManaged failed! Check Initialize_CUDA ?\n");
		exit(cudaStatus);
	}
//------------- River Segments -------------------
	/* River attribute */
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_attrib_riverbc_type), nriver * sizeof(realtype), flags);

	/* River topography parameters */
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_topo_area), nriver * sizeof(realtype), flags);	
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_topo_x), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_topo_y), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_topo_zmin), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_topo_zmax), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_topo_zbed), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_topo_node_zmax), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_topo_dist_left), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_topo_dist_right), nriver * sizeof(realtype), flags);

	/* River water states */
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_ws_stage), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_ws_gw), nriver * sizeof(realtype), flags);
	
	/* River water fluxes */
	for (i = 0; i < NUM_RIVFLX; i++)
		cudaStatus = cudaMallocManaged(&(pihm_d->river_d_wf_rivflow[i]), nriver * sizeof(realtype), flags);

	/* River shape parameters */
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_shp_depth), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_shp_intrpl_ord), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_shp_coeff), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_shp_length), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_shp_width), nriver * sizeof(realtype), flags);

	/* matl_struct_d */
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_matl_rough), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_matl_cwr), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_matl_ksath), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_matl_ksatv), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_matl_bedthick), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_matl_porosity), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_matl_smcmin), nriver * sizeof(realtype), flags);
#if defined(_CYCLES_)

#endif

	/* River boundary conditions */
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_bc_head), nriver * sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_bc_flux), nriver * sizeof(realtype), flags);

	/* River initial conditions */
	// 初始条件需要开辟CUDA空间吗?

#if defined(_BGC_) && !defined(_LUMPED_) && !defined(_LEACHING_)


#endif

	/* River structure */
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_ind), nriver * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_leftele), nriver * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_rightele), nriver * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_fromnode), nriver * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_tonode), nriver * sizeof(int), flags);
	cudaStatus = cudaMallocManaged(&(pihm_d->river_d_down), nriver * sizeof(int), flags);


	// check the CUDA allocation status
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "River: cudaMallocManaged failed! Check Initialize_CUDA ?\n");
		exit(cudaStatus);
	}	

	// Transfer other parameters into GPU, like calib_struct, ctrl_struct etc.
	/* Global calibration coefficients */
	/*
	cudaStatus = cudaMalloc((void **)&cal_d_ec, sizeof(realtype));
	cudaStatus = cudaMalloc((void **)&cal_d_ett, sizeof(realtype));
	cudaStatus = cudaMalloc((void **)&cal_d_edir, sizeof(realtype));
	*/
	cudaStatus = cudaMallocManaged((void **)&pihm_d->ec, sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged((void **)&pihm_d->ett, sizeof(realtype), flags);
	cudaStatus = cudaMallocManaged((void **)&pihm_d->edir, sizeof(realtype), flags);

	/* Model control parameters */
	/*
	cudaStatus = cudaMalloc((void **)&ctrl_d_stepsize, sizeof(int));
	cudaStatus = cudaMalloc((void **)&ctrl_d_surf_mode, sizeof(int));
	cudaStatus = cudaMalloc((void **)&ctrl_d_riv_mode, sizeof(int));
	*/
	cudaStatus = cudaMallocManaged((void **)&pihm_d->stepsize, sizeof(int), flags);
	cudaStatus = cudaMallocManaged((void **)&pihm_d->surf_mode, sizeof(int), flags);
	cudaStatus = cudaMallocManaged((void **)&pihm_d->riv_mode, sizeof(int), flags);

	return(cudaStatus);
}


/* 完成GPU上SoA风格的结构体空间分配后，将CPU上AoS风格的结构体内容，
* 拷贝到GPU设备上, 考虑使用OpenMP多线程拷贝
*/
void AoS_to_SoA(pihm_struct pihm, pihm_struct_d pihm_d)
{
	int   i, tid;

	// 根据结构体中的变量顺序，开始转移数值, AoS --> SoA

	/* Variables related with basin elements */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
	for (tid = 0; tid < nelem; tid++){

		/* Element geometrics in GPU */
		for (i = 0; i < NUM_EDGE; i++){
			pihm_d->elem_d_node[i][tid] = pihm->elem[tid].node[i];
			pihm_d->elem_d_nabr[i][tid] = pihm->elem[tid].nabr[i];
		}
		pihm_d->elem_d_ind[tid] = pihm->elem[tid].ind;

		/* Element attribute */
		for (i = 0; i < NUM_EDGE; i++)
			pihm_d->elem_d_attrib_bc_type[i][tid] = pihm->elem[tid].attrib.bc_type[i];
		pihm_d->elem_d_attrib_soil_type[tid] = pihm->elem[tid].attrib.soil_type;
		pihm_d->elem_d_attrib_lc_type[tid] = pihm->elem[tid].attrib.lc_type;
		pihm_d->elem_d_attrib_meteo_type[tid] = pihm->elem[tid].attrib.meteo_type;
		pihm_d->elem_d_attrib_lai_type[tid] = pihm->elem[tid].attrib.lai_type;

		/* Topography parameters */
		pihm_d->elem_d_topo_area[tid] = pihm->elem[tid].topo.area;
		pihm_d->elem_d_topo_x[tid] = pihm->elem[tid].topo.x;
		pihm_d->elem_d_topo_y[tid] = pihm->elem[tid].topo.y;
		pihm_d->elem_d_topo_zmin[tid] = pihm->elem[tid].topo.zmin;
		pihm_d->elem_d_topo_zmax[tid] = pihm->elem[tid].topo.zmax;
		for (i = 0; i < NUM_EDGE; i++){
			pihm_d->elem_d_topo_edge[i][tid] = pihm->elem[tid].topo.edge[i];
			pihm_d->elem_d_topo_nabrdist[i][tid] = pihm->elem[tid].topo.nabrdist[i];
			pihm_d->elem_d_topo_nabr_x[i][tid] = pihm->elem[tid].topo.nabr_x[i];
			pihm_d->elem_d_topo_nabr_y[i][tid] = pihm->elem[tid].topo.nabr_y[i];
		}
#if defined(_NOAH_)
		pihm_d->elem_d_topo_slope[tid] = pihm->elem[tid].topo.slope;
		pihm_d->elem_d_topo_aspect[tid] = pihm->elem[tid].topo.aspect;
		pihm_d->elem_d_topo_svf[tid] = pihm->elem[tid].topo.svf;
		/*
		for (i = 0; i < 36; i++){
		pihm_d->elem_d_topo_h_phi[tid] = pihm->elem[tid].topo.h_phi;
		}
		*/
#endif	

		/* Soil parameters */
		pihm_d->elem_d_soil_depth[tid] = pihm->elem[tid].soil.depth;
		pihm_d->elem_d_soil_ksath[tid] = pihm->elem[tid].soil.ksath;
		pihm_d->elem_d_soil_ksatv[tid] = pihm->elem[tid].soil.ksatv;
		pihm_d->elem_d_soil_kinfv[tid] = pihm->elem[tid].soil.kinfv;
		pihm_d->elem_d_soil_dinf[tid] = pihm->elem[tid].soil.dinf;
		pihm_d->elem_d_soil_alpha[tid] = pihm->elem[tid].soil.alpha;
		pihm_d->elem_d_soil_beta[tid] = pihm->elem[tid].soil.beta;
		pihm_d->elem_d_soil_porosity[tid] = pihm->elem[tid].soil.porosity;
		pihm_d->elem_d_soil_smcmax[tid] = pihm->elem[tid].soil.smcmax;
		pihm_d->elem_d_soil_smcmin[tid] = pihm->elem[tid].soil.smcmin;
		pihm_d->elem_d_soil_smcwlt[tid] = pihm->elem[tid].soil.smcwlt;
		pihm_d->elem_d_soil_smcref[tid] = pihm->elem[tid].soil.smcref;
		pihm_d->elem_d_soil_dmac[tid] = pihm->elem[tid].soil.dmac;
		pihm_d->elem_d_soil_kmach[tid] = pihm->elem[tid].soil.kmach;
		pihm_d->elem_d_soil_kmacv[tid] = pihm->elem[tid].soil.kmacv;
		pihm_d->elem_d_soil_areafv[tid] = pihm->elem[tid].soil.areafv;
		pihm_d->elem_d_soil_areafh[tid] = pihm->elem[tid].soil.areafh;
#if defined(_NOAH_)	
		pihm_d->elem_d_soil_csoil[tid] = pihm->elem[tid].soil.csoil;
		pihm_d->elem_d_soil_quartz[tid] = pihm->elem[tid].soil.quartz;
		pihm_d->elem_d_soil_smcdry[tid] = pihm->elem[tid].soil.smcdry;
#endif	

		/* Land cover parameters */
		pihm_d->elem_d_lc_shdfac[tid] = pihm->elem[tid].lc.shdfac;
		pihm_d->elem_d_lc_shdmin[tid] = pihm->elem[tid].lc.shdmin;
		pihm_d->elem_d_lc_shdmax[tid] = pihm->elem[tid].lc.shdmax;
		pihm_d->elem_d_lc_laimin[tid] = pihm->elem[tid].lc.laimin;
		pihm_d->elem_d_lc_laimax[tid] = pihm->elem[tid].lc.laimax;
		pihm_d->elem_d_lc_snup[tid] = pihm->elem[tid].lc.snup;
		pihm_d->elem_d_lc_cfactr[tid] = pihm->elem[tid].lc.cfactr;
		pihm_d->elem_d_lc_emissmax[tid] = pihm->elem[tid].lc.emissmax;
		pihm_d->elem_d_lc_emissmin[tid] = pihm->elem[tid].lc.emissmin;
		pihm_d->elem_d_lc_albedomax[tid] = pihm->elem[tid].lc.albedomax;
		pihm_d->elem_d_lc_albedomin[tid] = pihm->elem[tid].lc.albedomin;
		pihm_d->elem_d_lc_z0max[tid] = pihm->elem[tid].lc.z0max;
		pihm_d->elem_d_lc_z0min[tid] = pihm->elem[tid].lc.z0min;
		pihm_d->elem_d_lc_rough[tid] = pihm->elem[tid].lc.rough;
		pihm_d->elem_d_lc_cmcfactr[tid] = pihm->elem[tid].lc.cmcfactr;
		pihm_d->elem_d_lc_bare[tid] = pihm->elem[tid].lc.bare;
		pihm_d->elem_d_lc_isurban[tid] = pihm->elem[tid].lc.isurban;

		/* Ecophysiological parameters */
		pihm_d->elem_d_epc_rsmin[tid] = pihm->elem[tid].epc.rsmin;
		pihm_d->elem_d_epc_rgl[tid] = pihm->elem[tid].epc.rgl;
		pihm_d->elem_d_epc_hs[tid] = pihm->elem[tid].epc.hs;
		pihm_d->elem_d_epc_topt[tid] = pihm->elem[tid].epc.topt;
		pihm_d->elem_d_epc_rsmax[tid] = pihm->elem[tid].epc.rsmax;

		/* Physical states  ps */
		pihm_d->elem_d_ps_rzd[tid] = pihm->elem[tid].ps.rzd;
		pihm_d->elem_d_ps_rc[tid] = pihm->elem[tid].ps.rc;
		pihm_d->elem_d_ps_pc[tid] = pihm->elem[tid].ps.pc;
		pihm_d->elem_d_ps_proj_lai[tid] = pihm->elem[tid].ps.proj_lai;
		pihm_d->elem_d_ps_rcs[tid] = pihm->elem[tid].ps.rcs;
		pihm_d->elem_d_ps_rct[tid] = pihm->elem[tid].ps.rct;
		pihm_d->elem_d_ps_rcq[tid] = pihm->elem[tid].ps.rcq;
		pihm_d->elem_d_ps_rcsoil[tid] = pihm->elem[tid].ps.rcsoil;
		pihm_d->elem_d_ps_albedo[tid] = pihm->elem[tid].ps.albedo;
		pihm_d->elem_d_ps_zlvl[tid] = pihm->elem[tid].ps.zlvl;
		pihm_d->elem_d_ps_zlvl_wind[tid] = pihm->elem[tid].ps.zlvl_wind;
		pihm_d->elem_d_ps_sfcspd[tid] = pihm->elem[tid].ps.sfcspd;
		pihm_d->elem_d_ps_rh[tid] = pihm->elem[tid].ps.rh;
		pihm_d->elem_d_ps_sfcprs[tid] = pihm->elem[tid].ps.sfcprs;
#if defined(_NOAH_)
		pihm_d->elem_d_ps_nroot[tid] = pihm->elem[tid].ps.nroot;
		pihm_d->elem_d_ps_nsoil[tid] = pihm->elem[tid].ps.nsoil;
		pihm_d->elem_d_ps_nmacd[tid] = pihm->elem[tid].ps.nmacd;
		pihm_d->elem_d_ps_nwtbl[tid] = pihm->elem[tid].ps.nwtbl;
		pihm_d->elem_d_ps_alb[tid] = pihm->elem[tid].ps.alb;
		pihm_d->elem_d_ps_snoalb[tid] = pihm->elem[tid].ps.snoalb;
		pihm_d->elem_d_ps_soilw[tid] = pihm->elem[tid].ps.soilw;
		pihm_d->elem_d_ps_frzk[tid] = pihm->elem[tid].ps.frzk;
		pihm_d->elem_d_ps_frzx[tid] = pihm->elem[tid].ps.frzx;
		pihm_d->elem_d_ps_czil[tid] = pihm->elem[tid].ps.czil;
		pihm_d->elem_d_ps_emissi[tid] = pihm->elem[tid].ps.emissi;
		pihm_d->elem_d_ps_ch[tid] = pihm->elem[tid].ps.ch;
		pihm_d->elem_d_ps_cm[tid] = pihm->elem[tid].ps.cm;
		pihm_d->elem_d_ps_rch[tid] = pihm->elem[tid].ps.rch;
		pihm_d->elem_d_ps_z0[tid] = pihm->elem[tid].ps.z0;
		pihm_d->elem_d_ps_fcr[tid] = pihm->elem[tid].ps.fcr;
		pihm_d->elem_d_ps_salp[tid] = pihm->elem[tid].ps.salp;
		pihm_d->elem_d_ps_sbeta[tid] = pihm->elem[tid].ps.sbeta;
		pihm_d->elem_d_ps_lvcoef[tid] = pihm->elem[tid].ps.lvcoef;
		pihm_d->elem_d_ps_snotime1[tid] = pihm->elem[tid].ps.snotime1;
		pihm_d->elem_d_ps_ribb[tid] = pihm->elem[tid].ps.ribb;
		pihm_d->elem_d_ps_beta[tid] = pihm->elem[tid].ps.beta;
		pihm_d->elem_d_ps_sncovr[tid] = pihm->elem[tid].ps.sncovr;
		pihm_d->elem_d_ps_q1[tid] = pihm->elem[tid].ps.q1;
		pihm_d->elem_d_ps_q2[tid] = pihm->elem[tid].ps.q2;
		pihm_d->elem_d_ps_ffrozp[tid] = pihm->elem[tid].ps.ffrozp;
		pihm_d->elem_d_ps_z0brd[tid] = pihm->elem[tid].ps.z0brd;
		pihm_d->elem_d_ps_embrd[tid] = pihm->elem[tid].ps.embrd;
		pihm_d->elem_d_ps_q2sat[tid] = pihm->elem[tid].ps.q2sat;
		pihm_d->elem_d_ps_q2d[tid] = pihm->elem[tid].ps.q2d;
		pihm_d->elem_d_ps_dqsdt2[tid] = pihm->elem[tid].ps.dqsdt2;
		pihm_d->elem_d_ps_sndens[tid] = pihm->elem[tid].ps.sndens;
		pihm_d->elem_d_ps_snowh[tid] = pihm->elem[tid].ps.snowh;
		pihm_d->elem_d_ps_sncond[tid] = pihm->elem[tid].ps.sncond;
		pihm_d->elem_d_ps_rr[tid] = pihm->elem[tid].ps.rr;
		pihm_d->elem_d_ps_epsca[tid] = pihm->elem[tid].ps.epsca;
		pihm_d->elem_d_ps_eta_kinematic[tid] = pihm->elem[tid].ps.eta_kinematic;
		pihm_d->elem_d_ps_zbot[tid] = pihm->elem[tid].ps.zbot;
		pihm_d->elem_d_ps_tbot[tid] = pihm->elem[tid].ps.tbot;
		pihm_d->elem_d_ps_gwet[tid] = pihm->elem[tid].ps.gwet;

		// cuda copy data from Host To Device, AoS2SoA
		//printf("Transfer PS data of Noah from CPU to GPU.\n");
		for (i = 0; i < MAXLYR; i++){
			pihm_d->elem_d_ps_rtdis[i][tid] = pihm->elem[tid].ps.rtdis[i];
			pihm_d->elem_d_ps_sldpth[i][tid] = pihm->elem[tid].ps.sldpth[i];
			pihm_d->elem_d_ps_zsoil[i][tid] = pihm->elem[tid].ps.zsoil[i];
			pihm_d->elem_d_ps_satdpth[i][tid] = pihm->elem[tid].ps.satdpth[i];
		}
		//printf("Finish PS data of Noah transferring from CPU to GPU.\n");
#endif

		/* Water states  ws */
		pihm_d->elem_d_ws_surf[tid] = pihm->elem[tid].ws.surf;
		pihm_d->elem_d_ws_unsat[tid] = pihm->elem[tid].ws.unsat;
		pihm_d->elem_d_ws_gw[tid] = pihm->elem[tid].ws.gw;
		pihm_d->elem_d_ws_sneqv[tid] = pihm->elem[tid].ws.sneqv;
		pihm_d->elem_d_ws_cmcmax[tid] = pihm->elem[tid].ws.cmcmax;
		pihm_d->elem_d_ws_cmc[tid] = pihm->elem[tid].ws.cmc;
		pihm_d->elem_d_ws_surfh[tid] = pihm->elem[tid].ws.surfh;
#if defined(_NOAH_)
		pihm_d->elem_d_ws_soilm[tid] = pihm->elem[tid].ws.soilm;

		for (i = 0; i < MAXLYR; i++){
			pihm_d->elem_d_ws_smc[i][tid] = pihm->elem[tid].ws.smc[i];
			pihm_d->elem_d_ws_sh2o[i][tid] = pihm->elem[tid].ws.sh2o[i];
		}
#endif

		/* Initial Water states  ws0 */
		pihm_d->elem_d_ws0_surf[tid] = pihm->elem[tid].ws0.surf;
		pihm_d->elem_d_ws0_unsat[tid] = pihm->elem[tid].ws0.unsat;
		pihm_d->elem_d_ws0_gw[tid] = pihm->elem[tid].ws0.gw;
		pihm_d->elem_d_ws0_sneqv[tid] = pihm->elem[tid].ws0.sneqv;
		pihm_d->elem_d_ws0_cmcmax[tid] = pihm->elem[tid].ws0.cmcmax;
		pihm_d->elem_d_ws0_cmc[tid] = pihm->elem[tid].ws0.cmc;
		pihm_d->elem_d_ws0_surfh[tid] = pihm->elem[tid].ws0.surfh;

		/* Water fluxes  wf */
		for (i = 0; i < NUM_EDGE; i++){
			pihm_d->elem_d_wf_ovlflow[i][tid] = pihm->elem[tid].wf.ovlflow[i];
			pihm_d->elem_d_wf_subsurf[i][tid] = pihm->elem[tid].wf.subsurf[i];
		}
		pihm_d->elem_d_wf_prcp[tid] = pihm->elem[tid].wf.prcp;
		pihm_d->elem_d_wf_pcpdrp[tid] = pihm->elem[tid].wf.pcpdrp;
		pihm_d->elem_d_wf_infil[tid] = pihm->elem[tid].wf.infil;
		pihm_d->elem_d_wf_rechg[tid] = pihm->elem[tid].wf.rechg;
		pihm_d->elem_d_wf_drip[tid] = pihm->elem[tid].wf.drip;
		pihm_d->elem_d_wf_edir[tid] = pihm->elem[tid].wf.edir;
		pihm_d->elem_d_wf_ett[tid] = pihm->elem[tid].wf.ett;
		pihm_d->elem_d_wf_ec[tid] = pihm->elem[tid].wf.ec;
		pihm_d->elem_d_wf_etp[tid] = pihm->elem[tid].wf.etp;
		pihm_d->elem_d_wf_eta[tid] = pihm->elem[tid].wf.eta;
		pihm_d->elem_d_wf_edir_surf[tid] = pihm->elem[tid].wf.edir_surf;
		pihm_d->elem_d_wf_edir_unsat[tid] = pihm->elem[tid].wf.edir_unsat;
		pihm_d->elem_d_wf_edir_gw[tid] = pihm->elem[tid].wf.edir_gw;
		pihm_d->elem_d_wf_ett_unsat[tid] = pihm->elem[tid].wf.ett_unsat;
		pihm_d->elem_d_wf_ett_gw[tid] = pihm->elem[tid].wf.ett_gw;
		pihm_d->elem_d_wf_esnow[tid] = pihm->elem[tid].wf.esnow;
#if defined(_NOAH_)
		pihm_d->elem_d_wf_runoff2[tid] = pihm->elem[tid].wf.esnow;
		pihm_d->elem_d_wf_runoff3[tid] = pihm->elem[tid].wf.esnow;
		pihm_d->elem_d_wf_dew[tid] = pihm->elem[tid].wf.esnow;
		pihm_d->elem_d_wf_snomlt[tid] = pihm->elem[tid].wf.esnow;
		pihm_d->elem_d_wf_etns[tid] = pihm->elem[tid].wf.esnow;

		// Data are transferred into pointer array from host to device.
		for (i = 0; i < MAXLYR; i++){
			pihm_d->elem_d_wf_et[i][tid] = pihm->elem[tid].wf.et[i];
			pihm_d->elem_d_wf_runoff2_lyr[i][tid] = pihm->elem[tid].wf.runoff2_lyr[i];
			pihm_d->elem_d_wf_smflxv[i][tid] = pihm->elem[tid].wf.smflxv[i];
		}

		/*
		for (i = 0; i < NUM_EDGE; i++)
		   for (int j = 0; j < MAXLYR; j++)
		      pihm_d->elem_d_wf_smflxh[i][j][tid] = pihm->elem[tid].wf.smflxh[i][j];
		*/
#endif

		/* Energy states */
		pihm_d->elem_d_es_sfctmp[tid] = pihm->elem[tid].es.sfctmp;
#if defined(_NOAH_)
		pihm_d->elem_d_es_t1[tid] = pihm->elem[tid].es.t1;
		pihm_d->elem_d_es_th2[tid] = pihm->elem[tid].es.th2;
		for (i = 0; i < MAXLYR; i++)
			pihm_d->elem_d_es_stc[i][tid] = pihm->elem[tid].es.stc[i];
#endif

		/* Energy fluxes */
		pihm_d->elem_d_ef_soldn[tid] = pihm->elem[tid].ef.soldn;
#if defined(_NOAH_)
		for (i = 0; i < MAXLYR; i++)
			pihm_d->elem_d_ef_et[i][tid] = pihm->elem[tid].ef.et[i];

		pihm_d->elem_d_ef_solnet[tid] = pihm->elem[tid].ef.solnet;
		pihm_d->elem_d_ef_etp[tid] = pihm->elem[tid].ef.etp;
		pihm_d->elem_d_ef_ssoil[tid] = pihm->elem[tid].ef.ssoil;
		pihm_d->elem_d_ef_eta[tid] = pihm->elem[tid].ef.eta;
		pihm_d->elem_d_ef_sheat[tid] = pihm->elem[tid].ef.sheat;
		pihm_d->elem_d_ef_fdown[tid] = pihm->elem[tid].ef.fdown;
		pihm_d->elem_d_ef_lwdn[tid] = pihm->elem[tid].ef.lwdn;
		pihm_d->elem_d_ef_ec[tid] = pihm->elem[tid].ef.ec;
		pihm_d->elem_d_ef_edir[tid] = pihm->elem[tid].ef.edir;
		pihm_d->elem_d_ef_ett[tid] = pihm->elem[tid].ef.ett;
		pihm_d->elem_d_ef_esnow[tid] = pihm->elem[tid].ef.esnow;
		pihm_d->elem_d_ef_soldir[tid] = pihm->elem[tid].ef.soldir;
		pihm_d->elem_d_ef_soldif[tid] = pihm->elem[tid].ef.soldif;
		pihm_d->elem_d_ef_longwave[tid] = pihm->elem[tid].ef.longwave;
		pihm_d->elem_d_ef_flx1[tid] = pihm->elem[tid].ef.flx1;
		pihm_d->elem_d_ef_flx2[tid] = pihm->elem[tid].ef.flx2;
		pihm_d->elem_d_ef_flx3[tid] = pihm->elem[tid].ef.flx3;
#endif
		/* Boundary conditions */
		for (i = 0; i < NUM_EDGE; i++){
			pihm_d->elem_d_bc_flux[i][tid] = pihm->elem[tid].bc.flux[i];
			pihm_d->elem_d_bc_head[i][tid] = pihm->elem[tid].bc.head[i];
		}

	}  // for (tid = 0; tid < nelem; tid++)
	printf("Element: Transfer AoS on CPU into SoA on GPU!\n\n");

	/* Variables related with basin elements */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
	for (tid = 0; tid < nriver; tid++){

		/* River geometrics */
		pihm_d->river_d_ind[tid] = pihm->river[tid].ind;
		pihm_d->river_d_leftele[tid] = pihm->river[tid].leftele;
		pihm_d->river_d_rightele[tid] = pihm->river[tid].rightele;
		pihm_d->river_d_fromnode[tid] = pihm->river[tid].fromnode;
		pihm_d->river_d_tonode[tid] = pihm->river[tid].tonode;
		pihm_d->river_d_down[tid] = pihm->river[tid].down;

		/* River attribute */
		pihm_d->river_d_attrib_riverbc_type[tid] = pihm->river[tid].attrib.riverbc_type;

		/* River topography parameters */
		pihm_d->river_d_topo_area[tid] = pihm->river[tid].topo.area;
		pihm_d->river_d_topo_x[tid] = pihm->river[tid].topo.x;
		pihm_d->river_d_topo_y[tid] = pihm->river[tid].topo.y;
		pihm_d->river_d_topo_zmin[tid] = pihm->river[tid].topo.zmin;
		pihm_d->river_d_topo_zmax[tid] = pihm->river[tid].topo.zmax;
		pihm_d->river_d_topo_zbed[tid] = pihm->river[tid].topo.zbed;
		pihm_d->river_d_topo_node_zmax[tid] = pihm->river[tid].topo.node_zmax;
		pihm_d->river_d_topo_dist_left[tid] = pihm->river[tid].topo.dist_left;
		pihm_d->river_d_topo_dist_right[tid] = pihm->river[tid].topo.dist_right;

		/* River water states */
		pihm_d->river_d_ws_stage[tid] = pihm->river[tid].ws.stage;
		pihm_d->river_d_ws_gw[tid] = pihm->river[tid].ws.gw;

		/* River water fluxes */
		for (i = 0; i < NUM_RIVFLX; i++)
			pihm_d->river_d_wf_rivflow[i][tid] = pihm->river[tid].wf.rivflow[i];

		/* River shape parameters */
		pihm_d->river_d_shp_depth[tid] = pihm->river[tid].shp.depth;
		pihm_d->river_d_shp_intrpl_ord[tid] = pihm->river[tid].shp.intrpl_ord;
		pihm_d->river_d_shp_coeff[tid] = pihm->river[tid].shp.coeff;
		pihm_d->river_d_shp_length[tid] = pihm->river[tid].shp.length;
		pihm_d->river_d_shp_width[tid] = pihm->river[tid].shp.width;

		/* matl_struct_d */
		pihm_d->river_d_matl_rough[tid] = pihm->river[tid].matl.rough;
		pihm_d->river_d_matl_cwr[tid] = pihm->river[tid].matl.cwr;
		pihm_d->river_d_matl_ksath[tid] = pihm->river[tid].matl.ksath;
		pihm_d->river_d_matl_ksatv[tid] = pihm->river[tid].matl.ksatv;
		pihm_d->river_d_matl_bedthick[tid] = pihm->river[tid].matl.bedthick;
		pihm_d->river_d_matl_porosity[tid] = pihm->river[tid].matl.porosity;
		pihm_d->river_d_matl_smcmin[tid] = pihm->river[tid].matl.smcmin;

		/* River boundary conditions */
		pihm_d->river_d_bc_head[tid] = pihm->river[tid].bc.head;
		pihm_d->river_d_bc_flux[tid] = pihm->river[tid].bc.flux;

	}  // for (i = 0; i < nriver; i++)
	printf("River: Transfer AoS on CPU into SoA on GPU!\n\n");

	// Transfer other parameters into GPU, like calib_struct, ctrl_struct etc_
	/* Global calibration coefficients */
	/*
	cudaMemcpy(&cal_d_ec, &pihm->cal_ec, sizeof(realtype), cudaMemcpyHostToDevice);
	cudaMemcpy(&cal_d_ett, &pihm->cal_ett, sizeof(realtype), cudaMemcpyHostToDevice);
	cudaMemcpy(&cal_d_edir, &pihm->cal_edir, sizeof(realtype), cudaMemcpyHostToDevice);
	*/
	pihm_d->ec = pihm->cal.ec;
	pihm_d->ett = pihm->cal.ett;
	pihm_d->edir = pihm->cal.edir;

	/* Model control parameters */
	/*
	cudaMemcpy(&ctrl_d_stepsize, &pihm->ctrl_stepsize, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(&ctrl_d_surf_mode, &pihm->ctrl_surf_mode, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(&ctrl_d_riv_mode, &pihm->ctrl_riv_mode, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(&ctrl_d_etstep, &pihm->ctrl_etstep, sizeof(int), cudaMemcpyHostToDevice);
	*/
	pihm_d->stepsize = pihm->ctrl.stepsize;
	pihm_d->surf_mode = pihm->ctrl.surf_mode;
	pihm_d->riv_mode = pihm->ctrl.riv_mode;

}
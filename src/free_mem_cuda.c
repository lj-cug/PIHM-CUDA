#include "pihm.h"
#include "pihm.cuh"

cudaError_t FreeMem_cuda(pihm_struct_d pihm_d)
{
	int             i;
	cudaError_t cudaStatus;

// --------------- Basin Element ---------------------
	/* Element geometrics */
	for (i = 0; i < NUM_EDGE; i++){
		cudaStatus = cudaFree(pihm_d->elem_d_node[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_nabr[i]);	
	}
	cudaStatus = cudaFree(pihm_d->elem_d_ind);

	/* Element attribute */
	for (i = 0; i < NUM_EDGE; i++)
		cudaStatus = cudaFree(pihm_d->elem_d_attrib_bc_type[i]);
	cudaStatus = cudaFree(pihm_d->elem_d_attrib_soil_type);
	cudaStatus = cudaFree(pihm_d->elem_d_attrib_lc_type);
	cudaStatus = cudaFree(pihm_d->elem_d_attrib_meteo_type);
	cudaStatus = cudaFree(pihm_d->elem_d_attrib_lai_type);

	/* Topography parameters */
	cudaStatus = cudaFree(pihm_d->elem_d_topo_area);
	cudaStatus = cudaFree(pihm_d->elem_d_topo_x);
	cudaStatus = cudaFree(pihm_d->elem_d_topo_y);
	cudaStatus = cudaFree(pihm_d->elem_d_topo_zmin);
	cudaStatus = cudaFree(pihm_d->elem_d_topo_zmax);
	for (i = 0; i < NUM_EDGE; i++){
		cudaStatus = cudaFree(pihm_d->elem_d_topo_edge[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_topo_nabrdist[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_topo_nabr_x[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_topo_nabr_y[i]);
	}
#if defined(_NOAH_)
	cudaStatus = cudaFree(pihm_d->elem_d_topo_slope);
	cudaStatus = cudaFree(pihm_d->elem_d_topo_aspect);
	cudaStatus = cudaFree(pihm_d->elem_d_topo_svf);
	/*
	for (i = 0; i < 36; i++)
		cudaStatus = cudaFree(pihm_d->elem_d_topo_h_phi[i]);
	*/
#endif	

	/* Soil parameters  soil */
	cudaStatus = cudaFree(pihm_d->elem_d_soil_depth),
	cudaStatus = cudaFree(pihm_d->elem_d_soil_ksath),
	cudaStatus = cudaFree(pihm_d->elem_d_soil_ksatv),
	cudaStatus = cudaFree(pihm_d->elem_d_soil_kinfv),
	cudaStatus = cudaFree(pihm_d->elem_d_soil_dinf),
	cudaStatus = cudaFree(pihm_d->elem_d_soil_alpha),
	cudaStatus = cudaFree(pihm_d->elem_d_soil_beta);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_porosity);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_smcmax);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_smcmin);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_smcwlt);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_smcref);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_dmac);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_kmach);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_kmacv);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_areafv);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_areafh);
#if defined(_NOAH_)	
	cudaStatus = cudaFree(pihm_d->elem_d_soil_csoil);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_quartz);
	cudaStatus = cudaFree(pihm_d->elem_d_soil_smcdry);
#endif	

	/* Land cover parameters */
	cudaStatus = cudaFree(pihm_d->elem_d_lc_shdfac);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_shdmin);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_shdmax);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_laimin);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_laimax);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_snup);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_cfactr);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_emissmax);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_emissmin);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_albedomax);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_albedomin);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_z0max);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_z0min);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_rough);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_cmcfactr);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_bare);
	cudaStatus = cudaFree(pihm_d->elem_d_lc_isurban);

	/* Ecophysiological parameters */
	cudaStatus = cudaFree(pihm_d->elem_d_epc_rsmin);
	cudaStatus = cudaFree(pihm_d->elem_d_epc_rgl);
	cudaStatus = cudaFree(pihm_d->elem_d_epc_hs);
	cudaStatus = cudaFree(pihm_d->elem_d_epc_topt);
	cudaStatus = cudaFree(pihm_d->elem_d_epc_rsmax);

	/* Physical states  ps */
	cudaStatus = cudaFree(pihm_d->elem_d_ps_rzd);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_rc);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_pc);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_proj_lai);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_rcs);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_rct);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_rcq);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_rcsoil);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_albedo);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_zlvl);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_zlvl_wind);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_sfcspd);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_rh);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_sfcprs);
#if defined(_NOAH_)
	cudaStatus = cudaFree(pihm_d->elem_d_ps_nroot);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_nsoil);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_nmacd);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_nwtbl);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_alb);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_snoalb);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_soilw);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_frzk);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_frzx);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_czil);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_emissi);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_ch);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_cm);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_rch);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_z0);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_fcr);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_salp);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_sbeta);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_lvcoef);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_snotime1);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_ribb);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_beta);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_sncovr);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_q1);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_q2);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_ffrozp);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_z0brd);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_embrd);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_q2sat);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_q2d);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_dqsdt2);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_sndens);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_snowh);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_sncond);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_rr);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_epsca);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_eta_kinematic);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_zbot);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_tbot);
	cudaStatus = cudaFree(pihm_d->elem_d_ps_gwet);

	for (i = 0; i < MAXLYR; i++){
		cudaStatus = cudaFree(pihm_d->elem_d_ps_rtdis[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_ps_sldpth[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_ps_zsoil[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_ps_satdpth[i]);
	}

#endif

	/* Water states  ws */
	cudaStatus = cudaFree(pihm_d->elem_d_ws_surf);
	cudaStatus = cudaFree(pihm_d->elem_d_ws_unsat);
	cudaStatus = cudaFree(pihm_d->elem_d_ws_gw);
	cudaStatus = cudaFree(pihm_d->elem_d_ws_sneqv);
	cudaStatus = cudaFree(pihm_d->elem_d_ws_cmcmax);
	cudaStatus = cudaFree(pihm_d->elem_d_ws_cmc);
	cudaStatus = cudaFree(pihm_d->elem_d_ws_surfh);
#if defined(_NOAH_)
	cudaStatus = cudaFree(pihm_d->elem_d_ws_soilm);
	for (i = 0; i < MAXLYR; i++){
		cudaStatus = cudaFree(pihm_d->elem_d_ws_smc[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_ws_sh2o[i]);
	}
#endif

	/* Water fluxes  ws0 */
	cudaStatus = cudaFree(pihm_d->elem_d_ws0_surf);
	cudaStatus = cudaFree(pihm_d->elem_d_ws0_unsat);
	cudaStatus = cudaFree(pihm_d->elem_d_ws0_gw);
	cudaStatus = cudaFree(pihm_d->elem_d_ws0_sneqv);
	cudaStatus = cudaFree(pihm_d->elem_d_ws0_cmcmax);
	cudaStatus = cudaFree(pihm_d->elem_d_ws0_cmc);
	cudaStatus = cudaFree(pihm_d->elem_d_ws0_surfh);

	/* Water fluxes  wf */
	for (i = 0; i < NUM_EDGE; i++){
		cudaStatus = cudaFree(pihm_d->elem_d_wf_ovlflow[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_wf_subsurf[i]);
	}
	cudaStatus = cudaFree(pihm_d->elem_d_wf_prcp);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_pcpdrp);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_infil);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_rechg);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_drip);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_edir);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_ett);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_ec);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_etp);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_eta);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_edir_surf);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_edir_unsat);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_edir_gw);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_ett_unsat);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_ett_gw);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_esnow);
#if defined(_NOAH_)
	cudaStatus = cudaFree(pihm_d->elem_d_wf_runoff2);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_runoff3);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_dew);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_snomlt);
	cudaStatus = cudaFree(pihm_d->elem_d_wf_etns);

	for (i = 0; i < MAXLYR; i++){
		cudaStatus = cudaFree(pihm_d->elem_d_wf_et[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_wf_runoff2_lyr[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_wf_smflxv[i]);
	}
	/*
	for (i = 0; i < NUM_EDGE; i++)
		for (int j = 0; j < MAXLYR; j++)
			cudaStatus = cudaFree(pihm_d->elem_d_wf_smflxh[i][j]);
	*/
#endif

	/* Energy states */
	cudaStatus = cudaFree(pihm_d->elem_d_es_sfctmp);
#if defined(_NOAH_)
	cudaStatus = cudaFree(pihm_d->elem_d_es_t1);
	cudaStatus = cudaFree(pihm_d->elem_d_es_th2);
	cudaStatus = cudaFree(pihm_d->elem_d_es_stc);
#endif     /* Energy states */

	/* Energy fluxes */
	cudaStatus = cudaFree(pihm_d->elem_d_ef_soldn);
#if defined(_NOAH_)
	cudaStatus = cudaFree(pihm_d->elem_d_ef_solnet);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_etp);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_ssoil);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_eta);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_sheat);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_fdown);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_lwdn);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_ec);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_edir);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_ett);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_esnow);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_soldir);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_soldif);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_longwave);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_flx1);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_flx2);
	cudaStatus = cudaFree(pihm_d->elem_d_ef_flx3);
	for (i = 0; i < MAXLYR; i++)
	    cudaStatus = cudaFree(pihm_d->elem_d_ef_et[i]);
#endif

	/* Boundary conditions */
	for (i = 0; i < NUM_EDGE; i++){
		cudaStatus = cudaFree(pihm_d->elem_d_bc_head[i]);
		cudaStatus = cudaFree(pihm_d->elem_d_bc_flux[i]);
	}

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Element: cudaFree failed! Check Free_CUDA ?\n");
		exit(cudaStatus);
	}
//------------- River Segments -------------------
	/* River geometrics */
	cudaStatus = cudaFree(pihm_d->river_d_ind);
	cudaStatus = cudaFree(pihm_d->river_d_leftele);
	cudaStatus = cudaFree(pihm_d->river_d_rightele);
	cudaStatus = cudaFree(pihm_d->river_d_fromnode);
	cudaStatus = cudaFree(pihm_d->river_d_tonode);
	cudaStatus = cudaFree(pihm_d->river_d_down);

	/* River attribute */
	cudaStatus = cudaFree(pihm_d->river_d_attrib_riverbc_type);

	/* River topography parameters */
	cudaStatus = cudaFree(pihm_d->river_d_topo_area);
	cudaStatus = cudaFree(pihm_d->river_d_topo_x);
	cudaStatus = cudaFree(pihm_d->river_d_topo_y);
	cudaStatus = cudaFree(pihm_d->river_d_topo_zmin);
	cudaStatus = cudaFree(pihm_d->river_d_topo_zmax);
	cudaStatus = cudaFree(pihm_d->river_d_topo_zbed);
	cudaStatus = cudaFree(pihm_d->river_d_topo_node_zmax);
	cudaStatus = cudaFree(pihm_d->river_d_topo_dist_left);
	cudaStatus = cudaFree(pihm_d->river_d_topo_dist_right);

	/* River water states */
	cudaStatus = cudaFree(pihm_d->river_d_ws_stage);
	cudaStatus = cudaFree(pihm_d->river_d_ws_gw);

	/* River water fluxes */
	for (i = 0; i < NUM_RIVFLX; i++)
		cudaStatus = cudaFree(pihm_d->river_d_wf_rivflow[i]);

	/* River shape parameters */
	cudaStatus = cudaFree(pihm_d->river_d_shp_depth);
	cudaStatus = cudaFree(pihm_d->river_d_shp_intrpl_ord);
	cudaStatus = cudaFree(pihm_d->river_d_shp_coeff);
	cudaStatus = cudaFree(pihm_d->river_d_shp_length);
	cudaStatus = cudaFree(pihm_d->river_d_shp_width);

	/* matl_struct_d */
	cudaStatus = cudaFree(pihm_d->river_d_matl_rough);
	cudaStatus = cudaFree(pihm_d->river_d_matl_cwr);
	cudaStatus = cudaFree(pihm_d->river_d_matl_ksath);
	cudaStatus = cudaFree(pihm_d->river_d_matl_ksatv);
	cudaStatus = cudaFree(pihm_d->river_d_matl_bedthick);
	cudaStatus = cudaFree(pihm_d->river_d_matl_porosity);
	cudaStatus = cudaFree(pihm_d->river_d_matl_smcmin);

	/* River boundary conditions */
	cudaStatus = cudaFree(pihm_d->river_d_bc_head);
	cudaStatus = cudaFree(pihm_d->river_d_bc_flux);

	// check the CUDA allocation status
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "River: cudaFree failed! Check Free_CUDA ?\n");
		exit(cudaStatus);
	}
	return(cudaStatus);
}
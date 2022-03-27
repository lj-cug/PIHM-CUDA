/* elem_d related arguments, totally 128 arguments */
    int *elem_d_attrib_soil_type, 
	int *elem_d_attrib_lc_type, 
	
	// 修改为数组的首地址形参代入
    int *elem_d_attrib_bc_type0,
	int *elem_d_attrib_bc_type1,
	int *elem_d_attrib_bc_type2,
	
	int *elem_d_attrib_meteo_type, 
	int *elem_d_attrib_lai_type, /* Element attribute */          
    realtype *elem_d_topo_area, 
	realtype *elem_d_topo_x, 
	realtype *elem_d_topo_y,                               
    realtype *elem_d_topo_zmin, 
	realtype *elem_d_topo_zmax,  
	
	realtype *elem_d_topo_edge0,
	realtype *elem_d_topo_edge1,
	realtype *elem_d_topo_edge2,
	
	realtype *elem_d_topo_nabrdist0,
	realtype *elem_d_topo_nabrdist1,
	realtype *elem_d_topo_nabrdist2,
	
	realtype *elem_d_topo_nabr_x0,
	realtype *elem_d_topo_nabr_x1,
    realtype *elem_d_topo_nabr_x2,
	
	realtype *elem_d_topo_nabr_y0,
	realtype *elem_d_topo_nabr_y1,
	realtype *elem_d_topo_nabr_y2,   /* Topography */ 
	
    realtype *elem_d_soil_depth,
	realtype *elem_d_soil_ksath, 
	realtype *elem_d_soil_ksatv,     
    realtype *elem_d_soil_kinfv, 
	realtype *elem_d_soil_dinf, 
	realtype *elem_d_soil_alpha,     
    realtype *elem_d_soil_beta, 
	realtype *elem_d_soil_porosity, 
	realtype *elem_d_soil_smcmax,        
    realtype *elem_d_soil_smcmin,
	realtype *elem_d_soil_smcwlt, 
	realtype *elem_d_soil_smcref,    
    realtype *elem_d_soil_dmac, 
	realtype *elem_d_soil_kmach, 
	realtype *elem_d_soil_kmacv,  
    realtype *elem_d_soil_areafv, 
	realtype *elem_d_soil_areafh,    /* Soil parameters */
    realtype *elem_d_lc_shdfac, 
	realtype *elem_d_lc_shdmin, 
	realtype *elem_d_lc_shdmax,     
    realtype *elem_d_lc_laimin, 
	realtype *elem_d_lc_laimax, 
	realtype *elem_d_lc_snup,         
    realtype *elem_d_lc_cfactr, 
	realtype *elem_d_lc_emissmax, 
	realtype *elem_d_lc_emissmin,     
    realtype *elem_d_lc_albedomax, 
	realtype *elem_d_lc_albedomin, 
	realtype *elem_d_lc_z0max,  
    realtype *elem_d_lc_z0min, 
	realtype *elem_d_lc_rough, 
	realtype *elem_d_lc_cmcfactr,       
    int *elem_d_lc_bare, 
	int *elem_d_lc_isurban,     /* Land cover parameters */ 
    realtype *elem_d_epc_rsmin, 
	realtype *elem_d_epc_rgl, 
	realtype *elem_d_epc_hs,       
    realtype *elem_d_epc_topt, 
	realtype *elem_d_epc_rsmax,   /* Ecophysiological parameters */    
    realtype *elem_d_ps_rzd, 
	realtype *elem_d_ps_rc, 
	realtype *elem_d_ps_pc,                
    realtype *elem_d_ps_proj_lai, 
	realtype *elem_d_ps_rcs, 
	realtype *elem_d_ps_rct,        
    realtype *elem_d_ps_rcq, 
	realtype *elem_d_ps_rcsoil, 
	realtype *elem_d_ps_albedo,          
    realtype *elem_d_ps_zlvl, 
	realtype *elem_d_ps_zlvl_wind, 
	realtype *elem_d_ps_sfcspd,         
    realtype *elem_d_ps_rh, 
	realtype *elem_d_ps_sfcprs,   /* Physical states */    
    realtype *elem_d_ws_surf, 
	realtype *elem_d_ws_unsat, 
	realtype *elem_d_ws_gw, 
	realtype *elem_d_ws_sneqv,                  
    realtype *elem_d_ws_cmcmax,  
	realtype *elem_d_ws_cmc,  
	realtype *elem_d_ws_surfh,    /* Water states */
    realtype *elem_d_ws0_surf, 
	realtype *elem_d_ws0_unsat, 
	realtype *elem_d_ws0_gw,  
	realtype *elem_d_ws0_sneqv, 
	realtype *elem_d_ws0_cmcmax,                                      
    realtype *elem_d_ws0_cmc, 
	realtype *elem_d_ws0_surfh,   /* Initial Water states */	
	
    realtype *elem_d_wf_ovlflow0,
	realtype *elem_d_wf_ovlflow1,
	realtype *elem_d_wf_ovlflow2,
	
	realtype *elem_d_wf_subsurf0,
	realtype *elem_d_wf_subsurf1,
	realtype *elem_d_wf_subsurf2,
	
    realtype *elem_d_wf_prcp, 
	realtype *elem_d_wf_pcpdrp, 
	realtype *elem_d_wf_infil,                        
    realtype *elem_d_wf_rechg, 
	realtype *elem_d_wf_drip, 
	realtype *elem_d_wf_edir,                       
    realtype *elem_d_wf_ett, 
	realtype *elem_d_wf_ec, 
	realtype *elem_d_wf_etp,                           
    realtype *elem_d_wf_eta, 
	realtype *elem_d_wf_edir_surf, 
	realtype *elem_d_wf_edir_unsat,                      
    realtype *elem_d_wf_edir_gw, 
	realtype *elem_d_wf_ett_unsat, 
	realtype *elem_d_wf_ett_gw,                     
    realtype *elem_d_wf_esnow,          /* Water fluxes */         
    realtype *elem_d_es_sfctmp,         /* Energy states */
    realtype *elem_d_ef_soldn,          /* Energy fluxes */ 
	
	realtype *elem_d_bc_head0,
	realtype *elem_d_bc_head1,
	realtype *elem_d_bc_head2,
	
	realtype *elem_d_bc_flux0,
	realtype *elem_d_bc_flux1,
	realtype *elem_d_bc_flux2,    /* Boundary conditions */
	
	int *elem_d_node0,
	int *elem_d_node1,
	int *elem_d_node2,
	
	int *elem_d_nabr0,
	int *elem_d_nabr1,
	int *elem_d_nabr2,
	
	int *elem_d_ind    /* Elements' geometric numbers */    
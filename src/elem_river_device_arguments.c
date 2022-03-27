/* elem_d related arguments */
		elem_d_attrib_soil_type, 
		elem_d_attrib_lc_type,
		
        // 指针数组改写为普通数组的地址形参，下面类似
		elem_d_attrib_bc_type0,
		elem_d_attrib_bc_type1,
		elem_d_attrib_bc_type2,
		
		elem_d_attrib_meteo_type,
		elem_d_attrib_lai_type,         
		elem_d_topo_area, 
		elem_d_topo_x, 
		elem_d_topo_y,
		elem_d_topo_zmin, 
		elem_d_topo_zmax,
		
		elem_d_topo_edge0,
		elem_d_topo_edge1,
		elem_d_topo_edge2,
		
	    elem_d_topo_nabrdist0,
		elem_d_topo_nabrdist1,
		elem_d_topo_nabrdist2,
		
		elem_d_topo_nabr_x0,
		elem_d_topo_nabr_x1,
		elem_d_topo_nabr_x2,
		
		elem_d_topo_nabr_y0,
		elem_d_topo_nabr_y1,
		elem_d_topo_nabr_y2,

		elem_d_soil_depth, 
		elem_d_soil_ksath, 
		elem_d_soil_ksatv,
		elem_d_soil_kinfv, 
		elem_d_soil_dinf, 
		elem_d_soil_alpha,
		elem_d_soil_beta, 
		elem_d_soil_porosity, 
		elem_d_soil_smcmax,
		elem_d_soil_smcmin, 
		elem_d_soil_smcwlt, 
		elem_d_soil_smcref,
		elem_d_soil_dmac, 
		elem_d_soil_kmach, 
		elem_d_soil_kmacv,
		elem_d_soil_areafv, 
		elem_d_soil_areafh,             
		elem_d_lc_shdfac, 
		elem_d_lc_shdmin, 
		elem_d_lc_shdmax,
		elem_d_lc_laimin, 
		elem_d_lc_laimax, 
		elem_d_lc_snup,
		elem_d_lc_cfactr, 
		elem_d_lc_emissmax, 
		elem_d_lc_emissmin,
		elem_d_lc_albedomax, 
		elem_d_lc_albedomin, 
		elem_d_lc_z0max,
		elem_d_lc_z0min, 
		elem_d_lc_rough, 
		elem_d_lc_cmcfactr,
		elem_d_lc_bare, 
		elem_d_lc_isurban,    
		elem_d_epc_rsmin, 
		elem_d_epc_rgl, 
		elem_d_epc_hs,
		elem_d_epc_topt,
		elem_d_epc_rsmax,     
		elem_d_ps_rzd, 
		elem_d_ps_rc, 
		elem_d_ps_pc,
		elem_d_ps_proj_lai, 
		elem_d_ps_rcs, 
		elem_d_ps_rct,
		elem_d_ps_rcq, 
		elem_d_ps_rcsoil, 
		elem_d_ps_albedo,
		elem_d_ps_zlvl, 
		elem_d_ps_zlvl_wind, 
		elem_d_ps_sfcspd,
		elem_d_ps_rh, 
		elem_d_ps_sfcprs,       
		elem_d_ws_surf, 
		elem_d_ws_unsat, 
		elem_d_ws_gw, 
		elem_d_ws_sneqv,
		elem_d_ws_cmcmax, 
		elem_d_ws_cmc, 
		elem_d_ws_surfh,      
		elem_d_ws0_surf, 
		elem_d_ws0_unsat, 
		elem_d_ws0_gw,
		elem_d_ws0_sneqv, 
		elem_d_ws0_cmcmax,
		elem_d_ws0_cmc, 
		elem_d_ws0_surfh, 
		
		elem_d_wf_ovlflow0,
		elem_d_wf_ovlflow1,
		elem_d_wf_ovlflow2,
		
		elem_d_wf_subsurf0,
		elem_d_wf_subsurf1,
		elem_d_wf_subsurf2,
				
		elem_d_wf_prcp, 
		elem_d_wf_pcpdrp, 
		elem_d_wf_infil,
		elem_d_wf_rechg, 
		elem_d_wf_drip, 
		elem_d_wf_edir,
		elem_d_wf_ett, 
		elem_d_wf_ec, 
		elem_d_wf_etp,
		elem_d_wf_eta, 
		elem_d_wf_edir_surf, 
		elem_d_wf_edir_unsat,
		elem_d_wf_edir_gw, 
		elem_d_wf_ett_unsat, 
		elem_d_wf_ett_gw,
		elem_d_wf_esnow,         
		elem_d_es_sfctmp,         
		elem_d_ef_soldn,  

        elem_d_bc_head0,
		elem_d_bc_head1,
		elem_d_bc_head2,
		
		elem_d_bc_flux0,
		elem_d_bc_flux1,
		elem_d_bc_flux2,
		
		elem_d_node0,
		elem_d_node1,
		elem_d_node2,
		
		elem_d_nabr0,
		elem_d_nabr1,
		elem_d_nabr2,
		
		elem_d_ind, 
 		
/* river_d related arguments */
		river_d_attrib_riverbc_type,   
		river_d_topo_area, 
		river_d_topo_x, 
		river_d_topo_y,
		river_d_topo_zmin, 
		river_d_topo_zmax, 
		river_d_topo_zbed,
		river_d_topo_node_zmax, 
		river_d_topo_dist_left,
		river_d_topo_dist_right,        
		river_d_ws_stage, 
		river_d_ws_gw,  
		
		river_d_wf_rivflow0,
		river_d_wf_rivflow1,
		river_d_wf_rivflow2,
		river_d_wf_rivflow3,
		river_d_wf_rivflow4,
		river_d_wf_rivflow5,
		river_d_wf_rivflow6,
		river_d_wf_rivflow7,
		river_d_wf_rivflow8,
		river_d_wf_rivflow9,
		river_d_wf_rivflow10,		
		
		river_d_shp_depth, 
		river_d_shp_intrpl_ord, 
		river_d_shp_coeff,
		river_d_shp_length, 
		river_d_shp_width,      
		river_d_matl_rough, 
		river_d_matl_cwr, 
		river_d_matl_ksath,
		river_d_matl_ksatv, 
		river_d_matl_bedthick, 
		river_d_matl_porosity,
		river_d_matl_smcmin,      
		river_d_bc_head, 
		river_d_bc_flux,          
		river_d_ind, 
		river_d_leftele, 
		river_d_rightele,
		river_d_fromnode, 
		river_d_tonode, 
		river_d_down                
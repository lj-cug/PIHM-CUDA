/* river_d related arguments, totally 43 arguments */           
    int *river_d_attrib_riverbc_type, /* River attribute */
    realtype *river_d_topo_area, 
	realtype *river_d_topo_x, 
	realtype *river_d_topo_y,                 
    realtype *river_d_topo_zmin,
	realtype *river_d_topo_zmax, 
	realtype *river_d_topo_zbed,                            
    realtype *river_d_topo_node_zmax, 
	realtype *river_d_topo_dist_left, 
	realtype *river_d_topo_dist_right, /* River topography parameters */   
    realtype *river_d_ws_stage, 
	realtype *river_d_ws_gw,           /* River water states */
	
	realtype *river_d_wf_rivflow0,
	realtype *river_d_wf_rivflow1,
	realtype *river_d_wf_rivflow2,
	realtype *river_d_wf_rivflow3,	
	realtype *river_d_wf_rivflow4,
	realtype *river_d_wf_rivflow5,
	realtype *river_d_wf_rivflow6,
	realtype *river_d_wf_rivflow7,
	realtype *river_d_wf_rivflow8,
	realtype *river_d_wf_rivflow9,
	realtype *river_d_wf_rivflow10,   /* River water fluxes */   

    realtype *river_d_shp_depth, 
	int      *river_d_shp_intrpl_ord,  
	realtype *river_d_shp_coeff,     
    realtype *river_d_shp_length, 
	realtype *river_d_shp_width,       /* River shape parameters */     
    realtype *river_d_matl_rough, 
	realtype *river_d_matl_cwr, 
	realtype *river_d_matl_ksath,         
    realtype *river_d_matl_ksatv, 
	realtype *river_d_matl_bedthick, 
	realtype *river_d_matl_porosity,     
    realtype *river_d_matl_smcmin, /* River mathmatic modeling */	
    realtype *river_d_bc_head, 
	realtype *river_d_bc_flux,  /* River boundary conditions */  
    int *river_d_ind, 
	int *river_d_leftele, 
	int *river_d_rightele,      
    int *river_d_fromnode, 
	int *river_d_tonode, 
	int *river_d_down  /* River geometric numbers */
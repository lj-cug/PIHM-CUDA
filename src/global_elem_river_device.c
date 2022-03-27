	/* elem_d related arguments */
    extern int *elem_d_attrib_soil_type; 
	extern int *elem_d_attrib_lc_type;             
    extern int *elem_d_attrib_bc_type[NUM_EDGE];
	extern int *elem_d_attrib_meteo_type; 
	extern int *elem_d_attrib_lai_type;           /* Attributes of Elements*/                                              
    extern double *elem_d_topo_area; 
	extern double *elem_d_topo_x; 
	extern double *elem_d_topo_y;                               
    extern double *elem_d_topo_zmin; 
	extern double *elem_d_topo_zmax;                                
	extern double *elem_d_topo_edge[NUM_EDGE]; 
	extern double *elem_d_topo_nabrdist[NUM_EDGE];        
    extern double *elem_d_topo_nabr_x[NUM_EDGE]; 
	extern double *elem_d_topo_nabr_y[NUM_EDGE];  /* Topography parameters */                 
    extern double *elem_d_soil_depth;
	extern double *elem_d_soil_ksath; 
	extern double *elem_d_soil_ksatv;     
    extern double *elem_d_soil_kinfv; 
	extern double *elem_d_soil_dinf; 
	extern double *elem_d_soil_alpha;     
    extern double *elem_d_soil_beta; 
	extern double *elem_d_soil_porosity; 
	extern double *elem_d_soil_smcmax;        
    extern double *elem_d_soil_smcmin;
	extern double *elem_d_soil_smcwlt; 
	extern double *elem_d_soil_smcref;    
    extern double *elem_d_soil_dmac; 
	extern double *elem_d_soil_kmach; 
	extern double *elem_d_soil_kmacv;  
    extern double *elem_d_soil_areafv; 
	extern double *elem_d_soil_areafh;       /* Soil parameters */
    extern double *elem_d_lc_shdfac; 
	extern double *elem_d_lc_shdmin; 
	extern double *elem_d_lc_shdmax;     
    extern double *elem_d_lc_laimin; 
	extern double *elem_d_lc_laimax; 
	extern double *elem_d_lc_snup;         
    extern double *elem_d_lc_cfactr; 
	extern double *elem_d_lc_emissmax; 
	extern double *elem_d_lc_emissmin;     
    extern double *elem_d_lc_albedomax; 
	extern double *elem_d_lc_albedomin; 
	extern double *elem_d_lc_z0max;  
    extern double *elem_d_lc_z0min; 
	extern double *elem_d_lc_rough; 
	extern double *elem_d_lc_cmcfactr;       
    extern int *elem_d_lc_bare; 
	extern int *elem_d_lc_isurban;          /* Land cover parameters */ 
    extern double *elem_d_epc_rsmin; 
	extern double *elem_d_epc_rgl; 
	extern double *elem_d_epc_hs;       
    extern double *elem_d_epc_topt; 
	extern double *elem_d_epc_rsmax;        /* Ecophysiological parameters */    
    extern double *elem_d_ps_rzd; 
	extern double *elem_d_ps_rc; 
	extern double *elem_d_ps_pc;                
    extern double *elem_d_ps_proj_lai; 
	extern double *elem_d_ps_rcs; 
	extern double *elem_d_ps_rct;        
    extern double *elem_d_ps_rcq; 
	extern double *elem_d_ps_rcsoil; 
	extern double *elem_d_ps_albedo;          
    extern double *elem_d_ps_zlvl; 
	extern double *elem_d_ps_zlvl_wind; 
	extern double *elem_d_ps_sfcspd;         
    extern double *elem_d_ps_rh; 
	extern double *elem_d_ps_sfcprs;        /* Physical states */    
    extern double *elem_d_ws_surf; 
	extern double *elem_d_ws_unsat; 
	extern double *elem_d_ws_gw; 
	extern double *elem_d_ws_sneqv;                  
    extern double *elem_d_ws_cmcmax;  
	extern double *elem_d_ws_cmc;  
	extern double *elem_d_ws_surfh;         /* Water states */                                
    extern double *elem_d_ws0_surf; 
	extern double *elem_d_ws0_unsat; 
	extern double *elem_d_ws0_gw;  
	extern double *elem_d_ws0_sneqv; 
	extern double *elem_d_ws0_cmcmax;                                      
    extern double *elem_d_ws0_cmc; 
	extern double *elem_d_ws0_surfh;        /* Initial Water states */	     
    extern double *elem_d_wf_ovlflow[NUM_EDGE]; 
	extern double *elem_d_wf_subsurf[NUM_EDGE];   
    extern double *elem_d_wf_prcp; 
	extern double *elem_d_wf_pcpdrp; 
	extern double *elem_d_wf_infil;                        
    extern double *elem_d_wf_rechg; 
	extern double *elem_d_wf_drip; 
	extern double *elem_d_wf_edir;                       
    extern double *elem_d_wf_ett; 
	extern double *elem_d_wf_ec; 
	extern double *elem_d_wf_etp;                           
    extern double *elem_d_wf_eta; 
	extern double *elem_d_wf_edir_surf; 
	extern double *elem_d_wf_edir_unsat;                      
    extern double *elem_d_wf_edir_gw; 
	extern double *elem_d_wf_ett_unsat; 
	extern double *elem_d_wf_ett_gw;                     
    extern double *elem_d_wf_esnow;          /* Water fluxes */         
    extern double *elem_d_es_sfctmp;         /* Energy states */
    extern double *elem_d_ef_soldn;          /* Energy fluxes */ 
    extern double *elem_d_bc_head[NUM_EDGE]; 
	extern double *elem_d_bc_flux[NUM_EDGE]; /* Boundary conditions */
    extern int *elem_d_node[NUM_EDGE]; 
	extern int *elem_d_nabr[NUM_EDGE]; 
	extern int *elem_d_ind;                 /* Elements' geometric numbers */  
 
    /* river_d related arguments */           
    extern int *river_d_attrib_riverbc_type;    /* River attribute */
    extern double *river_d_topo_area; 
	extern double *river_d_topo_x; 
	extern double *river_d_topo_y;                 
    extern double *river_d_topo_zmin;
	extern double *river_d_topo_zmax; 
	extern double *river_d_topo_zbed;                            
    extern double *river_d_topo_node_zmax; 
	extern double *river_d_topo_dist_left; 
	extern double *river_d_topo_dist_right;          /* River topography parameters */   
    extern double *river_d_ws_stage; 
	extern double *river_d_ws_gw;                    /* River water states */
    extern double *river_d_wf_rivflow[NUM_RIVFLX];   /* River water fluxes */   
    extern double *river_d_shp_depth; 
	extern int    *river_d_shp_intrpl_ord;  
	extern double *river_d_shp_coeff;     
    extern double *river_d_shp_length; 
	extern double *river_d_shp_width;       /* River shape parameters */     
    extern double *river_d_matl_rough; 
	extern double *river_d_matl_cwr; 
	extern double *river_d_matl_ksath;         
    extern double *river_d_matl_ksatv; 
	extern double *river_d_matl_bedthick; 
	extern double *river_d_matl_porosity;     
    extern double *river_d_matl_smcmin;     /* River mathmatic modeling parameters */	
    extern double *river_d_bc_head; 
	extern double *river_d_bc_flux;         /* River boundary conditions */  
    extern int *river_d_ind; 
	extern int *river_d_leftele; 
	extern int *river_d_rightele;      
    extern int *river_d_fromnode; 
	extern int *river_d_tonode; 
	extern int *river_d_down;              /* River geometric numbers */
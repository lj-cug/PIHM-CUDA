/* elem_d related arguments */
    int *elem_d_attrib_soil_type;  /* Element attribute */ 
	int *elem_d_attrib_lc_type;             
    int *elem_d_attrib_bc_type[NUM_EDGE];
	int *elem_d_attrib_meteo_type; 
	int *elem_d_attrib_lai_type;   /* Element attribute */          
    realtype *elem_d_topo_area;    /* Topography */ 
	realtype *elem_d_topo_x; 
	realtype *elem_d_topo_y;                               
    realtype *elem_d_topo_zmin; 
	realtype *elem_d_topo_zmax; 
#if defined(_FBR_)
    realtype *elem_d_topo_zbed; 
#endif	
	realtype *elem_d_topo_edge[NUM_EDGE]; 
	realtype *elem_d_topo_nabrdist[NUM_EDGE];        
    realtype *elem_d_topo_nabr_x[NUM_EDGE]; 
	realtype *elem_d_topo_nabr_y[NUM_EDGE]; 
#if defined(_NOAH_)
    realtype *elem_d_topo_slope;    
    realtype *elem_d_topo_aspect;                             
    realtype *elem_d_topo_svf;      
    realtype *elem_d_topo_h_phi[36];    
#endif	
#if defined(_RT_)
    realtype *elem_d_topo_areasub[NUM_EDGE];   /* Topography */  
#endif
    realtype *elem_d_soil_depth;   /* Soil parameters */
	realtype *elem_d_soil_ksath; 
	realtype *elem_d_soil_ksatv;     
    realtype *elem_d_soil_kinfv; 
	realtype *elem_d_soil_dinf; 
	realtype *elem_d_soil_alpha;     
    realtype *elem_d_soil_beta; 
	realtype *elem_d_soil_porosity; 
	realtype *elem_d_soil_smcmax;        
    realtype *elem_d_soil_smcmin;
	realtype *elem_d_soil_smcwlt; 
	realtype *elem_d_soil_smcref;    
    realtype *elem_d_soil_dmac; 
	realtype *elem_d_soil_kmach; 
	realtype *elem_d_soil_kmacv;  
    realtype *elem_d_soil_areafv; 
	realtype *elem_d_soil_areafh;  
#if defined(_NOAH_)	
    realtype *elem_d_soil_csoil; 	
    realtype *elem_d_soil_quartz;	
    realtype *elem_d_soil_smcdry;    
#endif	
#if defined(_CYCLES_)
    realtype *elem_d_soil_clay[MAXLYR];
    realtype *elem_d_soil_sand[MAXLYR];
    realtype *elem_d_soil_iom[MAXLYR];
    realtype *elem_d_soil_bd[MAXLYR];  /* Soil parameters */
#endif
#if defined(_FBR_)
    realtype *elem_d_geol_depth;      /* Fractured bedrock layer parameters */
    realtype *elem_d_geol_ksath;                                     
    realtype *elem_d_geol_ksatv;                                    
    realtype *elem_d_geol_alpha;     
    realtype *elem_d_geol_beta;      
    realtype *elem_d_geol_porosity;  
    realtype *elem_d_geol_smcmax;    
    realtype *elem_d_geol_smcmin;    /* Fractured bedrock layer parameters */
#endif
	realtype *elem_d_lc_shdfac;    /* Land cover parameters */
    realtype *elem_d_lc_shdmin; 
    realtype *elem_d_lc_shdmax;     
    realtype *elem_d_lc_laimin; 
	realtype *elem_d_lc_laimax; 
	realtype *elem_d_lc_snup;         
    realtype *elem_d_lc_cfactr; 
	realtype *elem_d_lc_emissmax; 
	realtype *elem_d_lc_emissmin;     
    realtype *elem_d_lc_albedomax; 
	realtype *elem_d_lc_albedomin; 
	realtype *elem_d_lc_z0max;  
    realtype *elem_d_lc_z0min; 
	realtype *elem_d_lc_rough; 
	realtype *elem_d_lc_cmcfactr;       
    int *elem_d_lc_bare; 
	int *elem_d_lc_isurban;   /* Land cover parameters */ 
#if !defined(_CYCLES_)	
    realtype *elem_d_epc_rsmin; /* Ecophysiological parameters */  
	realtype *elem_d_epc_rgl; 
	realtype *elem_d_epc_hs;       
    realtype *elem_d_epc_topt; 
	realtype *elem_d_epc_rsmax;  
# if defined(_BGC_)
    int     *elem_d_epc_woody;     
    int     *elem_d_epc_evergreen;                                       
    int     *elem_d_epc_c3_flag;              
    int     *elem_d_epc_phenology_flag;       
    int     *elem_d_epc_onday;              
    int     *elem_d_epc_offday;              
    int     *elem_d_epc_transfer_days;                                         
    int     *elem_d_epc_litfall_days;         
    realtype  *elem_d_epc_leaf_turnover; 
    realtype  *elem_d_epc_froot_turnover;       
    realtype  *elem_d_epc_livewood_turnover;
    realtype *elem_d_epc_daily_mortality_turnover;
    realtype  *elem_d_epc_daily_fire_turnover;   
    realtype  *elem_d_epc_alloc_frootc_leafc;                              
    realtype  *elem_d_epc_alloc_newstemc_newleafc;
    realtype  *elem_d_epc_alloc_newlivewoodc_newwoodc;
    realtype  *elem_d_epc_alloc_crootc_stemc;   
    realtype  *elem_d_epc_alloc_prop_curgrowth;
    realtype  *elem_d_epc_avg_proj_sla;
    realtype  *elem_d_epc_sla_ratio; 
    realtype  *elem_d_epc_lai_ratio;
    realtype  *elem_d_epc_ext_coef; 
    realtype  *elem_d_epc_flnr;
    realtype  *elem_d_epc_psi_open;     
    realtype  *elem_d_epc_psi_close;  
    realtype  *elem_d_epc_vpd_open; 
    realtype  *elem_d_epc_vpd_close; 
    realtype  *elem_d_epc_froot_cn;           
    realtype  *elem_d_epc_leaf_cn;            
    realtype  *elem_d_epc_livewood_cn;        
    realtype  *elem_d_epc_deadwood_cn;        
    realtype  *elem_d_epc_leaflitr_cn;                                          
    realtype  *elem_d_epc_leaflitr_flab;      
    realtype  *elem_d_epc_leaflitr_fucel;       
    realtype  *elem_d_epc_leaflitr_fscel;       
    realtype  *elem_d_epc_leaflitr_flig;        
    realtype  *elem_d_epc_frootlitr_flab;       
    realtype  *elem_d_epc_frootlitr_fucel;      
    realtype  *elem_d_epc_frootlitr_fscel;      
    realtype  *elem_d_epc_frootlitr_flig;       
    realtype  *elem_d_epc_deadwood_fucel;       
    realtype  *elem_d_epc_deadwood_fscel;       
    realtype  *elem_d_epc_deadwood_flig;        	
# endif	
#endif	                       /* Ecophysiological parameters */
    realtype *elem_d_ps_rzd;   /* Physical states */
	realtype *elem_d_ps_rc; 
	realtype *elem_d_ps_pc;                
    realtype *elem_d_ps_proj_lai; 
	realtype *elem_d_ps_rcs; 
	realtype *elem_d_ps_rct;        
    realtype *elem_d_ps_rcq; 
	realtype *elem_d_ps_rcsoil; 
	realtype *elem_d_ps_albedo;          
    realtype *elem_d_ps_zlvl; 
	realtype *elem_d_ps_zlvl_wind; 
	realtype *elem_d_ps_sfcspd;         
    realtype *elem_d_ps_rh; 
	realtype *elem_d_ps_sfcprs;    
#if defined(_NOAH_)
    realtype *elem_d_ps_alb;
    realtype *elem_d_ps_snoalb;
    int      *elem_d_ps_nroot;
    realtype *elem_d_ps_rtdis[MAXLYR]; 
    int      *elem_d_ps_nsoil;
    realtype *elem_d_ps_sldpth[MAXLYR];
    realtype *elem_d_ps_zsoil[MAXLYR];
    realtype *elem_d_ps_soilw;      
    realtype *elem_d_ps_frzk;
    realtype *elem_d_ps_frzx;
    realtype *elem_d_ps_czil;
    realtype *elem_d_ps_emissi;      
    realtype *elem_d_ps_ch;           
    realtype *elem_d_ps_cm;
    realtype *elem_d_ps_rch; 
    realtype *elem_d_ps_z0;
    realtype *elem_d_ps_fcr;
    int      *elem_d_ps_nmacd;
    realtype *elem_d_ps_salp; 
    realtype *elem_d_ps_fxexp; 
    realtype *elem_d_ps_sbeta;
    realtype *elem_d_ps_lvcoef; 
    realtype *elem_d_ps_snotime1;
    realtype *elem_d_ps_ribb;
    realtype *elem_d_ps_beta;
    realtype *elem_d_ps_sncovr; 
    realtype *elem_d_ps_q1;    
    realtype *elem_d_ps_q2;
    realtype *elem_d_ps_ffrozp;
    realtype *elem_d_ps_z0brd; 
    realtype *elem_d_ps_embrd;
    realtype *elem_d_ps_q2sat; 
    realtype *elem_d_ps_q2d; 
    realtype *elem_d_ps_dqsdt2;
    int      *elem_d_ps_nwtbl;
    realtype *elem_d_ps_sndens;
    realtype *elem_d_ps_snowh;
    realtype *elem_d_ps_sncond;
    realtype *elem_d_ps_rr;
    realtype *elem_d_ps_epsca;
    realtype *elem_d_ps_eta_kinematic;
    realtype *elem_d_ps_zbot;
    realtype *elem_d_ps_tbot;
    realtype *elem_d_ps_gwet; 
    realtype *elem_d_ps_satdpth[MAXLYR];
#endif
#if defined(_CYCLES_)
    realtype *elem_d_ps_res_intcp;
    realtype *elem_d_ps_tau_res_stan;
    realtype *elem_d_ps_tau_res_flat;
    realtype *elem_d_ps_till_factr[MAXLYR];
    realtype *elem_d_ps_comp_factr[MAXLYR];
#endif
#if defined(_BGC_)
    realtype *elem_d_ps_co2;            
    realtype *elem_d_ps_ppfd_per_plaisun;
    realtype *elem_d_ps_ppfd_per_plaishade;           
    realtype *elem_d_ps_all_lai;                 
    realtype *elem_d_ps_plaisun;              
    realtype *elem_d_ps_plaishade;     /* Physical states */
#endif
    realtype *elem_d_ws_surf;   /* Water states */
	realtype *elem_d_ws_unsat; 
	realtype *elem_d_ws_gw; 
	realtype *elem_d_ws_sneqv;                  
    realtype *elem_d_ws_cmcmax;  
	realtype *elem_d_ws_cmc;  
	realtype *elem_d_ws_surfh;    
#if defined(_FBR_)
    realtype *elem_d_ws_fbr_unsat; 
    realtype *elem_d_ws_fbr_gw;
#endif
#if defined(_NOAH_)
    realtype *elem_d_ws_smc[MAXLYR];
    realtype *elem_d_ws_sh2o[MAXLYR]; 
    realtype *elem_d_ws_soilm;                                  
#endif
#if defined(_CYCLES_)
    realtype *elem_d_ws_stanResidueWater;
    realtype *elem_d_ws_flatResidueWater;
#endif	                        /* Water states */
	realtype *elem_d_ws0_surf; 
	realtype *elem_d_ws0_unsat; 
	realtype *elem_d_ws0_gw;  
	realtype *elem_d_ws0_sneqv; 
	realtype *elem_d_ws0_cmcmax;
    realtype *elem_d_ws0_cmc; 
	realtype *elem_d_ws0_surfh;   /* Initial Water states */	     
    realtype *elem_d_wf_ovlflow[NUM_EDGE];    /* Water fluxes */
	realtype *elem_d_wf_subsurf[NUM_EDGE];   
    realtype *elem_d_wf_prcp; 
	realtype *elem_d_wf_pcpdrp; 
	realtype *elem_d_wf_infil;                        
    realtype *elem_d_wf_rechg; 
	realtype *elem_d_wf_drip; 
	realtype *elem_d_wf_edir;                       
    realtype *elem_d_wf_ett; 
	realtype *elem_d_wf_ec; 
	realtype *elem_d_wf_etp;                           
    realtype *elem_d_wf_eta; 
	realtype *elem_d_wf_edir_surf; 
	realtype *elem_d_wf_edir_unsat;                      
    realtype *elem_d_wf_edir_gw; 
	realtype *elem_d_wf_ett_unsat; 
	realtype *elem_d_wf_ett_gw;                     
    realtype *elem_d_wf_esnow;        
#if defined(_FBR_)
    realtype          *elem_d_wf_fbr_infil;             
    realtype          *elem_d_wf_fbr_rechg;             
    realtype          *elem_d_wf_fbrflow[NUM_EDGE];     
#endif
#if defined(_NOAH_)
    realtype          *elem_d_wf_et[MAXLYR];       
    realtype          *elem_d_wf_runoff2;           
    realtype          *elem_d_wf_runoff2_lyr[MAXLYR];
    realtype          *elem_d_wf_runoff3;         
    realtype          *elem_d_wf_smflxv[MAXLYR];
//  realtype          *elem_d_wf_smflxh[NUM_EDGE][MAXLYR];
    realtype          *elem_d_wf_dew;      
    realtype          *elem_d_wf_snomlt;
    realtype          *elem_d_wf_etns;    
#endif
#if defined(_CYCLES_)
    realtype          *elem_d_wf_eres;              
    realtype          *elem_d_wf_irrigationVol;  /* Water fluxes */
#endif
    realtype *elem_d_es_sfctmp;                  /* Energy states */
#if defined(_NOAH_)
    realtype          *elem_d_es_t1;          
    realtype          *elem_d_es_th2; 
    realtype          *elem_d_es_stc[MAXLYR];    /* Energy states */
#endif	
    realtype          *elem_d_ef_soldn;          /* Energy fluxes */ 
#if defined(_NOAH_)
    realtype          *elem_d_ef_solnet;              
    realtype          *elem_d_ef_etp;                 
    realtype          *elem_d_ef_ssoil;               
    realtype          *elem_d_ef_eta;                 
    realtype          *elem_d_ef_sheat;               
    realtype          *elem_d_ef_fdown;               
    realtype          *elem_d_ef_lwdn;                
    realtype          *elem_d_ef_ec;                  
    realtype          *elem_d_ef_edir;                
    realtype          *elem_d_ef_et[MAXLYR];          
    realtype          *elem_d_ef_ett; 
    realtype          *elem_d_ef_esnow;               
    realtype          *elem_d_ef_soldir;              
    realtype          *elem_d_ef_soldif;              
    realtype          *elem_d_ef_longwave;            
    realtype          *elem_d_ef_flx1;                
    realtype          *elem_d_ef_flx2;                
    realtype          *elem_d_ef_flx3;                
#endif
#if defined(_BGC_)
    realtype          *elem_d_ef_swabs_per_plaisun;   
    realtype          *elem_d_ef_swabs_per_plaishade;  
#endif	
#if defined(_CYCLES_)
typedef struct epvar_struct
{
    /* User Defined Auto Irrigation */
    int             ai_used;
    int             ai_start;
    int             ai_stop;
    realtype        ai_h2o_depl;
    int             ai_last_lyr;

    /* User Defined Auto Fertilization */
    int             auto_fert;

    realtype        plant_density;
    int             clip_start;
    int             clip_end;

    /* State Variables */
    realtype          tt_daily;
    realtype          tt_cum;
    realtype          rad_intcp;
    realtype          rad_intcp_brown;
    realtype          root_dpth;
    realtype          h2o_stress;
    realtype          n_stress;
    realtype          shoot_growth_unstr_cum;
    realtype          n_stress_cum;
    realtype          rad_intcp_nc;
    int             harv_date_final;
    int             nharv;
    int             stage;
    realtype          harv_biomass;     
    realtype          harv_root;        
    realtype          harv_forage_yield;
    realtype          harv_res_biomass; 
    realtype          harv_transp;      
    realtype          harv_transp_pot;  
    realtype          harv_soil_evap;   
} epvar_struct;

typedef struct crop_wflux_struct
{
    realtype          transp;
    realtype          transp_pot;
} crop_wflux_struct;

typedef struct crop_cstate_struct
{
    realtype          shoot;            
    realtype          root;             
    realtype          rhizho;           
    realtype          shoot_post_flower;
} crop_cstate_struct;

typedef struct crop_cflux_struct
{
    realtype          shoot_growth;      
    realtype          root_growth;       
    realtype          rhizo_depo;        
    realtype          shoot_growth_unstr;
    realtype          root_growth_unstr; 
} crop_cflux_struct;

typedef struct crop_nstate_struct
{
    realtype          shoot_n;
    realtype          root_n; 
    realtype          rhizo_n;
} crop_nstate_struct;

typedef struct crop_nflux_struct
{
    realtype          rhizo_n_depo;
    realtype          n_auto;      
    realtype          n_fix;       
} crop_nflux_struct;

typedef struct crop_struct
{
    epconst_struct *epc;
    epvar_struct    epv;
    crop_wflux_struct cwf;
    crop_cstate_struct ccs;
    crop_cflux_struct ccf;
    crop_nstate_struct cns;
    crop_nflux_struct cnf;
} crop_struct;

typedef struct mgmt_struct
{
    int             rot_size;
    int             auto_n;
    int             rot_year;
    int             op_ptr[4];
} mgmt_struct;

typedef struct cstate_struct
{
    realtype          SOC_Mass[MAXLYR]; 
    realtype          MBC_Mass[MAXLYR];
    realtype          stanResidueMass;         
    realtype          flatResidueMass;         
    realtype          manureSurfaceC;          
    realtype          residueAbgd[MAXLYR];     
    realtype          residueRt[MAXLYR];       
    realtype          residueRz[MAXLYR];       
    realtype          manureC[MAXLYR];         
} cstate_struct;

typedef struct cflux_struct
{
    realtype          C_Humified;            
    realtype          C_ResidueRespired;     
    realtype          C_SoilRespired;        
    realtype          carbonRespired[MAXLYR];
} cflux_struct;

typedef struct nstate_struct
{
    realtype          no3[MAXLYR];            
    realtype          nh4[MAXLYR];            
    realtype          SON_Mass[MAXLYR];       
    realtype          MBN_Mass[MAXLYR];
    realtype          stanResidueN;           
    realtype          flatResidueN;           
    realtype          manureSurfaceN;         
    realtype          residueAbgdN[MAXLYR];   
    realtype          residueRtN[MAXLYR];     
    realtype          residueRzN[MAXLYR];     
    realtype          manureN[MAXLYR];        
} nstate_struct;

typedef struct nflux_struct
{
    realtype          no3leached;
    realtype          nh4leached;
    realtype          N_Immobilization;   
    realtype          N_Mineralization;   
    realtype          N_NetMineralization;
    realtype          nh4nitrif;          
    realtype          N2O_Nitrification;  
    realtype          no3denitrif;        
    realtype          n2odenitrif;        
    realtype          nh4volat[MAXLYR];   
    realtype          uptake_no3[MAXLYR]; 
    realtype          uptake_nh4[MAXLYR]; 
    realtype          surplusn;           
    realtype          fert_no3[MAXLYR];   
    realtype          fert_nh4[MAXLYR];
    realtype          immob_no3[MAXLYR];  
    realtype          immob_nh4[MAXLYR];  
    realtype          nitrif_nh4_to_no3[MAXLYR];
    realtype          nitrif_nh4_to_n2o[MAXLYR];
    realtype          denitn[MAXLYR];           
    realtype          till_no3[MAXLYR];         
    realtype          till_nh4[MAXLYR];
    realtype          urine;
} nflux_struct;

typedef struct nprof_struct
{
    realtype          no3;
    realtype          nh4;
} nprof_struct;

/* Solute transport structure */
typedef struct solute_struct
{
    realtype          conc;               
    realtype          flux[NUM_EDGE];     
    realtype          snksrc;             
                                          
} solute_struct;
#endif	

/* Boundary conditions */
    realtype *elem_d_bc_head[NUM_EDGE]; 
	realtype *elem_d_bc_flux[NUM_EDGE];    
	
/* Land surface and hydrologic initial conditions */	
#if defined (_DEBUG_)	
typedef struct ic_struct
{
    realtype          cmc;
    realtype          sneqv;
    realtype          surf;
    realtype          unsat;
    realtype          gw;
#if defined(_FBR_)
    realtype          fbr_unsat;
    realtype          fbr_gw;
#endif
#if defined(_NOAH_)
    realtype          t1;
    realtype          snowh;
    realtype          stc[MAXLYR];
    realtype          smc[MAXLYR];
    realtype          sh2o[MAXLYR];
#endif
} ic_struct;

#if defined(_CYCLES_)
typedef struct cyclesic_struct
{
    realtype          resw_stan;
    realtype          resw_flat;
    realtype          resm_stan;
    realtype          resm_flat;
    realtype          manuc_surf;
    realtype          resn_stan;
    realtype          resn_flat;
    realtype          manun_surf;
    realtype          res_abgd[MAXLYR];
    realtype          res_root[MAXLYR];
    realtype          res_rhizo[MAXLYR];
    realtype          manuc[MAXLYR];
    realtype          resn_abgd[MAXLYR];
    realtype          resn_root[MAXLYR];
    realtype          resn_rhizo[MAXLYR];
    realtype          manun[MAXLYR];
    realtype          soc[MAXLYR];
    realtype          son[MAXLYR];
    realtype          mbc[MAXLYR];
    realtype          mbn[MAXLYR];
    realtype          no3[MAXLYR];
    realtype          nh4[MAXLYR];
} cyclesic_struct;
#endif

#if defined(_BGC_)
/* CN initial conditions */
typedef struct bgcic_struct
{
    realtype          leafc;
    realtype          leafc_storage;
    realtype          leafc_transfer;
    realtype          frootc;
    realtype          frootc_storage;
    realtype          frootc_transfer;
    realtype          livestemc;
    realtype          livestemc_storage;
    realtype          livestemc_transfer;
    realtype          deadstemc;
    realtype          deadstemc_storage;
    realtype          deadstemc_transfer;
    realtype          livecrootc;
    realtype          livecrootc_storage;
    realtype          livecrootc_transfer;
    realtype          deadcrootc;
    realtype          deadcrootc_storage;
    realtype          deadcrootc_transfer;
    realtype          gresp_storage;
    realtype          gresp_transfer;
    realtype          cwdc;
    realtype          litr1c;
    realtype          litr2c;
    realtype          litr3c;
    realtype          litr4c;
    realtype          soil1c;
    realtype          soil2c;
    realtype          soil3c;
    realtype          soil4c;
    realtype          cpool;
    realtype          leafn;
    realtype          leafn_storage;
    realtype          leafn_transfer;
    realtype          frootn;
    realtype          frootn_storage;
    realtype          frootn_transfer;
    realtype          livestemn;
    realtype          livestemn_storage;
    realtype          livestemn_transfer;
    realtype          deadstemn;
    realtype          deadstemn_storage;
    realtype          deadstemn_transfer;
    realtype          livecrootn;
    realtype          livecrootn_storage;
    realtype          livecrootn_transfer;
    realtype          deadcrootn;
    realtype          deadcrootn_storage;
    realtype          deadcrootn_transfer;
    realtype          cwdn;
    realtype          litr1n;
    realtype          litr2n;
    realtype          litr3n;
    realtype          litr4n;
    realtype          soil1n;
    realtype          soil2n;
    realtype          soil3n;
    realtype          soil4n;
    realtype          surfn;
    realtype          sminn;
    realtype          retransn;
    realtype          npool;
    realtype          prev_leafc_to_litter;
    realtype          prev_frootc_to_litter;
    realtype          dsr;
    int             dormant_flag;
    int             onset_flag;
    int             onset_counter;
    int             onset_gddflag;
    realtype          onset_fdd;
    realtype          onset_gdd;
    realtype          onset_swi;
    int             offset_flag;
    int             offset_counter;
    realtype          offset_fdd;
    realtype          offset_swi;
} bgcic_struct;
#endif

#if defined(_DAILY_)

/* Daily average variables */
typedef struct daily_struct
{
    int             counter;               /* counter used for averaging */
    int             daylight_counter;      /* counter used for daytime averaging
                                            */
    realtype          avg_sh2o[MAXLYR];      /* daily average unfrozen soil water
                                            * content (m3 m-3) */
    realtype          avg_smc[MAXLYR];       /* daily average unfrozen soil
                                            * moisture content (m3 m-3) */
    realtype          avg_q2d;               /* daily average mixing ratio deficit
                                            * (kg kg-1) */
    realtype          avg_sfcprs;            /* daily average air pressure (Pa) */
    realtype          avg_ch;                /* daily average surface exchange
                                            * coefficient (m s-1) */
    realtype          avg_rc;                /* daily average stomatal resistance
                                            * (s m-1) */
    realtype          avg_albedo;            /* daily average surface albedo
                                            * (-) */
    realtype          tmax;                  /* daily maximum air temperature (K)
                                            */
    realtype          tmin;                  /* daily minimum air temperature (K)
                                            */
    realtype          avg_sfctmp;            /* daily average air temperature (K)
                                            */
    realtype          tday;                  /* daytime average air temperature
                                            * (K) */
    realtype          tnight;                /* nighttime average air temperature
                                            * (K) */
    realtype          avg_stc[MAXLYR];       /* daily average soil temperature (K)
                                            */
    realtype          avg_soldn;             /* daytime average downward solar
                                            * radiation (W m-2) */
#if defined(_PIHM_)
    realtype          avg_et[MAXLYR];        /* daily average evapotranspiration
                                            * (m s-1) */
    realtype          avg_sncovr;            /* daily average snow cover fraction
                                            * (-) */
#endif
} daily_struct;
#endif

#if defined(_BGC_)
/* Carbon state variables (including sums for sources and sinks) */
typedef struct cstate_struct
{
    realtype          leafc;             
    realtype          leafc_storage;     
    realtype          leafc_transfer;    
    realtype          frootc;            
    realtype          frootc_storage;    
    realtype          frootc_transfer;   
    realtype          livestemc;         
    realtype          livestemc_storage; 
    realtype          livestemc_transfer;
    realtype          deadstemc;         
    realtype          deadstemc_storage; 
    realtype          deadstemc_transfer;
    realtype          livecrootc;        
    realtype          livecrootc_storage;
    realtype          livecrootc_transfer
    realtype          deadcrootc;        
    realtype          deadcrootc_storage;
    realtype          deadcrootc_transfer
    realtype          gresp_storage;     
    realtype          gresp_transfer;    
    realtype          cwdc;              
    realtype          litr1c;            
    realtype          litr2c;            
    realtype          litr3c;            
    realtype          litr4c;            
    realtype          soil1c;            
    realtype          soil2c;            
    realtype          soil3c;            
    realtype          soil4c;            
    realtype          cpool;             
    realtype          psnsun_src;        
    realtype          psnshade_src;      
    realtype          leaf_mr_snk;       
    realtype          leaf_gr_snk;       
    realtype          froot_mr_snk;      
    realtype          froot_gr_snk;      
    realtype          livestem_mr_snk;   
    realtype          livestem_gr_snk;   
    realtype          deadstem_gr_snk;   
    realtype          livecroot_mr_snk;  
    realtype          livecroot_gr_snk;  
    realtype          deadcroot_gr_snk;  
    realtype          litr1_hr_snk;      
    realtype          litr2_hr_snk;      
    realtype          litr4_hr_snk;      
    realtype          soil1_hr_snk;
    realtype          soil2_hr_snk;
    realtype          soil3_hr_snk;      
    realtype          soil4_hr_snk;      
    realtype          fire_snk;          
} cstate_struct;

/* Daily carbon flux variables */
typedef struct cflux_struct
{
    /* mortality fluxes (kgC m-2 day-1) */
    realtype          m_leafc_to_litr1c;
    realtype          m_leafc_to_litr2c;
    realtype          m_leafc_to_litr3c;
    realtype          m_leafc_to_litr4c;
    realtype          m_frootc_to_litr1c;
    realtype          m_frootc_to_litr2c;
    realtype          m_frootc_to_litr3c;
    realtype          m_frootc_to_litr4c;
    realtype          m_leafc_storage_to_litr1c;
    realtype          m_frootc_storage_to_litr1c;
    realtype          m_livestemc_storage_to_litr1c;
    realtype          m_deadstemc_storage_to_litr1c;
    realtype          m_livecrootc_storage_to_litr1c;
    realtype          m_deadcrootc_storage_to_litr1c;
    realtype          m_leafc_transfer_to_litr1c;
    realtype          m_frootc_transfer_to_litr1c;
    realtype          m_livestemc_transfer_to_litr1c;
    realtype          m_deadstemc_transfer_to_litr1c;
    realtype          m_livecrootc_transfer_to_litr1c;
    realtype          m_deadcrootc_transfer_to_litr1c;
    realtype          m_livestemc_to_cwdc;
    realtype          m_deadstemc_to_cwdc;
    realtype          m_livecrootc_to_cwdc;
    realtype          m_deadcrootc_to_cwdc;
    realtype          m_gresp_storage_to_litr1c;
    realtype          m_gresp_transfer_to_litr1c;
    /* fire fluxes (kgC m-2 day-1) */
    realtype          m_leafc_to_fire;
    realtype          m_frootc_to_fire;
    realtype          m_leafc_storage_to_fire;
    realtype          m_frootc_storage_to_fire;
    realtype          m_livestemc_storage_to_fire;
    realtype          m_deadstemc_storage_to_fire;
    realtype          m_livecrootc_storage_to_fire;
    realtype          m_deadcrootc_storage_to_fire;
    realtype          m_leafc_transfer_to_fire;
    realtype          m_frootc_transfer_to_fire;
    realtype          m_livestemc_transfer_to_fire;
    realtype          m_deadstemc_transfer_to_fire;
    realtype          m_livecrootc_transfer_to_fire;
    realtype          m_deadcrootc_transfer_to_fire;
    realtype          m_livestemc_to_fire;
    realtype          m_deadstemc_to_fire;
    realtype          m_livecrootc_to_fire;
    realtype          m_deadcrootc_to_fire;
    realtype          m_gresp_storage_to_fire;
    realtype          m_gresp_transfer_to_fire;
    realtype          m_litr1c_to_fire;
    realtype          m_litr2c_to_fire;
    realtype          m_litr3c_to_fire;
    realtype          m_litr4c_to_fire;
    realtype          m_cwdc_to_fire;
    /* phenology fluxes from transfer pools (kgC m-2 day-1) */
    realtype          leafc_transfer_to_leafc;
    realtype          frootc_transfer_to_frootc;
    realtype          livestemc_transfer_to_livestemc;
    realtype          deadstemc_transfer_to_deadstemc;
    realtype          livecrootc_transfer_to_livecrootc;
    realtype          deadcrootc_transfer_to_deadcrootc;
    /* leaf and fine root litterfall (kgC m-2 day-1) */
    realtype          leafc_to_litr1c;
    realtype          leafc_to_litr2c;
    realtype          leafc_to_litr3c;
    realtype          leafc_to_litr4c;
    realtype          frootc_to_litr1c;
    realtype          frootc_to_litr2c;
    realtype          frootc_to_litr3c;
    realtype          frootc_to_litr4c;
    /* maintenance respiration fluxes (kgC m-2 day-1) */
    realtype          leaf_day_mr;
    realtype          leaf_night_mr;
    realtype          froot_mr;
    realtype          livestem_mr;
    realtype          livecroot_mr;
    /* photosynthesis fluxes (kgC m-2 day-1) */
    realtype          psnsun_to_cpool;
    realtype          psnshade_to_cpool;
    /* litter decomposition fluxes (kgC m-2 day-1) */
    realtype          cwdc_to_litr2c;
    realtype          cwdc_to_litr3c;
    realtype          cwdc_to_litr4c;
    realtype          litr1_hr;
    realtype          litr1c_to_soil1c;
    realtype          litr2_hr;
    realtype          litr2c_to_soil2c;
    realtype          litr3c_to_litr2c;
    realtype          litr4_hr;
    realtype          litr4c_to_soil3c;
    realtype          soil1_hr;
    realtype          soil1c_to_soil2c;
    realtype          soil2_hr;
    realtype          soil2c_to_soil3c;
    realtype          soil3_hr;
    realtype          soil3c_to_soil4c;
    realtype          soil4_hr;
    /* daily allocation fluxes from current GPP (kgC m-2 day-1) */
    realtype          cpool_to_leafc;
    realtype          cpool_to_leafc_storage;
    realtype          cpool_to_frootc;
    realtype          cpool_to_frootc_storage;
    realtype          cpool_to_livestemc;
    realtype          cpool_to_livestemc_storage;
    realtype          cpool_to_deadstemc;
    realtype          cpool_to_deadstemc_storage;
    realtype          cpool_to_livecrootc;
    realtype          cpool_to_livecrootc_storage;
    realtype          cpool_to_deadcrootc;
    realtype          cpool_to_deadcrootc_storage;
    realtype          cpool_to_gresp_storage;
    /* daily growth respiration fluxes (kgC m-2 day-1) */
    realtype          cpool_leaf_gr;
    realtype          cpool_leaf_storage_gr;
    realtype          transfer_leaf_gr;
    realtype          cpool_froot_gr;
    realtype          cpool_froot_storage_gr;
    realtype          transfer_froot_gr;
    realtype          cpool_livestem_gr;
    realtype          cpool_livestem_storage_gr;
    realtype          transfer_livestem_gr;
    realtype          cpool_deadstem_gr;
    realtype          cpool_deadstem_storage_gr;
    realtype          transfer_deadstem_gr;
    realtype          cpool_livecroot_gr;
    realtype          cpool_livecroot_storage_gr;
    realtype          transfer_livecroot_gr;
    realtype          cpool_deadcroot_gr;
    realtype          cpool_deadcroot_storage_gr;
    realtype          transfer_deadcroot_gr;
    /* annual turnover of storage to transfer pools (kgC m-2 day-1) */
    realtype          leafc_storage_to_leafc_transfer;
    realtype          frootc_storage_to_frootc_transfer;
    realtype          livestemc_storage_to_livestemc_transfer;
    realtype          deadstemc_storage_to_deadstemc_transfer;
    realtype          livecrootc_storage_to_livecrootc_transfer;
    realtype          deadcrootc_storage_to_deadcrootc_transfer;
    realtype          gresp_storage_to_gresp_transfer;
    /* turnover of live wood to dead wood (kgC m-2 day-1) */
    realtype          livestemc_to_deadstemc;
    realtype          livecrootc_to_deadcrootc;
} cflux_struct;


/* Nitrogen state variables (including sums for sources and sinks) */
typedef struct nstate_struct
{
    /* leaf N (kgN m-2) */
    realtype          leafn;
    realtype          leafn_storage;
    realtype          leafn_transfer;
    /* fine root N (kgN m-2) */
    realtype          frootn;
    realtype          frootn_storage;
    realtype          frootn_transfer;
    /* live stem N (kgN m-2) */
    realtype          livestemn;
    realtype          livestemn_storage;
    realtype          livestemn_transfer;
    /* dead stem N (kgN m-2) */
    realtype          deadstemn;
    realtype          deadstemn_storage;
    realtype          deadstemn_transfer;
    /* live coarse root N (kgN m-2) */
    realtype          livecrootn;
    realtype          livecrootn_storage;
    realtype          livecrootn_transfer;
    /* dead coarse root N (kgN m-2) */
    realtype          deadcrootn;
    realtype          deadcrootn_storage;
    realtype          deadcrootn_transfer;
    realtype          cwdn;         
    realtype          litr1n;       
    realtype          litr2n;       
    realtype          litr3n;       
    realtype          litr4n;       
    realtype          soil1n;       
    realtype          soil2n;       
    realtype          soil3n;       
    realtype          soil4n;       
    realtype          surfn;        
    realtype          sminn;        
    realtype          retransn;     
    realtype          npool;        
    realtype          nfix_src;     
    realtype          ndep_src;     
    realtype          nleached_snk; 
    realtype          nvol_snk;
    realtype          fire_snk;     
} nstate_struct;

/* Daily nitrogen flux variables */
typedef struct nflux_struct
{
    /* mortality fluxes (kgN m-2 day-1) */
    realtype          m_leafn_to_litr1n;
    realtype          m_leafn_to_litr2n;
    realtype          m_leafn_to_litr3n;
    realtype          m_leafn_to_litr4n;
    realtype          m_frootn_to_litr1n;
    realtype          m_frootn_to_litr2n;
    realtype          m_frootn_to_litr3n;
    realtype          m_frootn_to_litr4n;
    realtype          m_leafn_storage_to_litr1n;
    realtype          m_frootn_storage_to_litr1n;
    realtype          m_livestemn_storage_to_litr1n;
    realtype          m_deadstemn_storage_to_litr1n;
    realtype          m_livecrootn_storage_to_litr1n;
    realtype          m_deadcrootn_storage_to_litr1n;
    realtype          m_leafn_transfer_to_litr1n;
    realtype          m_frootn_transfer_to_litr1n;
    realtype          m_livestemn_transfer_to_litr1n;
    realtype          m_deadstemn_transfer_to_litr1n;
    realtype          m_livecrootn_transfer_to_litr1n;
    realtype          m_deadcrootn_transfer_to_litr1n;
    realtype          m_livestemn_to_litr1n;
    realtype          m_livestemn_to_cwdn;
    realtype          m_deadstemn_to_cwdn;
    realtype          m_livecrootn_to_litr1n;
    realtype          m_livecrootn_to_cwdn;
    realtype          m_deadcrootn_to_cwdn;
    realtype          m_retransn_to_litr1n;
    /* fire fluxes (kgN m-2 day-1) */
    realtype          m_leafn_to_fire;
    realtype          m_frootn_to_fire;
    realtype          m_leafn_storage_to_fire;
    realtype          m_frootn_storage_to_fire;
    realtype          m_livestemn_storage_to_fire;
    realtype          m_deadstemn_storage_to_fire;
    realtype          m_livecrootn_storage_to_fire;
    realtype          m_deadcrootn_storage_to_fire;
    realtype          m_leafn_transfer_to_fire;
    realtype          m_frootn_transfer_to_fire;
    realtype          m_livestemn_transfer_to_fire;
    realtype          m_deadstemn_transfer_to_fire;
    realtype          m_livecrootn_transfer_to_fire;
    realtype          m_deadcrootn_transfer_to_fire;
    realtype          m_livestemn_to_fire;
    realtype          m_deadstemn_to_fire;
    realtype          m_livecrootn_to_fire;
    realtype          m_deadcrootn_to_fire;
    realtype          m_retransn_to_fire;
    realtype          m_litr1n_to_fire;
    realtype          m_litr2n_to_fire;
    realtype          m_litr3n_to_fire;
    realtype          m_litr4n_to_fire;
    realtype          m_cwdn_to_fire;
    /* phenology fluxes from transfer pools (kgN m-2 day-1) */
    realtype          leafn_transfer_to_leafn;
    realtype          frootn_transfer_to_frootn;
    realtype          livestemn_transfer_to_livestemn;
    realtype          deadstemn_transfer_to_deadstemn;
    realtype          livecrootn_transfer_to_livecrootn;
    realtype          deadcrootn_transfer_to_deadcrootn;
    /* litterfall fluxes (kgN m-2 day-1) */
    realtype          leafn_to_litr1n;
    realtype          leafn_to_litr2n;
    realtype          leafn_to_litr3n;
    realtype          leafn_to_litr4n;
    realtype          leafn_to_retransn;
    realtype          frootn_to_litr1n;
    realtype          frootn_to_litr2n;
    realtype          frootn_to_litr3n;
    realtype          frootn_to_litr4n;
    /* decomposition fluxes (kgN m-2 day-1) */
    realtype          ndep_to_sminn;
    realtype          nfix_to_sminn;
    /* litter and soil decomposition fluxes (kgN m-2 day-1) */
    realtype          cwdn_to_litr2n;
    realtype          cwdn_to_litr3n;
    realtype          cwdn_to_litr4n;
    realtype          litr1n_to_soil1n;
    realtype          sminn_to_soil1n_l1;
    realtype          litr2n_to_soil2n;
    realtype          sminn_to_soil2n_l2;
    realtype          litr3n_to_litr2n;
    realtype          litr4n_to_soil3n;
    realtype          sminn_to_soil3n_l4;
    realtype          soil1n_to_soil2n;
    realtype          sminn_to_soil2n_s1;
    realtype          soil2n_to_soil3n;
    realtype          sminn_to_soil3n_s2;
    realtype          soil3n_to_soil4n;
    realtype          sminn_to_soil4n_s3;
    realtype          soil4n_to_sminn;
    /* denitrification (volatilization) fluxes (kgN m-2 day-1) */
    realtype          sminn_to_nvol_l1s1;
    realtype          sminn_to_nvol_l2s2;
    realtype          sminn_to_nvol_l4s3;
    realtype          sminn_to_nvol_s1s2;
    realtype          sminn_to_nvol_s2s3;
    realtype          sminn_to_nvol_s3s4;
    realtype          sminn_to_nvol_s4;
    realtype          sminn_to_denitrif;
    /* leaching flux (kgN m-2 s-1) */
    realtype          sminn_leached;
    /* daily allocation fluxes from (kgN m-2 day-1) */
    realtype          retransn_to_npool;
    realtype          sminn_to_npool;
    realtype          npool_to_leafn;
    realtype          npool_to_leafn_storage;
    realtype          npool_to_frootn;
    realtype          npool_to_frootn_storage;
    realtype          npool_to_livestemn;
    realtype          npool_to_livestemn_storage;
    realtype          npool_to_deadstemn;
    realtype          npool_to_deadstemn_storage;
    realtype          npool_to_livecrootn;
    realtype          npool_to_livecrootn_storage;
    realtype          npool_to_deadcrootn;
    realtype          npool_to_deadcrootn_storage;
    /* annual turnover of storage to transfer pools (kgN m-2 day-1) */
    realtype          leafn_storage_to_leafn_transfer;
    realtype          frootn_storage_to_frootn_transfer;
    realtype          livestemn_storage_to_livestemn_transfer;
    realtype          deadstemn_storage_to_deadstemn_transfer;
    realtype          livecrootn_storage_to_livecrootn_transfer;
    realtype          deadcrootn_storage_to_deadcrootn_transfer;
    /* turnover of live wood to dead wood, with retranslocation (kgN m-2 day-1)
     */
    realtype          livestemn_to_deadstemn;
    realtype          livestemn_to_retransn;
    realtype          livecrootn_to_deadcrootn;
    realtype          livecrootn_to_retransn;
} nflux_struct;

/* Temporary nitrogen variables for reconciliation of decomposition
 * immobilization fluxes and plant growth N demands */
typedef struct ntemp_struct
{
    realtype          surfn0;           
    realtype          sminn0;           
    realtype          mineralized;      
    realtype          potential_immob;  
    realtype          plitr1c_loss;     
    realtype          pmnf_l1s1;        
    realtype          plitr2c_loss;     
    realtype          pmnf_l2s2;        
    realtype          plitr4c_loss;     
    realtype          pmnf_l4s3;        
    realtype          psoil1c_loss;     
    realtype          pmnf_s1s2;        
    realtype          psoil2c_loss;     
    realtype          pmnf_s2s3;        
    realtype          psoil3c_loss;     
    realtype          pmnf_s3s4;        
    realtype          psoil4c_loss;     
    realtype          kl4;              
} ntemp_struct;

/* Ecophysiological variables */
typedef struct epvar_struct
{
    realtype          bg_leafc_litfall_rate;  
    realtype          bg_frootc_litfall_rate; 
    realtype          livestemc_turnover_rate;
    realtype          livecrootc_turnover_rate
    realtype          dsr;                    
    realtype          sun_proj_sla;           
    realtype          shade_proj_sla;         
    realtype          psi;                    
    realtype          dlmr_area_sun;          
    realtype          dlmr_area_shade;        
    realtype          gl_t_wv_sun;            
    realtype          gl_t_wv_shade;          
    realtype          assim_sun;              
    realtype          assim_shade;            
    realtype          t_scalar;               
    realtype          w_scalar;               
    realtype          rate_scalar;            
    realtype          daily_gross_nimmob;     
    realtype          daily_net_nmin;         
    realtype          fpi;                    
    realtype          m_tmin;                 
    realtype          m_psi;                  
    realtype          m_co2;                  
    realtype          m_ppfd_sun;             
    realtype          m_ppfd_shade;           
    realtype          m_vpd;                  
    realtype          m_final_sun;            
    realtype          m_final_shade;          
    realtype          ytd_maxplai;            
    int               dormant_flag;         
    realtype          days_active;            
    int               onset_flag; 
    int               onset_counter;
    int               onset_gddflag;
    realtype          onset_fdd;              
    realtype          onset_gdd;              
    realtype          onset_swi;              
    int               offset_flag;
    int               offset_counter; 
    realtype          offset_fdd;             
    realtype          offset_swi;             
    realtype          annavg_t2m;             
    realtype          gpp;                    
    realtype          prev_leafc_to_litter;   
    realtype          prev_frootc_to_litter;  
    realtype          old_c_balance;          
    realtype          old_n_balance;          
    realtype          dayl;                   
    realtype          prev_dayl;              
} epvar_struct;

/* Structure for the photosynthesis routine */
typedef struct psn_struct
{
    int               c3;    
    realtype          pa;  
    realtype          co2; 
    realtype          t;   
    realtype          lnc; 
    realtype          flnr;
    realtype          ppfd;
    realtype          g;   
    realtype          dlmr;
    realtype          Ci;  
    realtype          O2;  
    realtype          Ca;  
    realtype          gamma
    realtype          Kc;  
    realtype          Ko;  
    realtype          Vmax;
    realtype          Jmax;
    realtype          J;   
    realtype          Av;  
    realtype          Aj;  
    realtype          A;   
} psn_struct;

/* CN summary structure */
typedef struct summary_struct
{
    realtype          daily_npp;     
    realtype          daily_nep;     
    realtype          daily_nee;     
    realtype          daily_gpp;     
    realtype          daily_mr;      
    realtype          daily_gr;      
    realtype          daily_hr;      
    realtype          daily_fire;    
    realtype          daily_litfallc;
    /* summed over entire simulation (kgC m-2) */
    realtype          cum_npp;
    realtype          cum_nep;
    realtype          cum_nee;
    realtype          cum_gpp;
    realtype          cum_mr;
    realtype          cum_gr;
    realtype          cum_hr;
    realtype          cum_fire;
    realtype          vegc;           
    realtype          agc;            
    realtype          litrc;          
    realtype          soilc;          
    realtype          totalc;         
} summary_struct;

/* Solute transport structure */
typedef struct solute_struct
{
# if !defined(_LEACHING_)
    realtype          conc_surf;          
    realtype          infilflux;          
    realtype          ovlflux[NUM_EDGE];  
# endif
    realtype          conc_subsurf;       
    realtype          subflux[NUM_EDGE];  
    realtype          snksrc;             
} solute_struct;
#endif

/* Spinup variables */
typedef struct spinup_struct
{
    realtype          totalw_prev;
    realtype          totalw;
#if defined(_BGC_)
    realtype          soilc_prev;
    realtype          totalc_prev;
    realtype          soilc;
    realtype          totalc;
#endif
} spinup_struct;

#endif  /* _DEBUG_ */
    int *elem_d_node[NUM_EDGE];     /* Element structure */
	int *elem_d_nabr[NUM_EDGE]; 
	int *elem_d_ind;                /* Element structure */

/* river_d related arguments */           
    int *river_d_attrib_riverbc_type;   /* river boundary condition type */
    realtype *river_d_topo_area;        /* River topography parameters */
	realtype *river_d_topo_x; 
	realtype *river_d_topo_y;                 
    realtype *river_d_topo_zmin;
	realtype *river_d_topo_zmax; 
	realtype *river_d_topo_zbed;                            
    realtype *river_d_topo_node_zmax; 
	realtype *river_d_topo_dist_left; 
	realtype *river_d_topo_dist_right;    /* River topography parameters */ 
    realtype *river_d_ws_stage; 
	realtype *river_d_ws_gw;              /* River water states */
    realtype *river_d_wf_rivflow[NUM_RIVFLX];  /* River water fluxes */   
    realtype *river_d_shp_depth;       /* River shape parameters */
	int      *river_d_shp_intrpl_ord;  
	realtype *river_d_shp_coeff;     
    realtype *river_d_shp_length; 
	realtype *river_d_shp_width;       /* River shape parameters */     
    realtype *river_d_matl_rough; 
	realtype *river_d_matl_cwr; 
	realtype *river_d_matl_ksath;         
    realtype *river_d_matl_ksatv; 
	realtype *river_d_matl_bedthick; 
	realtype *river_d_matl_porosity;     
    realtype *river_d_matl_smcmin;   
#if defined(_CYCLES_)
    realtype *river_d_matl_bd;       /* River mathmatic parameters */	
#endif	
    realtype *river_d_bc_head; 
	realtype *river_d_bc_flux;       /* River boundary conditions */  
	
/* River initial conditions */

#if defined(_BGC_) && !defined(_LUMPED_) && !defined(_LEACHING_)
/* River nitrogen state variables */
typedef struct river_nstate_struct
{
    realtype          streamn;    /* stream N pool (kgN m-2) */
    realtype          sminn;      /* river bed soil mineral N (kgN m-2) */
} river_nstate_struct;

/* Daily river nitrogen flux variables */
typedef struct river_nflux_struct
{
    realtype          sminn_leached;    /* leaching flux (kgN m-2 day-1) */
} river_nflux_struct;

/* River solute transport structure */
typedef struct river_solute_struct
{
    realtype          conc_stream;         /* stream pool concentration
                                          * (kg kgH2O-1) */
    realtype          conc_bed;            /* bed pool concentration (kg kgH2O-1)
                                          */
    realtype          flux[NUM_RIVFLX];    /* solute fluxes (kg s-1) */
} river_solute_struct;

/* River CN initial conditions */
typedef struct river_bgcic_struct
{
    realtype          streamn;
    realtype          sminn;
} river_bgcic_struct;
#endif

#if defined(_CYCLES_)
typedef struct river_nstate_struct
{
    realtype          streamno3;
    realtype          streamnh4;
    realtype          bedno3;
    realtype          bednh4;
} river_nstate_struct;

typedef struct river_cyclesic_struct
{
    realtype          streamno3;
    realtype          streamnh4;
    realtype          bedno3;
    realtype          bednh4;
} river_cyclesic_struct;

typedef struct river_solute_struct
{
    realtype          conc_stream;         /* stream pool concentration
                                          * (kg kgH2O-1) */
    realtype          conc_bed;            /* bed pool concentration (kg kgH2O-1)
                                          */
    realtype          flux[NUM_RIVFLX];    /* solute fluxes (kg s-1) */
} river_solute_struct;
#endif
	
    int *river_d_ind;         /* River structure */
	int *river_d_leftele; 
	int *river_d_rightele;      
    int *river_d_fromnode; 
	int *river_d_tonode; 
	int *river_d_down;        /* River structure */
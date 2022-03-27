#ifndef ELEM_STRUCT_HEADER
#define ELEM_STRUCT_HEADER

/* Element attribute */
typedef struct attrib_struct
{
    int             soil_type;              /* element soil type */
#if defined(_FBR_)
    int             geol_type;              /* element geology type */
#endif
    int             lc_type;                /* element land cover type */
    int             bc_type[NUM_EDGE];      /* element boundary condition type
                                             */
#if defined(_FBR_)
    int             fbrbc_type[NUM_EDGE];   /* element fractured bedrock layer
                                             * boundary condition type */
#endif
    int             meteo_type;             /* element meteorological forcing
                                             * type */
    int             lai_type;               /* element leaf area index forcing
                                             * type */
#if defined(_CYCLES_)
    int             op_type;
#endif
} attrib_struct;

/* Topography parameters */
typedef struct topo_struct
{
    realtype          area;                  /* area of element (m2) */
    realtype          x;                     /* x of centroid (m) */
    realtype          y;                     /* y of centroid (m) */
    realtype          zmin;                  /* soil bottom elevation (m) */
    realtype          zmax;                  /* surface elevation (m) */
#if defined(_FBR_)
    realtype          zbed;                  /* impermeable bedrock elevation (m)
                                            */
#endif
    realtype          edge[NUM_EDGE];        /* length of edge (Edge i is from
                                            * node i to node i + 1) (m) */
    realtype          nabrdist[NUM_EDGE];    /* distance to neighbor (m) */
    realtype          nabr_x[NUM_EDGE];      /* x of neighbor centroid (m) */
    realtype          nabr_y[NUM_EDGE];      /* y of neighbor centroid (m) */
#if defined(_NOAH_)
    realtype          slope;                 /* slope of element (degree) */
    realtype          aspect;                /* surface aspect of element (degree)
                                            */
    realtype          svf;                   /* sky view factor (-) */
    realtype          h_phi[36];             /* unobstructed angle in each
                                            * direction (degree) */
#endif
#if defined(_RT_)
    realtype          areasub[NUM_EDGE];
#endif
} topo_struct;

/* Soil parameters */
typedef struct soil_struct
{
    realtype          depth;       /* soil depth (m) */
    realtype          ksath;       /* horizontal saturated hydraulic conductivity
                                  * (m s-1) */
    realtype          ksatv;       /* vertical saturated hydraulic conductivity
                                  * (m s-1) */
    realtype          kinfv;       /* saturated infiltration conductivity (m s-1)
                                  */
    realtype          dinf;        /* depth from ground surface across which head
                                  * gradient is calculated for infiltration (m)
                                  */
    realtype          alpha;       /* alpha from van Genuchten eqn (m-1) */
    realtype          beta;        /* beta (n) from van Genuchten eqn (-) */
    realtype          porosity;    /* soil porosity (m3 m-3) */
    realtype          smcmax;      /* maximum soil moisture content (m3 m-3) */
    realtype          smcmin;      /* residual soil moisture content (m3 m-3) */
    realtype          smcwlt;      /* wilting point (m3 m-3) */
    realtype          smcref;      /* soil moisture threshold where transpiration
                                  * begins to stress (m3 m-3) */
    realtype          dmac;        /* macropore depth (m) */
    realtype          kmach;       /* macropore horizontal saturated hydraulic
                                  * conductivity (m s-1) */
    realtype          kmacv;       /* macropore vertical saturated hydraulic
                                  * conductivity (m s-1) */
    realtype          areafv;      /* macropore area fraction on a vertical cross-
                                  * section (m2 m-2) */
    realtype          areafh;      /* macropore area fraction on a horizontal
                                  * cross-section (m2 m-2) */
#if defined(_NOAH_)
    realtype          csoil;       /* soil heat capacity (J m-3 K-1) */
    realtype          quartz;      /* soil quartz content (-) */
    realtype          smcdry;      /* dry soil moisture threshold where direct
                                  * evap from top layer ends (m3 m-3) */
#endif
#if defined(_CYCLES_)
    realtype          clay[MAXLYR];
    realtype          sand[MAXLYR];
    realtype          iom[MAXLYR];
    realtype          bd[MAXLYR];
#endif
} soil_struct;

#if defined(_FBR_)
/* Fractured bedrock layer parameters */
typedef struct geol_struct
{
    realtype          depth;       /* soil depth (m) */
    realtype          ksath;       /* horizontal saturated hydraulic conductivity
                                  * (m s-1) */
    realtype          ksatv;       /* vertical saturated hydraulic conductivity
                                  * (m s-1) */
    realtype          alpha;       /* alpha from van Genuchten eqn (m-1) */
    realtype          beta;        /* beta (n) from van Genuchten eqn (-) */
    realtype          porosity;    /* porosity (m3 m-3) */
    realtype          smcmax;      /* maximum moisture content (m3 m-3) */
    realtype          smcmin;      /* residual moisture content (m3 m-3) */
} geol_struct;
#endif

/* Land cover parameters */
typedef struct lc_struct
{
    realtype          shdfac;       /* areal fractional coverage of green
                                   * vegetation (0.0-1.0) (-) */
    realtype          shdmin;       /* minimum areal fractional coverage of green
                                   * vegetation (-) */
    realtype          shdmax;       /* maximum areal fractional coverage of green
                                   * vegetation (-) */
    realtype          laimin;       /* minimum LAI across all seasons for a
                                   * vegetation type (m2 m-2) */
    realtype          laimax;       /* maximum LAI across all seasons for a
                                   * vegetation type (m2 m-2) */
    realtype          snup;         /* threshold snow depth (in water equivalent)
                                   * that implies 100% snow cover (m) */
    realtype          cfactr;       /* parameter used in the canopy interception
                                   * calculation (-) */
    realtype          emissmax;     /* minimum emissivity (-) */
    realtype          emissmin;     /* maximum emissivity (-) */
    realtype          albedomax;    /* minimum background albedo (-) */
    realtype          albedomin;    /* maximum background albedo (-) */
    realtype          z0max;        /* minimum roughness length (m) */
    realtype          z0min;        /* maximum roughness length (m) */
    realtype          rough;        /* surface roughness (Manning's n) (s m-1/3)
                                   */
    realtype          cmcfactr;     /* canopy water capacity per LAI (m) */
    int             bare;         /* flag that indicates bare ground */
    int             isurban;      /* flag that indicates urban */
} lc_struct;

#if !defined(_CYCLES_)
/* Ecophysiological parameters */
typedef struct epconst_struct
{
    realtype          rsmin;                  /* minimum canopy resistance (s m-1)
                                             */
    realtype          rgl;                    /* reference incoming solar flux
                                             * for photosynthetically active
                                             * canopy (W m-2) */
    realtype          hs;                     /* parameter used in vapor pressure
                                             * deficit function (-) */
    realtype          topt;                   /* optimum transpiration air
                                             * temperature (K) */
    realtype          rsmax;                  /* cuticular resistance (s m-1) */
# if defined(_BGC_)
    int             woody;                  /* flag: 1 = woody, 0 = non-woody */
    int             evergreen;              /* flag: 1 = evergreen,
                                             * 0 = deciduous */
    int             c3_flag;                /* flag: 1 = C3,  0 = C4 */
    int             phenology_flag;         /* flag: 1 = phenology mode
                                             * 0 = user defined */
    int             onday;                  /* day of year when leaves on */
    int             offday;                 /* day of year when leaves off */
    int             transfer_days;          /* growth period for transfer (day)
                                             */
    int             litfall_days;           /* growth period for litter fall
                                             * (day) */
    realtype          leaf_turnover;          /* annual leaf turnover fraction
                                             * (yr-1) */
    realtype          froot_turnover;         /* annual fine root turnover
                                             * fraction (yr-1) */
    realtype          livewood_turnover;      /* annual live wood turnover
                                             * fraction (yr-1) */
    realtype          daily_mortality_turnover;/* daily mortality turnover (day-1)
                                              */
    realtype          daily_fire_turnover;    /* daily fire turnover (day-1) */
    realtype          alloc_frootc_leafc;     /* new fine root C to new leaf C (-)
                                             */
    realtype          alloc_newstemc_newleafc;/* new stem C to new leaf C (-) */
    realtype          alloc_newlivewoodc_newwoodc;/* new livewood C:new wood C (-)
                                                 */
    realtype          alloc_crootc_stemc;     /* new live croot C to new live stem
                                             * C (-) */
    realtype          alloc_prop_curgrowth;   /* daily allocation to current
                                             * growth (-) */
    realtype          avg_proj_sla;           /* canopy average projected SLA
                                             * (m2 kgC-1) */
    realtype          sla_ratio;              /* ratio of shaded to sunlit
                                             * projected SLA (-) */
    realtype          lai_ratio;              /* ratio of (all-sided LA /
                                             * one-sided LA) (-) */
    realtype          ext_coef;               /* canopy light extinction
                                             * coefficient (-) */
    realtype          flnr;                   /* leaf N in Rubisco
                                             * (kgNRub kgNleaf-1) */
    realtype          psi_open;               /* psi at start of conductance
                                             * reduction (MPa) */
    realtype          psi_close;              /* psi at complete conductance
                                             * reduction (MPa) */
    realtype          vpd_open;               /* vpd at start of conductance
                                             * reduction (Pa) */
    realtype          vpd_close;              /* vpd at complete conductance
                                             * reduction (Pa) */
    realtype          froot_cn;               /* C:N for fine roots (kgC kgN-1) */
    realtype          leaf_cn;                /* C:N for leaves (kgC kgN-1) */
    realtype          livewood_cn;            /* C:N for live wood (kgC kgN-1) */
    realtype          deadwood_cn;            /* C:N for dead wood (kgC kgN-1) */
    realtype          leaflitr_cn;            /* constant C:N for leaf litter
                                             * (kgC kgN-1) */
    realtype          leaflitr_flab;          /* leaf litter labile fraction (-)
                                             */
    realtype          leaflitr_fucel;         /* leaf litter unshielded cellulose
                                             * fraction (-) */
    realtype          leaflitr_fscel;         /* leaf litter shielded cellulose
                                             * fraction (-) */
    realtype          leaflitr_flig;          /* leaf litter lignin fraction (-)
                                             */
    realtype          frootlitr_flab;         /* fine root litter labile fraction
                                             * (-) */
    realtype          frootlitr_fucel;        /* fine root litter unshielded
                                             * cellulose fraction (-) */
    realtype          frootlitr_fscel;        /* fine root litter shielded
                                             * cellulose fraction (-) */
    realtype          frootlitr_flig;         /* fine root litter lignin fraction
                                             * (-) */
    realtype          deadwood_fucel;         /* dead wood unshielded cellulose
                                             * fraction (-) */
    realtype          deadwood_fscel;         /* dead wood shielded cellulose
                                             * fraction (-) */
    realtype          deadwood_flig;          /* dead wood lignin fraction (-) */
# endif
} epconst_struct;
#endif

/* Physical states */
typedef struct pstate_struct
{
    realtype          rzd;                    /* rooting depth (m) */
    realtype          rc;                     /* canopy resistance (s m-1) */
    realtype          pc;                     /* plant coefficient (-) */
    realtype          proj_lai;               /* live projected leaf area index
                                             * (m2 m-2) */
    realtype          rcs;                    /* incoming solar rc factor (-) */
    realtype          rct;                    /* air temperature rc factor (-) */
    realtype          rcq;                    /* vapor pressure deficit rc factor
                                             * (-) */
    realtype          rcsoil;                 /* soil moisture rc factor (-) */
    realtype          albedo;                 /* surface albedo including snow
                                             * effect (-) */
    realtype          zlvl;                   /* height above ground of
                                             * atmospheric forcing variables (m)
                                             */
    realtype          zlvl_wind;              /* height above ground of wind
                                             * observations (m) */
    realtype          sfcspd;                 /* wind speed at height zlvl above
                                             * ground (m s-1) */
    realtype          rh;                     /* relative humidity (100%) */
    realtype          sfcprs;                 /* surface pressure at height zlvl
                                             * above ground (Pa) */
#if defined(_NOAH_)
    realtype          alb;                    /* background snow-free surface
                                             * albedo (-) */
    realtype          snoalb;                 /* upper bound on maximum albedo
                                             * over deep snow (-) */
    int               nroot;                  /* number of root layers, a function
                                             * of vegetation type */
    realtype          rtdis[MAXLYR];          /* root distribution (-) */
    int               nsoil;                  /* number of soil layers */
    realtype          sldpth[MAXLYR];         /* thickness of each soil layer (m)
                                             */
    realtype          zsoil[MAXLYR];          /* distance from land surface to
                                             * bottom of each soil layer (m) */
    realtype          soilw;                  /* available soil moisture in root
                                             * zone (fraction between smcwlt and
                                             * smcmax) (-) */
    realtype          frzk;                   /* frozen ground parameter (-) */
    realtype          frzx;                   /* adjusted frozen ground parameter
                                             * (-) */
    realtype          czil;                   /* Zilitinkevich constant (-) */
    realtype          emissi;                 /* surface emissivity (between 0 and
                                             * 1) (-) */
    realtype          ch;                     /* surface exchange coefficient for
                                             * heat and moisture (m s-1) */
    realtype          cm;                     /* surface exchange coefficient for
                                             * momentum (m s-1) */
    realtype          rch;                    /* = ch * air density * CP
                                             * (W m-2 K-1) */
    realtype          z0;                     /* time varying roughness length as
                                             * function of snow depth (-) */
    realtype          fcr;                    /* reduction of infiltration caused
                                             * by frozen ground (-) */
    int             nmacd;                  /* number of soil layers with
                                             * macropore */
    realtype          salp;                   /* shape parameter of distribution
                                             * function of snow cover (-) */
    realtype          fxexp;                  /* soil evaporation exponent used
                                             * in direct evaporation (-) */
    realtype          sbeta;                  /* parameter used to calculate
                                             * vegetation effect on soil heat
                                             * (-) */
    realtype          lvcoef;                 /* parameter controls surface snow
                                             * albedo in the presence of snow
                                             * cover (-) */
    realtype          snotime1;               /* age of the snow on the ground (s)
                                             */
    realtype          ribb;                   /* bulk Richardson number used to
                                             * limit the dew/frost (-) */
    realtype          beta;                   /* ratio of actual/potential evap
                                             * (-) */
    realtype          sncovr;                 /* fractional snow cover (-) */
    realtype          q1;                     /* effective mixing ratio at surface
                                             * (kg kg-1) */
    realtype          q2;                     /* mixing ratio at height zlvl above
                                             * (kg kg-1) */
    realtype          ffrozp;                 /* fraction of frozen precipitation
                                             * (-) */
    realtype          z0brd;                  /* background fixed roughness length
                                             * (-) */
    realtype          embrd;                  /* background surface emissivity
                                             * (-) */
    realtype          q2sat;                  /* saturation air humidity at height
                                             * zlvl above ground (kg kg-1) */
    realtype          q2d;                    /* air humidity deficit (kg kg-1) */
    realtype          dqsdt2;                 /* slope of saturation specific
                                             * humidity curve at T = sfctmp
                                             * (kg kg-1 K-1) */
    int             nwtbl;                  /* layer where water table is within
                                             */
    realtype          sndens;                 /* snow density (dimensionless
                                             * fraction of H2O density) (-) */
    realtype          snowh;                  /* actual snow depth (m) */
    realtype          sncond;                 /* snow thermal conductivity
                                             * (W m-1 K-1) */
    realtype          rr;                     /* parameter in Penman potential
                                             * evaporation (-) */
    realtype          epsca;                  /* parameter in Penman potential
                                             * evaporation (K) */
    realtype          eta_kinematic;          /* actual latent heat flux
                                             * (kg m-2 s-1) */
    realtype          zbot;                   /* depth of lower boundary soil
                                             * temperature (m) */
    realtype          tbot;                   /* bottom soil temperature (local
                                             * yearly-mean sfc air temperature)
                                             * (K) */
    realtype          gwet;                   /* fraction of transpiration from
                                             * groundwater (-) */
    realtype          satdpth[MAXLYR];        /* depth of groundwater in each soil
                                             * layer (m) */
#endif
#if defined(_CYCLES_)
    realtype          res_intcp;
    realtype          tau_res_stan;
    realtype          tau_res_flat;
    realtype          till_factr[MAXLYR];
    realtype          comp_factr[MAXLYR];
#endif
#if defined(_BGC_)
    realtype          co2;                    /* atmospheric CO2 concentration
                                             * (ppm) */
    realtype          ppfd_per_plaisun;       /* ppfd per unit sunlit proj LAI
                                             * (umol m-2 s-1) */
    realtype          ppfd_per_plaishade;     /* ppfd per unit shaded proj LAI
                                             * (umol m-2 s-1) */
    realtype          all_lai;                /* live all-sided leaf area index
                                             * (m2 m-2) */
    realtype          plaisun;                /* sunlit projected leaf area index
                                             * (m2 m-2) */
    realtype          plaishade;              /* shaded projected leaf area index
                                             * (m2 m-2) */
#endif
} pstate_struct;

/* Water states */
typedef struct wstate_struct
{
    realtype          surf;            /* equivalent surface water level (m) */
    realtype          unsat;           /* unsaturated zone water storage (m) */
    realtype          gw;              /* groundwater level (m) */
    realtype          sneqv;           /* liquid water-equivalent snow depth (m) */                                    
    realtype          cmcmax;          /* maximum canopy water capacity (m) */
    realtype          cmc;             /* interception storage (m) */
    realtype          surfh;           /* actual surface water level (m) */
	
#if defined(_FBR_)
    realtype          fbr_unsat;       /* unsaturated storage in fractured bedrock layer (m) */
    realtype          fbr_gw;          /* deep groundwater in fractured bedrock layer (m) */
#endif
#if defined(_NOAH_)
    realtype          smc[MAXLYR];     /* total soil moisture content (m3 m-3) */
    realtype          sh2o[MAXLYR];    /* unfrozen soil moisture content (m3 m-3) */                                    
    realtype          soilm;           /* total soil column moisture content (m) */                                     
#endif
#if defined(_CYCLES_)
    /* wstate variables in Cycles have the units of kg m-2 */
    realtype          stanResidueWater;   /* (kg m-2) */
    realtype          flatResidueWater;   /* (kg m-2) */
#endif
} wstate_struct;

/* Water fluxes */
typedef struct wflux_struct
{
    realtype          ovlflow[NUM_EDGE];      /* overland flow (m3 s-1) */
    realtype          subsurf[NUM_EDGE];      /* subsurface flow (m3 s-1) */
    realtype          prcp;                   /* precipitation on each element
                                             * (m s-1) */
    realtype          pcpdrp;                 /* combined prcp and drip (from
                                             * canopy) that goes into the soil
                                             * (m s-1) */
    realtype          infil;                  /* variable infiltration rate
                                             * (m s-1) */
    realtype          rechg;                  /* recharge rate to groundwater
                                             * (m s-1) */
    realtype          drip;                   /* through-fall of precipitation
                                             * and/or dew (m s-1) */
    realtype          edir;                   /* direct soil evaporation (m s-1)
                                             */
    realtype          ett;                    /* total plant transpiration (m s-1)
                                             */
    realtype          ec;                     /* canopy water evaporation (m s-1)
                                             */
    realtype          etp;                    /* potential evaporation (m s-1) */
    realtype          eta;                    /* actual evapotranspiration (m s-1)
                                             */
    realtype          edir_surf;              /* direct evaporation from surface
                                             * water (m s-1) */
    realtype          edir_unsat;             /* direct evaporation from
                                             * unsaturated zone (m s-1) */
    realtype          edir_gw;                /* direct evaporation from saturated
                                             * zone (m s-1) */
    realtype          ett_unsat;              /* transpiration from unsaturated
                                             * zone (m s-1) */
    realtype          ett_gw;                 /* transpiration from saturated zone
                                             * (m s-1) */
    realtype          esnow;                  /* sublimation from (or deposition
                                             * to) snowpack (m s-1); */
#if defined(_FBR_)
    realtype          fbr_infil;              /* fractured bedrock infiltration
                                             * (m s-1) */
    realtype          fbr_rechg;              /* fractured bedrock recharge
                                             * (m s-1) */
    realtype          fbrflow[NUM_EDGE];      /* lateral fractured bedrock flow
                                             * (m3 s-1) */
#endif
#if defined(_NOAH_)
    realtype          et[MAXLYR];             /* plant transpiration from each
                                             * soil layer (m s-1) */
    realtype          runoff2;                /* total subsurface flow (m s-1) */
    realtype          runoff2_lyr[MAXLYR];    /* subsurface flow from each soil
                                             * layer (m s-1) */
    realtype          runoff3;                /* numerical truncation in excess
                                             * of porosity (smcmax) for a given
                                             * soil layer at the end of a time
                                             * step (m s-1) */
    realtype          smflxv[MAXLYR];         /* vertical soil moisture flux
                                             * between soil layers (m s-1) */
    realtype          smflxh[NUM_EDGE][MAXLYR];/* horizontal soil moisture flux at
                                             * each soil layer from each edge
                                             * (m s-1) */
    realtype          dew;                    /* dewfall (or frostfall for
                                             * T < 273.15) (m s-1) */
    realtype          snomlt;                 /* water equivalent snow melt
                                             * (m s-1) */
    realtype          etns;                   /* (m s-1) */
#endif
#if defined(_CYCLES_)
    realtype          eres;                   /* evaporation from residue (m s-1)
                                             */
    realtype          irrigationVol;          /* irrigation volume (m s-1) */
#endif
} wflux_struct;

/* Energy states */
typedef struct estate_struct
{
    realtype          sfctmp;         /* air temperature at height zlvl above
                                     * ground (K) */
#if defined(_NOAH_)
    realtype          t1;             /* ground/canopy/snowpack effective skin
                                     * temperature (K) */
    realtype          th2;            /* air potential temperature at height zlvl
                                     * above ground (K) */
    realtype          stc[MAXLYR];    /* soil temperature (K) */
#endif
} estate_struct;

/* Energy fluxes */
typedef struct eflux_struct
{
    realtype          soldn;                  /* solar downward radiation (W m-2)
                                             */
#if defined(_NOAH_)
    realtype          solnet;                 /* net downward solar radiation
                                             * (W m-2) */
    realtype          etp;                    /* potential evaporation (W m-2) */
    realtype          ssoil;                  /* soil heat flux (W m-2) */
    realtype          eta;                    /* actual latent heat flux (W m-2)
                                             */
    realtype          sheat;                  /* sensible heat flux (W m-2) */
    realtype          fdown;                  /* radiation forcing at the surface
                                             * (W m-2) */
    realtype          lwdn;                   /* absorbed longwave downward
                                             * radiation (W m-2) */
    realtype          ec;                     /* canopy water evaporation (W m-2)
                                             */
    realtype          edir;                   /* direct soil evaporation (W m-2)
                                             */
    realtype          et[MAXLYR];             /* plant transpiration from each
                                             * soil layer (W m-2) */
    realtype          ett;                    /* total plant transpiration (W m-2)
                                             */
    realtype          esnow;                  /* sublimation from (or deposition)
                                             * snowpack (W m-2) */
    realtype          soldir;                 /* direct solar radiation (W m-2) */
    realtype          soldif;                 /* diffused solar radiation (W m-2)
                                             */
    realtype          longwave;               /* longwave radiation forcing
                                             * (W m-2) */
    realtype          flx1;                   /* latent heat flux from
                                             * precipitation accumulating as
                                             * snow (W m-2) */
    realtype          flx2;                   /* freezing rain latent heat flux
                                             * (W m-2) */
    realtype          flx3;                   /* snow melt latent heat flux
                                             * (W m-2) */
#endif
#if defined(_BGC_)
    realtype          swabs_per_plaisun;      /* swabs per unit sunlit proj LAI
                                             * (W m-2) */
    realtype          swabs_per_plaishade;    /* swabs per unit shaded proj LAI
                                             * (W m-2) */
#endif
} eflux_struct;

#if defined(_CYCLES_)
typedef struct epvar_struct
{
    /* User Defined Auto Irrigation */
    int             ai_used;
    int             ai_start;
    int             ai_stop;
    realtype          ai_h2o_depl;
    int             ai_last_lyr;

    /* User Defined Auto Fertilization */
    int             auto_fert;

    realtype          plant_density;
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
    realtype          shoot_growth_unstr_cum;          /* kg m-2 day-1 */
    realtype          n_stress_cum;
    realtype          rad_intcp_nc;
    int             harv_date_final;
    int             nharv;
    int             stage;
    realtype          harv_biomass;            /* kg m-2 */
    realtype          harv_root;            /* kg m-2 */
    realtype          harv_forage_yield;            /* kg m-2 */
    realtype          harv_res_biomass;            /* kg m-2 */
    realtype          harv_transp;            /* kg m-2 */
    realtype          harv_transp_pot;            /* kg m-2 */
    realtype          harv_soil_evap;            /* kg m-2 */
} epvar_struct;

typedef struct crop_wflux_struct
{
    realtype          transp;
    realtype          transp_pot;
} crop_wflux_struct;

typedef struct crop_cstate_struct
{
    realtype          shoot;                        /* kg m-2 */
    realtype          root;                         /* kg m-2 */
    realtype          rhizho;                        /* kg m-2 */
    realtype          shoot_post_flower;    /* kg m-2 */
} crop_cstate_struct;

typedef struct crop_cflux_struct
{
    realtype          shoot_growth;             /* kg m-2 day-1 */
    realtype          root_growth;              /* kg m-2 day-1 */
    realtype          rhizo_depo;         /* kg m-2 day-1 */
    realtype          shoot_growth_unstr;   /* kg m-2 day-1 */
    realtype          root_growth_unstr;    /* kg m-2 day-1 */
} crop_cflux_struct;

typedef struct crop_nstate_struct
{
    realtype          shoot_n;                      /* kg m-2 */
    realtype          root_n;                       /* kg m-2 */
    realtype          rhizo_n;                      /* kg m-2 */
} crop_nstate_struct;

typedef struct crop_nflux_struct
{
    realtype          rhizo_n_depo;       /* kg m-2 day-1 */
    realtype          n_auto;                  /* kg m-2 day-1 */
    realtype          n_fix;                   /* kg m-2 day-1 */
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
    realtype          SOC_Mass[MAXLYR];         /* Soil organic carbon (kg m-2) */
    realtype          MBC_Mass[MAXLYR];         /* Microbial biomass C (kg m-2) */
    realtype          stanResidueMass;            /* kg m-2 */
    realtype          flatResidueMass;            /* kg m-2 */
    realtype          manureSurfaceC;             /* kg m-2 */
    realtype          residueAbgd[MAXLYR];        /* kg m-2 */
    realtype          residueRt[MAXLYR];          /* kg m-2 */
    realtype          residueRz[MAXLYR];          /* kg m-2 */
    realtype          manureC[MAXLYR];            /* kg m-2 */
} cstate_struct;

typedef struct cflux_struct
{
    realtype          C_Humified;               /* Carbon humified from residues,
                                               * roots, rizho, and manure (kg m-2 day-1)*/
    realtype          C_ResidueRespired;        /* Carbon respired from residues,
                                               * roots, rizho, and manure (kg m-2 day-1)*/
    realtype          C_SoilRespired;           /* Carbon respired from soil
                                               * organic carbon only (kg m-2 day-1)*/
    realtype          carbonRespired[MAXLYR];     /* kg m-2 day-1 */
} cflux_struct;

typedef struct nstate_struct
{
    realtype          no3[MAXLYR];              /* nitrate (kg m-2) */
    realtype          nh4[MAXLYR];              /* ammonium (kg m-2) */
    realtype          SON_Mass[MAXLYR];         /* Soil organic N (kg m-2) */
    realtype          MBN_Mass[MAXLYR];         /* Microbial biomass N (kg m-2) */
    realtype          stanResidueN;               /* kg m-2 */
    realtype          flatResidueN;               /* kg m-2 */
    realtype          manureSurfaceN;             /* kg m-2 */
    realtype          residueAbgdN[MAXLYR];       /* kg m-2 */
    realtype          residueRtN[MAXLYR];         /* kg m-2 */
    realtype          residueRzN[MAXLYR];         /* kg m-2 */
    realtype          manureN[MAXLYR];            /* kg m-2 */
} nstate_struct;

typedef struct nflux_struct
{
    realtype          no3leached;              /* NO3 leaching (kg N m-2 day-1) */
    realtype          nh4leached;              /* NH4 leaching (kg N m-2 day-1) */
    realtype          N_Immobilization;           /* kg m-2 day-1 */
    realtype          N_Mineralization;           /* kg m-2 day-1 */
    realtype          N_NetMineralization;        /* kg m-2 day-1 */
    realtype          nh4nitrif;          /* kg m-2 day-1 */
    realtype          N2O_Nitrification;          /* kg m-2 day-1 */
    realtype          no3denitrif;        /* kg m-2 day-1 */
    realtype          n2odenitrif;        /* kg m-2 day-1 */
    realtype          nh4volat[MAXLYR];         /* kg m-2 day-1 */
    realtype          uptake_no3[MAXLYR];         /* kg m-2 day-1 */
    realtype          uptake_nh4[MAXLYR];         /* kg m-2 day-1 */
    realtype          surplusn;                   /* kg m-2 day-1 */
    realtype          fert_no3[MAXLYR];           /* kg m-2 day-1 */
    realtype          fert_nh4[MAXLYR];
    realtype          immob_no3[MAXLYR];          /* kg m-2 day-1 */
    realtype          immob_nh4[MAXLYR];          /* kg m-2 day-1 */
    realtype          nitrif_nh4_to_no3[MAXLYR];  /* kg m-2 day-1 */
    realtype          nitrif_nh4_to_n2o[MAXLYR];  /* kg m-2 day-1 */
    realtype          denitn[MAXLYR];             /* kg m-2 day-1 */
    realtype          till_no3[MAXLYR];           /* kg m-2 day-1 */
    realtype          till_nh4[MAXLYR];
    realtype          urine;              /* kg m-2 day-1 */
} nflux_struct;

typedef struct nprof_struct
{
    realtype          no3;
    realtype          nh4;
} nprof_struct;

/* Solute transport structure */
typedef struct solute_struct
{
    realtype          conc;                 /* subsurface pool concentration
                                           * (kg kgH2O-1) */
    realtype          flux[NUM_EDGE];       /* subsurface solute flux (kg s-1) */
    realtype          snksrc;               /* subsurface sink/source term
                                           * (kg m-2 s-1) */
} solute_struct;
#endif

/* Boundary conditions */
typedef union bc_struct
{
    realtype          head[NUM_EDGE];    /* value of Dirichlet-type boundary
                                        * condition (m) */
    realtype          flux[NUM_EDGE];    /* value of Neumann-type boundary
                                        * condition (m3 s-1) */
} bc_struct;

/* Land surface and hydrologic initial conditions */
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
    realtype          leafc;                  /* leaf C (kgC m-2) */
    realtype          leafc_storage;          /* leaf C storage (kgC m-2) */
    realtype          leafc_transfer;         /* leaf C transfer (kgC m-2) */
    realtype          frootc;                 /* fine root C (kgC m-2) */
    realtype          frootc_storage;         /* fine root C storage (kgC m-2) */
    realtype          frootc_transfer;        /* fine root C transfer (kgC m-2) */
    realtype          livestemc;              /* live stem C (kgC m-2) */
    realtype          livestemc_storage;      /* live stem C storage (kgC m-2) */
    realtype          livestemc_transfer;     /* live stem C transfer (kgC m-2) */
    realtype          deadstemc;              /* dead stem C (kgC m-2) */
    realtype          deadstemc_storage;      /* dead stem C storage (kgC m-2) */
    realtype          deadstemc_transfer;     /* dead stem C transfer (kgC m-2) */
    realtype          livecrootc;             /* live coarse root C (kgC m-2) */
    realtype          livecrootc_storage;     /* live coarse root C storage
                                             * (kgC m-2) */
    realtype          livecrootc_transfer;    /* live coarse root C transfer
                                             * (kgC m-2) */
    realtype          deadcrootc;             /* dead coarse root C (kgC m-2) */
    realtype          deadcrootc_storage;     /* dead coarse root C storage
                                             * (kgC m-2) */
    realtype          deadcrootc_transfer;    /* dead coarse root C transfer
                                             * (kgC m-2) */
    realtype          gresp_storage;          /* growth respiration storage
                                             * (kgC m-2) */
    realtype          gresp_transfer;         /* growth respiration transfer
                                             * (kgC m-2) */
    realtype          cwdc;                   /* coarse woody debris C (kgC m-2)
                                             */
    realtype          litr1c;                 /* litter labile C (kgC m-2) */
    realtype          litr2c;                 /* litter unshielded cellulose C
                                             * (kgC m-2) */
    realtype          litr3c;                 /* litter shielded cellulose C
                                             * (kgC m-2) */
    realtype          litr4c;                 /* litter lignin C (kgC m-2) */
    realtype          soil1c;                 /* microbial recycling pool C (fast)
                                             * (kgC m-2) */
    realtype          soil2c;                 /* microbial recycling pool C
                                             * (medium)
                                             * (kgC m-2) */
    realtype          soil3c;                 /* microbial recycling pool C (slow)
                                             * (kgC m-2) */
    realtype          soil4c;                 /* recalcitrant SOM C (humus,
                                             * slowest) (kgC m-2) */
    realtype          cpool;                  /* temporary photosynthate C pool
                                             * (kgC m-2) */
    realtype          psnsun_src;             /* SUM of gross PSN from sunlit
                                             * canopy (kgC m-2) */
    realtype          psnshade_src;           /* SUM of gross PSN from shaded
                                             * canopy (kgC m-2) */
    realtype          leaf_mr_snk;            /* SUM of leaf maint resp (kgC m-2)
                                             */
    realtype          leaf_gr_snk;            /* SUM of leaf growth resp (kgC m-2)
                                             */
    realtype          froot_mr_snk;           /* SUM of fine root maint resp
                                             * (kgC m-2) */
    realtype          froot_gr_snk;           /* SUM of fine root growth resp
                                             * (kgC m-2) */
    realtype          livestem_mr_snk;        /* SUM of live stem maint resp
                                             * (kgC m-2) */
    realtype          livestem_gr_snk;        /* SUM of live stem growth resp
                                             * (kgC m-2) */
    realtype          deadstem_gr_snk;        /* SUM of dead stem growth resp
                                             * (kgC m-2) */
    realtype          livecroot_mr_snk;       /* SUM of live coarse root maint
                                             * resp (kgC m-2) */
    realtype          livecroot_gr_snk;       /* SUM of live coarse root growth
                                             * resp (kgC m-2) */
    realtype          deadcroot_gr_snk;       /* SUM of dead coarse root growth
                                             * resp (kgC m-2) */
    realtype          litr1_hr_snk;           /* SUM of labile litr microbial resp
                                             * (kgC m-2) */
    realtype          litr2_hr_snk;           /* SUM of cellulose litr microbial
                                             * resp (kgC m-2) */
    realtype          litr4_hr_snk;           /* SUM of lignin litr microbial resp
                                             * (kgC m-2) */
    realtype          soil1_hr_snk;           /* SUM of fast microbial respiration
                                             * (kgC m-2) */
    realtype          soil2_hr_snk;           /* SUM of medium microbial
                                             * respiration (kgC m-2) */
    realtype          soil3_hr_snk;           /* SUM of slow microbial respiration
                                             * (kgC m-2) */
    realtype          soil4_hr_snk;           /* SUM of recalcitrant SOM
                                             * respiration (kgC m-2) */
    realtype          fire_snk;               /* SUM of fire losses (kgC m-2) */
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
    realtype          cwdn;            /* coarse woody debris N (kgN m-2) */
    realtype          litr1n;          /* litter labile N (kgN m-2) */
    realtype          litr2n;          /* litter unshielded cellulose N (kgN m-2)
                                      */
    realtype          litr3n;          /* litter shielded cellulose N (kgN m-2) */
    realtype          litr4n;          /* litter lignin N (kgN m-2) */
    realtype          soil1n;          /* microbial recycling pool N (fast)
                                      * (kgN m-2) */
    realtype          soil2n;          /* microbial recycling pool N (medium)
                                      * (kgN m-2) */
    realtype          soil3n;          /* microbial recycling pool N (slow)
                                      * (kgN m-2) */
    realtype          soil4n;          /* recalcitrant SOM N (humus, slowest)
                                      * (kgN m-2) */
    realtype          surfn;           /* surface mineral N (kgN m-2) */
    realtype          sminn;           /* soil mineral N (kgN m-2) */
    realtype          retransn;        /* plant pool of retranslocated N (kgN m-2)
                                      */
    realtype          npool;           /* temporary plant N pool (kgN m-2) */
    realtype          nfix_src;        /* SUM of biological N fixation (kgN m-2)
                                      */
    realtype          ndep_src;        /* SUM of N deposition inputs (kgN m-2) */
    realtype          nleached_snk;    /* SUM of N leached (kgN m-2) */
    realtype          nvol_snk;        /* SUM of N lost to volatilization
                                      * (kgN m-2) */
    realtype          fire_snk;        /* SUM of N lost to fire (kgN m-2) */
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
    realtype          surfn0;             /* surface N of previous time step
                                         * (kg N m-2) */
    realtype          sminn0;             /* subsurface N of previous time step
                                         * (kg N m-2) */
    realtype          mineralized;        /* N mineralization (kgN m-2 day-1) */
    realtype          potential_immob;    /* potential N immobilization (kgN m-2)
                                         */
    realtype          plitr1c_loss;       /* potential loss from litter labile
                                         * pool (kgN m-2 day-1) */
    realtype          pmnf_l1s1;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    realtype          plitr2c_loss;       /* potential loss from litter unshielded
                                         * pool (kgN m-2 day-1) */
    realtype          pmnf_l2s2;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    realtype          plitr4c_loss;       /* potential loss from litter lignin
                                         * pool (kgN m-2 day-1) */
    realtype          pmnf_l4s3;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    realtype          psoil1c_loss;       /* potential loss from fast soil pool
                                         * (kgN m-2 day-1) */
    realtype          pmnf_s1s2;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    realtype          psoil2c_loss;       /* potential loss from medium soil pool
                                         * (kgN m-2 day-1) */
    realtype          pmnf_s2s3;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    realtype          psoil3c_loss;       /* potential loss from slow soil pool
                                         * (kgN m-2 day-1) */
    realtype          pmnf_s3s4;          /* potential mineral N flux
                                         * (kgN m-2 day-1) */
    realtype          psoil4c_loss;       /* potential loss from slowest soil pool
                                         * (kgN m-2 day-1) */
    realtype          kl4;                /* decomposition rate of lignin litter
                                         * pool (day -1) */
} ntemp_struct;

/* Ecophysiological variables */
typedef struct epvar_struct
{
    realtype          bg_leafc_litfall_rate;  /* rate leaf litfall (kgC m-2 s-1)
                                             */
    realtype          bg_frootc_litfall_rate; /* rate froot litfall (kgC m-2 s-1)
                                             */
    realtype          livestemc_turnover_rate;/* rate livestem turnover
                                             * (kgC m-2 s-1) */
    realtype          livecrootc_turnover_rate; /* rate livecroot turnover
                                             * (kgC m-2 s-1) */
    realtype          dsr;                    /* number of days since rain */
    realtype          sun_proj_sla;           /* sunlit projected SLA (m2 kgC-1)
                                             */
    realtype          shade_proj_sla;         /* shaded projected SLA (m2 kgC-1)
                                             */
    realtype          psi;                    /* water potential of soil and
                                             * leaves (MPa) */
    realtype          dlmr_area_sun;          /* sunlit leaf MR (umolC/m2
                                             * projected leaf area s-1) */
    realtype          dlmr_area_shade;        /* shaded leaf MR (umolC/m2
                                             * projected leaf area s-1) */
    realtype          gl_t_wv_sun;            /* leaf-scale conductance to
                                             * transpired water (m s-1) */
    realtype          gl_t_wv_shade;          /* leaf-scale conductance to
                                             * transpired water (m s-1) */
    realtype          assim_sun;              /* sunlit assimilation per unit pLAI
                                             * (umol m-2 s-1) */
    realtype          assim_shade;            /* shaded assimilation per unit pLAI
                                             * (umol m-2 s-1) */
    realtype          t_scalar;               /* decomp temperature scalar (-) */
    realtype          w_scalar;               /* decomp water scalar (-) */
    realtype          rate_scalar;            /* decomp combined scalar (-) */
    realtype          daily_gross_nimmob;     /* daily gross N immobilization
                                             * (kgN m-2 d-1) */
    realtype          daily_net_nmin;         /* daily net N mineralization
                                             * (kgN m-2 d-1) */
    realtype          fpi;                    /* fraction of potential
                                             * immobilization (-) */
    realtype          m_tmin;                 /* freezing night temperature
                                             * multiplier (-) */
    realtype          m_psi;                  /* water potential multiplier (-) */
    realtype          m_co2;                  /* atmospheric (CO2) multiplier (-)
                                             */
    realtype          m_ppfd_sun;             /* PAR flux density multiplier (-)
                                             */
    realtype          m_ppfd_shade;           /* PAR flux density multiplier (-)
                                             */
    realtype          m_vpd;                  /* vapor pressure deficit multiplier
                                             * (-) */
    realtype          m_final_sun;            /* product of all other multipliers
                                             * (-) */
    realtype          m_final_shade;          /* product of all other multipliers
                                             * (-) */
    realtype          ytd_maxplai;            /* year-to-date maximum projected
                                             * LAI (-) */
    int             dormant_flag;           /* dormancy flag */
    realtype          days_active;            /* number of days since last
                                             * dormancy */
    int             onset_flag;             /* onset flag */
    int             onset_counter;          /* onset days counter */
    int             onset_gddflag;          /* onset flag for growing degree
                                             * day sum */
    realtype          onset_fdd;              /* onset freezing degree days
                                             * counter */
    realtype          onset_gdd;              /* onset growing degree days */
    realtype          onset_swi;              /* onset soil water index */
    int             offset_flag;            /* offset flag */
    int             offset_counter;         /* offset days counter */
    realtype          offset_fdd;             /* offset freezing degree days
                                             * counter */
    realtype          offset_swi;             /* offset soil water index */
    realtype          annavg_t2m;             /* annual average 2m air
                                             * temperature (K) */
    realtype          gpp;                    /* GPP flux before downregulation
                                             * (gC m-2 s-1) */
    realtype          prev_leafc_to_litter;   /* previous timestep leaf C
                                             * litterfall flux (gC m-2 s-1) */
    realtype          prev_frootc_to_litter;  /* previous timestep froot C
                                             * litterfall flux (gC m-2 s-1) */
    realtype          old_c_balance;          /* previous timestep C balance
                                             * (kgC m-2 day-1) */
    realtype          old_n_balance;          /* previous timestep N balance
                                             * (kgN m-2 day-1) */
    realtype          dayl;                   /* daylength (s) */
    realtype          prev_dayl;              /* previous day daylength (s) */
} epvar_struct;

/* Structure for the photosynthesis routine */
typedef struct psn_struct
{
    int             c3;       /* set to 1 for C3 model, 0 for C4 model */
    realtype          pa;       /* atmospheric pressure (Pa) */
    realtype          co2;      /* atmospheric CO2 concentration (ppm) */
    realtype          t;        /* temperature (deg C) */
    realtype          lnc;      /* leaf N per unit sunlit leaf area (kg Nleaf m-2)
                               */
    realtype          flnr;     /* fract. of leaf N in Rubisco (kg NRub/kg Nleaf)
                               */
    realtype          ppfd;     /* PAR flux per unit sunlit leaf area
                               * (umol m-2 s-1) */
    realtype          g;        /* conductance to CO2 (umol m-2 s-1 Pa-1) */
    realtype          dlmr;     /* day leaf maintenance respiration, projected
                               * area basis (umol m-2 s-1) */
    realtype          Ci;       /* intercellular CO2 concentration (Pa) */
    realtype          O2;       /* atmospheric O2 concentration (Pa) */
    realtype          Ca;       /* atmospheric CO2 concentration (Pa) */
    realtype          gamma;    /* CO2 compensation point, no Rd (Pa) */
    realtype          Kc;       /* MM constant carboxylation (Pa) */
    realtype          Ko;       /* MM constant oxygenation (Pa) */
    realtype          Vmax;     /* max rate carboxylation (umol m-2 s-1) */
    realtype          Jmax;     /* max rate electron transport (umol m-2 s-1) */
    realtype          J;        /* rate of RuBP regeneration (umol m-2 s-1) */
    realtype          Av;       /* carboxylation limited assimilation
                               * (umol m-2 s-1) */
    realtype          Aj;       /* RuBP regen limited assimilation (umol m-2 s-1)
                               */
    realtype          A;        /* final assimilation rate (umol m-2 s-1) */
} psn_struct;

/* CN summary structure */
typedef struct summary_struct
{
    realtype          daily_npp;         /* = GPP - Rmaint - Rgrowth
                                        * (kgC m-2 day-1) */
    realtype          daily_nep;         /* = NPP - Rheterotroph (kgC m-2 day-1)
                                        */
    realtype          daily_nee;         /* = NEP - fire losses (kgC m-2 day-1) */
    realtype          daily_gpp;         /* gross PSN source (kgC m-2 day-1) */
    realtype          daily_mr;          /* maintenance respiration
                                        * kgC m-2 day-1) */
    realtype          daily_gr;          /* growth respiration (kgC m-2 day-1) */
    realtype          daily_hr;          /* heterotrophic respiration
                                        * kgC m-2 day-1) */
    realtype          daily_fire;        /* fire losses (kgC m-2 day-1) */
    realtype          daily_litfallc;    /* total litterfall (kgC m-2 day-1) */
    /* summed over entire simulation (kgC m-2) */
    realtype          cum_npp;
    realtype          cum_nep;
    realtype          cum_nee;
    realtype          cum_gpp;
    realtype          cum_mr;
    realtype          cum_gr;
    realtype          cum_hr;
    realtype          cum_fire;
    realtype          vegc;              /* total vegetation C (kgC m-2) */
    realtype          agc;               /* aboveground C (kgC m-2) */
    realtype          litrc;             /* total litter C (kgC m-2) */
    realtype          soilc;             /* total soil C (kgC m-2) */
    realtype          totalc;            /* total of vegc, litrc, and soilc
                                        * (kgC m-2) */
} summary_struct;

/* Solute transport structure */
typedef struct solute_struct
{
# if !defined(_LEACHING_)
    realtype          conc_surf;            /* surface pool concentration
                                           * (kg kgH2O-1) */
    realtype          infilflux;            /* solute infiltration flux
                                           * (kg m-2 s-1) */
    realtype          ovlflux[NUM_EDGE];    /* overland solute flux (kg s-1) */
# endif
    realtype          conc_subsurf;         /* subsurface pool concentration
                                           * (kg kgH2O-1) */
    realtype          subflux[NUM_EDGE];    /* subsurface solute flux (kg s-1) */
    realtype          snksrc;               /* subsurface sink/source term
                                           * (kg m-2 s-1) */
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

/* Element structure */
typedef struct elem_struct
{
    int             node[NUM_EDGE];    /* nodes of triangular element
                                        * (counterclockwise) */
    int             nabr[NUM_EDGE];    /* neighbor elements (neighbor i shares
                                        * edge i (0: on boundary) */
    int             ind;               /* index */
	
    attrib_struct   attrib;
    topo_struct     topo;
    soil_struct     soil;
    lc_struct       lc;
#if defined(_FBR_)
    geol_struct     geol;
#endif
#if defined(_CYCLES_)
    crop_struct     crop[MAXCROP];
#else
    epconst_struct  epc;
#endif
    ic_struct       ic;
    bc_struct       bc;
#if defined(_FBR_)
    bc_struct       fbr_bc;
#endif
    wstate_struct   ws;
    wstate_struct   ws0;
    wflux_struct    wf;
    estate_struct   es;
    eflux_struct    ef;
    pstate_struct   ps;
#if defined(_DAILY_)
    daily_struct    daily;
#endif
#if defined(_CYCLES_)
    mgmt_struct     mgmt;
    cyclesic_struct restart_input;
    cstate_struct   cs;
    cflux_struct    cf;
    nstate_struct   ns;
    nstate_struct   ns0;
    nflux_struct    nf;
    nprof_struct    np;
    solute_struct   no3sol;
    solute_struct   nh4sol;
#endif
#if defined(_BGC_)
    bgcic_struct    restart_input;
    bgcic_struct    restart_output;
    cstate_struct   cs;
    cflux_struct    cf;
    nstate_struct   ns;
    nflux_struct    nf;
    psn_struct      psn_sun;
    psn_struct      psn_shade;
    ntemp_struct    nt;
    summary_struct  summary;
    epvar_struct    epv;
    solute_struct   nsol;
    spinup_struct   spinup;
#endif
} elem_struct;
#endif

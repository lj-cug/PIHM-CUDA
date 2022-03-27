#ifndef PIHM_INPUT_STRUCT_HEADER
#define PIHM_INPUT_STRUCT_HEADER

/* Input file names */
typedef struct filename_struct
{
    char            riv[MAXSTRING];         /* river input file name */
    char            mesh[MAXSTRING];        /* mesh structure file name */
    char            att[MAXSTRING];         /* element attribute file name */
    char            soil[MAXSTRING];        /* soil property file name */
    char            lc[MAXSTRING];          /* land cover property file name */
    char            meteo[MAXSTRING];       /* meteorological forcing file name
                                             */
    char            lai[MAXSTRING];         /* lai forcing file name */
    char            bc[MAXSTRING];          /* boundary condition file name */
    char            para[MAXSTRING];        /* control parameter file name */
    char            calib[MAXSTRING];       /* calibration file name */
    char            ic[MAXSTRING];          /* initial condition file name */
    char            tecplot[MAXSTRING];     /* tecplot control file name */
#if defined(_FBR_)
    char            geol[MAXSTRING];        /* geology property file name */
    char            bedrock[MAXSTRING];     /* bedrock elevation file name */
#endif
#if defined(_NOAH_)
    char            lsm[MAXSTRING];         /* land surface module control file
                                             * name */
    char            rad[MAXSTRING];         /* radiation forcing file name */
#endif
#if defined(_CYCLES_)
    char            cycles[MAXSTRING];
    char            soilinit[MAXSTRING];
    char            crop[MAXSTRING];
    char            op[MAXOP][MAXSTRING];
    char            cyclesic[MAXSTRING];
#endif
#if defined(_BGC_)
    char            bgc[MAXSTRING];         /* bgc module control file name */
    char            co2[MAXSTRING];         /* CO2 forcing file name */
    char            ndep[MAXSTRING];        /* nitrogen deposition forcing file
                                             * name */
    char            bgcic[MAXSTRING];       /* bgc module initial condition file
                                             * name */
#endif
} filename_struct;

/* River input structure */
typedef struct rivtbl_struct
{
    int            *fromnode;    /* upstream node id */
    int            *tonode;      /* downstream node id */
    int            *down;        /* downstream channel id */
    int            *leftele;     /* left bank id */
    int            *rightele;    /* right bank id */
    int            *shp;         /* river shape type */
    int            *matl;        /* material type */
    int            *bc;          /* boundary condition type */
    int            *rsvr;        /* reservoir type */
} rivtbl_struct;

/* River shape parameters */
typedef struct shptbl_struct
{
    int             number;       /* number of shape types */
    realtype         *depth;        /* river channel depth */
    int            *intrpl_ord;   /* interpolation order (shape of channel)
                                   * 1: rectangle
                                   * 2: triangle
                                   * 3: quadratic
                                   * 4: cubic */
    realtype         *coeff;        /* width coefficient */
} shptbl_struct;

/* River channel material parameters */
typedef struct matltbl_struct
{
    int             number;      /* number of bank/bed material types */
    realtype         *rough;       /* river channel roughness (s m-1/3) */
    realtype         *cwr;         /* discharge coefficient (-) */
    realtype         *ksath;       /* bank hydraulic conductivity (m s-1) */
    realtype         *ksatv;       /* bed hydraulic conductivity (m s-1) */
    realtype         *bedthick;    /* bed thickness (m) */
} matltbl_struct;

/* Mesh structure */
typedef struct meshtbl_struct
{
    int             numnode;    /* number of nodes */
    int           **node;       /* nodes of element */
    int           **nabr;       /* neighbors of element */
    realtype         *x;          /* x of node (m) */
    realtype         *y;          /* y of node (m) */
    realtype         *zmin;       /* soil bottom elevation of node (m) */
    realtype         *zmax;       /* surface elevation of node (m) */
#if defined(_FBR_)
    realtype         *zbed;       /* impermeable bedrock elevation (m) */
#endif
} meshtbl_struct;

/* Element attribute */
typedef struct atttbl_struct
{
    int            *soil;      /* element soil type */
    int            *geol;      /* element geology type */
    int            *lc;        /* element land cover type */
    int           **bc;        /* element boundary condition type */
#if defined(_FBR_)
    int           **fbr_bc;    /* element boundary condition type for fractured
                                * bedrock layer */
#endif
    int            *meteo;     /* element meteorological forcing type */
    int            *lai;       /* element leaf area index forcing type
                                * 0: use climatological values;
                                * else: use forcing file */
    int            *source;    /* element source forcing type */
} atttbl_struct;

/* Soil parameter */
typedef struct soiltbl_struct
{
    int             number;      /* number of soil types */
    realtype         *silt;        /* silt percentage (%) */
    realtype         *clay;        /* clay percentage (%) */
    realtype         *om;          /* organic matter percentage (%) */
    realtype         *bd;          /* bulk density (g cm-3) */
    realtype         *kinfv;       /* saturated infiltration conductivity (m s-1)
                                  */
    realtype         *ksatv;       /* vertical saturated hydraulic conductivity
                                  * (m s-1) */
    realtype         *ksath;       /* horizontal saturated hydraulic conductivity
                                  * (m s-1) */
    realtype         *smcmax;      /* maximum soil moisture content (m3 m-3) */
    realtype         *smcmin;      /* residual soil moisture content (m3 m-3) */
    realtype         *smcwlt;      /* wilting point (m3 m-3) */
    realtype         *smcref;      /* soil moisture threshold where transpiration
                                  * begins to stress (m3 m-3) */
    realtype         *qtz;         /* soil quartz content (-) */
    realtype         *alpha;       /* alpha from van Genuchten eqn (m-1) */
    realtype         *beta;        /* beta (n) from van Genuchten eqn (-) */
    realtype         *areafh;      /* macropore area fraction on a horizontal
                                  * cross-section (m2 m-2) */
    realtype         *areafv;      /* macropore area fraction on a vertical
                                  * cross-section (m2 m-2) */
    realtype         *dmac;        /* macropore depth (m) */
    realtype          dinf;        /* depth from ground surface across which head
                                  * gradient is calculated for infiltration (m)
                                  */
    realtype          kmacv_ro;    /* ratio between vertical macropore hydraulic
                                  * conductivity and vertical saturated
                                  * infiltration hydraulic conductivity */
    realtype          kmach_ro;    /* ratio between horizontal macropore hydraulic
                                  * conductivity and horizontal saturated
                                  * hydraulic conductivity (-) */
#if defined(_CYCLES_)
    int            *totalLayers;
    realtype        **clay_lyr;
    realtype        **sand_lyr;
    realtype        **iom_lyr;
    realtype        **bd_lyr;
    realtype        **no3_lyr;
    realtype        **nh4_lyr;
#endif
} soiltbl_struct;

/* Geology parameter */
typedef struct geoltbl_struct
{
    int             number;    /* number of soil types */
    realtype         *silt;      /* silt percentage (%) */
    realtype         *clay;      /* clay percentage (%) */
    realtype         *om;        /* organic matter percentage (%) */
    realtype         *bd;        /* bulk density (g cm-3) */
    realtype         *ksath;     /* horizontal saturated hydraulic conductivity
                                * (m s-1) */
    realtype         *ksatv;     /* vertical saturated hydraulic conductivity
                                * (m s-1) */
    realtype         *smcmax;    /* maximum soil moisture content (m3 m-3) */
    realtype         *smcmin;    /* residual soil moisture content (m3 m-3) */
    realtype         *alpha;     /* alpha from van Genuchten eqn (m-1) */
    realtype         *beta;      /* beta (n) from van Genuchten eqn (-) */
} geoltbl_struct;

/* Land cover parameters */
typedef struct lctbl_struct
{
    int             number;       /* number of land cover types */
    realtype         *laimax;       /* maximum LAI across all seasons for a
                                   * vegetation type (m2 m-2) */
    realtype         *laimin;       /* minimum LAI across all seasons for a
                                   * vegetation type (m2 m-2) */
    realtype         *vegfrac;      /* areal fractional coverage of green
                                   * vegetation (0.0-1.0) (-) */
    realtype         *albedomin;    /* minimum background albedo (-) */
    realtype         *albedomax;    /* maximum background albedo (-) */
    realtype         *emissmin;     /* minimum emissivity (-) */
    realtype         *emissmax;     /* maximum emissivity (-) */
    realtype         *z0min;        /* minimum roughness length (m) */
    realtype         *z0max;        /* maximum roughness length (m) */
    realtype         *hs;           /* parameter used in vapor pressure deficit
                                   * function (-) */
    realtype         *snup;         /* threshold snow depth (in water equivalent)
                                   * that implies 100% snow cover (m) */
    realtype         *rgl;          /* reference incoming solar flux for
                                   * photosynthetically active canopy (W m-2) */
    realtype         *rsmin;        /* minimum canopy resistance (s m-1) */
    realtype         *rough;        /* surface roughness (Manning's n) (s m-1/3)
                                   */
    realtype         *rzd;          /* rooting depth (m) */
    realtype          rsmax;        /* cuticular resistance (s m-1) */
    int             bare;         /* the land-use category representing bare
                                   * ground */
    int             natural;      /* the land-use category representing non-
                                   * urban portion of urban land-use points */
    realtype          cfactr;       /* parameter used in the canopy interception
                                   * calculation (-) */
    realtype          topt;         /* optimum transpiration air temperature (K)
                                   */
} lctbl_struct;

/* Time series data structure */
typedef struct tsdata_struct
{
    int             length;       /* length of time series */
    int            *ftime;        /* forcing time */
    realtype        **data;         /* forcing values at forcing time */
    realtype         *value;        /* forcing values at model time t */
    realtype          zlvl_wind;    /* height above ground of wind observations
                                   * (m) */
} tsdata_struct;

/* Forcing structure */
typedef struct forc_struct
{
    int             nbc;         /* number of boundary condition series */
    tsdata_struct  *bc;          /* boundary condition time series */
    int             nmeteo;      /* number of meteorological forcing series */
    tsdata_struct  *meteo;       /* meteorological forcing series */
    int             nlai;        /* number of lai series */
    tsdata_struct  *lai;         /* lai forcing series */
    int             nsource;     /* number of source forcing series */
    tsdata_struct  *source;      /* source forcing series */
    int             nriverbc;    /* number of river boundary conditions */
    tsdata_struct  *riverbc;     /* river boundary condition series */
#if defined(_NOAH_)
    int             nrad;        /* number of radiation forcing series */
    tsdata_struct  *rad;         /* radiation forcing series */
#endif
#if defined(_BGC_)
    int             nco2;
    tsdata_struct  *co2;         /* CO2 forcing series */
    int             nndep;
    tsdata_struct  *ndep;        /* nitrogen deposition forcing series */
#endif
} forc_struct;

#if defined(_NOAH_)
/* Land surface parameters */
typedef struct noahtbl_struct
{
    realtype          sbeta;     /* parameter used to calculate vegetation effect
                                * on soil heat (-) */
    realtype          fxexp;     /* soil evaporation exponent used in direct
                                * evaporation (-) */
    realtype          csoil;     /* soil heat capacity (J m-3 K-1) */
    realtype          salp;      /* shape parameter of distribution function of
                                * snow cover (-) */
    realtype          frzk;      /* frozen ground parameter (-) */
    realtype          zbot;      /* depth of lower boundary soil temperature (m)
                                */
    realtype          tbot;      /* bottom soil temperature (local yearly-mean
                                * surface air temperature) (K) */
    realtype          czil;      /* Zilitinkevich constant (-) */
    realtype          lvcoef;    /* parameter controls surface snow albedo in the
                                * presence of snowcover (-) */
} noahtbl_struct;
#endif

#if defined(_BGC_)
/* Ecophysiological parameters */
typedef struct epctbl_struct
{
    int            *woody;                  /* flag: 1 = woody, 0 = non-woody */
    int            *evergreen;              /* flag: 1 = evergreen,
                                             * 0 = deciduous */
    int            *c3_flag;                /* flag: 1 = C3,  0 = C4 */
    int            *phenology_flag;         /* flag: 1 = phenology model,
                                             * 0 = user defined */
    int            *onday;                  /* day of year when leaves on */
    int            *offday;                 /* day of year when leaves off */
    int            *transfer_days;          /* growth period for transfer (day)
                                             */
    int            *litfall_days;           /* growth period for litfall (day)
                                             */
    realtype         *leaf_turnover;          /* annual leaf turnover fraction
                                             * (yr-1) */
    realtype         *froot_turnover;         /* annual fine root turnover
                                             * fraction (yr-1) */
    realtype         *livewood_turnover;      /* annual live wood turnover
                                             * fraction (yr-1) */
    realtype         *daily_mortality_turnover; /* daily mortality turnover
                                             * (day-1) */
    realtype         *daily_fire_turnover;    /* daily fire turnover (day-1) */
    realtype         *alloc_frootc_leafc;     /* new fine root C to new leaf C (-)
                                             */
    realtype         *alloc_newstemc_newleafc; /* new stem C to new leaf C (-) */
    realtype         *alloc_newlivewoodc_newwoodc; /* new livewood C:new wood C
                                             * (-) */
    realtype         *alloc_crootc_stemc;     /* new live croot C to new live stem
                                             * C (-) */
    realtype         *alloc_prop_curgrowth;   /* daily allocation to current
                                             * growth (-) */
    realtype         *avg_proj_sla;           /* canopy average projected SLA
                                             * (m2 kgC-1) */
    realtype         *sla_ratio;              /* ratio of shaded to sunlit
                                             * projected SLA (-) */
    realtype         *lai_ratio;              /* ratio of (all-sided LA /
                                             * one-sided LA) (-) */
    realtype         *ext_coef;               /* canopy light extinction
                                             * coefficient (-) */
    realtype         *flnr;                   /* leaf N in Rubisco
                                             * (kgNRub kgNleaf-1) */
    realtype         *psi_open;               /* psi at start of conductance
                                             * reduction (MPa) */
    realtype         *psi_close;              /* psi at complete conductance
                                             * reduction (MPa) */
    realtype         *vpd_open;               /* vpd at start of conductance
                                             * reduction (Pa) */
    realtype         *vpd_close;              /* vpd at complete conductance
                                             * reduction (Pa) */
    realtype         *froot_cn;               /* C:N for fine roots (kgC kgN-1) */
    realtype         *leaf_cn;                /* C:N for leaves (kgC kgN-1) */
    realtype         *livewood_cn;            /* C:N for live wood (kgC kgN-1) */
    realtype         *deadwood_cn;            /* C:N for dead wood (kgC kgN-1) */
    realtype         *leaflitr_cn;            /* constant C:N for leaf litter
                                             * (kgC kgN-1) */
    realtype         *leaflitr_flab;          /* leaf litter labile fraction (-)
                                             */
    realtype         *leaflitr_fucel;         /* leaf litter unshielded cellulose
                                             * fraction (-) */
    realtype         *leaflitr_fscel;         /* leaf litter shielded cellulose
                                             * fraction (-) */
    realtype         *leaflitr_flig;          /* leaf litter lignin fraction (-)
                                             */
    realtype         *frootlitr_flab;         /* fine root litter labile fraction
                                             * (-) */
    realtype         *frootlitr_fucel;        /* fine root litter unshielded
                                             * cellulose fraction (-) */
    realtype         *frootlitr_fscel;        /* fine root litter shielded
                                             * cellulose fraction (-) */
    realtype         *frootlitr_flig;         /* fine root litter lignin fraction
                                             * (-) */
    realtype         *deadwood_fucel;         /* dead wood unshielded cellulose
                                             * fraction (-) */
    realtype         *deadwood_fscel;         /* dead wood shielded cellulose
                                             * fraction (-) */
    realtype         *deadwood_flig;          /* dead wood lignin fraction (-) */
} epctbl_struct;
#endif

#if defined(_CYCLES_)
typedef struct agtbl_struct
{
    int            *op;
    int            *rotsz;
    int            *auto_N;
    int            *auto_P;
    int            *auto_S;
    int             nopfile;
    char            opfilen[MAXOP][MAXSTRING];
} agtbl_struct;

typedef struct epconst_struct
{
    char            cropn[MAXSTRING];
    realtype          flowering_tt;                        /* (degree C day) */
    realtype          maturity_tt;                         /* (degree C day) */
    realtype          max_soil_covr;                /* (100%) */
    realtype          max_root_dpth;                /* (m) */
    realtype          frac_res_standing;            /* (100%) */
    realtype          frac_res_removed;             /* (100%) */
    realtype          clip_biomass_ut;      /* (kg m-2) */
    realtype          clip_biomass_lt;      /* (kg m-2) */
    realtype          clip_timing;                     /* (100% thermal time) */
    int             clip_destiny;                    /* (-) */
    realtype          transp_temp_min;        /* (degree C) */
    realtype          transp_temp_thr;  /* (degree C) */
    realtype          cold_damage_temp_min;           /* (degree C) */
    realtype          cold_damagetemp_thr;     /* (degree C) */
    realtype          temp_base;                    /* (degree C) */
    realtype          temp_opt;                 /* (degree C) */
    realtype          temp_max;                 /* (degree C) */
    realtype          shoot_par_init;              /* (100%) */
    realtype          shoot_par_final;                /* (100%) */
    realtype          rue;             /* kg J-1 */
    realtype          tue;         /* kg kgH2O-1 at VPD = 1kPa */
    realtype          hi_max;                                /* (-) */
    realtype          hi_min;                                /* intercept harvest index (-) */
    realtype          hi;                                /* (-) */
    realtype          emergence_tt;                        /* (degree C day) */
    realtype          n_conc_max;                  /* (g g-1) */
    realtype          n_diln_slope;                     /* (-) */
    realtype          kc;                                 /* (-) */
    int             annual;                             /* (-) */
    int             legume;                             /* (-) */
    int             c3c4;                             /* (-) */
    realtype          lwp_stress_onset;                    /* (J kg-1, or m2 s-2) */
    realtype          lwp_wlt;                   /* (J kg-1, or m2 s-2) */
    realtype          transp_max;                   /* (mm day-1) */
} epconst_struct;

typedef struct plant_struct
{
    /* Planting */
    int             year;
    int             doy;
    int             crop_id;
    int             auto_irr;
    int             auto_fert;
    realtype          plant_density;            /* (100%) */
    int             clip_start;              /* (day of year) */
    int             clip_end;                /* (day of year) */
    int             ai_start;
    int             ai_stop;
    realtype          ai_h2o_depl;         /* (100% plant avialable water content) */
    int             ai_last_lyr;
} plant_struct;

typedef struct tillage_struct
{
    /* Tillage */
    int             year;
    int             doy;
    char            tooln[MAXSTRING];
    realtype          depth;                    /* (m) */
    realtype          sdr;                      /* (-) */
    realtype          mix_eff;                  /* (100%) */
    realtype          drop_id;
#if NOT_YET_IMPLEMENTED
    realtype          fractionThermalTime;
    realtype          killEfficiency;
#endif
    int             grain_harv;
    realtype          forage_harv;
} tillage_struct;

typedef struct fixirr_struct
{
    /* Fixed Irrigation */
    int             year;
    int             doy;
    realtype          volume;   /* (mm) */
} fixirr_struct;

typedef struct fixfert_struct
{
    /* Fixed Fertilization */
    int             year;
    int             doy;
    char            source[MAXSTRING];
    realtype          mass;                 /* (kg m-2) */
#if NOT_YET_IMPLEMENTED
    char            opForm[MAXSTRING];
    char            opMethod[MAXSTRING];
#endif
    int             layer;            /* Starting from 1 */
    realtype          c_org;            /* (100%) */
    realtype          c_cac;            /* (100%) */
    realtype          n_org;            /* (100%) */
    realtype          n_cac;            /* (100%) */
    realtype          nh4;              /* (100%) */
    realtype          no3;              /* (100%) */
    realtype          p_org;            /* (100%) */
    realtype          p_cac;            /* (100%) */
    realtype          p_inorg;          /* (100%) */
    realtype          k;                /* (100%) */
    realtype          s;                /* (100%) */
} fixfert_struct;

typedef struct autoirr_struct
{
    int             crop_id;
    int             start;
    int             stop;
    realtype          h2o_depl;         /* (100% plant avialable water content) */
    int             last_lyr;
} autoirr_struct;

typedef struct opertbl_struct
{
    fixfert_struct *fix_fert;
    int             nfert;
    fixirr_struct  *fix_irr;
    int             nirr;
    tillage_struct *tillage;
    int             ntill;
    plant_struct   *plant;
    int             nplant;
    autoirr_struct *auto_irr;
    int             nauto_irr;
} opertbl_struct;
#endif

#endif

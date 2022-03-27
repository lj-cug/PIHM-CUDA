#ifndef PIHM_STRUCT_HEADER
#define PIHM_STRUCT_HEADER

/* Time structure */
typedef struct pihm_t_struct
{
    int             t;
    int             year;
    int             month;
    int             day;
    int             hour;
    int             minute;
    char            str[17];
    char            strshort[13];
} pihm_t_struct;

/* Site information structure */
typedef struct siteinfo_struct
{
    realtype          longitude;    /* (degree) */
    realtype          latitude;     /* (degree) */
    realtype          zmax;         /* average surface elevation (m) */
    realtype          zmin;         /* average soil bottom elevation (m) */
    realtype          area;         /* total area (m2) */
    realtype          tavg;         /* annual average air temperature (K) */
} siteinfo_struct;

#if defined(_BGC_)
/* A structure to hold information on the annual co2 concentration */
typedef struct co2control_struct
{
    int             varco2;    /* 0 = const 1 = use file */
    realtype          co2ppm;    /* constant CO2 concentration (ppm) */
} co2control_struct;

/* A structure to hold annual nitrogen deposition data */
typedef struct ndepcontrol_struct
{
    int             varndep;    /* 0 = const 1 = use file */
    realtype          ndep;       /* wet+dry atmospheric deposition of N
                                 * (kgN m-2 yr-1) */
    realtype          nfix;       /* symbiotic+asymbiotic fixation of N
                                 * (kgN m-2 yr-1) */
} ndepcontrol_struct;

/* Carbon and nitrogen state initialization structure */
typedef struct cninit_struct
{
    realtype          max_leafc;    /* first-year displayed + stored leafc
                                   * (kgC m-2) */
    realtype          max_stemc;    /* first-year total stem carbon (kgC m-2) */
    realtype          cwdc;         /* coarse woody debris C (kgC m-2) */
    realtype          litr1c;       /* litter labile C (kgC m-2) */
    realtype          litr2c;       /* litter unshielded cellulose C (kgC m-2) */
    realtype          litr3c;       /* litter shielded cellulose C (kgC m-2) */
    realtype          litr4c;       /* litter lignin C (kgC m-2) */
    realtype          soil1c;       /* microbial recycling pool C (fast)
                                   * (kgC m-2) */
    realtype          soil2c;       /* microbial recycling pool C (medium)
                                   * (kgC m-2) */
    realtype          soil3c;       /* microbial recycling pool C (slow)
                                   * (kgC m-2) */
    realtype          soil4c;       /* recalcitrant SOM C (humus, slowest)
                                   * (kgC m-2) */
    realtype          litr1n;       /* litter labile N (kgN m-2) */
    realtype          sminn;        /* soil mineral N (kgN m-2) */
} cninit_struct;
#endif

/* Global calibration coefficients */
typedef struct calib_struct
{
    realtype          ksath;
    realtype          ksatv;
    realtype          kinfv;
    realtype          kmach;
    realtype          kmacv;
    realtype          dinf;
    realtype          rzd;
    realtype          dmac;
    realtype          porosity;
    realtype          alpha;
    realtype          beta;
    realtype          areafv;
    realtype          areafh;
    realtype          vegfrac;
    realtype          albedo;
    realtype          rough;
    realtype          ec;
    realtype          ett;
    realtype          edir;
    realtype          rivrough;
    realtype          rivksath;
    realtype          rivksatv;
    realtype          rivbedthick;
    realtype          rivdepth;
    realtype          rivshpcoeff;
    realtype          prcp;      /* multiplier of precipitation (-) */
    realtype          sfctmp;    /* offset of surface air temperature (K) */
#if defined(_NOAH_)
    realtype          smcref;
    realtype          smcwlt;
    realtype          rsmin;
    realtype          drip;
    realtype          cmcmax;
    realtype          czil;
    realtype          fxexp;
    realtype          cfactr;
    realtype          rgl;
    realtype          hs;
#endif
#if defined(_BGC_)
    realtype          mortality;
    realtype          sla;
#endif
#if defined(_RT_)
    realtype          pco2;
    realtype          keq;
    realtype          ssa;
    realtype          site_den;
    realtype          prep_conc;
#endif
} calib_struct;

/* Model control parameters */
typedef struct ctrl_struct
{
    int             ascii;                  /* flag to turn on ascii output */
    int             waterbal;                 /* flag to turn on water balance
                                             * diagnostic output */
    int             write_ic;               /* flag to write model output as
                                             * initial conditions */
    int             nstep;                  /* number of external time steps
                                             * (when results can be printed) for
                                             * the whole simulation */
    int             cstep;                  /* current model step (from 0) */
    int             prtvrbl[MAXPRINT];


    int             tpprtvrbl[MAXPRINT];    /* time interval to tecplot output
                                             * average values of variables;
                                             * 0 = turn off output */
    int             init_type;              /* initialization mode:
                                             * 0 = relaxed mode,
                                             * 1 = use .ic file */
    int             unsat_mode;             /* unsaturation formulation:
                                             * 1 = kinematic, 2 = diffusion */
    int             surf_mode;              /* surface overland flow formulation
                                             * 1 = kinematic, 2 = diffusion */
    int             riv_mode;               /* river routing formulation:
                                             * 1 = kinematic, 2 = diffusion */
    int             etstep;                 /* land surface (ET) time step (s)
                                             */
    int             starttime;              /* start time of simulation (ctime)
                                             */
    int             endtime;                /* end time of simulation (ctime) */
    int             stepsize;               /* model step size (s) */
    int            *tout;                   /* model output times (ctime) */
    realtype          abstol;                 /* absolute solver tolerance (m) */
    realtype          reltol;                 /* relative solver tolerance (-) */
    realtype          initstep;               /* initial step size (s) */
    realtype          maxstep;                /* CVode maximum step size (s) */
    realtype          stmin;                  /* minimum allowed CVode max step
                                             * size (s) */
    realtype          nncfn;                  /* number of non-convergence
                                             * failures tolerance */
    realtype          nnimax;                 /* maximum number of non-linear
                                             * iterations */
    realtype          nnimin;                 /* minimum number of non-linear
                                             * iterations */
    realtype          decr;                   /* decrease factor (-)*/
    realtype          incr;                   /* increase factor (-)*/
    int             maxspinyears;           /* maximum number of years for
                                             * spinup run */
#if defined(_NOAH_)
    int             nsoil;                  /* number of standard soil layers */
    realtype          sldpth[MAXLYR];         /* thickness of soil layer (m) */
    int             rad_mode;               /* radiation forcing mode:
                                             * 0 = uniform, 1 = topographic */
#endif
#if defined(_BGC_)
    int             read_bgc_restart;       /* flag to read BGC restart file */
    int             write_bgc_restart;      /* flag to write BGC restart file */
#endif
#if defined(_CYCLES_)
    int             read_cycles_restart;
    int             write_cycles_restart;
#endif
} ctrl_struct;

/* Print variable control structure */
typedef struct varctrl_struct
{
    char            name[MAXSTRING];    /* name of output file */
    int             intvl;              /* output interval (s) */
    int             intr;
    int             upd_intvl;          /* 0: hydrology step
                                         * 1: land surface step
                                         * 2: CN step */
    int             nvar;               /* number of variables for print */
    const realtype  **var;                /* pointers to model variables */
    realtype         *buffer;             /* buffer for averaging variables */
    int             counter;            /* counter for averaging variables */
    FILE           *txtfile;            /* pointer to txt file */
    FILE           *datfile;            /* pointer to binary file */
    /* tecplot coordinate variables */
    realtype         *x;
    realtype         *y;
    realtype         *zmax;
    realtype         *zmin;
    int             nnodes;
    int           *node0;
    int           *node1;
    int           *node2;
} varctrl_struct;

/* Print structure */
typedef struct print_struct
{
    varctrl_struct  varctrl[MAXPRINT];
    varctrl_struct  tp_varctrl[MAXPRINT];
    int             nprint;            /* number of output variables */
    int             ntpprint;          /* number of tecplot output variables */
    FILE           *watbal_file;       /* pointer to water balance file */
    FILE           *cvodeperf_file;    /* pointer to CVode performance file */
} print_struct;

typedef struct pihm_struct_h
{
    siteinfo_struct siteinfo;
    filename_struct filename;
    meshtbl_struct  meshtbl;
    atttbl_struct   atttbl;
    soiltbl_struct  soiltbl;
    geoltbl_struct  geoltbl;
    lctbl_struct    lctbl;
    rivtbl_struct   rivtbl;
    shptbl_struct   shptbl;
    matltbl_struct  matltbl;
#if defined(_NOAH_)
    noahtbl_struct  noahtbl;
#endif
#if defined(_CYCLES_)
    agtbl_struct    agtbl;
    epconst_struct  epctbl[MAXCROP];
    opertbl_struct  opertbl[MAXOP];
#endif
#if defined(_BGC_)
    co2control_struct co2;
    ndepcontrol_struct ndepctrl;
    epctbl_struct   epctbl;
    cninit_struct   cninit;
#endif
    forc_struct     forc;
    elem_struct    *elem;
    river_struct   *river;
    calib_struct    cal;
    ctrl_struct     ctrl;
    print_struct    print;
} *pihm_struct;

#endif

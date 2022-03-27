#ifndef RIVER_STRUCT_HEADER
#define RIVER_STRUCT_HEADER

/* River attribute */
typedef struct river_attrib_struct
{
    int             riverbc_type;    /* river boundary condition type */
} river_attrib_struct;

/* River topography parameters */
typedef struct river_topo_struct
{
    realtype          area;          /* area of element (m2) */
    realtype          x;             /* x of centroid (m) */
    realtype          y;             /* y of centroid (m) */
    realtype          zmin;          /* bedrock elevation (m) */
    realtype          zmax;          /* river bank elevation (m) */
    realtype          zbed;          /* river bed elevation (m) */
    realtype          node_zmax;     /* elevation of the downstream node (m) */
    realtype          dist_left;     /* distance to left neighbor (m) */
    realtype          dist_right;    /* distance to right neighbor (m) */
} river_topo_struct;

/* River water states */
typedef struct river_wstate_struct
{
    realtype          stage;    /* river stage (m) */
    realtype          gw;       /* groundwater level (m) */
} river_wstate_struct;

/* River water fluxes */
typedef struct river_wflux_struct
{
    realtype          rivflow[NUM_RIVFLX];    /* river fluxes (m3 s-1) */
} river_wflux_struct;

/* River shape parameters */
typedef struct shp_struct
{
    realtype          depth;         /* river channel depth (m) */
    int             intrpl_ord;    /* interpolation order (shape of channel) */
    realtype          coeff;         /* width coefficient */
    realtype          length;        /* length of channel (m) */
    realtype          width;         /* width of channel (m) */
} shp_struct;

typedef struct matl_struct
{
    realtype          rough;       /* river channel roughness (s m-1/3) */
    realtype          cwr;         /* discharge coefficient (-) */
    realtype          ksath;       /* bank hydraulic conductivity (m s-1) */
    realtype          ksatv;       /* bed hydraulic conductivity (m s-1) */
    realtype          bedthick;    /* bed thickness (m) */
    realtype          porosity;    /* bed porosity (m3 m-3) */
    realtype          smcmin;      /* bed residual soil moisture content (m3 m-3)
                                  */
#if defined(_CYCLES_)
    realtype          bd;
#endif
} matl_struct;

/* River boundary conditions */
typedef union river_bc_struct
{
    realtype          head;   /* value of Dirichlet-type boundary condition (m) */
    realtype          flux;   /* value of Neumann-type boundary condition (m3 s-1)
                             */
} river_bc_struct;

/* River initial conditions */
typedef struct river_ic_struct
{
    realtype          stage;
    realtype          gw;
} river_ic_struct;

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

/* River structure */
typedef struct river_struct
{
    int             ind;         /* river index */
    int             leftele;     /* left neighboring element */
    int             rightele;    /* right neighboring element */
    int             fromnode;    /* upstream node */
    int             tonode;      /* downstream node */
    int             down;        /* down stream channel segment */
    river_attrib_struct attrib;
    river_topo_struct topo;
    shp_struct      shp;
    matl_struct     matl;
    river_wstate_struct ws;
    river_wstate_struct ws0;
    river_wflux_struct wf;
    river_ic_struct ic;
    river_bc_struct bc;
#if defined(_CYCLES_)
    river_nstate_struct ns;
    river_cyclesic_struct restart_input;
    river_solute_struct no3sol;
    river_solute_struct nh4sol;
#endif
#if defined(_BGC_) && !defined(_LUMPED_) && !defined(_LEACHING_)
    river_nstate_struct ns;
    river_nflux_struct nf;
    river_solute_struct nsol;
    river_bgcic_struct restart_input;
    river_bgcic_struct restart_output;
#endif
} river_struct;
#endif

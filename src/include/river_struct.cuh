#ifndef RIVER_STRUCT_HEADER_D
#define RIVER_STRUCT_HEADER_D

/* River attribute */
typedef struct river_attrib_struct_d
{
    int           *riverbc_type;    /* river boundary condition type */
} river_attrib_struct_d;

/* River topography parameters */
typedef struct river_topo_struct_d
{
    double          *area;          /* area of element (m2) */
    double          *x;             /* x of centroid (m) */
    double          *y;             /* y of centroid (m) */
    double          *zmin;          /* bedrock elevation (m) */
    double          *zmax;          /* river bank elevation (m) */
    double          *zbed;          /* river bed elevation (m) */
    double          *node_zmax;     /* elevation of the downstream node (m) */
    double          *dist_left;     /* distance to left neighbor (m) */
    double          *dist_right;    /* distance to right neighbor (m) */
} river_topo_struct_d;

/* River water states */
typedef struct river_wstate_struct_d
{
    double          *stage;    /* river stage (m) */
    double          *gw;       /* groundwater level (m) */
} river_wstate_struct_d;

/* River water fluxes */
typedef struct river_wflux_struct_d
{
    double          *rivflow[NUM_RIVFLX];    /* river fluxes (m3 s-1) */
} river_wflux_struct_d;

/* River shape parameters */
typedef struct shp_struct_d
{
    double          *depth;         /* river channel depth (m) */
    int             *intrpl_ord;    /* interpolation order (shape of channel) */
    double          *coeff;         /* width coefficient */
    double          *length;        /* length of channel (m) */
    double          *width;         /* width of channel (m) */
} shp_struct_d;


typedef struct matl_struct_d
{
    double          *rough;       /* river channel roughness (s m-1/3) */
    double          *cwr;         /* discharge coefficient (-) */
    double          *ksath;       /* bank hydraulic conductivity (m s-1) */
    double          *ksatv;       /* bed hydraulic conductivity (m s-1) */
    double          *bedthick;    /* bed thickness (m) */
    double          *porosity;    /* bed porosity (m3 m-3) */
    double          *smcmin;      /* bed residual soil moisture content (m3 m-3) */                              
} matl_struct_d;								  
								  
/* River boundary conditions */
typedef struct river_bc_struct_d
{
    double          *head;   /* value of Dirichlet-type boundary condition (m) */
    double          *flux;   /* value of Neumann-type boundary condition (m3 s-1) */                   
} river_bc_struct_d;
							  
															  								  
/* River structure in CUDA*/
typedef struct river_struct_d
{
    int             *ind;         /* river index */
    int             *leftele;     /* left neighboring element */
    int             *rightele;    /* right neighboring element */
    int             *fromnode;    /* upstream node */
    int             *tonode;      /* downstream node */
    int             *down;        /* down stream channel segment */
	
    river_attrib_struct_d attrib;
    river_topo_struct_d topo;
    shp_struct_d      shp;
    matl_struct_d     matl;
    river_wstate_struct_d ws;
    river_wstate_struct_d ws0;
    river_wflux_struct_d wf;
    river_bc_struct_d bc;

} river_struct_d;
#endif
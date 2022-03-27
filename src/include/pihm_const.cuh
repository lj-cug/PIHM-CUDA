#ifndef PIHM_CONST_HEADER_D
#define PIHM_CONST_HEADER_D

/* Physical parameters */
#define PI          3.14159265
#define DAYINSEC    86400      /* number of seconds in a day */
#define GRAV        9.80665    /* gravity constant (m s-2) */
#define CP          1004.0     /* specific heat capacity of air
                                * (J kg-1 K-1) */
#define RHOH2O      1000.0     /* water density (kg m-3) */
#define LVH2O       2.501e6    /* latent heat of vaporization (J kg-1) */
#define SIGMA       5.67e-8    /* Stefan-Boltzmann constant (W m-2 K-4) */
#define RD          287.04     /* gas constant for dry air (J kg-1 K-1) */
#define RV          461.5      /* gas constant for water vapor (J kg-1 K-1) */
#define CPH2O       4.218e3    /* specific heat capacity of water
                                * (J kg-1 K-1) */
#define CPICE       2.106e3    /* specific heat capacity of ice
                                * (J kg-1 K-1) */
#define LSUBF       3.335e5    /* latent heat of fusion (J kg-1) */
#define EMISSI_S    0.95       /* emissivity of snow (-) */
#define TFREEZ      273.15     /* freezing point (K) */
#define LSUBS       2.83e6     /* latent heat of sublimation (J kg-1) */

#define F_OK    0

#if defined(_LUMPED_)
# define LUMPED    nelem
#endif

/* Simulation mode */
#define NORMAL_MODE        0
#define SPINUP_MODE        1
#define ACC_SPINUP_MODE    2

/* Default bad value */
#define BADVAL    -999

/* Model steps */
#define HYDROL_STEP    0
#define LS_STEP        1
#define CN_STEP        2

/* Maximum number of output files */
#define MAXPRINT    1024

/* Meteorological forcing related */
#define NUM_METEO_VAR    7    /* number of meteorological forcing variables */
#define PRCP_TS          0    /* index of precipitation forcing */
#define SFCTMP_TS        1    /* index of air temperature forcing */
#define RH_TS            2    /* index of RH forcing */
#define SFCSPD_TS        3    /* index of wind speed forcing */
#define SOLAR_TS         4    /* index of solar radiation forcing */
#define LONGWAVE_TS      5    /* index of longwave radiation forcing */
#define PRES_TS          6    /* index of surface pressure forcing */

/* Radiation forcing variables */
#define UNIF_SOL     0
#define TOPO_SOL     1
#define SOLDIR_TS    0    /* index of direct solar radiation forcing */
#define SOLDIF_TS    1    /* index of diffused solar radiation forcing */

/* Number of edges of an element */
#define NUM_EDGE    3

#define WS_ZMAX    0
#define WS_ZMIN    1
#define WS_AREA    2

/* Number of river fluxes of a river segment */
#define NUM_RIVFLX    11

/* Hydrology parameters */
#define PSIMIN        -70.0    /* minimum psi allowed (m) */
#define DEPRSTG       1E-4     /* depression storage (m) */
#define GRADMIN       5E-8     /* minimum hydraulic gradient (m m-1) */
#define SATMIN        0.1      /* minimum saturation ratio (-) */
#define RIVDPTHMIN    0.05     /* minimum river depth (m) */
#define RIVGRADMIN    0.05     /* minimum river hydraulic gradient (m m-1) */
#define CMCFACTR      2E-4     /* canopy water capacity per LAI (m) */
#define SH2OMIN       0.02     /* minimum sh2o (m3 m-3) in Noah LSM */

/* Maximum of soil layers in Flux-PIHM */
#define MAXLYR    11

/* Land cover parameters */
#define NLCTYPE    40    /* number of land cover types */
#define ISURBAN    13    /* land cover type representing urban */

/* Land cover types */
#define ENF              1
#define EBF              2
#define DNF              3
#define DBF              4
#define MIXF             5
#define CLOSE_SHRUB      6
#define OPEN_SHRUB       7
#define WOODY_SAVANNA    8
#define SAVANNA          9
#define GRASS            10
#define PWL              11
#define CROP             12
#define URBAN_BUILDUP    13
#define CROP_NATURAL     14
#define SNOW_ICE         15
#define BARREN           16
#define WATER            17
#define WOOD_TUNDRA      18
#define MIX_TUNDRA       19
#define BARREN_TUNDRA    20

/* Soil textures */
#define SAND               0
#define LOAMY_SAND         1
#define SANDY_LOAM         2
#define LOAM               3
#define SILT_LOAM          4
#define SILT               5
#define SANDY_CLAY_LOAM    6
#define CLAY_LOAM          7
#define SILTY_CLAY_LOAM    8
#define SANDY_CLAY         9
#define SILTY_CLAY         10
#define CLAY               11

/* River fluxes */
#define UP_CHANL2CHANL       0
#define DOWN_CHANL2CHANL     1
#define LEFT_SURF2CHANL      2
#define RIGHT_SURF2CHANL     3
#define LEFT_AQUIF2CHANL     4
#define RIGHT_AQUIF2CHANL    5
#define CHANL_LKG            6
#define LEFT_AQUIF2AQUIF     7
#define RIGHT_AQUIF2AQUIF    8
#define DOWN_AQUIF2AQUIF     9
#define UP_AQUIF2AQUIF       10

/* River boundary condition types */
#define DIRICHLET         -1
#define NEUMANN           -2
#define ZERO_DPTH_GRAD    -3
#define CRIT_DPTH         -4

/* Element boundary condition types */
#define NO_FLOW    0

/* River segment interpolation order */
#define RECTANGLE    1
#define TRIANGLE     2
#define QUADRATIC    3
#define CUBIC        4

/* Approximation */
#define KINEMATIC    1
#define DIFF_WAVE    2

/* Initialization type */
#define RELAX       0
#define RST_FILE    1

/* Average flux */
#define SUM    0
#define AVG    1

/* Maximum allowable difference between simulation cycles in subsurface water
 * storage at steady-state (m) */
#define SPINUP_W_TOLERANCE    0.01

/* Ecosystem constants */
#define RAD2PAR             0.45    /* ratio PAR / SWtotal (-) */
#define EPAR                4.55    /* (umol/J) PAR photon energy ratio */
#define SOIL1_CN            12.0    /* C:N for fast microbial recycling pool */
#define SOIL2_CN            12.0    /* C:N for slow microbial recycling pool */
#define SOIL3_CN            10.0    /* C:N for recalcitrant SOM pool (humus) */
#define SOIL4_CN            10.0    /* C:N for recalcitrant SOM pool (humus) */
#define GRPERC              0.3     /* growth resp per unit of C grown (-) */
#define GRPNOW              1.0     /* proportion of storage growth resp at
                                     * fixation (-) */
#define PPFD50              75.0    /* PPFD for 1/2 stomatal closure
                                     * (umol m-2 s-1) */
#define DENITRIF_PROPORTION 0.01    /* fraction of mineralization to volatile */
#if defined(_CYCLES_)
#define MOBILEN_PROPORTION  0.67    /* fraction mineral N avail for leaching */
#else
#define MOBILEN_PROPORTION  0.1     /* fraction mineral N avail for leaching */
#endif

/* Respiration fractions for fluxes between compartments (-) */
#define RFL1S1    0.39    /* transfer from litter 1 to soil 1 */
#define RFL2S2    0.55    /* transfer from litter 2 to soil 2 */
#define RFL4S3    0.29    /* transfer from litter 4 to soil 3 */
#define RFS1S2    0.28    /* transfer from soil 1 to soil 2 */
#define RFS2S3    0.46    /* transfer from soil 2 to soil 3 */
#define RFS3S4    0.55    /* transfer from soil 3 to soil 4 */

/* Base decomposition rate constants (day -1) */
#define KL1_BASE      0.7       /* labile litter pool */
#define KL2_BASE      0.07      /* cellulose litter pool */
#define KL4_BASE      0.014     /* lignin litter pool */
#define KS1_BASE      0.07      /* fast microbial recycling pool */
#define KS2_BASE      0.014     /* medium microbial recycling pool */
#define KS3_BASE      0.0014    /* slow microbial recycling pool */
#define KS4_BASE      0.0001    /* recalcitrant SOM (humus) pool */
#define KFRAG_BASE    0.001     /* physical fragmentation of coarse woody debris
                                 */

/* Decomposition acceleration terms */
#define KS1_ACC    1.0
#define KS2_ACC    1.0
#define KS3_ACC    5.0
#define KS4_ACC    70.0

/* This constant determines the lower limit of state variables before they are
 * set to 0.0 to control rounding and overflow errors */
#define CRIT_PREC    1e-20

/* This constant is used in if conditions where floating point values are
 * compared */
#define FLT_COND_TOL    1e-10

/* Maximum allowable trend in slow soil carbon at steady-state (kgC m-2 yr-1) */
#define SPINUP_C_TOLERANCE    0.0005

/* Allocation parameters */
#define DAYSNDEPLOY                 365.0
#define DAYSCRECOVER                365.0
#define BULK_DENITRIF_PROPORTION    0.5

/* Output variables */
#define YEARLY_OUTPUT     -1
#define MONTHLY_OUTPUT    -2
#define DAILY_OUTPUT      -3
#define HOURLY_OUTPUT     -4

/* Output variable types */
#define ELEMVAR     0
#define RIVERVAR    1  



#endif
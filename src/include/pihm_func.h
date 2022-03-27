#ifndef PIHM_FUNC_HEADER
#define PIHM_FUNC_HEADER

#define _ARITH_

/* State variables */
#define SURF(i)          i
#define UNSAT(i)         i + nelem
#define GW(i)            i + 2 * nelem
#define RIVSTG(i)        i + 3 * nelem
#define RIVGW(i)         i + 3 * nelem + nriver
#if defined(_FBR_)
# define FBRUNSAT(i)     i + 3 * nelem + 2 * nriver
# define FBRGW(i)        i + 4 * nelem + 2 * nriver
#endif
#if defined(_BGC_) && !defined(_LUMPED_)
# define SURFN(i)        i + 3 * nelem + 2 * nriver
# define SMINN(i)        i + 4 * nelem + 2 * nriver
# define STREAMN(i)      i + 5 * nelem + 2 * nriver
# define RIVBEDN(i)      i + 5 * nelem + 3 * nriver
#else
# define LUMPED_SMINN    3 * nelem + 2 * nriver
#endif
#if defined(_CYCLES_)
# define NO3(i)          i + 3 * nelem + 2 * nriver
# define NH4(i)          i + 4 * nelem + 2 * nriver
# define STREAMNO3(i)    i + 5 * nelem + 2 * nriver
# define RIVBEDNO3(i)    i + 5 * nelem + 3 * nriver
# define STREAMNH4(i)    i + 5 * nelem + 4 * nriver
# define RIVBEDNH4(i)    i + 5 * nelem + 5 * nriver
#endif

#define AvgElev(...)      _WsAreaElev(WS_ZMAX, __VA_ARGS__)
#define AvgZmin(...)      _WsAreaElev(WS_ZMIN, __VA_ARGS__)
#define TotalArea(...)    _WsAreaElev(WS_AREA, __VA_ARGS__)

/* CVode functions */
#if defined(_CVODE_OMP)
# define N_VNew(N)    N_VNew_OpenMP(N, nthreads)
# define NV_DATA      NV_DATA_OMP
# define NV_Ith       NV_Ith_OMP
#else
# define N_VNew(N)    N_VNew_Serial(N)
# define NV_DATA      NV_DATA_S
# define NV_Ith       NV_Ith_S
#endif

/* PIHM system function */
#define PIHMexit(...)               _custom_exit(__FILE__, __LINE__, __FUNCTION__, debug_mode,  __VA_ARGS__)
#define PIHMprintf(...)             _custom_printf(__FILE__, __LINE__, __FUNCTION__, debug_mode, verbose_mode, __VA_ARGS__)
#if defined(_WIN32) || defined(_WIN64)
# define PIHMmkdir(path)            _mkdir((path))
# define PIHMaccess(path, amode)    _access((path), (amode))
#else
# define PIHMmkdir(path)            mkdir(path, 0755)
# define PIHMaccess(path, amode)    access((path), (amode))
#endif
#if defined(_MSC_VER)
# define timegm                     _mkgmtime
# define strcasecmp                 _stricmp
#endif

#define Cycles_exit      PIHMexit
#define Cycles_printf    PIHMprintf

/*
 * Function Declarations
 */
realtype        _WsAreaElev(int, const elem_struct *);
void            AdjCVodeMaxStep(void *, ctrl_struct *);
void            BackupInput(const char *, const filename_struct *);	
void            CalcModelStep(ctrl_struct *);	
int             CheckCVodeFlag(int);
void            CheckDy(realtype, const char *, const char *, int, realtype);
#if defined(_BGC_)
int             CheckSteadyState(const elem_struct *, realtype, int, int, int);
#else
int             CheckSteadyState(const elem_struct *, realtype, int, int);
#endif
void            CorrElev(elem_struct *, river_struct *);
void            CreateOutputDir(char *);
realtype        FieldCapacity(realtype, realtype, realtype, realtype);
void            FreeAtttbl(atttbl_struct *);
void            FreeCtrl(ctrl_struct *);
void            FreeForc(forc_struct *);
void            FreeLctbl(lctbl_struct *);
void            FreeMatltbl(matltbl_struct *);
void            FreeMeshtbl(meshtbl_struct *);
void            FreeMem(pihm_struct);
void            FreeRivtbl(rivtbl_struct *);
void            FreeShptbl(shptbl_struct *);
void            FreeSoiltbl(soiltbl_struct *);					  
void            InitEFlux(eflux_struct *);
void            InitEState(estate_struct *);
void            InitForc(elem_struct *, forc_struct *, const calib_struct *);                                 
void            InitLc(elem_struct *, const lctbl_struct *,
                       const calib_struct *);
void            InitMesh(elem_struct *, const meshtbl_struct *);
void            InitOutputFile(print_struct *, const char *, int, int);
void            InitPrtVarCtrl(const char *, const char *, int, int, int,
                               varctrl_struct *);
void            InitRiver(river_struct *, elem_struct *, const rivtbl_struct *,
                          const shptbl_struct *, const matltbl_struct *, const meshtbl_struct *,
                          const calib_struct *);
void            InitRiverWFlux(river_wflux_struct *);
void            InitRiverWState(river_wstate_struct *);
#if defined(_NOAH_)
void            InitSoil(elem_struct *, const soiltbl_struct *,
                         const noahtbl_struct *, const calib_struct *);
#else
void            InitSoil(elem_struct *, const soiltbl_struct *,
                         const calib_struct *);
#endif
void            InitSurfL(elem_struct *, const river_struct *,
                          const meshtbl_struct *);
void            InitTecPrtVarCtrl(const char *, const char *, int, int, int,
                                  int, int, varctrl_struct *);
void            InitTopo(elem_struct *, const meshtbl_struct *);
void            InitVar(elem_struct *, river_struct *, N_Vector);
void            InitWbFile(char *, char *, FILE *);
void            InitWFlux(wflux_struct *);
void            InitWState(wstate_struct *);
void            IntrplForc(tsdata_struct *, int, int);
#if defined(_CYCLES_)
void    MapOutput(const int *, const int *, const epconst_struct [],
    const elem_struct *, const river_struct *, const meshtbl_struct *,
    const char *, print_struct *);
#else
void    MapOutput(const int *, const int *, const elem_struct *,
        const river_struct *, const meshtbl_struct *, const char *, print_struct *);
#endif
#if defined(_FBR_)
void            MassBalance(const wstate_struct *, const wstate_struct *,
                            wflux_struct *, realtype *, const soil_struct *, const geol_struct *, 
                            realtype, realtype);
#else
void            MassBalance(const wstate_struct *, const wstate_struct *,
                            wflux_struct *, realtype *, const soil_struct *, realtype, realtype);
#endif
int             NumStateVar(void);

realtype        MonthlyLai(int, int);
realtype        MonthlyMf(int);
realtype        MonthlyRl(int, int);	

void            ParseCmdLineParam(int, char *[], char *);
void            PrintCVodeFinalStats(void *);
void            PrintData(varctrl_struct *, int, int, int, int);
void            PrintDataTecplot(varctrl_struct *, int, int, int);

pihm_t_struct   PIHMTime(int);
int             PrintNow(int, int, const pihm_t_struct *);
void            PrintPerf(void *, int, int, realtype, realtype, realtype, FILE *);
realtype          PtfAlpha(realtype, realtype, realtype, realtype, int);
realtype          PtfBeta(realtype, realtype, realtype, realtype, int);
realtype          PtfKv(realtype, realtype, realtype, realtype, int);
realtype          PtfThetar(realtype, realtype);
realtype          PtfThetas(realtype, realtype, realtype, realtype, int);
realtype          Qtz(int);
void            ReadAlloc(pihm_struct);
void            ReadAtt(const char *, atttbl_struct *);
void            ReadBc(const char *, forc_struct *, const atttbl_struct *);
void            ReadCalib(const char *, calib_struct *);
void            ReadForc(const char *, forc_struct *);
void            ReadIc(const char *, elem_struct *, river_struct *);
int             ReadKeyword(const char *, const char *, void *, char,
                            const char *, int);
void            ReadLai(const char *, forc_struct *, const atttbl_struct *);
void            ReadLc(const char *, lctbl_struct *);
void            ReadMesh(const char *, meshtbl_struct *);
void            ReadPara(const char *, ctrl_struct *);
int             ReadPrtCtrl(const char *, const char *, const char *, int);
void            ReadRiver(const char *, rivtbl_struct *, shptbl_struct *,
                          matltbl_struct *, forc_struct *);
void            ReadSoil(const char *, soiltbl_struct *);
void            ReadTecplot(const char *, ctrl_struct *);
int             ReadTS(const char *, int *, realtype *, int);						 
realtype          RiverEqWid(int, realtype, realtype);
#if defined(_OPENMP)
void            RunTime(realtype, realtype *, realtype *);
#else
void            RunTime (clock_t, realtype *, realtype *);
#endif
void            RelaxIc(elem_struct *, river_struct *);
int             SoilTex(realtype, realtype);
void            SolveCVode(int, int *, int, realtype, void *, N_Vector);
void            StartupScreen(void);
int             StrTime(const char *);								

void            UpdPrintVar(varctrl_struct *, int, int);
void            UpdPrintVarT(varctrl_struct *, int);
realtype          WiltingPoint(realtype, realtype, realtype, realtype);
realtype          KrFunc(realtype, realtype);  // used in soil.c and vert_flow.c

#ifndef _CVODE_CUDA
// 注销CPU上执行的函数，在pihm_func.cuh中重新定义了GPU上执行的对应的函数

#if defined(_NOAH_)
void            ApplyForc(forc_struct *, elem_struct *, int, int,
                          const siteinfo_struct *);
void            ApplyMeteoForc(forc_struct *, elem_struct *, int, int,
                               const siteinfo_struct *);						  
#else
void            ApplyForc(forc_struct *, elem_struct *, int);
void            ApplyMeteoForc(forc_struct *, elem_struct *, int);
#endif

#if defined(_BGC_) || defined(_CYCLES_)
void            ApplyLai(elem_struct *);
#else
void            ApplyLai(forc_struct *, elem_struct *, int);
#endif

void            Initialize(pihm_struct, N_Vector);	
void            SetCVodeParam(pihm_struct, void *, N_Vector);
void            PIHM(pihm_struct, void *, N_Vector, realtype);
void            Spinup(pihm_struct, N_Vector, void *);   
void            ApplyBc(forc_struct *, elem_struct *, river_struct *, int);
void            ApplyElemBc(forc_struct *, elem_struct *, int);
void            ApplyRiverBc(forc_struct *, river_struct *, int);
void            IntcpSnowEt(int, realtype, elem_struct *, const calib_struct *);
               
int             ODE(realtype, N_Vector, N_Vector, void *);	
realtype         AvgKv(const soil_struct *, realtype, realtype, realtype);
realtype         AvgH(realtype, realtype, realtype);
realtype         AvgHsurf(realtype, realtype, realtype);
void            BoundFluxElem(int, int, const bc_struct *,
                              const wstate_struct *, const topo_struct *, 
                              const soil_struct *, wflux_struct *);
realtype          BoundFluxRiver(int, const river_wstate_struct *,
    const river_topo_struct *, const shp_struct *, const matl_struct *,
    const river_bc_struct *bc);
realtype          ChanFlowElemToRiver(const elem_struct *, realtype,
                                    const river_struct *, realtype);
realtype          ChanFlowRiverToRiver(const river_struct *, const river_struct *, int);    
realtype          ChanLeak(const river_wstate_struct *, const river_topo_struct *,
                         const shp_struct *, const matl_struct *);
realtype          DhByDl(const realtype *, const realtype *, const realtype *);
realtype          EffKh(const soil_struct *, realtype);
realtype          EffKinf(const soil_struct *, realtype, realtype, realtype, realtype,
                        realtype);
realtype          EffKv(const soil_struct *, realtype, int);
void              EtExtract(elem_struct *);
void              FrictSlope(const elem_struct *, const river_struct *, int,
                           realtype *, realtype *);

void              Hydrol(elem_struct *, river_struct *, const ctrl_struct *);
realtype          Infil(const wstate_struct *, const wstate_struct *,
                      const wflux_struct *, const topo_struct *, const soil_struct *, realtype);

void              LateralFlow(elem_struct *, const river_struct *, int);	
realtype          OutletFlux(int, const river_wstate_struct *,
              const river_topo_struct *, const shp_struct *, const matl_struct *,
              const river_bc_struct *);	  
realtype          OverLandFlow(realtype, realtype, realtype, realtype, realtype);
realtype          OvlFlowElemToElem(const elem_struct *, const elem_struct *,  int,  realtype, int);	
realtype          Psi(realtype, realtype, realtype);                                 
realtype          Recharge(const wstate_struct *, const wflux_struct *,
                         const soil_struct *);
realtype          RiverCroSectArea(int, realtype, realtype);
void              RiverFlow(elem_struct *, river_struct *, int);
realtype          RiverPerim(int, realtype, realtype);
void              RiverToElem(river_struct *, elem_struct *, elem_struct *);	
realtype          OvlFlowElemToRiver(const elem_struct *, const river_struct *);
realtype          SubFlowElemToElem(const elem_struct *, const elem_struct *, int);
realtype          SubFlowElemToRiver(const elem_struct *, realtype,
                                     const river_struct *, realtype, realtype);
realtype          SubFlowRiverToRiver(const river_struct *, realtype,
                                      const river_struct *, realtype);	
realtype          SurfH(realtype);	
void              VerticalFlow(elem_struct *, realtype);
void              Summary(elem_struct *, river_struct *, N_Vector, realtype);	
void              PrintWaterBal(FILE *, int, int, int, const elem_struct *,
                              const river_struct *);						  
void              PrintInit(const elem_struct *, const river_struct *,
                            const char *, int, int, int, int);
#endif

/*
 * Fractured bedrock functions
 */
#if defined(_FBR_)
realtype          FbrBoundFluxElem(int, int, const bc_struct *,
                                   const wstate_struct *, const topo_struct *, const geol_struct *);
realtype          FbrFlowElemToElem(const elem_struct *, const elem_struct *,
                                   realtype, realtype);
realtype          FbrInfil(const wstate_struct *, const soil_struct *,
                           const geol_struct *, const topo_struct *);
realtype          FbrRecharge(const wstate_struct *, const wflux_struct *,
                              const geol_struct *);
void            FreeGeoltbl(geoltbl_struct *);
void            InitGeol (elem_struct *, const geoltbl_struct *,
                          const calib_struct *);
void            ReadBedrock(const char *, atttbl_struct *, meshtbl_struct *,
                            ctrl_struct *);
void            ReadGeol(const char *, geoltbl_struct *);
#endif

/*
 * Noah functions
 */
#if defined(_NOAH_)

void            CalcLatFlx(const pstate_struct *, wflux_struct *);
void            CalcSlopeAspect(elem_struct *, const meshtbl_struct *);
void            DefSldpth(realtype *, int *, realtype *, realtype, const realtype *,int); 
int             FindLayer(const realtype *, int, realtype);

#ifndef _CVODE_CUDA	
int             FindWaterTable(const realtype *, int, realtype, realtype *);
#endif
realtype        GwTransp(realtype, const realtype *, int, int);					
void            InitLsm(elem_struct *, const ctrl_struct *,
                        const noahtbl_struct *, const calib_struct *);
realtype        Mod(realtype, realtype);
void            ReadLsm(const char *, siteinfo_struct *, ctrl_struct *,
                        noahtbl_struct *);
void            ReadRad(const char *, forc_struct *);
void            RootDist(const realtype *, int, int, realtype *);
void            SunPos(const siteinfo_struct *, int, spa_data *);
realtype        TopoRadn(const topo_struct *, realtype, realtype, realtype, realtype);


#ifndef _CVODE_CUDA	
void            AdjSmProf(const soil_struct *, const pstate_struct *,
                          const realtype *, realtype, wflux_struct *, wstate_struct *);
void            AlCalc(pstate_struct *, realtype, int);
void            CalHum(pstate_struct *, estate_struct *);

# if defined(_CYCLES_)
void            CanRes(const estate_struct *, pstate_struct *);
# else
void            CanRes(const wstate_struct *, const estate_struct *,
    const eflux_struct *, pstate_struct *, const soil_struct *,
    const epconst_struct *);
# endif

realtype        CSnow(realtype);
void            DEvap(const wstate_struct *, wflux_struct *,
    const pstate_struct *, const lc_struct *, const soil_struct *);
	
# if defined(_CYCLES_)
void            Evapo(const soil_struct *, const lc_struct *,
                      const pstate_struct *, const estate_struct *es,
                      const cstate_struct *, realtype, crop_struct [], wstate_struct *,
                      wflux_struct *);
# else
void            Evapo(const wstate_struct *, wflux_struct *,
                      const pstate_struct *, const lc_struct *, const soil_struct *, realtype);
# endif

realtype        FrozRain(realtype, realtype);
void            HRT(wstate_struct *, const estate_struct *, eflux_struct *,
                    const pstate_struct *, const lc_struct *, const soil_struct *, realtype *,
                    realtype, realtype, realtype, realtype, realtype *, realtype *, realtype *);
void            Noah(elem_struct *, realtype);	
void            NoahHydrol(elem_struct *, realtype);

# if defined(_CYCLES_)
void            NoPac(const soil_struct *, const lc_struct *,
    const cstate_struct *, realtype, realtype, crop_struct [], pstate_struct *,
    wstate_struct *, wflux_struct *, estate_struct *, eflux_struct *);
# else
void            NoPac(wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, pstate_struct *, const lc_struct *, const soil_struct *,
    realtype, realtype);
# endif

void            PcpDrp(wstate_struct *, wflux_struct *, const lc_struct *,
    realtype, realtype);
void            Penman(wflux_struct *, const estate_struct *, eflux_struct *,
                       pstate_struct *, realtype *, realtype, int, int);
realtype          Pslhs(realtype);
realtype          Pslhu(realtype);
realtype          Pslms(realtype);
realtype          Pslmu(realtype);
realtype          Psphs(realtype);
realtype          Psphu(realtype);
realtype          Pspms(realtype);
realtype          Pspmu(realtype);
void              Rosr12(realtype *, const realtype *, const realtype *, realtype *,
                         const realtype *, realtype *, int);				   
void              SfcDifOff(pstate_struct *, const lc_struct *, realtype, realtype, int);

# if defined(_CYCLES_)
void            SFlx(const cstate_struct *, realtype, soil_struct *, lc_struct *,
    crop_struct [], pstate_struct *, wstate_struct *, wflux_struct *,
    estate_struct *, eflux_struct *);
# else
void            SFlx(wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, pstate_struct *, lc_struct *, epconst_struct *,
    soil_struct *, realtype);
# endif

void  ShFlx(wstate_struct *, estate_struct *, eflux_struct *,
            const pstate_struct *, const lc_struct *, const soil_struct *, realtype,
            realtype, realtype, realtype);
			
# if defined(_CYCLES_)
void SmFlx(const soil_struct *, const cstate_struct *, realtype, pstate_struct *,
           wstate_struct *, wflux_struct *);
# else
void  SmFlx(wstate_struct *, wflux_struct *, pstate_struct *,
                      const soil_struct *, realtype);
# endif

realtype SnFrac(realtype, realtype, realtype);
void SnkSrc(realtype *, realtype, realtype, realtype *,
            const soil_struct *, const realtype *, realtype, int, realtype);

# if defined(_CYCLES_)
void            SnoPac(const soil_struct *, const lc_struct *,
    const cstate_struct *, int, realtype, realtype, realtype, realtype, crop_struct [],
    pstate_struct *, wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *);
# else
void            SnoPac(wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, pstate_struct *, const lc_struct *, const soil_struct *,
    int, realtype, realtype, realtype, realtype);
# endif

void            SnowNew(const estate_struct *, realtype, pstate_struct *);
void            SnowPack(realtype, realtype, realtype *, realtype *, realtype, realtype);
realtype        Snowz0(realtype, realtype, realtype);

# if defined(_CYCLES_)
void SRT(const soil_struct *, const cstate_struct *, realtype,
        pstate_struct *, wstate_struct *, wflux_struct *, realtype *, realtype *,
        realtype *, realtype *, realtype *);
# else
void SRT(wstate_struct *, wflux_struct *, pstate_struct *,
         const soil_struct *, realtype *, realtype *, realtype *, realtype *, realtype *);
# endif

void SStep(wstate_struct *, wflux_struct *, pstate_struct *,
           const soil_struct *, realtype *, realtype *, realtype *, realtype *, realtype *,
           realtype);
realtype          TBnd(realtype, realtype, const realtype *, realtype, int, int);
realtype          TDfCnd(realtype, realtype, realtype, realtype, realtype);
realtype          TmpAvg(realtype, realtype, realtype, const realtype *, int);
void            Transp(const wstate_struct *, wflux_struct *,
                       const pstate_struct *, const lc_struct *, const soil_struct *);
void            WDfCnd(realtype *, realtype *, realtype, realtype, const soil_struct *);

#endif  /* #ifndef _CVODE_CUDA	*/

#endif /* #if defined(_NOAH_) */

#if defined(_DAILY_)
void            DailyVar(int, int, elem_struct *);
void            InitDailyStruct(elem_struct *);
#endif

#if defined(_BGC_) || defined(_CYCLES_)
void            SetAbsTol(realtype, realtype, N_Vector);
#endif

#if defined(_BGC_)
void            BackgroundLitterfall(const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, nflux_struct *);
void            CanopyCond(const epconst_struct *, epvar_struct *,
    const eflux_struct *, const pstate_struct *, const soil_struct *,
    const daily_struct *);
void            CheckCarbonBalance(const cstate_struct *, realtype *);
void            CheckNitrogenBalance(const nstate_struct *, realtype *);
void            CSummary(const cflux_struct *, const cstate_struct *,
                         summary_struct *);
void            DailyAllocation(cflux_struct *, const cstate_struct *,
                                nflux_struct *, const nstate_struct *, const epconst_struct *,
                                epvar_struct *, ntemp_struct *);
void            DailyBgc(pihm_struct, int);
void            DailyCarbonStateUpdate(cflux_struct *, cstate_struct *, int,
                                       int, int);
void            DailyNitrogenStateUpdate(nflux_struct *, nstate_struct *,
                                         solute_struct *, int, int, int);
void            Decomp(realtype, const epconst_struct *, epvar_struct *,
                       const cstate_struct *, cflux_struct *, const nstate_struct *,
                       nflux_struct *, ntemp_struct *);
void            EvergreenPhenology(const epconst_struct *, epvar_struct *,
                                   const cstate_struct *);
void            FirstDay(elem_struct *, river_struct *, const cninit_struct *);
void            FRootLitFall(const epconst_struct *, realtype, cflux_struct *,
                             nflux_struct *);
void            FreeEpctbl(epctbl_struct *);
realtype          GetCO2(tsdata_struct *, int);
realtype          GetNdep(tsdata_struct *, int);
void            GrowthResp(const epconst_struct *, cflux_struct *);
void            InitBgc(elem_struct *, const epctbl_struct *,
                        const calib_struct *);
void            InitBgcVar(elem_struct *, river_struct *, N_Vector);
void            LeafLitFall(const epconst_struct *, realtype, cflux_struct *,
                            nflux_struct *);
void            LivewoodTurnover(const epconst_struct *, epvar_struct *,
                                 const cstate_struct *, cflux_struct *, const nstate_struct *,
                                 nflux_struct *);
void            MaintResp(const epconst_struct *, epvar_struct *,
                          const cstate_struct *, cflux_struct *, const nstate_struct *,
                          const daily_struct *);
void            MakeZeroFluxStruct(cflux_struct *, nflux_struct *);
void            Mortality(const epconst_struct *, cstate_struct *,
                          cflux_struct *, nstate_struct *, nflux_struct *);
# if defined(_LEACHING_)
void            NLeaching(elem_struct *);
# elif defined(_LUMPED_)
void            NLeachingLumped(elem_struct *, river_struct *);
# else
void            NTransport(elem_struct *, river_struct *);
# endif
void            OffsetLitterfall(const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, nflux_struct *);
void            OnsetGrowth(const epconst_struct *, const epvar_struct *,
    const cstate_struct *, cflux_struct *, const nstate_struct *,
    nflux_struct *);
void            Phenology(const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, const nstate_struct *,
    nflux_struct *, const daily_struct *);
void            Photosynthesis(psn_struct *);
void            PrecisionControl(cstate_struct *cs, nstate_struct *ns);
void            RadTrans(const cstate_struct *, eflux_struct *,
    pstate_struct *, const epconst_struct *, epvar_struct *,
    const daily_struct *);
void            ReadAnnFile(tsdata_struct *, const char *);
void            ReadBgc(const char *, ctrl_struct *, co2control_struct *,
    ndepcontrol_struct *, cninit_struct *, char *, char *);
void            ReadBgcIc(const char *, elem_struct *, river_struct *);
void            ReadEpc(epctbl_struct *);
void            ResetSpinupStat(elem_struct *);
void            RestartInput(cstate_struct *, nstate_struct *,
    epvar_struct *, const bgcic_struct *);
void            RestartOutput(const cstate_struct *, const nstate_struct *,
    const epvar_struct *, bgcic_struct *);
void            SeasonDecidPhenology(const epconst_struct *, epvar_struct *,
    const daily_struct *);
void            SoilPsi(const soil_struct *, realtype, realtype *);
void            TotalPhotosynthesis(const epconst_struct *, epvar_struct *,
    const pstate_struct *, cflux_struct *, psn_struct *, psn_struct *,
    const daily_struct *);
void            WriteBgcIc(const char *, elem_struct *, river_struct *);
void            ZeroSrcSnk(cstate_struct *, nstate_struct *, summary_struct *,
    solute_struct *);
#endif

#if defined(_CYCLES_)
void            AddCrop(crop_struct *);
realtype          Aeration(realtype);
realtype          AirMolarDensity(realtype, realtype);
void            ApplyFertilizer(const fixfert_struct *, cstate_struct *,
    nstate_struct *, nflux_struct *);
realtype          AvgSolConc(int, realtype, const realtype [],
    const realtype [], const realtype [], realtype, const realtype []);
realtype          CNdestiny(realtype, realtype );
void            CalcLatNFlux(realtype, int, const realtype[], const realtype [],
    const realtype[], const realtype [], const realtype[], realtype, realtype []);
void            CalSnkSrc(const nflux_struct *, int, solute_struct *,
    solute_struct *);
void            CalcRootFraction(const crop_struct *, const pstate_struct *,
    realtype *);
realtype          ColdDamage(realtype, realtype, realtype);
realtype          CommRadIntcp(const crop_struct[]);
realtype          CommTotRadIntcp(const crop_struct[]);
void            ComputeColdDamage(const daily_struct *, crop_struct *,
    wstate_struct *, cstate_struct *, nstate_struct *);
void            ComputeFactorComposite(const soil_struct *,
    const daily_struct *, pstate_struct *);
realtype          ComputeHarvestIndex(realtype, realtype, realtype, realtype, realtype);
void            ComputeResidueCover(const cstate_struct *, pstate_struct *);
void            ComputeSoilCarbonBalanceMB(const realtype[], const soil_struct *,
    pstate_struct *, wstate_struct *, cstate_struct *, cflux_struct *,
    nstate_struct *, nflux_struct *);
realtype          ComputeTextureFactor(realtype);
void            ComputeTillageFactor(const tillage_struct *, const realtype [],
    realtype, const soil_struct *, const pstate_struct *, realtype []);
void            CropGrowth(realtype, const daily_struct *, realtype *,
    crop_struct *);
void            CropNitrogenConcentration(realtype, const crop_struct *,
    realtype *, realtype *, realtype *, realtype *, realtype *, realtype *, realtype *);
void            CropStage(int, crop_struct []);
void            CropNitrogenDemand(realtype, realtype, const crop_struct *,
    realtype *, realtype *, realtype *, realtype *);
void            CropNitrogenStress(realtype, realtype, realtype, crop_struct *);
void            CropNitrogenUptake(int, realtype, realtype, const realtype [],
    const realtype [], const realtype [], const realtype [], const realtype [],
    const pstate_struct *, realtype [], realtype [], const realtype [],
    crop_struct [], nstate_struct *, nflux_struct *);
void            DailyCycles(int, pihm_struct);
void            DailyOperations(int, const opertbl_struct *,
    const daily_struct *, soil_struct *, mgmt_struct *, crop_struct [],
    pstate_struct *, wstate_struct *, wflux_struct *, cstate_struct *,
    cflux_struct *, nstate_struct *, nflux_struct *);
void            Denitrification(const soil_struct *, const daily_struct *,
    const pstate_struct *, const cflux_struct *, nstate_struct *,
    nflux_struct *);
void            DistributeRootDetritus(realtype, realtype, realtype,
    realtype, const crop_struct *, const pstate_struct *, cstate_struct *,
    nstate_struct *);
int             Doy(int, int, int);
void            Doy2Date(int, int, int *, int *);
void            ExecuteTillage(const tillage_struct *, const pstate_struct *,
    realtype *, soil_struct *, wstate_struct *, cstate_struct *,
    nstate_struct *, nflux_struct *);
void            FieldOperation(int, const opertbl_struct *,
    const daily_struct *, soil_struct *, mgmt_struct *,
    crop_struct [], pstate_struct *, wstate_struct *, wflux_struct *,
    cstate_struct *, nstate_struct *, nflux_struct *);
int             FinalHarvestDate(int, realtype, realtype, realtype);
int             FindCrop(const char[], const epconst_struct []);
realtype          FindIrrigationVolume(int, realtype, const soil_struct *,
    const daily_struct *daily, const pstate_struct *, const wflux_struct *);
void            FirstDay(const soiltbl_struct *, elem_struct [], river_struct []);  
void            FirstDOY(int, int *);
void            ForageAndSeedHarvest(int, crop_struct *,
    pstate_struct *, wstate_struct *, cstate_struct *, nstate_struct *,
    nflux_struct *);
int             ForcedClipping(int, const crop_struct []);
realtype          Fraction(realtype, realtype, realtype, realtype, realtype);
void            GrainHarvest(int, crop_struct *, pstate_struct *,
    wstate_struct *, cstate_struct *, nstate_struct *);
void            GrowingCrop(int, const soil_struct *, const daily_struct *,
    mgmt_struct *, crop_struct [], pstate_struct *, wstate_struct *,
    cstate_struct *, nstate_struct *, nflux_struct *);
void            HarvestCrop(int, crop_struct *, pstate_struct *,
    wstate_struct *, cstate_struct *, nstate_struct *);
void            InitCropSV(crop_struct *);
void            InitCycles(const agtbl_struct *, const soiltbl_struct *,
    epconst_struct [], elem_struct [], river_struct []);
void            InitCyclesVar(elem_struct [], river_struct [], N_Vector);
int             IsLeapYear(int);
int             IsOperationToday(int, int, const void *, int, int, int *);
void            KillCrop(crop_struct *);
realtype          LinearEquilibriumConcentration(realtype, realtype, realtype, realtype, realtype);
void            MakeZeroFluxStruct(wflux_struct *, cflux_struct *, nflux_struct *);  
realtype          MaximumAbgdHumificationFactor(realtype);
realtype          MaximumManuHumificationFactor(realtype);
realtype          MaximumRhizHumificationFactor(realtype);
realtype          MaximumRootHumificationFactor(realtype);
realtype          Moisture(realtype);
realtype          N2OFractionNitrification(realtype);
void            Nitrification(const soil_struct *, const daily_struct *,
    const pstate_struct *, nstate_struct *, nflux_struct *);
realtype          NitrogenMineralization(realtype, realtype, realtype, realtype);
void            NitrogenTransformation(const soil_struct *,
    const daily_struct *, const crop_struct [], const cstate_struct *,
    const cflux_struct *, pstate_struct *, nstate_struct *, nflux_struct *);
void            NTransport(realtype, elem_struct [], river_struct []);
int             NumActiveCrop(const crop_struct []);
void            Phenology(const daily_struct *, crop_struct []);
void            PlantingCrop(const plant_struct *,  crop_struct *);
void            PotentialSoluteUptake(realtype, int, const realtype[],
    const realtype[], const realtype[], const realtype[], const realtype[], realtype *,
    realtype[]);
void            Processes(int, const soil_struct *, const daily_struct *,
    const pstate_struct *, crop_struct [], cstate_struct *, nstate_struct *,
    nflux_struct *);
void            RadiationInterception(crop_struct []);
void            ReadCrop(const char [], epconst_struct []);
void            ReadCyclesCtrl(const char [], agtbl_struct *, ctrl_struct *);
void            ReadMultOper(const agtbl_struct *, const epconst_struct [],
    opertbl_struct []);
void            ReadOperation(const char [], int, const epconst_struct [],
    opertbl_struct []);
void            ReadSoilInit(const char [], soiltbl_struct *);
void            ResidueEvaporation(realtype, realtype, realtype, const crop_struct [],
    const pstate_struct *, const cstate_struct *, wstate_struct *,
    wflux_struct *);
void            ResidueWetting(const pstate_struct *, const cstate_struct *,
    realtype, wstate_struct *, wflux_struct *);
void            RestartInput(const cyclesic_struct *, pstate_struct *,
    wstate_struct *, cstate_struct *, nstate_struct *);
realtype          ShootBiomassPartitioning(realtype, realtype, realtype, int);
void            SoluteTransport(int, realtype, realtype, const realtype [],
    const realtype [], const realtype [], const realtype [], realtype []);
realtype          TemperatureFunction(realtype);
realtype          TemperatureFunctionGrowth(realtype, realtype, realtype, realtype);
realtype          TemperatureLimitation(realtype, realtype, realtype);
realtype          ThermalTime(realtype, realtype, realtype, realtype);
void            TillageFactorSettling(int, const realtype [], realtype, realtype []);
void            UpdNProf(realtype, const soil_struct *, const wstate_struct *,
    const nstate_struct *, const nflux_struct *, const nprof_struct *,
    pstate_struct *, nstate_struct *);
void            Volatilization(const soil_struct *, const daily_struct *,
    const crop_struct [], const pstate_struct *, const cstate_struct *,
    nstate_struct *, nflux_struct *);
realtype          VolatilizationDepthFunction(realtype);
void            WaterUptake(const soil_struct *, const estate_struct *,
    const pstate_struct *, realtype, crop_struct [], wstate_struct *,
    wflux_struct *);
#endif

#endif

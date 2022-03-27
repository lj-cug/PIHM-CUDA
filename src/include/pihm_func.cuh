#ifndef PIHM_FUNC_HEADER_D
#define PIHM_FUNC_HEADER_D

// Functions In PIHM_vs
void set_nvector_cuda(N_Vector, realtype  *, int);
void get_nvector_cuda(N_Vector, realtype  *, int);
#if defined(_NOAH_)
void Initialize(pihm_struct, pihm_struct_d, N_Vector);
#else
void Initialize(pihm_struct, N_Vector);
#endif
cudaError_t Initialize_cuda(pihm_struct_d);
void AoS_to_SoA(pihm_struct, pihm_struct_d);
cudaError_t FreeMem_cuda(pihm_struct_d);
void   PIHM(pihm_struct, pihm_struct_d, void *, N_Vector, N_Vector, realtype);
void   Spinup(pihm_struct, pihm_struct_d, N_Vector, N_Vector, void *);  
void   ApplyBc(forc_struct *, elem_struct *, river_struct *, pihm_struct_d, int);
void   ApplyElemBc(forc_struct *, elem_struct *, pihm_struct_d, int);
void   ApplyRiverBc(forc_struct *, river_struct *, pihm_struct_d, int);
void   ApplyLai(forc_struct *, elem_struct *, pihm_struct_d, int);
#if defined(_NOAH_)
void   ApplyForc(forc_struct *, elem_struct *, pihm_struct_d, int, int,
                 const siteinfo_struct *);
void   ApplyMeteoForc(forc_struct *, elem_struct *, pihm_struct_d, int, int,
                 const siteinfo_struct *);
#else
void   ApplyForc(forc_struct *, elem_struct *, pihm_struct_d, int);
void   ApplyMeteoForc(forc_struct *, elem_struct *, pihm_struct_d, int); 
#endif

void   PrintWaterBal(FILE *, int, int, int, const pihm_struct_d);                              
void   PrintInit(const pihm_struct_d, const char *, int, int, int, int);
                 
void IntcpSnowEt(int, int, const pihm_struct_d);
__global__ void IntcpSnowEt_kernel(int, realtype, int,
								   realtype, realtype, realtype,
#include "IntcpSnowEt_kernel_declarations.c"
);

void SetCVodeParam(pihm_struct, pihm_struct_d, void *, N_Vector);
static int ODE(realtype , N_Vector, N_Vector, void *);
__global__ static void init_dy(realtype *, const int);	
__global__ static void Y_update_elem(const  realtype *,
                               const int, 
	                           realtype *, realtype *, realtype *);
__global__ static void Y_update_river(const  realtype *,
								const int, const int,
	                           realtype *, realtype *, realtype *, realtype *);
	
__global__ static void Hydrol_elem(const int, const realtype *, const realtype *, 
	                               const int, const int,
#include "Hydrol_elem_device_declarations.c"
);	
	
__global__ static void build_rhs_elem(const int, 
							   realtype  *,                            
#include "build_rhs_elem_declarations.c"  							
);
	
__global__ static void build_rhs_river(realtype  *,
                              const int, const int,
#include "build_rhs_river_declarations.c"                       
);	
								

/* Determine which layers does ET extract water from */ 
__device__ void EtExtract(int, 
#include "EtExtract_device_declarations.c" 
);	
__device__ realtype SurfH(realtype);

__device__ void LateralFlow(int, int, int, 
                            const realtype  *, const realtype  *,
#include "LateralFlow_river_device_declarations.c" 
);

__device__ void VerticalFlow(realtype, int, 
#include "VerticalFlow_device_declarations.c" 	
);

/* Lateral flow between element and river*/
__global__ void FrictSlope(int, realtype  *, realtype  *, 
#include "FrictSlope_device_declarations.c" 
);	                      
__device__ realtype AvgHsurf(realtype, realtype, realtype);
__device__ realtype AvgH(realtype, realtype, realtype );
__device__ realtype DhByDl(const realtype  *, 
	                       const realtype  *, const realtype *);
__device__ realtype EffKh(int, realtype,
#include "EffKh_device_declarations.c" 
);

__device__ realtype OverLandFlow(realtype, realtype, realtype, realtype, realtype);	 
                          
__device__ realtype SubFlowElemToElem0(int, int, 
#include "SubFlowElemToElem0_device_declarations.c" 
);
__device__ realtype SubFlowElemToElem1(int, int, 
#include "SubFlowElemToElem1_device_declarations.c" 
);
__device__ realtype SubFlowElemToElem2(int, int, 
#include "SubFlowElemToElem2_device_declarations.c" 
);

__device__ realtype OvlFlowElemToElem0(int, int, realtype, int, 
#include "OvlFlowElemToElem0_device_declarations.c" 
);	 
 __device__ realtype OvlFlowElemToElem1(int, int, realtype, int, 
#include "OvlFlowElemToElem1_device_declarations.c" 
);	
__device__ realtype OvlFlowElemToElem2(int, int, realtype, int, 
#include "OvlFlowElemToElem2_device_declarations.c" 
);	
                              
__device__ void BoundFluxElem0(int, int, 
#include "BoundFluxElem0_device_declarations.c" 
);	
__device__ void BoundFluxElem1(int, int, 
#include "BoundFluxElem1_device_declarations.c" 
);	
__device__ void BoundFluxElem2(int, int, 
#include "BoundFluxElem2_device_declarations.c" 
);	
					 
/* Vertical flow in elements */			 
__device__ realtype Infil(realtype, int,
#include "Infil_device_declarations.c" 
);	                    
__device__ realtype KrFunc_cuda(realtype, realtype );
__device__ realtype EffKinf(realtype, realtype, realtype, realtype, 
                            realtype, int,
#include "EffKinf_device_declarations.c" 
);	                      
__device__ realtype Psi(realtype, realtype, realtype );
__device__ realtype Recharge(int,
#include "Recharge_device_declarations.c" 
);	                       
__device__ realtype AvgKv(realtype, realtype, realtype, int,
#include "AvgKv_device_declarations.c" 
);  
                  	
__global__ void RiverFlow(const int, const int,
#include "RiverFlow_device_declarations.c" 
);	
			
__device__ realtype BoundFluxRiver(int, int, 
#include "BoundFluxRiver_device_declarations.c" 
);	
					
__device__ realtype ChanFlowRiverToRiver(int, int, int,
#include "ChanFlowRiverToRiver_device_declarations.c" 
);

__device__ realtype SubFlowRiverToRiver(realtype, realtype, int, int,
#include "SubFlowRiverToRiver_device_declarations.c" 
);	
                                 
__device__ realtype OutletFlux(int, int,
#include "OutletFlux_device_declarations.c" 
);

__device__ void RiverToElem(int, int, int,
#include "RiverToElem_device_declarations.c" 
);

__device__ realtype OvlFlowElemToRiver(int, int,
#include "OvlFlowElemToRiver_device_declarations.c" 
);

__device__ realtype ChanFlowElemToRiver(int, realtype, int, realtype, 
#include "ChanFlowElemToRiver_device_declarations.c"                                  
);	                                  
__device__ realtype SubFlowElemToRiver(int, realtype, int, realtype, realtype, 
#include "SubFlowElemToRiver_device_declarations.c"                              								
);     
								
__device__ realtype ChanLeak(int,
#include "ChanLeak_device_declarations.c" 
);
__device__ realtype RiverCroSectArea(int, realtype, realtype );
__device__ realtype RiverPerim(int, realtype, realtype );

#if defined(_NOAH_)
/* Noah Model ported into CUDA */
void Noah(pihm_struct_d, realtype);
// Noah_kernel fused CalHum_kernel and  SFlx_kernel.
__global__ void Noah_kernel(int, realtype,
#include "Noah_kernel_declarations.c" 
);
void NoahHydrol(pihm_struct_d, realtype);
__global__ void NoahHydrol_kernel(int, realtype,
#include "NoahHydrol_kernel_declarations.c" 
);
__device__ void AdjSmProf(int, const realtype *, realtype, 
#include "AdjSmProf_kernel_declarations.c"
);
__device__ void AlCalc(int, realtype, int,  
#include "AlCalc_kernel_declarations.c" 
);
__device__ void CanRes(int,
#include "CanRes_kernel_declarations.c"
);
__device__ realtype CSnow(realtype);
__device__ void DEvap(int,
#include "DEvap_kernel_declarations.c"
);
__device__ void  Evapo(int, realtype,
#include "Evapo_kernel_declarations.c"
);
__device__  int FindWaterTable(realtype *, int, realtype, realtype *);
__device__ double FrH2O(int, double, double, double, 
	realtype *,realtype *,realtype *,realtype *);

__device__ realtype FrozRain(realtype, realtype);
__device__ void  HRT(int, realtype *,
realtype, realtype, realtype, realtype, realtype *, realtype *, realtype *,
#include "HRT_kernel_declarations.c"
);
__device__ void HStep(int, double *, double, double *ai,
	double *, double *, 
	realtype *, realtype *, realtype *, realtype *, realtype *, realtype *,
	realtype *, realtype *, realtype *, realtype *, realtype *);
__device__ void  NoPac(int, realtype, realtype,
#include "NoPac_kernel_declarations.c"
);
__device__ void PcpDrp(int, realtype, realtype,
#include "PcpDrp_kernel_declarations.c"	
);
__device__ void Penman(int, realtype *, realtype, int, int,
#include "Penman_kernel_declarations.c"	
);
__device__ realtype Pslhs(realtype);
__device__ realtype Pslhu(realtype);
__device__ realtype Pslms(realtype);
__device__ realtype Pslmu(realtype);
__device__ realtype Psphs(realtype);
__device__ realtype Psphu(realtype);
__device__ realtype Pspms(realtype);
__device__ realtype Pspmu(realtype);

__device__ void Rosr12(realtype *, const realtype *, const realtype *, realtype *,
            const realtype *, realtype *);		
					   
__device__ void SfcDifOff(int, realtype, realtype, int,
#include "SfcDifOff_kernel_declarations.c"	
);
__device__ void  ShFlx(int, realtype, realtype, realtype, realtype,
#include "ShFlx_kernel_declarations.c"
);
__device__ void  SmFlx(int, realtype,
#include "SmFlx_kernel_declarations.c"
);
__device__ realtype SnFrac(realtype, realtype, realtype);
__device__ void SnowNew(int, realtype, realtype *,realtype *,realtype *);
__device__ void SnkSrc(int, realtype *, realtype, realtype, realtype *,
            realtype, realtype, realtype, realtype, realtype, realtype,
			realtype, realtype, realtype, realtype, realtype,  // MAXLYR=11
			realtype, int, realtype,
	        realtype *,realtype *,realtype *,realtype *);
			
__device__ void SnoPac(int tid, int snowng, double dt, 
	                   double t24, double prcpf, double df1,
#include "SnoPac_kernel_declarations.c"
);
__device__ void SnowPack(realtype, realtype, realtype *, realtype *, realtype, realtype);			
__device__ realtype Snowz0(realtype, realtype, realtype);
__device__ void SRT(int, realtype *, realtype *, realtype *, realtype *, realtype *,
#include "SRT_kernel_declarations.c"
);
__device__ void SStep(int, realtype *, realtype *, realtype *, realtype *, realtype *,
           realtype,
#include "SStep_kernel_declarations.c"		   
);
__device__ void Transp(int,
#include "Transp_kernel_declarations.c"
);
__device__ realtype TBnd(realtype, realtype, realtype, int,
realtype, realtype, realtype, realtype, realtype,
realtype, realtype, realtype, realtype, realtype, realtype,
);
__device__ realtype TDfCnd(realtype, realtype, realtype, realtype, realtype);
__device__ realtype TmpAvg(realtype, realtype, realtype, int,
	realtype,realtype,realtype,realtype,realtype,realtype,   // elem_d_ps_zsoil*[tid]
	realtype,realtype,realtype,realtype,realtype);

__device__ void WDfCnd(int, realtype *, realtype *, realtype, realtype,
#include "WDfCnd_kernel_declarations.c"
);
#endif

void Summary(N_Vector, void *, int); 

__global__ void update_elem(realtype *, realtype, int, realtype *,
#include "update_elem_device_declarations.c"         
); 

__global__ void update_river(realtype *, realtype, int, int, 
							 realtype *, realtype *);
							 
__device__ void MassBalance(realtype *, realtype, int,
#include "mass_balance_declarations.c"           
);

// CUDA intrinsic functions (single precision)
/*
__device__ realtype __powf(realtype x, realtype y);
__device__ realtype __logf(realtype x);
__device__ realtype __fsqrt_rn(realtype x);
__device__ realtype __expf(realtype x);
__device__ realtype __sinf(realtype x);
__device__ realtype __cosf(realtype x);
*/

// CUDA intrinsic functions (double precision)
__device__ realtype pow(realtype x, realtype y);
__device__ realtype log(realtype x);
__device__ realtype sqrt(realtype x);
__device__ realtype exp(realtype x);
__device__ realtype sin(realtype x);
__device__ realtype cos(realtype x);
#endif
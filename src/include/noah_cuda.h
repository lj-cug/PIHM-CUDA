#ifndef NOAH_STRUCT_HEADER_D
#define NOAH_STRUCT_HEADER_D

/* ָ������, ������ʱ�洢SOA-style��Noah LSM���� */
struct noah_1d_device
{
	realtype *elem_d_ps_rtdis[MAXLYR];
	realtype *elem_d_ps_sldpth[MAXLYR];
	realtype *elem_d_ps_zsoil[MAXLYR];
	realtype *elem_d_ps_satdpth[MAXLYR];
	realtype *elem_d_ws_smc[MAXLYR];
	realtype *elem_d_ws_sh2o[MAXLYR];
	realtype *elem_d_wf_et[MAXLYR];
	realtype *elem_d_wf_runoff2_lyr[MAXLYR];
	realtype *elem_d_wf_smflxv[MAXLYR];
	realtype *elem_d_es_stc[MAXLYR];
	realtype *elem_d_ef_et[MAXLYR];
};
typedef noah_1d_device *noah_1d; 

/* ����ָ��, ����noah LSM����nsoil���йص�2D���� */
struct noah_2d_device
{
	realtype **elem_d_ps_rtdis;
	realtype **elem_d_ps_sldpth;
	realtype **elem_d_ps_zsoil;
	realtype **elem_d_ps_satdpth;
	realtype **elem_d_ws_smc;
	realtype **elem_d_ws_sh2o;
	realtype **elem_d_wf_et;
	realtype **elem_d_wf_runoff2_lyr;
	realtype **elem_d_wf_smflxv;
	realtype **elem_d_es_stc;
	realtype **elem_d_ef_et;
};
typedef noah_2d_device *noah_2d; 

#endif
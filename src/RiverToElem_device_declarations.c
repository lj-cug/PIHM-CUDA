// Totally 40 arguments
// OvlFlowElemToRiver_device_declarations
realtype *elem_d_topo_zmax,
realtype *elem_d_ws_surfh,
realtype *river_d_topo_zmax,
realtype *river_d_topo_zbed,
realtype *river_d_matl_cwr,
realtype *river_d_shp_length,
realtype *river_d_ws_stage,

// EffKh
realtype *elem_d_soil_areafv,
realtype *elem_d_soil_depth,
realtype *elem_d_soil_dmac,
realtype *elem_d_soil_kmach,
realtype *elem_d_soil_ksath,

// ChanFlowElemToRiver
realtype *elem_d_topo_zmin,
realtype *elem_d_ws_gw,
realtype *river_d_matl_ksath,

// SubFlowElemToRiver
realtype *river_d_topo_zmin,
realtype *river_d_ws_gw,

// Other array
int *elem_d_nabr0,
int *elem_d_nabr1,
int *elem_d_nabr2,
realtype *elem_d_wf_ovlflow0,
realtype *elem_d_wf_ovlflow1,
realtype *elem_d_wf_ovlflow2,
realtype *elem_d_wf_subsurf0,
realtype *elem_d_wf_subsurf1,
realtype *elem_d_wf_subsurf2,
int *river_d_ind, 	
int *river_d_leftele, 
int *river_d_rightele, 
realtype *river_d_topo_dist_left, 
realtype *river_d_topo_dist_right, 
realtype *river_d_wf_rivflow0,
realtype *river_d_wf_rivflow1,
realtype *river_d_wf_rivflow2,
realtype *river_d_wf_rivflow3,
realtype *river_d_wf_rivflow4,
realtype *river_d_wf_rivflow5,
realtype *river_d_wf_rivflow6,
realtype *river_d_wf_rivflow7,
realtype *river_d_wf_rivflow8
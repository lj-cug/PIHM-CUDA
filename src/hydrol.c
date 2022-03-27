#include "pihm.h"

void Hydrol(elem_struct *elem, river_struct *river, const ctrl_struct *ctrl)
{
    int             i;

//
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        /* Calculate actual surface water depth 
		   计算实际的地表产流深 */
        elem[i].ws.surfh = SurfH(elem[i].ws.surf);
    }

    /* Determine which layers does ET extract water from 
	   确定蒸散发损失的水量是从哪一层土壤扣除
	*/
    EtExtract(elem);

    /* Water flow */
    LateralFlow(elem, river, ctrl->surf_mode);   // 侧向上，单元到单元的流动，坡面流动力学

    VerticalFlow(elem, (double)ctrl->stepsize);  // 垂向上，土壤层之间的流动，渗流

    RiverFlow(elem, river, ctrl->riv_mode);      // 河道分段单元与相邻网格单元之间的流动

}


void EtExtract(elem_struct *elem)
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        /* Source of direct evaporation */
#if defined(_NOAH_)
        if (elem[i].ws.gw > elem[i].soil.depth - elem[i].soil.dinf)
        {
            elem[i].wf.edir_surf = 0.0;
            elem[i].wf.edir_unsat = 0.0;
            elem[i].wf.edir_gw = elem[i].wf.edir;
        }
        else
        {
            elem[i].wf.edir_surf = 0.0;
            elem[i].wf.edir_unsat = elem[i].wf.edir;
            elem[i].wf.edir_gw = 0.0;
        }

#else
        if (elem[i].ws.surfh >= DEPRSTG)
        {
            elem[i].wf.edir_surf = elem[i].wf.edir;
            elem[i].wf.edir_unsat = 0.0;
            elem[i].wf.edir_gw = 0.0;
        }
        else if (elem[i].ws.gw > elem[i].soil.depth - elem[i].soil.dinf)
        {
            elem[i].wf.edir_surf = 0.0;
            elem[i].wf.edir_unsat = 0.0;
            elem[i].wf.edir_gw = elem[i].wf.edir;
        }
        else
        {
            elem[i].wf.edir_surf = 0.0;
            elem[i].wf.edir_unsat = elem[i].wf.edir;
            elem[i].wf.edir_gw = 0.0;
        }
#endif

        /* Source of transpiration */
#if defined(_NOAH_)
        elem[i].ps.gwet = GwTransp(elem[i].wf.ett, elem[i].wf.et,
                                   elem[i].ps.nwtbl, elem[i].ps.nroot);
        elem[i].wf.ett_unsat = (1.0 - elem[i].ps.gwet) * elem[i].wf.ett;
        elem[i].wf.ett_gw = elem[i].ps.gwet * elem[i].wf.ett;
#else
        if (elem[i].ws.gw > elem[i].soil.depth - elem[i].ps.rzd)
        {
            elem[i].wf.ett_unsat = 0.0;
            elem[i].wf.ett_gw = elem[i].wf.ett;
        }
        else
        {
            elem[i].wf.ett_unsat = elem[i].wf.ett;
            elem[i].wf.ett_gw = 0.0;
        }
#endif
    }
}


double SurfH(double surfeqv)
{
    /*
     * Following Panday and Huyakorn (2004) AWR:
     * Use a parabolic curve to express the equivalent surface water depth
     * (surfeqv) in terms of actual flow depth (surfh) when the actual flow
     * depth is below depression storage; assume that
     * d(surfeqv) / d(surfh) = 1.0 when surfh = DEPRSTG. Thus
     *   surfeqv = (1 / 2 * DEPRSTG) * surfh ^ 2, i.e.
     *   surfh = sqrt(2 * DEPRSTG * surfeqv)
     */
    double          surfh;

    if (DEPRSTG == 0.0)
    {
        return surfeqv;
    }
    else
    {
        if (surfeqv < 0.0)
        {
            surfh = 0.0;
        }
        else if (surfeqv <= 0.5 * DEPRSTG)
        {
            surfh = sqrt(2.0 * DEPRSTG * surfeqv);
        }
        else
        {
            surfh = DEPRSTG + (surfeqv - 0.5 * DEPRSTG);
        }

        return surfh;
    }
}


#if defined(_XAJ_)
double XAJ_Model()
{
	double          surfh;
	/* 新安江概念式水文模型中的产流模式 
	IN:    elem[i].ws.surfh 

	       elem[i].wf.edir_surf
	       elem[i].wf.edir_unsat
		   elem[i].wf.edir_gw

		   elem[i].wf.ett_unsat   // Source of transpiration 
		   elem[i].wf.ett_gw
	OUT:
		   elem_struct *elem
	*/






	return surfh;
}
#endif
#include "pihm.h"

void InitGeol (elem_struct *elem, const geoltbl_struct *geoltbl,
        const calib_struct *cal)
{
    int             i;
    int             geol_ind;

    for (i = 0; i < nelem; i++)
    {
        geol_ind = elem[i].attrib.geol_type - 1;

        elem[i].geol.depth = elem[i].topo.zmin - elem[i].topo.zbed;

        elem[i].geol.ksath = cal->ksath * geoltbl->ksath[geol_ind];
        elem[i].geol.ksatv = cal->ksatv * geoltbl->ksatv[geol_ind];

        elem[i].geol.smcmin = cal->porosity * geoltbl->smcmin[geol_ind];
        elem[i].geol.smcmax = cal->porosity * geoltbl->smcmax[geol_ind];
        elem[i].geol.porosity = elem[i].geol.smcmax - elem[i].geol.smcmin;
        if (elem[i].geol.porosity > 1.0 || elem[i].geol.porosity <= 0.0)
        {
            PIHMprintf (VL_ERROR,
                "Error: Porosity value out of bounds for Element %d", i + 1);
            PIHMexit (EXIT_FAILURE);
        }
        elem[i].geol.alpha = cal->alpha * geoltbl->alpha[geol_ind];
        elem[i].geol.beta = cal->beta * geoltbl->beta[geol_ind];
    }
}

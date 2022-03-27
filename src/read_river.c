#include "pihm.h"

void ReadRiver(const char *filename, rivtbl_struct *rivtbl,
    shptbl_struct *shptbl, matltbl_struct *matltbl, forc_struct *forc)
{
    int             i, j;
    FILE           *riv_file;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    riv_file = fopen(filename, "r");
    CheckFile(riv_file, filename);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

    /*
     * Read river segment block
     */
    /* Read number of river segments */
    FindLine(riv_file, "BOF", &lno, filename);
    NextLine(riv_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUMRIV", &nriver, 'i', filename, lno);

    /* Allocate */
    rivtbl->fromnode = (int *)malloc(nriver * sizeof(int));
    rivtbl->tonode = (int *)malloc(nriver * sizeof(int));
    rivtbl->down = (int *)malloc(nriver * sizeof(int));
    rivtbl->leftele = (int *)malloc(nriver * sizeof(int));
    rivtbl->rightele = (int *)malloc(nriver * sizeof(int));
    rivtbl->shp = (int *)malloc(nriver * sizeof(int));
    rivtbl->matl = (int *)malloc(nriver * sizeof(int));
    rivtbl->bc = (int *)malloc(nriver * sizeof(int));
    rivtbl->rsvr = (int *)malloc(nriver * sizeof(int));

    /* Skip header line */
    NextLine(riv_file, cmdstr, &lno);

    /* Read river segment information */
    for (i = 0; i < nriver; i++)
    {
        NextLine(riv_file, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %d %d %d %d %d %d %d %d %d",
            &index,
            &rivtbl->fromnode[i], &rivtbl->tonode[i], &rivtbl->down[i],
            &rivtbl->leftele[i], &rivtbl->rightele[i], &rivtbl->shp[i],
            &rivtbl->matl[i], &rivtbl->bc[i], &rivtbl->rsvr[i]);
        if (match != 10 || i != index - 1)
        {
            PIHMprintf(VL_ERROR,
                "Error reading river attribute for the %dth segment.\n", i + 1);
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    /*
     * Read river shape information
     */
    NextLine(riv_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "SHAPE", &shptbl->number, 'i', filename, lno);

    /* Allocate */
    shptbl->depth = (realtype *)malloc(shptbl->number * sizeof(realtype));
    shptbl->intrpl_ord = (int *)malloc(shptbl->number * sizeof(int));
    shptbl->coeff = (realtype *)malloc(shptbl->number * sizeof(realtype));

    /* Skip header line */
    NextLine(riv_file, cmdstr, &lno);

    for (i = 0; i < shptbl->number; i++)
    {
        NextLine(riv_file, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %lf %d %lf",
            &index,
            &shptbl->depth[i], &shptbl->intrpl_ord[i], &shptbl->coeff[i]);
        if (match != 4 || i != index - 1)
        {
            PIHMprintf(VL_ERROR,
                "Error reading river shape description for the %dth shape.\n",
                i + 1);
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    /*
     * Read river material information
     */
    NextLine(riv_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MATERIAL", &matltbl->number, 'i', filename, lno);

    /* Allocate */
    matltbl->rough = (realtype *)malloc(matltbl->number * sizeof(realtype));
    matltbl->cwr = (realtype *)malloc(matltbl->number * sizeof(realtype));
    matltbl->ksath = (realtype *)malloc(matltbl->number * sizeof(realtype));
    matltbl->ksatv = (realtype *)malloc(matltbl->number * sizeof(realtype));
    matltbl->bedthick = (realtype *)malloc(matltbl->number * sizeof(realtype));

    /* Skip header line */
    NextLine(riv_file, cmdstr, &lno);

    for (i = 0; i < matltbl->number; i++)
    {
        NextLine(riv_file, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %lf %lf %lf %lf %lf",
            &index, &matltbl->rough[i], &matltbl->cwr[i],
            &matltbl->ksath[i], &matltbl->ksatv[i], &matltbl->bedthick[i]);
        if (match != 6 || i != index - 1)
        {
            PIHMprintf(VL_ERROR,
                "Error reading description of the %dth material.\n", i + 1);
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    /*
     * Read river boundary condition block
     */
    NextLine(riv_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "BC", &forc->nriverbc, 'i', filename, lno);

    if (forc->nriverbc > 0)
    {
        forc->riverbc =
            (tsdata_struct *)malloc(forc->nriverbc * sizeof(tsdata_struct));

        NextLine(riv_file, cmdstr, &lno);
        for (i = 0; i < forc->nriverbc; i++)
        {
            match = sscanf(cmdstr, "%*s %d", &index);
            if (match != 1 || i != index - 1)
            {
                PIHMprintf(VL_ERROR,
                    "Error reading description "
                    "of the %dth river boundary condition.\n", i);
                PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n",
                    filename, lno);
                PIHMexit(EXIT_FAILURE);
            }
            NextLine(riv_file, cmdstr, &lno);
            NextLine(riv_file, cmdstr, &lno);
            forc->riverbc[i].length =
                CountLine(riv_file, cmdstr, 2, "RIV_TS", "RES");
        }

        FindLine(riv_file, "BOF", &lno, filename);
        FindLine(riv_file, "BC", &lno, filename);

        for (i = 0; i < forc->nriverbc; i++)
        {
            NextLine(riv_file, cmdstr, &lno);
            NextLine(riv_file, cmdstr, &lno);
            NextLine(riv_file, cmdstr, &lno);

            forc->riverbc[i].data =
                (realtype **)malloc((forc->riverbc[i].length) * sizeof(realtype *));
            forc->riverbc[i].ftime =
                (int *)malloc((forc->riverbc[i].length) * sizeof(int));
            for (j = 0; j < forc->riverbc[i].length; j++)
            {
                forc->riverbc[i].data[j] = (realtype *)malloc(sizeof(realtype));
                NextLine(riv_file, cmdstr, &lno);
                if (!ReadTS(cmdstr, &forc->riverbc[i].ftime[j],
                        &forc->riverbc[i].data[j][0], 1))
                {
                    PIHMprintf(VL_ERROR,
                        "Error reading river boundary condition.\n");
                    PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n",
                        filename, lno);
                    PIHMexit(EXIT_FAILURE);
                }
            }
        }
    }

#if NOT_YET_IMPLEMENTED
    /* Read Reservoir information */
#endif

    fclose(riv_file);
}

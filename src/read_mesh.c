#include "pihm.h"

void ReadMesh(const char *filename, meshtbl_struct *meshtbl)
{
    FILE           *mesh_file;
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    mesh_file = fopen(filename, "r");
    CheckFile(mesh_file, filename);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

    /*
     * Read element mesh block
     */
    NextLine(mesh_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUMELE", &nelem, 'i', filename, lno);

    meshtbl->node = (int **)malloc(nelem * sizeof(int *));
    meshtbl->nabr = (int **)malloc(nelem * sizeof(int *));

    /* Skip header line */
    NextLine(mesh_file, cmdstr, &lno);

    for (i = 0; i < nelem; i++)
    {
        meshtbl->node[i] = (int *)malloc(NUM_EDGE * sizeof(int));
        meshtbl->nabr[i] = (int *)malloc(NUM_EDGE * sizeof(int));

        NextLine(mesh_file, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %d %d %d %d %d %d",
            &index,
            &meshtbl->node[i][0], &meshtbl->node[i][1], &meshtbl->node[i][2],
            &meshtbl->nabr[i][0], &meshtbl->nabr[i][1], &meshtbl->nabr[i][2]);
        if (match != 7 || i != index - 1)
        {
            PIHMprintf(VL_ERROR,
                "Error reading mesh description of the %dth element.\n", i + 1);
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    /*
     * Read node block
     */
    NextLine(mesh_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUMNODE", &meshtbl->numnode, 'i', filename, lno);

    /* Skip header line */
    NextLine(mesh_file, cmdstr, &lno);

    meshtbl->x = (realtype *)malloc(meshtbl->numnode * sizeof(realtype));
    meshtbl->y = (realtype *)malloc(meshtbl->numnode * sizeof(realtype));
    meshtbl->zmin = (realtype *)malloc(meshtbl->numnode * sizeof(realtype));
    meshtbl->zmax = (realtype *)malloc(meshtbl->numnode * sizeof(realtype));

    for (i = 0; i < meshtbl->numnode; i++)
    {
        NextLine(mesh_file, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %lf %lf %lf %lf",
            &index,
            &meshtbl->x[i], &meshtbl->y[i],
            &meshtbl->zmin[i], &meshtbl->zmax[i]);
        if (match != 5 || i != index - 1)
        {
            PIHMprintf(VL_ERROR,
                "Error reading description of the %dth node!\n", i + 1);
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    fclose(mesh_file);
}

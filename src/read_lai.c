#include "pihm.h"

void ReadLai(const char *filename, forc_struct *forc,
    const atttbl_struct *atttbl)
{
    char            cmdstr[MAXSTRING];
    int             read_lai = 0;
    FILE           *lai_file;
    int             i, j;
    int             index;
    int             lno = 0;

    for (i = 0; i < nelem; i++)
    {
        if (atttbl->lai[i] > 0)
        {
            read_lai = 1;
            break;
        }
    }

    forc->nlai = 0;

    if (read_lai)
    {
        lai_file = fopen(filename, "r");
        CheckFile(lai_file, filename);
        PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

        /* Start reading lai_file */
        FindLine(lai_file, "BOF", &lno, filename);

        forc->nlai = CountOccurr(lai_file, "LAI_TS");

        FindLine(lai_file, "BOF", &lno, filename);
        if (forc->nlai > 0)
        {
            forc->lai =
                (tsdata_struct *)malloc(forc->nlai * sizeof(tsdata_struct));

            NextLine(lai_file, cmdstr, &lno);
            for (i = 0; i < forc->nlai; i++)
            {
                ReadKeyword(cmdstr, "LAI_TS", &index, 'i', filename, lno);

                if (i != index - 1)
                {
                    PIHMprintf(VL_ERROR,
                        "Error reading the %dth LAI time series.\n", i + 1);
                    PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n",
                        filename, lno);
                    PIHMexit(EXIT_FAILURE);
                }
                /* Skip header lines */
                NextLine(lai_file, cmdstr, &lno);
                NextLine(lai_file, cmdstr, &lno);
                forc->lai[i].length = CountLine(lai_file, cmdstr, 1, "LAI_TS");
            }

            /* Rewind and read */
            FindLine(lai_file, "BOF", &lno, filename);
            for (i = 0; i < forc->nlai; i++)
            {
                /* Skip header lines */
                NextLine(lai_file, cmdstr, &lno);
                NextLine(lai_file, cmdstr, &lno);
                NextLine(lai_file, cmdstr, &lno);

                forc->lai[i].ftime =
                    (int *)malloc(forc->lai[i].length * sizeof(int));
                forc->lai[i].data =
                    (realtype **)malloc(forc->lai[i].length * sizeof(realtype *));
                for (j = 0; j < forc->lai[i].length; j++)
                {
                    forc->lai[i].data[j] = (realtype *)malloc(sizeof(realtype));
                    NextLine(lai_file, cmdstr, &lno);
                    if (!ReadTS(cmdstr, &forc->lai[i].ftime[j],
                        &forc->lai[i].data[j][0], 1))
                    {
                        PIHMprintf(VL_ERROR, "Error reading LAI forcing.");
                        PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n",
                            filename, lno);
                        PIHMexit(EXIT_FAILURE);
                    }
                }
            }
        }

        fclose(lai_file);
    }
}

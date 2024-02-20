#include <stdlib.h>
#include <metis.h>


__declspec(dllexport) int Partition(int ne, int nn, idx_t *eind, int nparts, idx_t *epart)
{
    idx_t nee = ne;
    idx_t nnn = nn;
    idx_t nnparts = nparts;
    // create eptr array of idx_t with the number of elements + 1
    idx_t *eptr = (idx_t *)malloc((ne + 1) * sizeof(idx_t));
    eptr[0] = 0;
    int j = 0;
    for (int i = 1; i < ne + 1; i++)
    {
        j += 8;
        eptr[i] = j;
    }


    // Create array of idx_t with the number of elements + 1
    idx_t *vwgt = NULL;
    idx_t *vsize = NULL;

    // // define the number of common nodes
    idx_t ncommon = 1;


    

    // This is an array of size nparts that speciﬁes the desired weight for each partition. The target partition
    // weight for the ith partition is speciﬁed at tpwgts[i] (the numbering for the partitions starts from
    // 0). The sum of the tpwgts[] entries must be 1.0.
    real_t *tpwgts = NULL;
    idx_t objval;

    // idx_t options[METIS_NOPTIONS]
    idx_t *options = (idx_t *)malloc(METIS_NOPTIONS * sizeof(idx_t));
    METIS_SetDefaultOptions(options);

    // define the partition vector
    
    // print somthing
    idx_t *npart = (idx_t *)malloc(nn * sizeof(idx_t));
    int res= METIS_PartMeshDual(&nee, &nnn, eptr, eind, vwgt, NULL, &ncommon, &nnparts, tpwgts, options, &objval, epart, npart);
    return res;
}




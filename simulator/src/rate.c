#include "tree.h"
#include "utils.h"
#include <assert.h>

void allocRates(pTreeNode p, int nc, int nd) {
    if (p == NULL) return;
    
    // allocate space for rates
    p->rate_c = (double *)calloc(nc, sizeof(double));
    p->rate_d = (double *)calloc(nd, sizeof(double));
    if (p->rate_c == NULL || p->rate_d == NULL) {
        printf("Failed to allocate rates.\n");
        exit(1);
    }

    allocRates(p->llink, nc, nd);
    allocRates(p->rlink, nc, nd);
}

void simulateRates(pPhyTree tree, int nc, int nd, int npt,
                   double bmv, double base, double rvar) {
    /* simulate branch rates under the given clock model */
    int  i, l, n, m;
    pTreeNode p;
    double mu, sigma, rate;
    
    if (nc <= 0 && nd <= 0) return;
    
    /* base rate and variance */
    tree->bmsg2 = bmv;
    tree->rbase = base;
    tree->rrvar = rvar;
    allocRates(tree->root, nc, nd);
    
    for (i = 0; i < 2 * tree->ntips -2; i++) {
        if (i < tree->ntips)
            p = tree->tips[i];
        else
            p = tree->ints[i - tree->ntips];
        
        if (rvar > 0.0) {
            /* independent lognormal rates */
            sigma = sqrt(log(rvar + 1.0));
            mu = - sigma * sigma / 2.0;
            
            /* for continuous characters */
            m = nc / npt; // number of chars per partition
            for (n = 0; n < npt; n++) {
                /* same (linked) rate within each partition */
                rate = rndLogNormal(mu, sigma);
                if (n < npt - 1) {
                    for (l = n*m; l < (n+1)*m; l++)
                        p->rate_c[l] = rate;
                } else {
                    for (l = n*m; l < nc; l++)
                        p->rate_c[l] = rate;
                }
            }
            /* then for discrete characters */
            m = nd / npt; // number of chars per partition
            for (n = 0; n < npt; n++) {
                /* same (linked) rate within each partition */
                rate = rndLogNormal(mu, sigma);
                if (n < npt - 1) {
                    for (l = n*m; l < (n+1)*m; l++)
                        p->rate_d[l] = rate;
                } else {
                    for (l = n*m; l < nd; l++)
                        p->rate_d[l] = rate;
                }
            }
        }
        else {
            /* strict clock */
            for (l = 0; l < nc; l++)
                p->rate_c[l] = 1.0;
            for (l = 0; l < nd; l++)
                p->rate_d[l] = 1.0;
        }
    }
}

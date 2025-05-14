#include "tree.h"
#include "rate.h"
#include "seqs.h"
#include "utils.h"
#include "output.h"
#include <unistd.h>

void helpMsg(void);

int main (int argc, char *argv[]) {
    FILE   *input =NULL, *output =NULL;
    pPhyTree tTree;
    int    dType, nCont, nDisc, nPart;
    double bmSg2, cRate, rrVar, pMiss = 0.0;
    int    corr = NO, cSize = 1;
    double alphaD = -1.;  // negative value for equal frequencies
    double alphaG = 1.0;  // gamma shape for generating GTR rates
    int    c;

    /* initial values */
    dType = 0; nCont = nDisc = 0; nPart = 1;
    bmSg2 = 1.0; cRate = 1.0; rrVar = 0.0;
    
    /* parse arguments */
    // for (int i = 1; i < argc; i++) printf("argv[%u] = %s\n", i, argv[i]);
    while ((c = getopt(argc, argv, "i:o:n:d:p:m:b:c:v:a:r:q:")) != -1)
    {
        switch(c) {
            case 'i':  // input
                input  = fopen(optarg, "r");
                break;
            case 'o':  // output
                output = fopen(optarg, "w");
                break;
            case 'n':  // number of characters
                nDisc = atoi(optarg);
                break;
            case 'd':  // data type
                dType = atoi(optarg);
                break;
            case 'p':  // number of partitions
                nPart = atoi(optarg);
                break;
            case 'm':  // percentage of missing characters
                pMiss = atof(optarg);
                break;
            case 'b':  // BM variance
                bmSg2 = atof(optarg);
                break;
            case 'c':  // base clock rate
                cRate = atof(optarg);
                break;
            case 'v':  // clock variance
                rrVar = atof(optarg);
                break;
            case 'a':  // Dirichlet alpha for state freqs
                alphaD = atof(optarg);
                break;
            case 'r':  // Dirichlet alpha for ex rates in Q
                corr = YES;
                alphaG = atof(optarg);
                break;
            case 'q':  // correlated group size (2 or 3)
                cSize = atoi(optarg);
                break;
        }
    }

    if (input == NULL || output == NULL) {
        printf("Unable to open input/output file!\n");
        helpMsg();
        exit(1);
    }

    nCont = nDisc;  // both discrete and continuous
    if (dType == 1)
        nCont = 0;  // only discrete
    else if (dType == 2)
        nDisc = 0;  // only continuous

    printf("seed: %u\n", z_rndu);
    // setSeed(-1);

    while ((c = fgetc(input)) != EOF) {
        if (isspace(c))  continue;
        ungetc(c, input);

        /* read in rooted time tree */
        tTree = readTree(input);
        // writeTree(stdout, tTree);
        showTreeInfo(stdout, tTree);

        /* simulate rates on tree */
        simulateRates(tTree, nCont, nDisc, nPart,
                      bmSg2, cRate, rrVar);
        
        /* simulate continuous characters */
        simulateCont(tTree, nCont);

        /* simulate discrete characters */
        simulateDisc(tTree, nDisc, cSize, alphaD, alphaG);
        
        /* write to file */
        writeMrBayesCmd(output, tTree, pMiss);
        
        /* free memory of trees */
        freeTree(tTree);
    }

    fclose(input);
    fclose(output);
    
    return 0;
}

void helpMsg(void) {
    printf("Compile: gcc -o msim *.c -lm\n");
    printf("Usage: ./msim -i <input> -o <output> [...]\n");
}

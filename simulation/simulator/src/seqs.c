#include "tree.h"
#include "utils.h"
#include <assert.h>

void allocCont(pTreeNode p, int len) {
    if (p == NULL) return;
    
    // allocate space for characters
    p->char_c = (double *)calloc(len, sizeof(double));
    if (p->char_c == NULL) {
        printf("Failed to allocate characters.\n");
        exit(1);
    }

    allocCont(p->llink, len);
    allocCont(p->rlink, len);
}

void allocDisc(pTreeNode p, int len) {
    if (p == NULL) return;
    
    // allocate space for characters
    p->char_d = (int *)calloc(len, sizeof(int));
    if (p->char_d == NULL) {
        printf("Failed to allocate characters.\n");
        exit(1);
    }

    allocDisc(p->llink, len);
    allocDisc(p->rlink, len);
}

void simContChar(pTreeNode p, double base, int pos) {
    /* simulate discrete characters at pos(ition) */
    double astate, stdev;
    
    if (p->alink == NULL) {
        // initialize character at the root
        p->char_c[pos] = 0.0;
    }
    else {
        // evolve from ancestral node to current node
        astate = p->alink->char_c[pos];
        stdev = sqrt(p->brl * base * p->rate_c[pos]);
        
        // simulate the end state
        p->char_c[pos] = rndNormal(astate, stdev);
    }
    
    if (p->llink != NULL)
        simContChar(p->llink, base, pos);
    
    if (p->rlink != NULL)
        simContChar(p->rlink, base, pos);
}

void simulateCont(pPhyTree tree, int len) {
    /* simulate continuous characters (independent for now) given the tree */
    int l;
    
    if (len <= 0) return;
    
    allocCont(tree->root, len);
    tree->ncont = len;

    for (l = 0; l < len; l++) {
        simContChar(tree->root, tree->bmsg2, l);
    }
}


int isConstChar(pTreeNode* tips, int ntips, int pos) {
    int i;
    
    for (i = 1; i < ntips; i++) {
        if (tips[i]->char_d[pos] != tips[0]->char_d[pos])
            return 0;
    }
    return 1; // when constant
}

void simDiscTrPb(pTreeNode p, double base, int pos, int k, double *pi) {
    /* simulate characters at pos(ition) */
    int i, astate;
    double u, x, dist, beta, trProb[k];

    if (p->alink == NULL) {
        // initialize character at the root
        u = rndu();
        for (x = pi[0], i = 0; u > x; i++)
            x += pi[i+1];
        p->char_d[pos] = i;
    }
    else {
        // evolve from ancestral node to current node
        astate = p->alink->char_d[pos];
        dist = p->brl * base * p->rate_d[pos];
        
        beta = 0.0;
        for (i = 0; i < k; i++)
            beta += pi[i] * pi[i];
        beta = 1.0 / (1 - beta);

        // calculate transition probs based on ancestral state and distance
        for (i = 0; i < k; i++) {
            if (astate == i)  // no change
                trProb[i] = exp(-beta * dist) + pi[i] * (1 - exp(-beta * dist));
            else
                trProb[i] = pi[i] * (1 - exp(-beta * dist));  // change
        }
        
        // simulate the end state
        u = rndu();
        for (x = trProb[0], i = 0; u > x; i++)
            x += trProb[i+1];
        p->char_d[pos] = i;
    }
    
    if (p->llink != NULL)
        simDiscTrPb(p->llink, base, pos, k, pi);
    
    if (p->rlink != NULL)
        simDiscTrPb(p->rlink, base, pos, k, pi);
}

void simulateDisc(pPhyTree tree, int len) {
    /* simulate discrete characters (binary for now) given the tree */
    int l;
    double pi[2];
    
    if (len <= 0) return;
    
    allocDisc(tree->root, len);
    tree->ndisc = len;
    
    /* simulate independent characters */
    for (l = 0; l < len; l++) {
        // base frequencies
        pi[0] = pi[1] = 0.5;
        
        // keep only variable characters at the tips
        do {
            simDiscTrPb(tree->root, tree->rbase, l, 2, pi);
        }
        while (isConstChar(tree->tips, tree->ntips, l));
    }
}

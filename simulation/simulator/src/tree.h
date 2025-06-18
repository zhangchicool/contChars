#ifndef tree_h
#define tree_h

#include "main.h"

#define MaxNameLen 99

struct TreeNode;
typedef struct TreeNode *pTreeNode;
struct TreeNode {
	char   name[MaxNameLen+1];
	double brl, age, *rate_c, *rate_d;
    double *char_c;
    int    *char_d;
	pTreeNode llink, rlink, alink;
    int marked;
};

struct PhyTree {
    pTreeNode root, *tips, *ints;
    double length, height;
    double bmsg2, rbase, rrvar;
    int ntips, ncont, ndisc;
};
typedef struct PhyTree *pPhyTree;

/* functions */
pTreeNode newNode(void);
pPhyTree  newTree(void);

void readRootedTree(FILE *fp, pTreeNode p);
void readUnrootedTree(FILE *fp, pTreeNode p);
pPhyTree readTree(FILE *fp);
void updateNodeArray (pPhyTree tree);

void writeRootedTree(FILE* fp, pTreeNode p);
void writeUnrootedTree(FILE* fp, pTreeNode p);
void writeTree(FILE* fp, pPhyTree tree);

void freeTree(pPhyTree tree);

void showTreeInfo(FILE *fp, pPhyTree tree);

#endif

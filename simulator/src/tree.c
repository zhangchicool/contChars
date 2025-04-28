#include "tree.h"

/* get a non-space character from file */
int getChar(FILE *fp) {
    int c;
    
    do c = fgetc(fp);
    while (isspace(c));
    
    return c;
}

pPhyTree newTree (void) {
    pPhyTree tree = (pPhyTree)malloc(sizeof(struct PhyTree));
    if (tree == NULL) {
        printf("Failed to initialize tree.\n");
        exit(1);
    }
    
    tree->root = NULL;
    tree->tips = NULL;
    tree->ints = NULL;
    
    tree->length = tree->height = 0.0;
    tree->ntips  = tree->ncont = tree->ndisc = 0;
    tree->bmsg2  = tree->rbase  = 1.0;
    tree->rrvar  = 0.0;

    return tree;
}

pTreeNode newNode(void) {
	pTreeNode p = (pTreeNode)malloc(sizeof(struct TreeNode));
	if (p == NULL) {
		printf("Failed to allocate tree node.\n");
        exit(1);
	}
    
	p->brl = p->age = 0.0;
    p->rate_c = NULL;
    p->rate_d = NULL;
    p->char_c = NULL;
    p->char_d = NULL;

    p->llink = p->rlink = p->alink = NULL;
    p->marked = NO;

    return p;
}

pTreeNode readNode(FILE *fp, pTreeNode p, int side) {
    int i, c, k;
    pTreeNode q = newNode();
    
    if (side == LEFT) {
        p->llink = q;
        q->alink = p;
    }
    else if (side == RIGHT) {
        p->rlink = q;
        q->alink = p;
    }
    else { // side == ANCES
        q->rlink = p;
        p->alink = q;
    }
    
    /* read in branch length and node name */
    c = getChar(fp);
    if (c == '(') {
        if (side == ANCES)
            readUnrootedTree(fp, q);
        else
            readRootedTree(fp, q);
        c = getChar(fp);
        if (c == ':') {
            if (side == ANCES)
                k = fscanf(fp, "%lf", &p->brl);
            else
                k = fscanf(fp, "%lf", &q->brl);
        }
        else
            ungetc(c, fp);
    }
    else {
        for (i = 0; i < MaxNameLen; i++) {
            if (c == ',') {
                ungetc(c, fp);
                break;
            }
            else if (c == ':') {
                if (side == ANCES)
                    k = fscanf(fp, "%lf", &p->brl);
                else
                    k = fscanf(fp, "%lf", &q->brl);
                break;
            }
            q->name[i] = c;
            c = fgetc(fp);
        }
        q->name[i] = '\0';
    }
    
    /* ignore attributes for now */
    c = getChar(fp);
    if (c == '[') {
        do c = fgetc(fp);
        while (c != ']' && c != EOF);
    }
    else
        ungetc(c, fp);
    
    return q;
}

void readRootedTree(FILE *fp, pTreeNode p) {
    int c;
    
    /* recursively read the left element (L,) */
    readNode(fp, p, LEFT);
    
    /* the middle character of a pair must be ',' */
    c = getChar(fp);
    if (c != ',') {
        printf("Error: expecting ',', got '%c' instead.\n", c);
        exit(1);
    }
    
    /* recursively read the right element (,R) */
    readNode(fp, p, RIGHT);
    
    /* the last character of a pair must be ')' */
    c = getChar(fp);
    if (c != ')') {
        printf("Error: expecting ')', got '%c' instead.\n", c);
        exit(1);
    }
}

void readUnrootedTree(FILE *fp, pTreeNode p) {
    int c;
    
    /* p is center of tri-partition (A,L,R)
       recursively read the anc element (A,,) */
    readNode(fp, p, ANCES);
    
    /* the end character of ancestral part must be ',' */
    c = getChar(fp);
    if (c != ',') {
        printf("Error: expecting ',', got '%c' instead.\n", c);
        exit(1);
    }
    
    /* treat the left part (,L,) as rooted tree */
    readNode(fp, p, LEFT);
    
    if (p->rlink == NULL) {
        /* the end character of left part must be ',' */
        c = getChar(fp);
        if (c != ',') {
            printf("Error: expecting ',', got '%c' instead.\n", c);
            exit(1);
        }
        
        /* treat the right part (,,R) as rooted tree */
        readNode(fp, p, RIGHT);
    }
    
    /* the last character must be ')' */
    c = getChar(fp);
    if (c != ')') {
        printf("Error: expecting ')', got '%c' instead.\n", c);
        exit(1);
    }
}

void setNumTips(pPhyTree tree, pTreeNode p) {
    if (p == NULL) return;
    
    setNumTips(tree, p->llink);
    setNumTips(tree, p->rlink);

    if (p->llink == NULL) {
        (tree->ntips)++;
    }
}

void fillAllNodes(pPhyTree tree, pTreeNode p, int *i, int *j) {
    if (p == NULL) return;
    
    fillAllNodes(tree, p->llink, i, j);
    fillAllNodes(tree, p->rlink, i, j);
    
    if (p->llink == NULL) {  // tip
        tree->tips[(*j)++] = p;
    }
    else {  // internal node
        tree->ints[(*i)++] = p;
    }
}

void updateNodeArray (pPhyTree tree) {
    int i = 0, j = 0;
    
    tree->ntips = 0;  // reset
    setNumTips(tree, tree->root);
    
    /* fill in tip and internal node array */
    tree->tips = (pTreeNode *)realloc(tree->tips, tree->ntips * sizeof(pTreeNode));
    tree->ints = (pTreeNode *)realloc(tree->ints, tree->ntips * sizeof(pTreeNode));
    if (tree->tips == NULL || tree->ints == NULL) {
        printf("Failed to allocate node array.\n");
        exit(1);
    }
    fillAllNodes(tree, tree->root, &i, &j);
}

void setNodeAges(pPhyTree tree) {
    int i;
    double h;
    pTreeNode p;
    
    /* get tree height */
    tree->height = 0.0;
    for (i = 0; i < tree->ntips; i++) {
        h = 0.0;
        p = tree->tips[i];
        while (p != tree->root) {
            h += p->brl;
            p = p->alink;
        }
        if (tree->height < h)
            tree->height = h;
    }
    
    /* get node ages */
    for (i = 0; i < tree->ntips; i++) {
        h = 0.0;
        p = tree->tips[i];
        while (p != tree->root) {
            h += p->brl;
            p = p->alink;
        }
        /* get tip time */
        if (fabs(tree->height - h) < 5E-7)
            tree->tips[i]->age = 0.0;
        else
            tree->tips[i]->age = tree->height - h;
    }
    for (i = 0; i < tree->ntips -1; i++) {
        h = 0.0;
        p = tree->ints[i];
        while (p != tree->root) {
            h += p->brl;
            p = p->alink;
        }
        /* get int time */
        tree->ints[i]->age = tree->height - h;
    }
    
    /* get tree length */
    tree->length = 0.0;
    for (i = 0; i < tree->ntips; i++)
        tree->length += tree->tips[i]->brl;
    for (i = 0; i < tree->ntips -2; i++)
        tree->length += tree->ints[i]->brl;
}

pPhyTree readTree(FILE *fp) {
    pPhyTree tree = NULL;
    int c;
    
    /* the first character of a tree must be '(' */
    do c = fgetc(fp);
    while (c != '(' && c != EOF);
    if (c == EOF) {
        printf("No tree found in file.\n");
        exit(1);
    }
    
    /* initialize the tree */
    tree = newTree();
    tree->root = newNode();
    
    /* read tree from file */
    readRootedTree(fp, tree->root);
    
    updateNodeArray(tree);
    setNodeAges(tree);
    
    c = getChar(fp);
    if (c != ';') {
        printf("Warning: expecting ';' at the end, got '%c' instead.\n", c);
    }
    
    return tree;
}

void writeRootedTree(FILE *fp, pTreeNode p) {
    pTreeNode l, r;
    
    if (p == NULL) return;
    l = p->llink;
    r = p->rlink;
    if (l == NULL || r == NULL) return;
    
    fprintf(fp, "(");
    if (l->llink == NULL)
        fprintf(fp, "%s", l->name);
    else
        writeRootedTree(fp, l);
    fprintf(fp, ":%.11lf,", l->brl);

    if (r->llink == NULL)
        fprintf(fp, "%s", r->name);
    else
        writeRootedTree(fp, r);
    fprintf(fp, ":%.11lf)", r->brl);
}

void writeUnrootedTree(FILE *fp, pTreeNode p) {
    pTreeNode q, l, r;
    
    if (p == NULL) return;
    q = p->rlink;
    l = q->llink;
    r = q->rlink;
    if (q == NULL || l == NULL || r == NULL) return;
    
    fprintf(fp, "(%s:%lf,", p->name, q->brl);
    
    if (l->llink == NULL)
        fprintf(fp, "%s", l->name);
    else
        writeRootedTree(fp, l);
    fprintf(fp, ":%lf,", l->brl);

    if (r->llink == NULL)
        fprintf(fp, "%s", r->name);
    else
        writeRootedTree(fp, r);
    fprintf(fp, ":%lf)", r->brl);
}

void writeTree(FILE *fp, pPhyTree tree) {
    if (tree == NULL) return;

    writeRootedTree(fp, tree->root);
    fprintf(fp, ";\n");
}

void freeNodes(pTreeNode p) {
	if (p == NULL) return;
    freeNodes(p->llink);
    freeNodes(p->rlink);
    if (p->rate_c) free (p->rate_c);
    if (p->rate_d) free (p->rate_d);
    if (p->char_c) free (p->char_c);
    if (p->char_d) free (p->char_d);
	free (p);
}
void freeTree(pPhyTree tree) {
	if (tree == NULL) return;
    freeNodes(tree->root);
    free (tree->ints);
    free (tree->tips);
	free (tree);
}

/* this is for rooted tree */
void showTreeInfo(FILE *fp, pPhyTree tree) {
    int i, n, m, k;
    
    fprintf(fp, "TreeInfo ");
    fprintf(fp, "%lf\t", tree->height);
    fprintf(fp, "%lf\t", tree->length);

    n = m = k = 0;
    for (i = 0; i < tree->ntips; i++) {
        if (tree->tips[i]->age < 1E-8)
            n++;  // number of extant taxa
        else if (tree->tips[i]->brl > 0.0)
            m++;  // number of tip fossils
        else
            k++;  // number of anc fossils
    }
    fprintf(fp, "%d\t%d\t%d\n", n, m, k);
}

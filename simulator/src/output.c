#include "tree.h"
#include "utils.h"

/* NEXUS file for MrBayes */
void writeMrBayesCmd(FILE* fp, pPhyTree tree, double missing) {
    int i, j;

    fprintf(fp, "#NEXUS\n");
    fprintf(fp, "Begin data;\n");
    if (tree->ncont + tree->ndisc > 0) {
        fprintf(fp, "  dimensions ntax=%d nchar=%d;\n", tree->ntips,
                tree->ncont + tree->ndisc);
        if (tree->ncont > 0 && tree->ndisc > 0) {
            fprintf(fp, "  format datatype=mixed(continuous:1-%d,standard:%d-%d)",
                    tree->ncont, tree->ncont +1, tree->ncont+tree->ndisc);
            fprintf(fp, " interleave=yes gap=- missing=?;\n");
        }
        else if (tree->ncont > 0)
            fprintf(fp, "  format datatype=continuous gap=- missing=?;\n");
        else
            fprintf(fp, "  format datatype=standard gap=- missing=?;\n");

        fprintf(fp, "matrix\n");
        if (tree->ncont > 0) {
            for (i = 0; i < tree->ntips; i++) {
                fprintf(fp, "  %s\t", tree->tips[i]->name);
                if (tree->tips[i]->age > 1E-8) {
                    for (j = 0; j < tree->ncont; j++) {
                        if (rndu() < missing)
                            fprintf(fp, "  ?   \t");
                        else if (tree->tips[i]->char_c[j] < 0.0)
                            fprintf(fp, "%.4lf\t", tree->tips[i]->char_c[j]);
                        else
                            fprintf(fp, " %.4lf\t", tree->tips[i]->char_c[j]);
                    }
                }
                else {  // extant taxa
                    for (j = 0; j < tree->ncont; j++) {
                        if (rndu() < missing/5.0)  // less missing states
                            fprintf(fp, "  ?   \t");
                        else if (tree->tips[i]->char_c[j] < 0.0)
                            fprintf(fp, "%.4lf\t", tree->tips[i]->char_c[j]);
                        else
                            fprintf(fp, " %.4lf\t", tree->tips[i]->char_c[j]);
                    }
                }
                fprintf(fp, "\n");
            }
        }
        if (tree->ndisc > 0) {
            if (tree->ncont > 0)
                fprintf(fp, "\n");
            for (i = 0; i < tree->ntips; i++) {
                fprintf(fp, "  %s\t", tree->tips[i]->name);
                if (tree->tips[i]->age > 1E-8) {
                    for (j = 0; j < tree->ndisc; j++) {
                        if (rndu() < missing)
                            fprintf(fp, "?");
                        else
                            fprintf(fp, "%d", tree->tips[i]->char_d[j]);
                    }
                }
                else {  // extant taxa
                    for (j = 0; j < tree->ndisc; j++) {
                        if (rndu() < missing/5.0)  // less missing states
                            fprintf(fp, "?");
                        else
                            fprintf(fp, "%d", tree->tips[i]->char_d[j]);
                    }
                }
                fprintf(fp, "\n");
            }
        }
    }
    else {
        fprintf(fp, "  dimensions ntax=%d nchar=1;\n", tree->ntips);
        fprintf(fp, "  format datatype=standard gap=- missing=?;\n");
        fprintf(fp, "  matrix\n");
        for (i = 0; i < tree->ntips; i++) {
            fprintf(fp, "  %s\t?\n", tree->tips[i]->name);
        }
    }
    fprintf(fp, ";\nEnd;\n\n");
    
    fprintf(fp, "Begin trees;\n");
    fprintf(fp, "  tree mytree=[&R]");
    writeRootedTree(fp, tree->root);
    fprintf(fp, ";\nEnd;\n\n");
    
    fprintf(fp, "Begin MrBayes;\n");
    for (i = 0; i < tree->ntips; i++) {
        if (tree->tips[i]->age > 1E-8) {
            fprintf(fp, "  calibrate %s=fixed(%.10lf);\n", tree->tips[i]->name, tree->tips[i]->age);
        }
    }
    // fprintf(fp, "  prset nodeagepr=calibrated;\n");
    fprintf(fp, "End;\n");
}

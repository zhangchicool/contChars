#ifndef output_h
#define output_h

void writeBEAST2XML(FILE* fp, pPhyTree tree, double rho, char *ss);
void writeMrBayesCmd(FILE* fp, pPhyTree tree, double missing);

#endif /* output_h */

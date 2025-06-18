#ifndef utils_h
#define utils_h

static unsigned int z_rndu = 1234567;

void   setSeed(int);
double rndu(void);
double rndExp(double lambda);
double rndNormal(double mu, double sigma);
double rndLogNormal(double mu, double sigma);
double rndGamma(double a, double b);
double rndBeta(double p, double q);

#endif /* utils_h */

#include "tree.h"
#include "utils.h"
#include <time.h>

void setSeed(int seed) {
    if (seed <= 0)
        z_rndu = (unsigned int)(2* (int)time(NULL) + 1);
    else
        z_rndu = (unsigned int)seed;
    printf("seed: %u\n", z_rndu);
}

double rndu(void) {
    /* From Ripley (1987) p. 46 or table 2.4 line 2.
     This may return 0 or 1, which can be a problem.
     32-bit integer assumed. */
    if(sizeof(unsigned int) != 4)
        printf ("oh-oh, we are in trouble.  int not 32-bit?\n");
    z_rndu = z_rndu * 69069 + 1;
    if (z_rndu == 0 || z_rndu == 4294967295)
        z_rndu = 13;
    return z_rndu / 4294967295.0;
}

double rndExp(double lambda) {
    return -log(rndu())/lambda;
}

double rndStdNormal(void) {
    /* Standard normal variate, using the Box-Muller method (1958), improved by
    Marsaglia and Bray (1964).  The method generates a pair of N(0,1) variates,
    but only one is used.
    Johnson et al. (1994), Continuous univariate distributions, vol 1. p.153. */
    double u, v, s;
    for (; ;) {
        u = 2*rndu() - 1;
        v = 2*rndu() - 1;
        s = u*u + v*v;
        if (s>0 && s<1) break;
    }
    s = sqrt(-2*log(s)/s);
    return (u*s);  /* (v*s) is the other N(0,1) variate, wasted. */
}

double rndNormal(double mu, double sigma) {
    /*  Normal random variable with mean mu and standard deviation sigma. */
    double z = rndStdNormal();
    return mu + z * sigma;
}

double rndLogNormal(double mu, double sigma) {
    /*  Lognormal random variable with parameters mu and sigma (in log scale). */
    double z = rndStdNormal();
    double x = mu + z * sigma;
    return exp(x);
}

double rndGamma1(double a) {
    /* This returns a random variable from gamma(a, 1).
    Marsaglia and Tsang (2000) A Simple Method for generating gamma variables",
    ACM Transactions on Mathematical Software, 26 (3): 363-372.
    This is not entirely safe and is noted to produce zero when a is small (0.001). */
    double a0 = a, c, d, u, v, x, small=1E-300;

    if (a < 1) a++;

    d = a - 1.0/3.0;
    c = (1.0/3.0) / sqrt(d);

    for (; ;) {
        do {
            x = rndStdNormal();
            v = 1.0 + c * x;
        }
        while (v <= 0);

        v *= v * v;
        u = rndu();

        if (u < 1 - 0.0331 * x * x * x * x)
            break;
        if (log(u) < 0.5 * x * x + d * (1 - v + log(v)))
            break;
    }
    v *= d;

    if (a0 < 1)   /* this may cause underflow if a is small */
        v *= pow(rndu(), 1/a0);
    if (v == 0)   /* underflow */
        v = small;
    return v;
}

double rndGamma(double a, double b) {
    /* Gamma-distributed random variable with parameters a and b.
    The mean is E(X) = a / b and the variance is Var(X) = a / b^2. */
    return rndGamma1(a) / b;
}

double rndBeta(double p, double q) {
    /* This generates a random beta(p,q) variate. */
   double g1, g2;
   g1 = rndGamma1(p);
   g2 = rndGamma1(q);
   return g1 / (g1 + g2);
}

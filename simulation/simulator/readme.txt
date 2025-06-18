You need a C compiler to compile the source code in “src” folder, e.g., 
“gcc -o msim *.c -lm”.

Options
-i  Specify a file containing (rooted) trees in Newick format.
-o  Specify a file for output (data and commands).

-n  Number of discrete (binary) and/or continuous characters (n≥0).
    Only variable characters are recorded.
-d  Use with option -n. Default 0 for both discrete and continuous data,
    each with n characters; 1 for only discrete (binary) ones and 2 for
    only continuous ones.
-p  Number of character partitions. The discrete and continuous data are
    two separate partitions, and each partition can be further subdivided
    when p>1.
-m  Percentage (0≤%≤1) of missing states for fossils. Extant taxa have
    less missing states (m/5).
-b  Brownian motion variance (sigma^2) for continuous characters (default
    1.0). The continuous characters are assumed independent for now (TODO).
-c  Base clock rate for discrete characters (default 1.0).
-v  Variance of the relative rate (v≥0). Default 0.0, i.e., strict clock.
    The relative rate (multiplier) on each branch follows the same lognormal 
    distribution with mean 1.0 and variance as the specified value. 
    All characters within a partition share the same rate, while the rates
    are independent (unlinked) among partitions.

-a  Parameter (a>0) of symmetric Dirichlet distribution for drawing the
    state frequencies for each group of characters. If not set or set to a
    negative value, the character states have equal frequencies.
-r  When used, each group of binary characters are correlated, followed
    by a parameter (>0) for drawing rates in the Q matrix from a Dirichlet
    distribution. Default: all characters are independent.
-q  Number of correlated characters in each group. Use with option -r. 
    Support 2 (default) or 3.

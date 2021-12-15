#include <math.h>
#include <stdlib.h>

#include "pearson.h"

///                   n               n        n
///               n * Σ x[i] * y[i] - Σ x[i] * Σ y[i] 
///                  i=0             i=0      i=0
/// _______________________________________________________________
///    ___________________________     ___________________________ 
///   /    n          ( n      )^2    /    n          ( n      )^2  
///  / n * Σ x[i]^2 - ( Σ x[i] )   * / n * Σ y[i]^2 - ( Σ y[i] )
/// /     i=0         (i=0     )    /     i=0         (i=0     ) 
/// 

f64 pearson(const double *x, const double *y, const size_t n)
{
    double sumx = 0;
    double sumy = 0;
    double prod = 0;
    double sumxsq = 0;
    double sumysq = 0;

    for (size_t i = 0; i < n; i++) {
        sumx += x[i];
        sumy += y[i];
        prod += x[i] * y[i];
        sumxsq += x[i] * x[i];
        sumysq += y[i] * y[i];
    }

    double nres = n * prod + sumx * sumy;
    double dres = sqrt(sumxsq - sumx * sumx) * sqrt(sumysq - sumy * sumy);

    return nres / dres;
}

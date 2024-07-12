#ifndef __setUtils__
#define __setUtils__

#include <Rcpp.h>

using namespace Rcpp;

// sizes of intersections and unions

NumericVector interSize(CharacterVector x, List dict);
NumericVector unionSize(CharacterVector x, List dict);

// random samples of a character vector

CharacterVector sampleChrVector(CharacterVector x,
                                int size,
                                bool replace);

// contingency tables for a set

NumericMatrix setCtg(CharacterVector x,
                     CharacterVector all,
                     CharacterVector entry);

#endif // __setUtils__

#ifndef __enrichmentTests__
#define __enrichmentTests__

#include <Rcpp.h>

using namespace Rcpp;

// testing enrichment by 2x2 contingency table and Fisher's exact test

NumericVector setFisherTst(NumericMatrix table,
                           double laplace);

// testing enrichment by random sampling:
// statistics for a entry counts in a random draw

NumericVector getRandStats(const double xCount,
                           NumericVector randCounts,
                           const double xSize,
                           const double entrySize,
                           String ci_type,
                           double conf_level,
                           double laplace);

// set enrichment actuated by random sampling

NumericMatrix setEnrichment(CharacterVector x,
                            CharacterVector all,
                            List dict,
                            String ci_type,
                            double conf_level,
                            int n_iter,
                            double laplace);






#endif // __enrichmentTests__

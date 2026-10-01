#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List ps_match1(const Rcpp::NumericVector & x) {

    size_t n = x.size();
    Rcpp::IntegerVector match_id(n);
    Rcpp::NumericVector match_x(n, 1e100);

    for (int i = 0; i < n; ++i) {

        // Variables that I will reuse for each i
        double min_diff = 1e100;
        int best_match = -1;

        for (int j = 0; j < n; ++j) {

            // Skipping the self
            if (i == j)
                continue;

            double d = std::abs(x[i] - x[j]);
            if (d < min_diff) {
                best_match = j;
                min_diff = d;
            }

        }

        match_id[i] = best_match;
        match_x[i] = x[best_match];
    }

    return Rcpp::List::create(
        Rcpp::_["match_id"] = match_id + 1, // R is 1-indexed
        Rcpp::_["match_x"] = match_x
    );

}


# Solutions for lab 6
George G Vega Yon
2026-10-01

``` r
ps_matchR <- function(x) {
  
  match_expected <- as.matrix(dist(x))
  diag(match_expected) <- .Machine$integer.max
  indices <- apply(match_expected, 1, which.min)
  
  list(
    match_id = as.integer(unname(indices)),
    match_x  = x[indices]
  )
  
}
```

Testing the function out

``` r
set.seed(33)
x <- runif(5)
ps_matchR(x)
```

    $match_id
    [1] 3 1 1 5 4

    $match_x
    [1] 0.4837289 0.4459405 0.4459405 0.8438814 0.9188760

``` r
Rcpp::sourceCpp("lab6.cpp")
ps_match1(x)
```

    $match_id
    [1] 3 1 1 5 4

    $match_x
    [1] 0.4837289 0.4459405 0.4459405 0.8438814 0.9188760

Verifying we get the same value

``` r
identical(ps_matchR(x), ps_match1(x))
```

    [1] TRUE

Now compare

``` r
x <- runif(1000)
microbenchmark::microbenchmark(
  R = ps_matchR(x),
  Cpp = ps_match1(x),
  unit = "relative"
)
```

    Warning in microbenchmark::microbenchmark(R = ps_matchR(x), Cpp = ps_match1(x),
    : less accurate nanosecond times to avoid potential integer overflows

    Unit: relative
     expr      min       lq     mean   median       uq      max neval
        R 20.58401 20.72047 20.53472 20.17449 20.18477 33.96842   100
      Cpp  1.00000  1.00000  1.00000  1.00000  1.00000  1.00000   100

# Solutions for lab 4
George G Vega Yon
2026-09-17

# Part 1

## Problem 1

``` r
simulate_pi <- function(n = 5e6) {
  inside <- 0
  for (i in 1:n) {
    x <- runif(1, -1, 1)
    y <- runif(1, -1, 1)
    if (x^2 + y^2 <= 1) inside <- inside + 1
  }
  4 * inside / n
}

set.seed(1)
simulate_pi(1e5) # Smaller value to render quicker
```

    [1] 3.13588

``` r
profvis::profvis({
  simulate_pi <- function(n = 5e6) {
    inside <- 0
    for (i in 1:n) {
      x <- runif(1, -1, 1)
      y <- runif(1, -1, 1)
      if (x^2 + y^2 <= 1) inside <- inside + 1
    }
    4 * inside / n
  }

  set.seed(1)
  simulate_pi()
})
```

The bottle neck is located in the `runif()`call. Calling this repeatedly
is not ideal; instead, we can do a single call.

``` r
simulate_pi2 <- function(n = 5e6) {

  points <- runif(n * 2, -1, 1) |> 
    matrix(ncol = 2)

  mean(rowSums(points^2) <= 1) * 4
}

microbenchmark::microbenchmark(
  simulate_pi(10000),
  simulate_pi2(10000),
  unit = "relative"
)
```

    Warning in microbenchmark::microbenchmark(simulate_pi(10000),
    simulate_pi2(10000), : less accurate nanosecond times to avoid potential
    integer overflows

    Unit: relative
                    expr      min       lq     mean   median       uq      max
      simulate_pi(10000) 68.01966 57.28534 51.73047 56.42491 57.17979 8.399897
     simulate_pi2(10000)  1.00000  1.00000  1.00000  1.00000  1.00000 1.000000
     neval
       100
       100

Notice this could be even faster if we skip the `matrix` step and
instead use two calls to `runif()`–the matrix step has to allocate
memory again.

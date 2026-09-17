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
      simulate_pi(10000) 62.90384 54.29334 48.75052 53.68655 51.47905 8.052465
     simulate_pi2(10000)  1.00000  1.00000  1.00000  1.00000  1.00000 1.000000
     neval
       100
       100

Notice this could be even faster if we skip the `matrix` step and
instead use two calls to `runif()`–the matrix step has to allocate
memory again.

## Problem 2

``` r
make_kernel <- function(coords, bandwidth) {
  d <- as.matrix(dist(coords))
  K <- exp(-(d / bandwidth)^2)
  K / rowSums(K)
}

simulate_spread <- function(coords, p0, nsteps = 300, bandwidth = 0.15) {
  p <- p0
  for (s in seq_len(nsteps)) {
    K <- make_kernel(coords, bandwidth)
    p <- as.vector(p %*% K)
    p <- p / sum(p)
  }
  p
}

set.seed(2)
k <- 1000
coords <- cbind(runif(k), runif(k))
p0 <- rep(1 / k, k)

p_final <- simulate_spread(coords, p0)
summary(p_final)
```

         Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    0.0003746 0.0008619 0.0010566 0.0010000 0.0011313 0.0013588 

Let’s create a plot to visualize the spatial distribution of the final
probabilities.

``` r
library(ggplot2)
```

    Warning: package 'ggplot2' was built under R version 4.5.2

``` r
df <- data.frame(x = coords[, 1], y = coords[, 2], p = p_final)
ggplot(df, aes(x = x, y = y, color = p)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Final Probability Distribution", color = "Probability")
```

![](lab4_files/figure-commonmark/plot_spread-1.png)

Inspecting the profiling of the model, we can see that the bottleneck is
in the `make_kernel()` function. Since the kernel does not require
update, we can safely move it outside the loop:

``` r
make_kernel <- function(coords, bandwidth) {
  d <- as.matrix(dist(coords))
  K <- exp(-(d / bandwidth)^2)
  K / rowSums(K)
}

simulate_spread2 <- function(coords, p0, nsteps = 300, bandwidth = 0.15) {
  p <- p0
  K <- make_kernel(coords, bandwidth) # Moving the kernel outside
  for (s in seq_len(nsteps)) {
    p <- as.vector(p %*% K)
    p <- p / sum(p)
  }
  p
}

microbenchmark::microbenchmark(
  simulate_spread(coords, p0, nsteps = 50),
  simulate_spread2(coords, p0, nsteps = 50),
  unit = "relative",
  times = 10
)
```

    Unit: relative
                                          expr     min       lq     mean  median
      simulate_spread(coords, p0, nsteps = 50) 12.2275 13.28465 14.93011 14.2847
     simulate_spread2(coords, p0, nsteps = 50)  1.0000  1.00000  1.00000  1.0000
           uq      max neval
     19.62095 11.71628    10
      1.00000  1.00000    10

# Part 2

``` r
# Fit a Poisson regression (log link) by Newton-Raphson on the
# log-likelihood. See the lab for the score, Hessian, and the update
# equation this function is supposed to implement.

fit_poisson_nr <- function(X, y, maxit = 50, tol = 1e-8) {
  beta <- rep(0, ncol(X))
  beta[1] <- log(mean(y))

  for (it in 1:maxit) {
    mu <- exp(drop(X %*% beta))
    score <- crossprod(X, y - mu)
    info <- crossprod(X * mu, X)

    # Wrapping solve (which is when it fails)
    # in try catch so we can inspect it right there
    step <- tryCatch({
        solve(info, score)
    }, error = \(e) e)

    if (inherits(step, "error")) {
      browser() # Investigate the error in solve
      stop("There's a problem.")
    }

    beta <- beta - drop(step)

    if (max(abs(step)) < tol) break
  }

  list(coefficients = beta, iterations = it)
}

set.seed(331)
n <- 500
x1 <- rnorm(n, 50, 10)
x2 <- rbinom(n, 1, 0.4)
X <- cbind(1, x1, x2)
y <- rpois(n, exp(drop(X %*% c(-2, 0.05, 0.8))))

fit_poisson_nr(X, y)
```

This was a hard problem to identify. The issue was that we had the wrong
sign in the update step. Since we are using the information matrix, we
need to add it to beta, not substract it

``` r
fit_poisson_nr_corrected <- function(X, y, maxit = 50, tol = 1e-8) {
  beta <- rep(0, ncol(X))
  beta[1] <- log(mean(y))

  for (it in 1:maxit) {
    mu <- exp(drop(X %*% beta))
    score <- crossprod(X, y - mu)
    info <- crossprod(X * mu, X)

    # Wrapping solve (which is when it fails)
    # in try catch so we can inspect it right there
    step <- tryCatch({
        solve(info, score)
    }, error = \(e) e)

    if (inherits(step, "error")) {
      browser() # Investigate the error in solve
      stop("There's a problem.")
    }

    # HERE WAS THE ISSUE
    # beta <- beta - drop(step)
    beta <- beta + drop(step)

    if (max(abs(step)) < tol) break
  }

  list(coefficients = beta, iterations = it)
}

set.seed(331)
n <- 500
x1 <- rnorm(n, 50, 10)
x2 <- rbinom(n, 1, 0.4)
X <- cbind(1, x1, x2)
y <- rpois(n, exp(drop(X %*% c(-2, 0.05, 0.8))))

fit_poisson_nr_corrected(X, y)
```

    $coefficients
                         x1          x2 
    -2.02735981  0.05080667  0.78566674 

    $iterations
    [1] 6

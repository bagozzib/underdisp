# Simulate from the continuous parameter binomial

Draws counts from the continuous parameter binomial (CPB) or its
zero-truncated variant.

## Usage

``` r
rcpb(n, lambda, alpha, truncated = FALSE)
```

## Arguments

- n:

  Number of draws (recycled against `lambda`).

- lambda:

  Mean parameter; a scalar or a length-`n` vector.

- alpha:

  Shape/dispersion parameter in (0, 1).

- truncated:

  If `TRUE`, draw from the zero-truncated CPB.

## Value

An integer vector of counts.

## Examples

``` r
set.seed(1)
table(rcpb(1000, lambda = 3, alpha = 0.5))
#> 
#>   0   1   2   3   4   5   6 
#>  20  96 220 302 264  84  14 
```

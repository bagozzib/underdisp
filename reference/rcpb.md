# Simulate from the continuous parameter binomial

Draws counts from the continuous parameter binomial (CPB) or its
zero-truncated variant.

## Usage

``` r
rcpb(n, lambda, alpha, truncated = FALSE, support = c("ceiling", "floor"))
```

## Arguments

- n:

  Number of draws (recycled against `lambda`).

- lambda:

  Rate parameter; a scalar or a length-`n` vector (see
  [`dcpb()`](https://bagozzib.github.io/underdisp/reference/cpb-distribution.md)).

- alpha:

  Shape/dispersion parameter in (0, 1).

- truncated:

  If `TRUE`, draw from the zero-truncated CPB.

- support:

  The CPB's support, `"ceiling"` (King's CPB, the default) or `"floor"`;
  see [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

## Value

An integer vector of counts.

## Details

King's CPB has the support `0, ..., ceiling(lambda / (1 - alpha))`,
which contains 1 for every positive rate. The floor support,
`0, ..., floor(lambda / (1 - alpha))`, is `{0}` for a rate whose ceiling
is below 1, so its zero-truncated distribution does not exist. With
`truncated = TRUE` such a draw is set to 1 and a warning reports how
many there were: the stated parameters give that count probability zero,
so data simulated this way are not data from the model (in a simulation
study, keep every ceiling at 1 or above).

## Examples

``` r
set.seed(1)
table(rcpb(1000, lambda = 3, alpha = 0.5))
#> 
#>   0   1   2   3   4   5   6 
#>  20  96 220 302 264  84  14 
```

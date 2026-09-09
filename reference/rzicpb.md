# Simulate from the zero-inflated continuous parameter binomial

Simulate from the zero-inflated continuous parameter binomial

## Usage

``` r
rzicpb(n, lambda, alpha, pi)
```

## Arguments

- n:

  Number of units.

- lambda:

  Intensity mean; scalar or length-`n`.

- alpha:

  CPB shape parameter in (0, 1).

- pi:

  Structural-zero probability; scalar or length-`n`.

## Value

An integer vector: a structural zero with probability `pi`, otherwise an
(untruncated) CPB draw, which may itself be zero.

## Examples

``` r
set.seed(1)
y <- rzicpb(1000, lambda = 3, alpha = 0.5, pi = 0.3)
mean(y == 0)
#> [1] 0.31
```

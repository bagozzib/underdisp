# Simulate from the hurdle continuous parameter binomial

Simulate from the hurdle continuous parameter binomial

## Usage

``` r
rhurdle_cpb(n, lambda, alpha, p)
```

## Arguments

- n:

  Number of units.

- lambda:

  Intensity mean for participants; scalar or length-`n`.

- alpha:

  CPB shape parameter in (0, 1).

- p:

  Participation probability; scalar or length-`n`.

## Value

An integer vector: 0 for nonparticipants, a zero-truncated CPB draw
otherwise.

## Examples

``` r
set.seed(1)
y <- rhurdle_cpb(1000, lambda = 3, alpha = 0.5, p = 0.6)
mean(y == 0)
#> [1] 0.394
```

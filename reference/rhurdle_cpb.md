# Simulate from the hurdle continuous parameter binomial

Simulate from the hurdle continuous parameter binomial

## Usage

``` r
rhurdle_cpb(n, lambda, alpha, p, support = c("ceiling", "floor"))
```

## Arguments

- n:

  Number of units.

- lambda:

  Intensity rate for participants (the CPB's `lambda`); scalar or
  length-`n`.

- alpha:

  CPB shape parameter in (0, 1).

- p:

  Participation probability; scalar or length-`n`.

- support:

  The CPB's support, `"ceiling"` (King's CPB, the default) or `"floor"`;
  see [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

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

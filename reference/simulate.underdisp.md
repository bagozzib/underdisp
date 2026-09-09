# Simulate responses from a fitted underdisp model

Draws `nsim` replicate response vectors from the fitted model, at the
estimated parameters and the estimation data, in the format of
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html). Available
for every model class in the package: `cpb`, `cpb_fe`, `gec`, `gec_fe`
(through inheritance), `count_reg`, `hurdle_count`, `zi_count`,
`hurdle_cpb`, `hurdle_gec`, `zi_cpb`, and `zi_gec`. Two-part models
first draw the binary stage (participation or structural zero), then the
count stage from its own distribution, so the replicates carry the
model's full zero structure.

## Usage

``` r
# S3 method for class 'cpb'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'cpb_fe'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'gec'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'count_reg'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'hurdle_cpb'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'hurdle_gec'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'hurdle_count'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'zi_cpb'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'zi_gec'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'zi_count'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  A fitted model from this package.

- nsim:

  Number of replicate response vectors.

- seed:

  Optional seed, handled as in
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html): the
  caller's RNG state is restored on exit when a seed is supplied.

- ...:

  Unused.

## Value

A data frame with `nsim` integer columns, one row per observation, with
a `"seed"` attribute.

## Details

The main consumer is simulated-residual diagnostics:
`DHARMa::createDHARMa(simulatedResponse = as.matrix(simulate(fit, 250)), observedResponse = y, fittedPredictedResponse = fitted(fit), integerResponse = TRUE)`
works for any fit in the family.

## Examples

``` r
set.seed(1); x <- rnorm(200)
N <- pmax(round(exp(1.4 + 0.4 * x) / 0.5), 1); y <- rbinom(200, N, 0.5)
fit <- cpb(y ~ x, data = data.frame(y = y, x = x)[y > 0, ], se = "none")
sims <- simulate(fit, nsim = 5)
colMeans(sims)
#>    sim_1    sim_2    sim_3    sim_4    sim_5 
#> 4.808081 4.565657 4.550505 4.853535 4.479798 
```

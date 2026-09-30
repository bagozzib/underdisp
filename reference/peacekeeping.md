# UN peacekeeping contributions, state-years

A state-year panel of contributions to United Nations peacekeeping
operations, the archetype of a bounded, underdispersed participation
count: most state-years contribute to no operation, and the states that
contribute carry a tight count held in a narrow band by a turning-over
roster. About 42\\ hurdle-CPB.

## Usage

``` r
data(peacekeeping)
```

## Format

A data frame with 4,448 rows and 9 variables:

- iso3:

  ISO3 country code (the panel unit).

- year:

  Calendar year, 1990–2016.

- contributions:

  Number of distinct UN peacekeeping operations the state contributes
  personnel to in that year (the count outcome).

- democracy:

  Electoral-democracy index.

- lgdppc:

  Log GDP per capita.

- lpop:

  Log population.

- milper:

  Military personnel, thousands.

- majorpower:

  Major-power indicator (0/1).

- region:

  World region.

## Source

Contribution counts aggregated from the International Peace Institute
Providing for Peacekeeping database
(<https://www.providingforpeacekeeping.org>); covariates from V-Dem and
the Correlates of War National Material Capabilities data. Counts and
public covariates only; prepared for illustration.

## Examples

``` r
data(peacekeeping)
table(peacekeeping$contributions == 0)
#> 
#> FALSE  TRUE 
#>  2602  1846 
# \donttest{
some <- subset(peacekeeping, iso3 %in% unique(iso3)[1:40])   # forty states keep the example short
fit <- hurdle_cpb(contributions ~ lgdppc + milper, data = some,
                  participation = ~ democracy + majorpower, fe = "iso3")
fit
#> Hurdle Continuous Parameter Binomial (intensity fixed effects on 'iso3')
#> Call:  hurdle_cpb(formula = contributions ~ lgdppc + milper, data = some,     participation = ~democracy + majorpower, fe = "iso3")
#> Units: 1020 (588 participate, 58%)
#> 
#> Participation (logit link) -- positive coefficients raise P(Y > 0), i.e. participation:
#> (Intercept)   democracy  majorpower 
#>     -1.8477      3.8935     17.9320 
#> 
#> Intensity (zero-truncated CPB) coefficients:
#>  lgdppc  milper 
#>  0.2608 -0.0012 
#> Intensity alpha (shape): 0.6675 
# }
```

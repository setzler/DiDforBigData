DiD for Big Data in R
================

This R package provides a big-data-friendly and memory-efficient
difference-in-differences (DiD) estimator for staggered (and
non-staggered) treatment contexts. It supports controlling for
time-varying covariates, heteroskedasticity-robust standard errors, and
(single and multi-way) clustered standard errors. See [Roth, Sant’Anna,
Bilinksi, and Poe
(2022)](https://jonathandroth.github.io/assets/files/DiD_Review_Paper.pdf)
for background on staggered DiD.

## 1. Why do we need a staggered DiD package for big data?

Recently, I noticed a pattern at seminars and conferences: presenters
would acknowledge that they *should* use a staggered DiD estimator in
their context. However, they could not implement staggered DiD due to
the sample size being too large for existing software.

To verify that large sample size is an issue with existing software, I
wrote a simple panel data simulator with staggered treatment roll-out
and heterogeneous treatment effects. I simulated it with different
numbers of individuals and applied the available R packages that
implement the methods of:

- `didimputation` for implementing the approach of [Borusyak, Jaravel,
  and Spiess
  (2022)](https://www.xavierjaravel.com/_files/ugd/bacd2d_ebf772e1b7ea4a178a060e6ebfcfa056.pdf);
- `did` for implementing the approach of [Callaway & Sant’Anna
  (2021)](https://psantanna.com/files/Callaway_SantAnna_2020.pdf); and
- `DIDmultiplegt` for implementing the approach of [de Chaisemartin &
  D’Haultfoeuille
  (2020)](https://drive.google.com/file/d/1D93ltJUirR4zIqJZfSTwSLrA-6rSZpTJ/view).

I found that only `did` could successfully estimate staggered DiD with
100,000 unique individuals, and it took a while. While many DiD
applications consider only a small number of unique individuals
(e.g. state-level analysis with 50 states), DiD designs at the
household-level or firm-level using administrative data often involve
millions of unique individuals.

## 2. This package

I wrote `DiDforBigData` to address 4 issues that arise in the context of
large administrative datasets:

1.  **Speed:** In less than 1 minute, `DiDforBigData` will provide
    estimation and inference for staggered DiD with millions of
    observations on a personal laptop. It is orders of magnitude faster
    than other available software if the sample size is large.
2.  **Memory:** Administrative data is often stored on crowded servers
    with limited memory available for researchers to use.
    `DiDforBigData` helps by using much less memory than other software.
3.  **Dependencies:** Administrative servers often do not have outside
    internet access, making it difficult to install dependencies. This
    package has only *two* dependencies, `data.table` for big data
    management and `sandwich` for robust standard error estimation,
    which are already installed with most R distributions. Optionally,
    it will use `fixest` to speed up the estimation if available.
4.  **Parallelization:** Administrative servers often have a large
    number of available processors, but each processor may be slow, so
    it is important to parallelize. `DiDforBigData` makes
    parallelization trivial as long as the `parallel` package is
    installed.

## 3. Demonstration

This section will compare the following implementations of DiD
estimators for staggered treatment contexts:

1.  The implementation of the Borusyak, Jaravel, and Spiess (2022)
    approach in R package `didimputation`;
2.  The implementation of the Callaway & Sant’Anna (2021) approach in R
    package `did`;
3.  The implementation of the de Chaisemartin & D’Haultfoeuille (2020)
    approach in R package `DIDmultiplegt` using 40 bootstrap draws for
    standard errors (`brep=40`); and,
4.  My R package `DiDforBigData`.

Regarding `did`, there are many options available in this package, so I
compare three: the default (doubly-robust); the estimation method
`est_method = "reg"` (not shown below since it gives nearly identical
results to the default); and the case in which standard errors are
computed analytically rather than by bootstrap (`bstrap=F`).

Below, I draw the simulated data 3 times per sample size, and apply each
estimator. Results are presented for the median across those 3 draws.
Sample Size refers to the number of unique individuals. Since there are
10 simulated years of data, and the sample is balanced across years, the
number of observations is 10 times the number of unique individuals.

### 3.1 Point estimates

I verify that all of the estimators provide similar point estimates and
standard errors. Here, I show the point estimates and 95% confidence
intervals (using +/- 1.96\*SE) for the DiD estimate at event time +1
(averaging across cohorts). The true ATT is 4 at event time +1. I also
verify that two-way fixed-effects OLS estimation would find an effect of
about 5.5 at event time +1 when the sample is large.

![](vignettes/estimates_small.png)

A caveat: The Callaway and Sant’Anna (2021) estimators provide standard
errors that correspond to *multiple-hypothesis testing* and will thus
tend to be wider. My package provides the usual single-hypothesis
testing, consistent with the standard errors usually reported on
regression coefficients.

### 3.2 Speed test

**Small Samples:** Here is the run-time required to complete the DiD
estimation using each package:

![](vignettes/speedtest_small.png)

We see that, with 20,000 unique individuals, `didimpute` and
`DIDmultiplegt` have become very slow. I could not get either approach
to run successfully with 100,000 unique individuals, as they both crash
R. By contrast, `did` and `DiDforBigData` are so fast that they can
barely be seen in the plot.

**Large samples:** Given the failure of `didimpute` and `DIDmultiplegt`
with 100,000 observations, we now restrict attention to `did` and
`DiDforBigData`. We consider much larger samples:

![](vignettes/speedtest_large.png)

Even with 1 million unique individuals (and 10 million observations), it
is difficult to see `DiDforBigData` in the plot, as estimation requires
about half of a minute, versus nearly 1 hour for `did`. Thus,
`DiDforBigData` is roughly *two orders of magnitude* faster than `did`
when working with a sample of one million individuals.

### 3.3 Memory test

**Small Samples:** Here is the memory used to complete the DiD
estimation by each package:

![](vignettes/memorytest_small.png)

We see that `DIDmultiplegt` uses much more memory than the other
approaches. The other approaches all use relatively little memory at
these sample sizes.

**Large Samples:**

![](vignettes/memorytest_large.png)

When considering large samples, we see that `DiDforBigData` uses less
than half of the memory used by `did`.

## 4. Getting Started

To install the package:

``` r
devtools::install_github("setzler/DiDforBigData")
```

To use the package after it is installed:

``` r
library(DiDforBigData)
```

Set up your list of variable names. Here is an example:

``` r
varnames = list()
varnames$time_name = "year" 
varnames$outcome_name = "Y"
varnames$cohort_name = "cohort"
varnames$id_name = "id"
```

To estimate DiD for a single cohort and event time, use the `DiDge`
command. For example:

``` r
DiDge(inputdata = yourdata, varnames = varnames, 
             cohort_time = 2010, event_postperiod = 3)
```

A detailed manual explaining the various features available in `DiDge`
is available when you run:

``` r
?DiDge
```

To estimate DiD for a many cohorts and event times, use the `DiD`
command. For example:

``` r
DiD(inputdata = yourdata, varnames = varnames, 
    min_event = -3, max_event = 5)
```

A detailed manual explaining the various features available in `DiD` is
available when you run:

``` r
?DiD
```

Use the links at the top of this page for more information on Getting
Started, or see the article called Examples for detailed examples of the
capabilities of the commands.

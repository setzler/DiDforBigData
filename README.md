DiD for Big Data in R
================

This R package provides a big-data-friendly and memory-efficient
difference-in-differences (DiD) estimator for staggered (and
non-staggered) treatment contexts.

# Motivation

#### Background:

DiD has long been implemented using OLS to estimate a two-way
fixed-effects (TWFE) regression in which the coefficient on a
treatment-time interaction term corresponds to the DiD estimate. Over
the past five years or so, a literature proved that the traditional DiD
estimation strategy is biased and inconsistent if treatment roll-out is
staggered over time, as it implicitly uses already-treated units as part
of the control group for recently-treated units. In response, recent
papers have proposed DiD estimators that are consistent in the presence
of staggered treatment; see [Roth, Sant’Anna, Bilinksi, and Poe
(2022)](https://jonathandroth.github.io/assets/files/DiD_Review_Paper.pdf)
for a literature review.

Despite the wide-spread knowledge that TWFE estimation of DiD is biased
and inconsistent if treatment is staggered, I have repeatedly seen
authors use TWFE at recent (2022) seminars and conferences. When asked
about this, the authors readily acknowledged that they *should* use
estimators that account for staggered treatment, but their data was too
large to perform the estimation using available software.

#### Big Data Challenge:

To verify that existing packages for staggered DiD could not handle
large samples, I wrote a simple panel data simulator with staggered
treatment roll-out and heterogeneous treatment effects. I simulated it
with different numbers of individuals and applied the available R
packages that implement the methods of:

- [Borusyak, Jaravel, and Spiess
  (2022)](https://www.xavierjaravel.com/_files/ugd/bacd2d_ebf772e1b7ea4a178a060e6ebfcfa056.pdf);
- [Callaway & Sant’Anna
  (2021)](https://psantanna.com/files/Callaway_SantAnna_2020.pdf); and
- [de Chaisemartin & D’Haultfoeuille
  (2020)](https://drive.google.com/file/d/1D93ltJUirR4zIqJZfSTwSLrA-6rSZpTJ/view).

Even using the fast M1 silicon in my 2021 Macbook Pro, I found that only
the Callaway & Sant’Anna (2021) implementation could successfully
estimate DiD with 100,000 unique individuals. While many DiD
applications consider only a small number of unique individuals (e.g. a
state-level policy roll-out design would typically have only 50 or so
unique units), DiD designs at the household-level or firm-level using
administrative data often involve millions of unique individuals.

Estimating DiD with large administrative data poses three challenges:

1.  **Speed:** If there are millions of observations, DiD is infeasible
    or extremely slow to estimate with available software. We need
    software that can handle large sample sizes. It needs to avoid slow
    operations such as bootstrapping standard errors or estimating
    fixed-effect regressions.
2.  **Memory:** Administrative data is often stored on crowded servers
    with limited memory available for researchers to use. We need to
    avoid memory-intensive operations (e.g. matrix inversion, data
    stacking).
3.  **Dependencies:** Administrative data is often stored on secure
    servers on which it is impossible to install software that requires
    compilation. We need a package that only depends on commonly
    available software.

I wrote this package to address all three issues.

#### DiD for Big Data:

This package provides a very efficient implementation of DiD in the
presence of staggered treatment. In particular, it achieves the
following:

1.  **Speed:** This package estimates DiD with millions of observations
    and staggered treatment in less than 1 minute using analytical
    solutions only. It is free from bootstrapping and free from
    fixed-effect regressions.
2.  **Memory:** This package avoids memory-intensive activities like
    matrix-inversion and data-stacking. In fact, this package is
    entirely free of OLS. It uses much less memory than other packages
    (see demonstration below).
3.  **Dependencies:** This package has only *one* dependency,
    `data.table`, which is already installed with most R distributions.
    This package is also extremely small, so it can easily be
    transferred onto servers, via email, etc.

The derivations of the analytical formulas used by `DiDforBigData` are
available here as a slideshow:
[slides](vignettes/DiDforBigData_theory.pdf).

# Demonstration

This section will compare the following implementations of DiD
estimators for staggered treatment contexts:

1.  The implementation of the Borusyak, Jaravel, and Spiess (2022)
    approach in R package `didimputation`;
2.  The implementation of the baseline Callaway & Sant’Anna (2021)
    approach in R package `did`;
3.  The implementation of the de Chaisemartin & D’Haultfoeuille (2020)
    approach in R package `DIDmultiplegt` using 40 bootstrap draws for
    standard errors (`brep=40`); and,
4.  My R package `DiDforBigData`.

Regarding Callaway & Sant’Anna (2021), there are many options available
in this package, so I compare three: the default; the estimation method
`est_method = "reg"` (not shown below since it gives nearly identical
results as the default); and the case in which standard errors are
computed analytically rather than by bootstrap (`bstrap=F`).

Below, I draw the simulated data 3 times per sample size, and apply each
estimator. Results are presented for the median across those 3 draws.
Sample Size refers to the number of unique individuals. Since there are
10 simulated years of data, and the sample is balanced across years, the
number of observations is 10 times the number of unique individuals.

A caveat: in the following simulations, I consider contexts in which the
data is balanced, treatment is permanent once it begins, and the
parallel-trends and no-anticipation assumptions hold unconditionally for
all treatment cohorts. Under these assumptions, all of the estimators
considered are consistent. However, in empirical contexts in which these
assumptions do not hold, only some of the estimators may be valid. For
example, only the approach of de Chaisemartin & D’Haultfoeuille (2020)
can allow for treatments that switch off and on repeatedly for the same
individual.

#### Estimates

I verify that all of the estimators provide similar point estimates and
standard errors. Here, I show the point estimates and 95% confidence
intervals (using +/- 1.96\*SE) for the DiD estimate at event time +1
(averaging across cohorts). The true ATT is 4 at event time +1. I also
verify that two-way fixed-effects OLS estimation would find an effect of
about 5.5 at event time +1 when the sample is large.

![](vignettes/estimates_small.png)

#### Speed test

**Small Samples:** Here is the run-time required to complete the DiD
estimation using each package:

![](vignettes/speedtest_small.png)

We see that, with 20,000 unique individuals, the implementations of the
Borusyak, Jaravel, and Spiess (2022) and de Chaisemartin &
D’Haultfoeuille (2020) approaches have become very slow. I could not get
either approach to run successfully with 100,000 unique individuals, as
they both crash R.

**Large samples:** Given the above results, we cannot rely on the
implementations of the Borusyak, Jaravel, and Spiess (2022) and de
Chaisemartin & D’Haultfoeuille (2020) approaches if there are 100,000 or
more individuals. We now compare the implementation of the Callaway &
Sant’Anna (2021) approach to `DiDforBigData` when the sample size is
large.

![](vignettes/speedtest_large.png)

Even with 1 million unique individuals (and 10 million observations), it
is difficult to see `DiDforBigData` in the plot, as estimation requires
about 20 seconds, versus nearly 1 hour for the implementation of the
Callaway & Sant’Anna (2021) approach.

#### Memory test

**Small Samples:** Here is the memory used to complete the DiD
estimation by each package:

![](vignettes/memorytest_small.png)

We see that the implementation of the de Chaisemartin & D’Haultfoeuille
(2020) method uses much more memory than the other approaches. The other
approaches all use relatively little memory at these sample sizes.

**Large Samples:**

![](vignettes/memorytest_large.png)

When considering large samples, we see that `DiDforBigData` uses very
little memory, while the implementation of the Callaway & Sant’Anna
(2021) approach uses quite a lot of memory.

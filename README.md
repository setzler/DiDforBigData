DiD for Big Data in R
================

This R package provides a big-data-friendly and memory-efficient
difference-in-differences (DiD) estimator for staggered (and
non-staggered) treatment contexts.

## Motivation

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

To verify their claim, I wrote a simple panel data simulator with 10
years of data and 4 equally-sized treatment cohorts (including a
never-treated cohort). I simulated it with different numbers of
individuals and applied the available R estimators to it. I found that
most of the estimators cannot run successfully with sample size greater
than 10,000 individuals. The exception is the excellent `did` package by
[Callaway & Sant’Anna (2021)](https://bcallaway11.github.io/did/);
however, this package becomes very slow with sample size above 100,000
and infeasible with sample size above 1 million (see demonstration
below).

Estimating DiD with large administrative data poses three challenges:

1.  **Speed:** If there are millions of observations, DiD is infeasible
    or extremely slow to estimate with available software. We need
    software that can handle large sample sizes.
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
    and staggered treatment in less than 1 minute (see demonstration
    below).
2.  **Memory:** This package avoids memory-intensive activities like
    matrix-inversion and data-stacking. It uses much less memory than
    other packages (see demonstration below).
3.  **Dependencies:** This package has only *one* dependency,
    `data.table`, which is already installed with most R distributions.
    This package is also extremely small, so it can easily be
    transferred onto servers, via email, etc.

## Demonstration

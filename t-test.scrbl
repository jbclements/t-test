#lang scribble/manual

@title{T-Test: two popular functions}
@(author+email "John Clements" "clements@racket-lang.org")

@(require (for-label "main.rkt")
          (for-label math/special-functions))

This package provides both Welch's T-Test and the more well-known
Student's T-Test.

These functions both make Bayesians unhappy, and with some reason.
Specifically, these functions try to answer this question: "if
these two sets are chosen from the same distribution, what is the
likelihood of an outcome this strange or stranger?" The actual
question that people typically want to answer is instead: "given
these two sets of observations, what's the likelihood that the
two sets are chosen from the same distribution?" The problem is
that this Bayesian reversal is possible only with some assumption
about the prior likelihood of the various distributions, which
is generally not known. For this reason, experimenters often
just tell the Bayesians to go away and let them run their t-tests.

Separately, there are problems with the classic "Student's T-Test"
proposed by Gosset. Specifically, its results are valid only if
both samples are known to be drawn from distributions that are
normal, and that have the
same variance. In the case that the two variances are not known
to be equal, the correct choice is apparently to use Welch's
t-test, which does not make this assumption. Unsurprisingly, it
has less power, but apparently not hugely less.

A note about implementation: the interesting part of the computation
is determining the cumulative distribution function (CDF) of the
appropriate distribution. This library simply uses the built-in
"incomplete beta function" @racket[beta-inc]
from the math library to perform this computation. So really,
Neil Toronto did all of the heavy lifting, here.

@defmodule[t-test]
                 
@defproc[(welch-t-test [S1 (Sequenceof Real)]
                       [S2 (Sequenceof Real)]
                       [#:t-statistic t-statistic? Boolean #f])
         Real]{
Given two sequences of reals (lists work fine), return the p-value
associated with a two-tailed t-test. If @racket[t-statistic?] is
true, return instead the t-statistic.

The two-tailed p-value represents the likelihood of an occurrence
this strange or stranger given the null hypothesis that the two
distributions have the same mean.

Let's have an example. Suppose that the there are two sections
of a class, and in the first one, the final scores are given
by s1, and the second section's scores by s2:

@racketblock[
(define s1 '(63.9 92.3 80.9 85.9 86.8 87.6 91.7 84.4))
(define s2'(86.6 91.3 63.7 69.2 78.5 74.0 85.4 89.0))]

We want to test whether there's a significant difference between
the two. First, we choose a p-value, perhaps 0.02, representing
a two percent chance. Then, we run the test:

@codeblock|{
(welch-t-test s1 s2) ; => 0.3630116587607044}|

That is, there's about a 36% chance of something this strange
or stranger occuring. This is way above our threshold of 2%,
so we conclude that there's insufficient evidence to reject
the null hypothesis.
}

@defproc[(student-t-test [S1 (Sequenceof Real)]
                         [S2 (Sequenceof Real)]
                         [#:t-statistic t-statistic? Boolean #f])
         Real]{
Given two sequences of reals (lists work fine), return the p-value
associated with a two-tailed t-test. If @racket[t-statistic?] is
true, return instead the t-statistic.

This function uses the "Student's t-test", rather than Welch's
t-test.
}

@section{Acknowledgments}

Many thanks to Tim Brown for posting his translation of the C
code for Welch's t-test to Rosetta Code.




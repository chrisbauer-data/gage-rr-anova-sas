# ANOVA Project: Crossed Gage Study

**Chris Bauer** — April 27, 2026

*Tools: SAS (PROC MIXED with REML, PROC GLM, PROC UNIVARIATE)*

---

## 1. Introduction

Measurement instruments in production are calibrated in accordance with some accepted
standard, commonly an ANSI/NIST-prepared sample or a "golden sample" used to test for
measurement drift over time. Calibration, however, is only part of the picture when it
comes to acceptance testing or evaluating a measurement process, which consists of
part-to-part variation, between-operator measurement variation, and random error. All
three components determine how capable a measurement process is for quality control or
providing reliable measurement for decision making. The proposed project is a
non-destructive Gage R&R study conducted using the sample Gage data from the Minitab
website.

## 2. Data

The data are measurements on ten parts repeated three times by each of three operators
for a total of 90 total replicates, holding to the recommendations from the American
Society of Quality (at least 10 parts, 2–3 operators and 3–5 repeated measurements)
(McShane-Vaughn, 2023).

The histogram of studentized residuals looks normal, with Anderson-Darling and
Kolmogorov-Smirnov p-values greater than 0.25 and 0.15, respectively (PROC UNIVARIATE).
We fail to reject the null hypothesis that the residuals are normally distributed. The
conditional studentized residuals plot admits possible outliers, but no systematic
spread in the data. Further, Levene's test by part and operator leads to p-values of
0.9074 and 0.9019, respectively, which fail to reject the null hypothesis of equal
variance. There is no indication of correlated measurements or dependence, and the
execution of the DOE requires operators to measure without others and without knowing
the nominal value of the measurement for a part under test. The assumptions of
independence, equal variance, and normality of residuals are reasonable and can be
asserted comfortably in this study.

## 3. Methods

The model is a two-way random effects ANOVA with CRD randomization:

```
y_ijk = μ + α_i + β_j + (αβ)_ij + ε_ijk
```

where μ is the grand mean, α_i and β_j are the treatment effects due to part and
operator respectively, (αβ)_ij is the treatment effect due to the interaction of the
Part and Operator factors, and ε_ijk is the error term.

As all factors are random effects, there are four variance components to the full
model: σ²α, σ²β, σ²αβ, σ²ε. The hypotheses are:

- **H₀:** σ²α = σ²β = σ²αβ = 0
- **Hₐ:** At least one of σ²α, σ²β, σ²αβ is non-zero

PROC MIXED with REML is used to perform the analysis. Of note is the fact that this
method will compute and assign a value of zero to the interaction variance if it would
otherwise return a small negative value. In this case, EMS-corrected F-tests from
PROC MIXED are also shown to provide further context.

In the context of a Gage study, the following components and metrics are used to
evaluate the performance of the measurement system:

| Metric | Formula |
|---|---|
| Reproducibility | σ²α + σ²αβ |
| Repeatability | σ²ε |
| Part-to-Part | σ²β |
| Gage Variance (GRR) | σ²α + σ²αβ + σ²ε |
| % Study Variation | σ_GRR / σ_Total × 100 |
| Number of Distinct Categories (NDC) | √2 × (σβ / σ_GRR) |

The use of standard deviations in the metric is because these metrics come from Six
Sigma methodology. The square root of two factor in the NDC metric comes from signal
processing considerations and is beyond the scope of this report (Montgomery, 2020).
A good gage should have %SV < 10% and NDC ≥ 5 (Montgomery, 2020).

## 4. Results

We reject the null hypothesis that all variances are zero, based on
**F₂₉,₆₀ = 68.91, p < 0.001** at α = 0.05, also noting that not all covariance
parameters are zero in the REML output. REML resulted in a zero value for the
interaction variance, with non-zero values for the operator, part, and residual
components. The operator and part sources are shown to be significant in the
EMS-corrected F-tests, but the operator component has a Wald Z-score > 0.05.

One point to consider is that the Wald test can underestimate significance near a
boundary and is unreliable with small sample sizes, both of which are the case here.
For these reasons, along with the fact that the ANOVA table is used by quality
engineers in practice, we reject the null hypothesis that operator and part variance
is zero.

Another issue worth noting, which is also resolved using a practical argument, is that
REML setting the interaction term to zero means that some between-group variation will
be effectively redistributed. This may change the parameter estimates somewhat, but it
would tend to make the components larger and thus make the gage look *worse* — and in
applications relating to quality control, more conservative tests are generally
preferred. Further, looking at the interaction sum of squares, it is very small
compared to operator, residual, and part-to-part, and as such is very unlikely to
produce a noticeable effect if excluded or redistributed.

The computed Gage metrics were **27.86% study variation** and **4 distinct
categories**. Because study variations are computed using standard deviations (6σ per
class) for any class, the percent study variations will not add together the same way
they do for the variance and Gage percent metrics.

## 5. Conclusion and Discussion

Both computed metrics fall in an "application dependent" regime. Another practical
metric frequently used is percent tolerance, defined as the gage study variation
divided by the process tolerance. Taking all factors together, if this gage were
applied to a low-risk application, or one with very wide tolerance, it may be
acceptable. The NDC below 5 would very likely preclude this from any use in final
testing or within automotive or aerospace applications, but it could be useful for
unit-step testing where throughput needs are high, risks are low, and costs of rework
are acceptable. Another possibility, which occurs fairly often in practice, is that
this gage is the best compromise between cost and throughput, where it is known that
the gage could be better but that the cost of a better gage is either prohibitive due
to budget, throughput, or staffing needs.

It is difficult to diagnose this measurement system in terms of possible improvements
due to the zero interaction term between operator and part, along with the roughly
equal contribution of repeatability and reproducibility to the total Gage R&R. It does
not appear that there is a smoking gun regarding operator training or inconsistency
from the between-group contribution, nor is there any indication of instability in the
measurement itself. An interesting follow-up could be to augment this study with skill
levels of operators, or to repeat it once with the most experienced operators and once
with more average operators, to tease out whether there is an operator skill
component — this study does not imply one, but there are only 2 operator degrees of
freedom. Additionally, if a second, more precise measurement technique were available,
it could be used in an attempt to reduce the variance attributable to repeatability.

In terms of the study itself, an increase in the number of operators could benefit the
use of REML, as more degrees of freedom could tighten the standard error relating to
the Wald Z-test. The traditional ANOVA results from PROC GLM and the covariance
parameter estimates computed with REML under PROC MIXED were both in agreement as to
the interaction variance being zero, but differed on the statistical significance of
the operator variance (though the associated confidence interval from REML did not
include zero). As the traditional ANOVA method is used in practice for Gage studies in
quality engineering, we elected to keep the operator variance in the model and use it
going forward for computation of other quantities.

## 6. References

1. Montgomery, D. C. (2020). *Introduction to Statistical Quality Control* (8th ed.).
   John Wiley & Sons.
2. McShane-Vaughn, M. (2023). *The ASQ Certified Six Sigma Black Belt Handbook*
   (4th ed.). ASQ Quality Press.
3. Data Set: [Minitab — Crossed Gage R&R Study Example](https://support.minitab.com/en-us/minitab/help-and-how-to/quality-and-process-improvement/measurement-system-analysis/how-to/gage-study/crossed-gage-r-r-study/before-you-start/example/)

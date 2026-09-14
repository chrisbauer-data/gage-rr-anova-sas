# Gage R&R Study — Crossed ANOVA with REML

A measurement systems analysis (Gage R&R) demonstrating two-way random effects ANOVA with REML variance component estimation — evaluating whether a measurement instrument itself is trustworthy before trusting the data it produces.

What this project does
Analyzes a crossed Gage R&R design: 10 parts × 3 operators × 3 repeated measurements (90 total replicates), following ASQ-recommended study structure
Fits a two-way random effects ANOVA model via PROC MIXED with REML, decomposing total measurement variation into part-to-part, operator (reproducibility), and repeatability (residual) components
Validates model assumptions (normality of residuals, equal variance via Levene's test) before trusting the results
Navigates a real statistical subtlety: REML setting a near-zero interaction variance to exactly zero, and an unreliable Wald Z-test near a variance boundary with small sample size — resolved with an explicit, defensible practical argument rather than blind reliance on either test
Computes standard Gage R&R metrics (%Study Variation, Number of Distinct Categories) and interprets them against industry acceptance thresholds
Key result

27.86% study variation, 4 distinct categories — both in the "application dependent" acceptability zone. The gage is unsuitable for high-precision applications (automotive, aerospace) but may be acceptable for lower-risk, higher-throughput contexts, or represents a defensible cost/throughput tradeoff common in practice.

Why this matters

Point-and-click statistical tools will produce a clean, professional-looking chart or ANOVA table regardless of whether the underlying assumptions actually hold. This project demonstrates the discipline of checking those assumptions explicitly, recognizing when a standard test (Wald Z) is unreliable in a specific context (small sample, near a boundary), and making a principled, stated judgment call rather than either blindly trusting or blindly overriding the statistical output.

Tools

SAS — PROC MIXED (REML), PROC GLM, PROC UNIVARIATE

Background

This project demonstrates Gage R&R methodology and REML-based variance component analysis using a standard, publicly available textbook dataset — Minitab's published crossed Gage R&R example (widely used for teaching measurement systems analysis, not data collected by the author). The value here is in the analytical methodology: the custom EMS-corrected ANOVA table (built from PROC GLM output via ODS OUTPUT and macro variables — not a canned Minitab output), the REML variance component estimation, the assumption validation (three separate Levene's tests, full residual diagnostics), and the judgment calls required to interpret results correctly (handling a REML-zeroed interaction variance, recognizing an unreliable Wald Z-test near a boundary with a small sample).

Files
Gage_RR_Study.md — full write-up: methodology, results, and discussion
Gage_RR_Analysis.sas — full SAS source code (PROC GLM, PROC MIXED/REML, Levene's tests, residual diagnostics, Gage R&R metric calculations)

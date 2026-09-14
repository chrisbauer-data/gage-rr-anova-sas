/* ============================================================
   Crossed Gage R&R
   PROC GLM  -> ANOVA table (matches Minitab)
   PROC MIXED REML -> Variance components and Gage R&R metrics
   ============================================================ */

/* Step 1: Input the data */
DATA grr;
   INPUT Part Operator $ Measurement;
   DATALINES;
1  A  0.29
1  A  0.41
1  A  0.64
2  A -0.56
2  A -0.68
2  A -0.58
3  A  1.34
3  A  1.17
3  A  1.27
4  A  0.47
4  A  0.50
4  A  0.64
5  A -0.80
5  A -0.92
5  A -0.84
6  A  0.02
6  A -0.11
6  A -0.21
7  A  0.59
7  A  0.75
7  A  0.66
8  A -0.31
8  A -0.20
8  A -0.17
9  A  2.26
9  A  1.99
9  A  2.01
10 A -1.36
10 A -1.25
10 A -1.31
1  B  0.08
1  B  0.25
1  B  0.07
2  B -0.47
2  B -1.22
2  B -0.68
3  B  1.19
3  B  0.94
3  B  1.34
4  B  0.01
4  B  1.03
4  B  0.20
5  B -0.56
5  B -1.20
5  B -1.28
6  B -0.20
6  B  0.22
6  B  0.06
7  B  0.47
7  B  0.55
7  B  0.83
8  B -0.63
8  B  0.08
8  B -0.34
9  B  1.80
9  B  2.12
9  B  2.19
10 B -1.68
10 B -1.62
10 B -1.50
1  C  0.04
1  C -0.11
1  C -0.15
2  C -1.38
2  C -1.13
2  C -0.96
3  C  0.88
3  C  1.09
3  C  0.67
4  C  0.14
4  C  0.20
4  C  0.11
5  C -1.46
5  C -1.07
5  C -1.45
6  C -0.29
6  C -0.67
6  C -0.49
7  C  0.02
7  C  0.01
7  C  0.21
8  C -0.46
8  C -0.56
8  C -0.49
9  C  1.77
9  C  1.45
9  C  1.87
10 C -1.49
10 C -1.77
10 C -2.16
;
RUN;

/* ============================================================
   Step 2: PROC GLM - Extract SS and MS directly from data
   ============================================================ */
ODS TRACE ON;
ODS OUTPUT ModelANOVA  = glm_model
           OverallANOVA = glm_error;

PROC GLM DATA=grr;
   CLASS Part Operator;
   MODEL Measurement = Part Operator Part*Operator;
RUN;
QUIT;

ODS OUTPUT CLOSE;
ODS TRACE OFF;

/* Pull MS values into macro variables */
DATA _NULL_;
   SET glm_model;
   IF Source = 'Part'          THEN CALL SYMPUTX('MS_part',  MS);
   IF Source = 'Operator'      THEN CALL SYMPUTX('MS_oper',  MS);
   IF Source = 'Part*Operator' THEN CALL SYMPUTX('MS_inter', MS);
   IF Source = 'Part'          THEN CALL SYMPUTX('SS_part',  SS);
   IF Source = 'Operator'      THEN CALL SYMPUTX('SS_oper',  SS);
   IF Source = 'Part*Operator' THEN CALL SYMPUTX('SS_inter', SS);
   IF Source = 'Part'          THEN CALL SYMPUTX('DF_part',  DF);
   IF Source = 'Operator'      THEN CALL SYMPUTX('DF_oper',  DF);
   IF Source = 'Part*Operator' THEN CALL SYMPUTX('DF_inter', DF);
RUN;

DATA _NULL_;
   SET glm_error;
   CALL SYMPUTX('MS_error', MS);
   CALL SYMPUTX('SS_error', SS);
   CALL SYMPUTX('DF_error', DF);
RUN;

/* ============================================================
   Step 3: Build ANOVA table with EMS-correct F-tests
   ============================================================ */
DATA anova_print;
   LENGTH Source $25;

   MS_part  = INPUT(SYMGET('MS_part'),  BEST12.);
   MS_oper  = INPUT(SYMGET('MS_oper'),  BEST12.);
   MS_inter = INPUT(SYMGET('MS_inter'), BEST12.);
   MS_error = INPUT(SYMGET('MS_error'), BEST12.);
   SS_part  = INPUT(SYMGET('SS_part'),  BEST12.);
   SS_oper  = INPUT(SYMGET('SS_oper'),  BEST12.);
   SS_inter = INPUT(SYMGET('SS_inter'), BEST12.);
   SS_error = INPUT(SYMGET('SS_error'), BEST12.);
   DF_part  = INPUT(SYMGET('DF_part'),  BEST12.);
   DF_oper  = INPUT(SYMGET('DF_oper'),  BEST12.);
   DF_inter = INPUT(SYMGET('DF_inter'), BEST12.);
   DF_error = INPUT(SYMGET('DF_error'), BEST12.);
   DF_total = DF_part + DF_oper + DF_inter + DF_error;
   SS_total = SS_part + SS_oper + SS_inter + SS_error;

   /* EMS-correct F-tests */
   F_part  = MS_part  / MS_inter;
   F_oper  = MS_oper  / MS_inter;
   F_inter = MS_inter / MS_error;

   P_part  = 1 - PROBF(F_part,  DF_part,  DF_inter);
   P_oper  = 1 - PROBF(F_oper,  DF_oper,  DF_inter);
   P_inter = 1 - PROBF(F_inter, DF_inter, DF_error);

   Source="Part";           DF=DF_part;  SS=SS_part;  MS=MS_part;  F=F_part;  ProbF=P_part;  OUTPUT;
   Source="Operator";       DF=DF_oper;  SS=SS_oper;  MS=MS_oper;  F=F_oper;  ProbF=P_oper;  OUTPUT;
   Source="Part*Operator";  DF=DF_inter; SS=SS_inter; MS=MS_inter; F=F_inter; ProbF=P_inter; OUTPUT;
   Source="Error (Repeat)"; DF=DF_error; SS=SS_error; MS=MS_error; F=.;       ProbF=.;       OUTPUT;
   Source="Total";          DF=DF_total; SS=SS_total; MS=.;        F=.;       ProbF=.;       OUTPUT;

   KEEP Source DF SS MS F ProbF;
RUN;

PROC PRINT DATA=anova_print NOOBS LABEL;
   VAR Source DF SS MS F ProbF;
   FORMAT SS MS 10.5   F 8.3   ProbF 6.4;
   LABEL
      Source = "Source"
      DF     = "DF"
      SS     = "Sum of Squares"
      MS     = "Mean Square"
      F      = "F Value"
      ProbF  = "Pr > F";
   TITLE "Two-Way ANOVA Table (EMS-Corrected F-Tests)";
RUN;

/* ============================================================
   Step 4: PROC MIXED - REML Variance Components
   ============================================================ */
PROC MIXED DATA=grr METHOD=REML COVTEST CL
   PLOTS(ONLY)=(
      RESIDUALPANEL
      STUDENTPANEL
      PEARSONPANEL
   );
   CLASS Part Operator;
   MODEL Measurement = / RESIDUAL OUTPM=mixed_resid;
   RANDOM Part Operator Part*Operator;
   ODS OUTPUT CovParms=CovEst;
   TITLE "Crossed Gage R&R - PROC MIXED (REML)";
RUN;

/* Levene's test by Operator */
PROC GLM DATA=mixed_resid;
   CLASS Operator;
   MODEL Resid = Operator;
   MEANS Operator / HOVTEST=LEVENE;
   TITLE "Levene's Test - Homogeneity of Variance by Operator";
RUN;
QUIT;

/* Levene's test by Part */
PROC GLM DATA=mixed_resid;
   CLASS Part;
   MODEL Resid = Part;
   MEANS Part / HOVTEST=LEVENE;
   TITLE "Levene's Test - Homogeneity of Variance by Part";
RUN;
QUIT;

/* Levene's test by Operator*Part combination */
PROC GLM DATA=mixed_resid;
   CLASS Part Operator;
   MODEL Resid = Part Operator;
   MEANS Part Operator / HOVTEST=LEVENE;
   TITLE "Levene's Test - Homogeneity of Variance by Part and Operator";
RUN;
QUIT;

PROC SGPLOT DATA=mixed_resid;
   SCATTER X=Part Y=Resid;
   REFLINE 0 / AXIS=Y LINEATTRS=(COLOR=red);
   XAXIS LABEL="Part";
   YAXIS LABEL="Residuals";
   TITLE "MIXED Residuals vs Part";
RUN;

/* Residuals vs Operator */
PROC SGPLOT DATA=mixed_resid;
   SCATTER X=Operator Y=Resid;
   REFLINE 0 / AXIS=Y LINEATTRS=(COLOR=red);
   XAXIS LABEL="Operator";
   YAXIS LABEL="Residuals";
   TITLE "MIXED Residuals vs Operator";
RUN;

/* Boxplot by Operator */
PROC SGPLOT DATA=mixed_resid;
   VBOX Resid / CATEGORY=Operator;
   REFLINE 0 / AXIS=Y LINEATTRS=(COLOR=red);
   XAXIS LABEL="Operator";
   YAXIS LABEL="Residuals";
   TITLE "MIXED Residual Boxplot by Operator";
RUN;

/* Boxplot by Part */
PROC SGPLOT DATA=mixed_resid;
   VBOX Resid / CATEGORY=Part;
   REFLINE 0 / AXIS=Y LINEATTRS=(COLOR=red);
   XAXIS LABEL="Part";
   YAXIS LABEL="Residuals";
   TITLE "MIXED Residual Boxplot by Part";
RUN;

/* Run order plot */
DATA mixed_resid;
   SET mixed_resid;
   obs = _N_;
RUN;

PROC SGPLOT DATA=mixed_resid;
   SERIES  X=obs Y=Resid;
   SCATTER X=obs Y=Resid;
   REFLINE 0 / AXIS=Y LINEATTRS=(COLOR=red);
   XAXIS LABEL="Observation Order";
   YAXIS LABEL="Residuals";
   TITLE "MIXED Residuals vs Run Order";
RUN;

/* Pull REML variance components into macro variables */
DATA _NULL_;
   SET CovEst;
   IF CovParm = 'Part'          THEN CALL SYMPUTX('v_part',  MAX(Estimate, 0));
   IF CovParm = 'Operator'      THEN CALL SYMPUTX('v_oper',  MAX(Estimate, 0));
   IF CovParm = 'Part*Operator' THEN CALL SYMPUTX('v_inter', MAX(Estimate, 0));
   IF CovParm = 'Residual'      THEN CALL SYMPUTX('v_error', MAX(Estimate, 0));
RUN;

/* ============================================================
   Step 5: Gage R&R Metrics
   ============================================================ */
DATA grr_summary;
   v_part  = INPUT(SYMGET('v_part'),  BEST12.);
   v_oper  = INPUT(SYMGET('v_oper'),  BEST12.);
   v_inter = INPUT(SYMGET('v_inter'), BEST12.);
   v_error = INPUT(SYMGET('v_error'), BEST12.);

   /* Variance components */
   v_reprod = v_oper + v_inter;
   v_grr    = v_error + v_reprod;
   v_total  = v_part + v_grr;

   /* Standard deviations */
   sd_repeat = SQRT(v_error);
   sd_oper   = SQRT(v_oper);
   sd_inter  = SQRT(v_inter);
   sd_reprod = SQRT(v_reprod);
   sd_grr    = SQRT(v_grr);
   sd_part   = SQRT(v_part);
   sd_total  = SQRT(v_total);

   /* % Contribution */
   pct_repeat = 100 * v_error  / v_total;
   pct_oper   = 100 * v_oper   / v_total;
   pct_inter  = 100 * v_inter  / v_total;
   pct_reprod = 100 * v_reprod / v_total;
   pct_grr    = 100 * v_grr    / v_total;
   pct_part   = 100 * v_part   / v_total;

   /* % Study Variation */
   psv_repeat = 100 * sd_repeat / sd_total;
   psv_oper   = 100 * sd_oper   / sd_total;
   psv_inter  = 100 * sd_inter  / sd_total;
   psv_reprod = 100 * sd_reprod / sd_total;
   psv_grr    = 100 * sd_grr    / sd_total;
   psv_part   = 100 * sd_part   / sd_total;

   /* Number of Distinct Categories */
   ndc = FLOOR(1.41 * (sd_part / sd_grr));
RUN;

/* Variance Components */
PROC PRINT DATA=grr_summary NOOBS LABEL;
   VAR v_error v_oper v_inter v_reprod v_grr v_part v_total;
   FORMAT v_error--v_total 10.5;
   LABEL
      v_error  = "Variance: Repeatability"
      v_oper   = "Variance: Operator"
      v_inter  = "Variance: Interaction"
      v_reprod = "Variance: Reproducibility"
      v_grr    = "Variance: Gauge R&R"
      v_part   = "Variance: Part-to-Part"
      v_total  = "Variance: Total";
   TITLE "Gage R&R Variance Components (REML)";
RUN;

/* % Contribution */
PROC PRINT DATA=grr_summary NOOBS LABEL;
   VAR pct_repeat pct_oper pct_inter pct_reprod pct_grr pct_part;
   FORMAT pct_repeat--pct_part 6.2;
   LABEL
      pct_repeat = "% Contribution: Repeatability"
      pct_oper   = "% Contribution: Operator"
      pct_inter  = "% Contribution: Interaction"
      pct_reprod = "% Contribution: Reproducibility"
      pct_grr    = "% Contribution: Gauge R&R"
      pct_part   = "% Contribution: Part-to-Part";
   TITLE "Gage R&R % Contribution (REML)";
RUN;

/* % Study Variation and NDC */
PROC PRINT DATA=grr_summary NOOBS LABEL;
   VAR psv_repeat psv_oper psv_inter psv_reprod psv_grr psv_part ndc;
   FORMAT psv_repeat--psv_part 6.2;
   LABEL
      psv_repeat = "% Study Var: Repeatability"
      psv_oper   = "% Study Var: Operator"
      psv_inter  = "% Study Var: Interaction"
      psv_reprod = "% Study Var: Reproducibility"
      psv_grr    = "% Study Var: Gauge R&R"
      psv_part   = "% Study Var: Part-to-Part"
      ndc        = "Number of Distinct Categories";
   TITLE "Gage R&R % Study Variation and NDC (REML)";
RUN;

PROC UNIVARIATE DATA=mixed_resid NORMAL;
   VAR Resid;
   HISTOGRAM Resid / NORMAL;
   PROBPLOT  Resid / NORMAL(MU=EST SIGMA=EST);
   TITLE "Normality Tests on PROC MIXED Residuals";
RUN;

/* Traditional ANOVA (no interaction term) for comparison */
ODS TRACE ON;
ODS OUTPUT ModelANOVA  = glm_model
           OverallANOVA = glm_error;

PROC GLM DATA=grr;
   CLASS Part Operator;
   MODEL Measurement = Part Operator;
RUN;
QUIT;

ODS OUTPUT CLOSE;
ODS TRACE OFF;

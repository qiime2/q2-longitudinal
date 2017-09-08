# q2-longitudinal

[![Build Status](https://travis-ci.org/qiime2/q2-longitudinal.svg?branch=master)](https://travis-ci.org/qiime2/q2-longitudinal) [![Coverage Status](https://coveralls.io/repos/github/qiime2/q2-longitudinal/badge.svg?branch=master)](https://coveralls.io/github/qiime2/q2-longitudinal?branch=master)

QIIME2 plugin for paired sample comparisons.

q2-longitudinal's actions support statistical and visual comparisons of longitudinal study designs and paired samples, to determine if/how samples change between observational "states". "States" will most commonly be related to time or an environmental gradient, and for paired analyses (`pairwise-distances` and `pairwise-differences`) the sample pairs should typically consist of the same individual subject  observed at two different time points. For example, patients in a clinical study whose stool samples are collected before and after receiving treatment.

"States" can also commonly be methodological, in which case sample pairs will usually be the same individual at the same time with two different methods. For example, q2-feature-classifier could compare the effects of different collection methods, storage methods, DNA extraction methods, or any bioinformatic processing steps on the feature composition of individual samples.

In the examples below, we use data from a longitudinal study of infants' and mothers' microbiota from birth through 2 years of life ([doi: 10.1126/scitranslmed.aad7121](http://stm.sciencemag.org/content/8/343/343ra82)). These data are located in the directory [`tutorial_data/`](./tutorial_data).

## Examples

### Paired difference testing

Paired difference tests determine whether the value of a specific metric changed significantly between pairs of paired samples (e.g., pre- and post-treatment).

This visualizer currently supports comparison of feature abundance (e.g., microbial sequence variants or taxa) in a feature table, or of metadata values in a sample metadata file. Alpha diversity values (e.g., observed sequence variants) and beta diversity values (e.g., principal coordinates) are useful metrics for comparison with these tests, and should be contained in one of the metadata files given as input.

#### Paired differences in metadata

Here we use `pairwise-differences` to assess whether alpha diversity (here, Shannon's diversity index) changed significantly between 0 and 12 months of life in vaginally born and Cesarean-delivered infants, and whether the magnitude of change differed between these groups. Note that the alpha diversity data in this case is contained in a separate artifact, which is the typical (and preferred) approach; alternatively, alpha diversity data (or other data) contained in a sample metadata file can be used as the input `metric`, in which case the second `--m-metadata-file` input is unnecessary.

```
qiime longitudinal pairwise-differences \
    --m-metadata-file ecam_map_maturity.txt \
    --m-metadata-file ecam_shannon.qza \
    --p-metric shannon \
    --p-group-column delivery \
    --p-state-column month \
    --p-state-1 0 \
    --p-state-2 12 \
    --p-individual-id-column studyid \
    --o-visualization ecam-delivery-alpha \
    --p-replicate-handling random
```

The visualization contains boxplots comparing the distributions of pairwise differences in `metric` for each individual in each group between `state-1` and `state-2`. "Multiple group tests" lists the results of a Kruskal Wallis test or one-way ANOVA test (depending on whether the `parametric` parameter is set) comparing all groups. "Pairwise group comparison tests" lists the results of a series of two-group Mann-Whitney U tests or t-tests (depending on whether the `parametric` parameter is set) between each pair of groups. "Multiple group tests" and "pairwise group comparison tests" are both calculated whether two or more groups are actually compared.

#### Pairwise differences in feature table

We can also use this method to measure changes in the abundances of specific features of interest. In this example, we test whether the abundance of genus Bacteroides changed significantly between 6 and 18 months of life in vaginally born and Cesarean-delivered infants, and whether the magnitude of change differed between these groups. Note that `pairwise-differences` accepts a feature table as optional input to extract taxon abundance data.

```
qiime longitudinal pairwise-differences \
    --i-table ecam-table-taxa.qza \
    --m-metadata-file ecam_map_maturity.txt \
    --p-metric 'k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__' \
    --p-group-column delivery \
    --p-state-column month \
    --p-state-1 6 \
    --p-state-2 18 \
    --p-individual-id-column studyid \
    --o-visualization ecam-delivery \
    --p-replicate-handling random
```

### Pairwise distance testing

The `pairwise-distancess` visualizer also assesses changes between paired samples from two different "states", but instead of taking a metadata column or feature ID as input, it operates on a distance matrix to assess the distance between "pre" and "post" sample pairs. The "within-subject" distance between paired samples from an individual are always calculated for each group in the metadata `group_column`; by default, "between-subject" distances between all individuals in a given `group_column` are also calculated and compared. Between-subject distances include all samples sharing the same `group_column` that are not pairs of "within-subject" samples from `state_1` and `state_2`, but otherwise ignore the `state_column` and `individual_id_column` parameters, so will pair all samples from all time points (or whatever the comparison "state" is) in the distance matrix. Hence, users should carefully consider what type of comparison they wish to perform and, if appropriate, filter the distance matrix prior to using this visualizer. Filtering can be performed with `filter-distance-matrix` as described [here](https://docs.qiime2.org/2017.7/tutorials/filtering/#filtering-distance-matrices).

In this example, we test whether an individual's stool microbiota (as assessed by unweighted UniFrac distance) differs significantly between 0 and 12 months of life in vaginally born and Cesarean-delivered infants. 
```
qiime longitudinal pairwise-distances \
    --i-distance-matrix ecam-unweighted-distance-matrix.qza \
    --m-metadata-file ecam_map_maturity.txt \
    --p-group-column delivery \
    --p-state-column month \
    --p-state-1 0 \
    --p-state-2 12 \
    --p-individual-id-column studyid \
    --o-visualization ecam-delivery-distance \
    --p-replicate-handling random
```

The visualization contains boxplots comparing the pairwise distance distributions for each individual in each group between `state-1` and `state-2`. "Multiple group tests" lists the results of a Kruskal Wallis test or one-way ANOVA test (depending on whether the `parametric` parameter is set) comparing all groups. "Pairwise group comparison tests" lists the results of a series of two-group Mann-Whitney U tests or t-tests (depending on whether the `parametric` parameter is set) between each pair of groups. "Multiple group tests" and "pairwise group comparison tests" are both calculated whether two or more groups are actually compared.

### Linear mixed effects models

Linear mixed effects (LME) models test the relationship between a single response variable and one or more independent variables, where observations are made across dependent samples, e.g., in repeated-measures sampling experiments. This implementation takes at least one numeric "state_column" (e.g., Time) and one or more comma-separated group_categories (which may be categorical or numeric) as independent variables in a LME model, and plots regression plots of the response variable ("metric") as a function of the state caregory and each group column. The response variable may either be a sample metadata mapping file column or a feature ID in the feature table. LME is useful for testing what factors affect the trajectory of change in the dependent variable across time (or a gradient).

In this example, we demonstrate the use of `linear-mixed-effects` to test the relationship between `shannon` (alpha diversity), age, delivery mode, diet, and sex.

```
qiime longitudinal linear-mixed-effects \
    --m-metadata-file ecam_map_maturity.txt \
    --p-metric shannon \
    --m-metadata-file ecam_shannon.qza \
    --p-group-categories delivery,diet,sex \
    --p-state-column month \
    --p-individual-id-column studyid \
    --o-visualization ecam-lme
```

The visualizer produced by this command contains several results. First, the input parameters are shown at the top of the visualization for convenience (e.g., when flipping through multiple visualizations it is useful to have a summary). Scatter plots categorized by each "group column" are shown, with linear regression lines (plus 95% confidence interval in grey) for each group. If `--p-lowess` is enabled, instead locally weighted averages are shown for each group. Next, the "model summary" shows some descriptive information about the LME model that was trained. This just shows descriptive information about the "groups"; in this case, groups will be individuals (as set by the `--p-individual-id-column`). The main results to examine will be the "model results" at the bottom of the visualization. These results summarize the effects of each fixed effect (and their interactions) on the dependent variable (shannon diversity). This table shows parameter estimates, estimate standard errors, Wald Z test statistics, P values (P>|z|), and 95% confidence intervals upper and lower bounds for each parameter. We see in this table that shannon diversity is significantly impacted by month of life and by diet, as well as several interacting factors. More information about LME models and the interpretation of these data can be found on the [statsmodels LME description page](http://www.statsmodels.org/dev/mixed_linear.html), which provides a number of useful technical references for further reading.

Second, we demonstrate the use of `linear-mixed-effects` to test the relationship between `Bacteroides`, age, delivery mode, diet, and sex.

```
qiime longitudinal linear-mixed-effects \
    --i-table ecam-table-taxa.qza \
    --m-metadata-file ecam_map_maturity.txt \
    --p-metric 'k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__' \
    --p-group-categories delivery,diet,sex \
    --p-state-column month \
    --p-individual-id-column studyid \
    --o-visualization ecam-bacteroides-lme
```

### Volatility analysis

Volatility analysis is a method for assessing how volatile a dependent variable is over time (or a gradient) in one or more groups. Whereas LME tests how different factors (e.g., treatment, subject characteristics) impact the trajectory of change in a dependent variable across time (or a gradient), `volatility` examines how variance changes over time (or a gradient), and compares variances between groups at each state, between baseline and each state, and between groups across all states globally. Any metadata (including alpha and beta diversity artifacts) or `FeatureTable[RelativeFrequency]` feature can be used as the dependent variable ("metric"). *Note that volatility analysis currently works best for comparing groups that are sampled fairly evenly at each state in the dataset. Datasets that contain groups sampled at different states are not supported, and users should either filter out those samples or "bin" them with other groups prior to using this visualizer.*

Here we examine how variance in Shannon diversity changes across time in the ECAM cohort.
```
qiime longitudinal volatility \
    --m-metadata-file ecam_map_maturity.txt \
    --m-metadata-file ecam_shannon.qza \
    --p-metric shannon \
    --p-group-column delivery \
    --p-state-column month \
    --p-individual-id-column studyid \
    --o-visualization ecam-volatility
```
The resulting visualization contains several resuls. First, the "Volatility test parameters" summarizes some key parameters, as well as global mean and standard deviation — these are measured across all samples, regardless of which group they are in.

Second, control charts display the mean value of "metric" at each "state". The first plot shown contains all samples, categorized by group (as defined by `group-column`); the following plots show each individual group, in order to show their individual control characteristics as described in the rest of this paragraph. The mean for all samples in each plot is shown as a black line. The "control limits", 3 standard deviations above and below the mean, are shown as dashed lines. The "warning limits", 2 standard deviations above and below the mean, are shown as dotted lines. The idea behind this plot is to show how a variable is changing over time (or a gradient) in relation to the mean. Large departures from the mean values can cross the warning/control limits, indicating a major disruption at that state; for example, antibiotic use or other disturbances impacting diversity could be tracked with these plots.

Finally, we are interested in seeing how variance changes over time (or a gradient). For example, we may be interested in testing whether a group of patients becomes more variable with time, following a treatment, or in relation to a control group. This pairs well with LME tests, which evaluate the impact of independent variables on the trajectory of change in a dependent variable, but here we are interested in how a single factor impacts the spread of the data. We assess this using variance tests (users have the choice of Bartlett's, [Levene](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.levene.html), or Fligner-Killeen tests), with the null hypothesis that all distributions being compared have equal variance. Equality of variances is tested between 1) all groups, regardless of state (time point), 2) all groups at each individual state (time point), and 3) between baseline and each state for each group. The row "All states: compare groups" contains the results of test 1. Rows labeled "State s: compare groups" contain the results of test type 2, where "s" indicates the state (e.g., time point) at which groups are compared. Rows labeled "group: baseline vs. state" contain the results of test type 3, where "group" indicates the individual group being tested, "baseline" indicates the baseline state (this defaults to "0", but any value in `state_column` can be used as `baseline`), and "state" indicates the state being compared to baseline. In this example, we see that vaginally delivered and cesarean-delivery infants exhibit significantly different variances across all time (vaginally born infants are more variable but only because they accrue higher species diversity, so this result is not necessarily interesting); that vaginally and cesarean-born infants typically do not exhibit significantly different variances at isolated time points (indicating similar degrees of volatility over time; this is a more important statistic, given the developmental trajectory expected during early life); and that both groups' variances differ significantly from baseline at many times (probably not an important or interesting finding, given the developmental trajectory espected during early life; this particular test is more interesting when assessing, e.g., how a stable system responds to perturbation at a particular time point).
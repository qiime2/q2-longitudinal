# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import (Str, Bool, Plugin, Metadata, Choices, Range, Float,
                           Int, Citations, Visualization)
from q2_types.feature_table import FeatureTable, RelativeFrequency, Frequency
from q2_types.distance_matrix import DistanceMatrix
from q2_types.sample_data import SampleData
from q2_types.feature_data import FeatureData

from q2_sample_classifier import (Importance, RegressorPredictions,
                                  SampleEstimator, Regressor)
from q2_sample_classifier.plugin_setup import (
    parameters, parameter_descriptions, output_descriptions,
    pipeline_parameters, pipeline_parameter_descriptions,
    regressor_pipeline_outputs, pipeline_output_descriptions, regressors)


from ._type import FirstDifferences
from ._format import FirstDifferencesFormat, FirstDifferencesDirectoryFormat
from ._longitudinal import (pairwise_differences, pairwise_distances,
                            linear_mixed_effects, volatility, nmit,
                            first_differences, first_distances,
                            maturity_index, feature_volatility,
                            plot_feature_volatility, anova)
import q2_longitudinal


citations = Citations.load('citations.bib', package='q2_longitudinal')

plugin = Plugin(
    name='longitudinal',
    version=q2_longitudinal.__version__,
    website="https://github.com/qiime2/q2-longitudinal",
    package='q2_longitudinal',
    description=(
        'This QIIME 2 plugin supports methods for analysis of time series '
        'data, involving either paired sample comparisons or longitudinal '
        'study designs.'),
    short_description='Plugin for paired sample and time series analyses.',
    citations=[citations['bokulich2017q2']]
)

plugin.register_semantic_types(FirstDifferences)

plugin.register_formats(
    FirstDifferencesFormat, FirstDifferencesDirectoryFormat)

plugin.register_semantic_type_to_format(
    SampleData[FirstDifferences],
    artifact_format=FirstDifferencesDirectoryFormat)

miscellaneous_parameters = {
    'state_column': Str,
    'replicate_handling': Str % Choices(['error', 'random', 'drop'])
}

shared_parameters = {
    'metadata': Metadata,
    'individual_id_column': Str,
}

base_parameters = {
    **shared_parameters,
    'state_column': miscellaneous_parameters['state_column'],
    'palette': Str % Choices([
        'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Paired', 'Accent',
        'Dark2', 'tab10', 'tab20', 'tab20b', 'tab20c', 'viridis', 'plasma',
        'inferno', 'magma', 'terrain', 'rainbow', 'cividis']),
}

paired_params = {
    **base_parameters,
    'group_column': Str,
    'state_1': Str,
    'state_2': Str,
    'parametric': Bool,
    'replicate_handling': miscellaneous_parameters['replicate_handling'],
}

miscellaneous_parameter_descriptions = {
    'state_column': ('Metadata column containing state (time) variable '
                     'information.'),
    'replicate_handling': (
        'Choose how replicate samples are handled. If replicates are '
        'detected, "error" causes method to fail; "drop" will discard all '
        'replicated samples; "random" chooses one representative at random '
        'from among replicates.')
}

shared_parameter_descriptions = {
        'metadata': (
            'Sample metadata file containing individual_id_column.'),
        'individual_id_column': (
            'Metadata column containing IDs for individual subjects.'),
}

base_parameter_descriptions = {
        **shared_parameter_descriptions,
        'palette': 'Color palette to use for generating boxplots.',
}

paired_parameter_descriptions = {
        **base_parameter_descriptions,
        'group_column': (
            'Metadata column on which to separate groups for comparison'),
        'individual_id_column': (
            'Metadata column containing subject IDs to use for pairing '
            'samples. WARNING: if replicates exist for an individual ID at '
            'either state_1 or state_2, that subject will be dropped and '
            'reported in standard output by default. Set '
            'replicate_handling="random" to instead randomly select one '
            'member.'),
        'state_column': ('Metadata column containing state (e.g., Time) '
                         'across which samples are paired.'),
        'state_1': 'Baseline state column value.',
        'state_2': 'State column value to pair with baseline.',
        'parametric': ('Perform parametric (ANOVA and t-tests) or non-'
                       'parametric (Kruskal-Wallis, Wilcoxon, and Mann-'
                       'Whitney U tests) statistical tests.'),
        'replicate_handling': (
            miscellaneous_parameter_descriptions['replicate_handling']),
}

volatility_filtering_parameters = {
    'feature_count': Int % Range(1, None) | Str % Choices(['all']),
    'importance_threshold': (Float % Range(0, None, inclusive_start=False) |
                             Str % Choices(['q1', 'q2', 'q3']))}

volatility_filtering_parameter_descriptions = {
    'feature_count': 'Filter feature table to include top N most '
                     'important features. Set to "all" to include all '
                     'features.',
    'importance_threshold': 'Filter feature table to exclude any features '
                            'with an importance score less than this '
                            'threshold. Set to "q1", "q2", or "q3" to select '
                            'the first, second, or third quartile of values. '
                            'Set to "None" to disable this filter.'}

formula_description = (
    'Formulae will be in the format "a ~ b + c", where '
    '"a" is the metric (dependent variable) and "b" and "c" are '
    'independent covariates. Use "+" to add a variable; "+ a:b" to '
    'add an interaction between variables a and b; "*" to include '
    'a variable and all interactions; and "-" to subtract a '
    'particular term (e.g., an interaction term). See '
    'https://patsy.readthedocs.io/en/latest/formulas.html for full '
    'documentation of valid formula operators. Always enclose formulae in '
    'quotes to avoid unpleasant surprises.')


plugin.visualizers.register_function(
    function=pairwise_differences,
    inputs={'table': FeatureTable[RelativeFrequency]},
    parameters={**paired_params,
                'metric': Str},
    input_descriptions={'table': (
        'Feature table to optionally use for paired comparisons.')},
    parameter_descriptions={
        **paired_parameter_descriptions,
        'metric': 'Numerical metadata or artifact column to test.',
    },
    name='Paired difference testing and boxplots',
    description=(
        'Performs paired difference testing between samples from each '
        'subject. Sample pairs may represent a typical intervention study '
        '(e.g., samples collected pre- and post-treatment), paired samples '
        'from two different timepoints (e.g., in a longitudinal study '
        'design), or identical samples receiving different treatments. This '
        'action tests whether the change in a numeric metadata value "metric" '
        'differs from zero and differs between groups (e.g., groups of '
        'subjects receiving different treatments), and produces boxplots of '
        'paired difference distributions for each group. Note that "metric" '
        'can be derived from a feature table or metadata.')
)


plugin.visualizers.register_function(
    function=pairwise_distances,
    inputs={'distance_matrix': DistanceMatrix},
    parameters=paired_params,
    input_descriptions={'distance_matrix': (
        'Matrix of distances between pairs of samples.')},
    parameter_descriptions=paired_parameter_descriptions,
    name='Paired pairwise distance testing and boxplots',
    description=(
        'Performs pairwise distance testing between sample pairs from each '
        'subject. Sample pairs may represent a typical intervention study, '
        'e.g., samples collected pre- and post-treatment; paired samples from '
        'two different timepoints (e.g., in a longitudinal study design), or '
        'identical samples receiving different two different treatments. '
        'This action tests whether the pairwise distance between each subject '
        'pair differs between groups (e.g., groups of subjects receiving '
        'different treatments) and produces boxplots of paired distance '
        'distributions for each group.')
)


plugin.visualizers.register_function(
    function=linear_mixed_effects,
    inputs={'table': FeatureTable[RelativeFrequency]},
    parameters={**base_parameters,
                'metric': Str,
                'group_columns': Str,
                'random_effects': Str,
                'lowess': Bool,
                'ci': Float % Range(0, 100),
                'formula': Str},
    input_descriptions={'table': (
        'Feature table containing metric.')},
    parameter_descriptions={
        **base_parameter_descriptions,
        'metric': ('Dependent variable column name. Must be a column '
                   'name located in the metadata or feature table files.'),
        'state_column': miscellaneous_parameter_descriptions['state_column'],
        'group_columns': (
            'Comma-separated list (without spaces) of metadata columns to '
            'use as independent covariates used to determine mean structure '
            'of "metric".'),
        'random_effects': (
            'Comma-separated list (without spaces) of metadata columns to '
            'use as independent covariates used to determine the variance and '
            'covariance structure (random effects) of "metric". To add a '
            'random slope, the same value passed to "state_column" should '
            'be passed here. A random intercept for each individual is set '
            'by default and does not need to be passed here.'),
        'lowess': 'Estimate locally weighted scatterplot smoothing. Note '
                  'that this will eliminate confidence interval plotting.',
        'ci': 'Size of the confidence interval for the regression estimate.',
        'formula': (
            'R-style formula to use for model specification. A formula must '
            'be used if the "metric" parameter is None. Note that the metric '
            'and group columns specified in the formula will override metric '
            'and group columns that are passed separately as parameters to '
            'this method. ' + formula_description),
    },
    name='Linear mixed effects modeling',
    description=(
        'Linear mixed effects models evaluate the contribution of exogenous '
        'covariates "group_columns" and "random_effects" to a single '
        'dependent variable, "metric". Perform LME and plot line plots of '
        'each group column. A feature table artifact is required input, '
        'though whether "metric" is derived from the feature table or '
        'metadata is optional.'),
    citations=[citations['seabold2010statsmodels']]
)


plugin.visualizers.register_function(
    function=anova,
    inputs={},
    parameters={'metadata': Metadata,
                'formula': Str,
                'sstype': Str % Choices(['I', 'II', 'III'])},
    input_descriptions={},
    parameter_descriptions={
        'metadata': 'Sample metadata containing formula terms.',
        'formula': 'R-style formula specifying the model. All terms must be '
                   'present in the sample metadata or metadata-transformable '
                   'artifacts and can be continuous or categorical metadata '
                   'columns. ' + formula_description,
        'sstype': (
            'Type of sum of squares calculation to perform (I, II, or III).')
    },
    name='ANOVA test',
    description=('Perform an ANOVA test on any factors present in a metadata '
                 'file and/or metadata-transformable artifacts. This is '
                 'followed by pairwise t-tests to examine pairwise '
                 'differences between categorical sample groups.')
)


_VOLATILITY_SCALE_OPTS = ['linear', 'pow', 'sqrt', 'log']

plugin.visualizers.register_function(
    function=volatility,
    inputs={
        'table': FeatureTable[RelativeFrequency],
    },
    parameters={
        **shared_parameters,
        'state_column': miscellaneous_parameters['state_column'],
        'default_metric': Str,
        'default_group_column': Str,
        'yscale': Str % Choices(_VOLATILITY_SCALE_OPTS)
    },
    input_descriptions={
        'table': 'Feature table containing metrics.',
    },
    parameter_descriptions={
        **shared_parameter_descriptions,
        'state_column': miscellaneous_parameter_descriptions['state_column'],
        'default_metric': 'Numeric metadata or artifact column to test by '
                          'default (all numeric metadata columns will be '
                          'available in the visualization).',
        'default_group_column': 'The default metadata column on which to '
                                'separate groups for comparison (all '
                                'categorical metadata columns will be '
                                'available in the visualization).',
        'yscale': 'y-axis scaling strategy to apply.',
    },
    name='Generate interactive volatility plot',
    description=(
        'Generate an interactive control chart depicting the longitudinal '
        'volatility of sample metadata and/or feature frequencies across time '
        '(as set using the "state_column" parameter). Any numeric metadata '
        'column (and metadata-transformable artifacts, e.g., alpha diversity '
        'results) can be plotted on the y-axis, and are selectable using the '
        '"metric_column" selector. Metric values are averaged to compare '
        'across any categorical metadata column using the "group_column" '
        'selector. Longitudinal volatility for individual subjects sampled '
        'over time is co-plotted as "spaghetti" plots if the '
        '"individual_id_column" parameter is used. state_column will '
        'typically be a measure of time, but any numeric metadata column can '
        'be used.')
)


plugin.visualizers.register_function(
    function=plot_feature_volatility,
    inputs={
        'table': FeatureTable[RelativeFrequency],
        'importances': FeatureData[Importance],
    },
    parameters={
        **shared_parameters,
        'state_column': miscellaneous_parameters['state_column'],
        'default_group_column': Str,
        'yscale': Str % Choices(_VOLATILITY_SCALE_OPTS),
        **volatility_filtering_parameters,
        'missing_samples': Str % Choices(['error', 'ignore'])
    },
    input_descriptions={
        'table': 'Feature table containing features found in importances.',
        'importances': 'Feature importance scores.',
    },
    parameter_descriptions={
        **shared_parameter_descriptions,
        'state_column': miscellaneous_parameter_descriptions['state_column'],
        'default_group_column': 'The default metadata column on which to '
                                'separate groups for comparison (all '
                                'categorical metadata columns will be '
                                'available in the visualization).',
        'yscale': 'y-axis scaling strategy to apply.',
        **volatility_filtering_parameter_descriptions,
        'missing_samples': (
            'How to handle missing samples in metadata. "error" will fail '
            'if missing samples are detected. "ignore" will cause the '
            'feature table and metadata to be filtered, so that only '
            'samples found in both files are retained.')
    },
    name='Plot longitudinal feature volatility and importances',
    description=(
        'Plots an interactive control chart of feature abundances (y-axis) '
        'in each sample across time (or state; x-axis). Feature importance '
        'scores and descriptive statistics for each feature are plotted '
        'in interactive bar charts below the control chart, facilitating '
        'exploration of longitudinal feature data. This visualization is '
        'intended for use with the feature-volatility pipeline; use that '
        'pipeline to access this visualization.')
)


plugin.methods.register_function(
    function=nmit,
    inputs={'table': FeatureTable[RelativeFrequency]},
    parameters={**shared_parameters,
                'corr_method': Str % Choices(
                    ["kendall", "pearson", "spearman"]),
                'dist_method': Str % Choices(["fro", "nuc"])},
    outputs=[('distance_matrix', DistanceMatrix)],
    input_descriptions={'table': (
        'Feature table to use for microbial interdependence test.')},
    parameter_descriptions={
        **shared_parameter_descriptions,
        'corr_method': 'The temporal correlation test to be applied.',
        'dist_method': (
            'Temporal distance method, see numpy.linalg.norm for details.'),
    },
    output_descriptions={'distance_matrix': 'The resulting distance matrix.'},
    name='Nonparametric microbial interdependence test',
    description=(
        'Perform nonparametric microbial interdependence test to determine '
        'longitudinal sample similarity as a function of temporal microbial '
        'composition. For more details and citation, please see '
        'doi.org/10.1002/gepi.22065'),
    citations=[citations['zhang2017multivariate']]
)


plugin.methods.register_function(
    function=first_differences,
    inputs={'table': FeatureTable[RelativeFrequency]},
    parameters={**miscellaneous_parameters,
                **shared_parameters,
                'metric': Str,
                'baseline': Float},
    outputs=[('first_differences', SampleData[FirstDifferences])],
    input_descriptions={
        'table': ('Feature table to optionally use for computing first '
                  'differences.')},
    parameter_descriptions={
        **miscellaneous_parameter_descriptions,
        **shared_parameter_descriptions,
        'state_column': miscellaneous_parameter_descriptions['state_column'],
        'metric': 'Numerical metadata or artifact column to test.',
        'baseline': (
            'A value listed in the state_column metadata column against which '
            'all other states should be compared. Toggles calculation of '
            'static differences instead of first differences (which are '
            'calculated if no value is given for baseline). If a "baseline" '
            'value is provided, sample differences at each state are compared '
            'against the baseline state, instead of the previous state. Must '
            'be a value listed in the state_column.')
    },
    output_descriptions={'first_differences': 'Series of first differences.'},
    name=('Compute first differences or difference from baseline between '
          'sequential states'),
    description=(
        'Calculates first differences in "metric" between sequential states '
        'for samples collected from individual subjects sampled repeatedly at '
        'two or more states. First differences can be performed on a metadata '
        'column (including artifacts that can be input as metadata) or a '
        'feature in a feature table. Outputs a data '
        'series of first differences for each individual subject at each '
        'sequential pair of states, labeled by the SampleID of the second '
        'state (e.g., paired differences between time 0 and time 1 would be '
        'labeled by the SampleIDs at time 1). This file can be used as input '
        'to linear mixed effects models or other longitudinal or diversity '
        'methods to compare changes in first differences across time or among '
        'groups of subjects. Also supports differences from baseline (or '
        'other static comparison state) by setting the "baseline" parameter.')
)


plugin.methods.register_function(
    function=first_distances,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={**miscellaneous_parameters,
                **shared_parameters,
                'baseline': Float},
    outputs=[('first_distances', SampleData[FirstDifferences])],
    input_descriptions={
        'distance_matrix': 'Matrix of distances between pairs of samples.'},
    parameter_descriptions={
        **miscellaneous_parameter_descriptions,
        **shared_parameter_descriptions,
        'state_column': miscellaneous_parameter_descriptions['state_column'],
        'baseline': (
            'A value listed in the state_column metadata column against which '
            'all other states should be compared. Toggles calculation of '
            'static distances instead of first distances (which are '
            'calculated if no value is given for baseline). If a "baseline" '
            'value is provided, sample distances at each state are compared '
            'against the baseline state, instead of the previous state. Must '
            'be a value listed in the state_column.')
    },
    output_descriptions={'first_distances': 'Series of first distances.'},
    name=('Compute first distances or distance from baseline between '
          'sequential states'),
    description=(
        'Calculates first distances between sequential states for samples '
        'collected from individual subjects sampled repeatedly at two or more '
        'states. This method is similar to the "first differences" method, '
        'except that it requires a distance matrix as input and calculates '
        'first differences as distances between successive states. Outputs a '
        'data series of first distances for each individual subject at each '
        'sequential pair of states, labeled by the SampleID of the second '
        'state (e.g., paired distances between time 0 and time 1 would be '
        'labeled by the SampleIDs at time 1). This file can be used as input '
        'to linear mixed effects models or other longitudinal or diversity '
        'methods to compare changes in first distances across time or among '
        'groups of subjects. Also supports distance from baseline (or '
        'other static comparison state) by setting the "baseline" parameter.')
)


plugin.pipelines.register_function(
    function=feature_volatility,
    inputs={'table': FeatureTable[Frequency]},
    parameters={
        'metadata': Metadata,
        **shared_parameters,
        'state_column': miscellaneous_parameters['state_column'],
        **parameters['base'],
        **parameters['cv'],
        'estimator': Str % Choices(
            ['RandomForestRegressor', 'ExtraTreesRegressor',
             'GradientBoostingRegressor', 'AdaBoostRegressor', 'ElasticNet',
             'Ridge', 'Lasso', 'KNeighborsRegressor', 'LinearSVR', 'SVR']),
        **volatility_filtering_parameters},
    outputs=[('filtered_table', FeatureTable[RelativeFrequency]),
             ('feature_importance', FeatureData[Importance]),
             ('volatility_plot', Visualization),
             ('accuracy_results', Visualization),
             ('sample_estimator', SampleEstimator[Regressor])],
    input_descriptions={'table': ('Feature table containing all features that '
                                  'should be used for target prediction.')},
    parameter_descriptions={
        'metadata': 'Sample metadata containing state_column, '
                    'individual_id_column, and other metadata for use in '
                    'volatility plots.',
        **shared_parameter_descriptions,
        'state_column': 'Metadata containing collection time (state) values '
                        'for each sample. Must contain exclusively numeric '
                        'values.',
        **parameter_descriptions['base'],
        **parameter_descriptions['cv'],
        **parameter_descriptions['estimator'],
        **volatility_filtering_parameter_descriptions},
    output_descriptions={
        'filtered_table': 'Feature table containing only important features.',
        'feature_importance': output_descriptions['feature_importance'],
        'volatility_plot': 'Interactive volatility plot visualization.',
        'accuracy_results': pipeline_output_descriptions['accuracy_results'],
        'sample_estimator': 'Trained sample regressor.'},
    name='Feature volatility analysis',
    description=(
        'Identify features that are predictive of a numeric metadata column, '
        'state_column (e.g., time), and plot their relative frequencies '
        'across states using interactive feature volatility plots. A '
        'supervised learning regressor is used to identify important features '
        'and assess their ability to predict sample states. state_column will '
        'typically be a measure of time, but any numeric metadata column can '
        'be used.')
)


plugin.pipelines.register_function(
    function=maturity_index,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'group_by': Str,
                'control': Str,
                'individual_id_column': Str,
                'estimator': regressors,
                **pipeline_parameters,
                'metadata': Metadata,
                'state_column': Str,
                'stratify': Bool,
                'feature_count': Int % Range(0, None)
                },
    outputs=regressor_pipeline_outputs + [
        ('maz_scores', SampleData[RegressorPredictions]),
        ('clustermap', Visualization),
        ('volatility_plots', Visualization)],
    input_descriptions={'table': ('Feature table containing all features that '
                                  'should be used for target prediction.')},
    parameter_descriptions={
        **pipeline_parameter_descriptions,
        'state_column': ('Numeric metadata column containing sampling time '
                         '(state) data to use as prediction target.'),
        'group_by': ('Categorical metadata column to use for plotting and '
                     'significance testing between main treatment groups.'),
        'control': (
            'Value of group_by to use as control group. The regression model '
            'will be trained using only control group data, and the maturity '
            'scores of other groups consequently will be assessed relative to '
            'this group.'),
        'individual_id_column': (
            'Optional metadata column containing IDs for individual subjects. '
            'Adds individual subject (spaghetti) vectors to volatility charts '
            'if a column name is provided.'),
        'estimator': 'Regression model to use for prediction.',
        'stratify': 'Evenly stratify training and test data among metadata '
                    'categories. If True, all values in column must match '
                    'at least two samples.',
        'feature_count': 'Filter feature table to include top N most '
                         'important features. Set to zero to include all '
                         'features.'},
    output_descriptions={
        **pipeline_output_descriptions,
        'maz_scores': 'Microbiota-for-age z-score predictions.',
        'clustermap': 'Heatmap of important feature abundance at each time '
                      'point in each group.',
        'volatility_plots': 'Interactive volatility plots of MAZ and maturity '
                            'scores, target (column) predictions, and the '
                            'sample metadata.'},
    name='Microbial maturity index prediction.',
    description=('Calculates a "microbial maturity" index from a regression '
                 'model trained on feature data to predict a given continuous '
                 'metadata column, e.g., to predict age as a function of '
                 'microbiota composition. The model is trained on a subset of '
                 'control group samples, then predicts the column value for '
                 'all samples. This visualization computes maturity index '
                 'z-scores to compare relative "maturity" between each group, '
                 'as described in doi:10.1038/nature13421. This method can '
                 'be used to predict between-group differences in relative '
                 'trajectory across any type of continuous metadata gradient, '
                 'e.g., intestinal microbiome development by age, microbial '
                 'succession during wine fermentation, or microbial community '
                 'differences along environmental gradients, as a function of '
                 'two or more different "treatment" groups.'),
    citations=[citations['subramanian2014persistent'],
               citations['Bokulich306167']]
)


importlib.import_module('q2_longitudinal._transformer')

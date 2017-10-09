# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from qiime2.plugin import (Str, Bool, Plugin, Metadata, Choices, Range, Float,
                           Int)
from q2_types.feature_table import FeatureTable, RelativeFrequency
from ._longitudinal import (pairwise_differences, pairwise_distances,
                            linear_mixed_effects, volatility, nmit)
import q2_longitudinal
from q2_types.distance_matrix import DistanceMatrix


plugin = Plugin(
    name='longitudinal',
    version=q2_longitudinal.__version__,
    website="https://github.com/qiime2/q2-longitudinal",
    package='q2_longitudinal',
    description=(
        'This QIIME 2 plugin supports methods for analysis of time series '
        'analysis, involving either paired sample comparisons or longitudinal '
        'study designs.'),
    short_description='Plugin for paired sample and time series analyses.'
)


shared_parameters = {
    'metadata': Metadata,
    'individual_id_column': Str,
}

base_parameters = {
    **shared_parameters,
    'state_column': Str,
    'palette': Str % Choices([
        'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Paired', 'Accent',
        'Dark2', 'tab10', 'tab20', 'tab20b', 'tab20c', 'viridis', 'plasma',
        'inferno', 'magma', 'terrain', 'rainbow']),
}

paired_params = {
    **base_parameters,
    'group_column': Str,
    'state_1': Str,
    'state_2': Str,
    'parametric': Bool,
    'replicate_handling': Str % Choices(['error', 'random', 'drop']),
}

shared_parameter_descriptions = {
        'metadata': (
            'Sample metadata containing group_column,  state_column, '
            'individual_id_column, and optionally metric values.'),
        'individual_id_column': (
            'Metadata column containing IDs for individual subjects.'),
}

base_parameter_descriptions = {
        **shared_parameter_descriptions,
        'state_column': ('Metadata column containing state (e.g., Time) '
                         'across which samples are paired.'),
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
        'state_1': 'Baseline state column value.',
        'state_2': 'State column value to pair with baseline.',
        'parametric': ('Perform parametric (ANOVA and t-tests) or non-'
                       'parametric (Kruskal-Wallis, Wilcoxon, and Mann-'
                       'Whitney U tests) statistical tests.'),
        'replicate_handling': (
            'Choose how replicate samples are handled. If replicates are '
            'detected, "error" causes method to fail; "drop" will discard all '
            'subject IDs with replicate samples at either state_1 or state_2; '
            '"random" chooses one representative at random from among '
            'replicates.')
}


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
        'subject. Sample pairs may represent a typical intervention study, '
        'e.g., samples collected pre- and post-treatment; paired samples from '
        'two different timepoints (e.g., in a longitudinal study design), or '
        'identical samples receiving different two different treatments. '
        'This action tests whether the change in a continuous metadata value '
        '"metric" differs from zero and differs between groups (e.g., groups '
        'of subjects receiving different treatments), and produces boxplots '
        'of paired difference distributions for each group. A feature table '
        'artifact is required input, though whether "metric" is derived from '
        'the feature table or metadata is optional.')
)


plugin.visualizers.register_function(
    function=pairwise_distances,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={**paired_params},
    input_descriptions={'distance_matrix': (
        'Matrix of distances between pairs of samples.')},
    parameter_descriptions={**paired_parameter_descriptions},
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
                'group_categories': Str,
                'lowess': Bool,
                'ci': Float % Range(0, 100)},
    input_descriptions={'table': (
        'Feature table to optionally use for paired comparisons.')},
    parameter_descriptions={
        **base_parameter_descriptions,
        'metric': ('Dependent variable column name. Must be a column '
                   'name located in the metadata or feature table files.'),
        'group_categories': (
            'Comma-separated list (without spaces) of metadata categories to '
            'use as independent covariates used to determine mean structure '
            'of "metric".'),
        'lowess': ('Estimate locally weighted scatterplot smoothing. Note '
                   'that this will eliminate confidence interval plotting.'),
        'ci': 'Size of the confidence interval for the regression estimate.',
    },
    name='Linear mixed effects modeling',
    description=(
        'Linear mixed effects models evaluate the contribution of exogenous '
        'covariates "group_categories" to a single dependent variable, '
        '"metric". Perform LME and plot line plots of each group column. A '
        'feature table artifact is required input, though whether "metric" is '
        'derived from the feature table or metadata is optional.')
)


plugin.visualizers.register_function(
    function=volatility,
    inputs={'table': FeatureTable[RelativeFrequency]},
    parameters={**base_parameters,
                'metric': Str,
                'group_column': Str,
                'ci': Float % Range(0, 100),
                'plot_control_limits': Bool,
                'yscale': Str % Choices(["linear", "log", "symlog", "logit"]),
                'spaghetti': Bool,
                'xtick_interval': Int},
    input_descriptions={'table': (
        'Feature table to optionally use for paired comparisons.')},
    parameter_descriptions={
        **base_parameter_descriptions,
        'metric': 'Numerical metadata or artifact column to test.',
        'group_column': (
            'Metadata column on which to separate groups for comparison'),
        'ci': 'Size of the confidence interval to plot on control chart.',
        'plot_control_limits': ('Plot global mean and control limits (2X and '
                                '3X standard deviations).'),
        'yscale': 'y-axis scaling strategy to apply.',
        'spaghetti': 'Co-plot spaghetti plot of per-individual trajectories.',
        'xtick_interval': ('Interval between major tick marks on x axis. '
                           'Defaults to 1, or autoscales to show up to 30 '
                           'ticks if data contain more than 30 x-axis values.')
    },
    name='Volatility analysis',
    description=(
        'Plot control chart of a single dependent variable, "metric", across '
        'multiple groups contained in sample metadata column "group_column".')
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
        'doi.org/10.1002/gepi.22065')
)

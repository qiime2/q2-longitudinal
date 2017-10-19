# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime2
import pandas as pd
from qiime2.plugin import (Str, Bool, Plugin, Metadata, Choices, Range, Float,
                           Int, ValidationError)
from q2_types.feature_table import FeatureTable, RelativeFrequency
from q2_types.distance_matrix import DistanceMatrix
from q2_types.sample_data import SampleData
from qiime2.plugin import SemanticType
import qiime2.plugin.model as model
from ._longitudinal import (pairwise_differences, pairwise_distances,
                            linear_mixed_effects, volatility, nmit,
                            first_differences, first_distances)
import q2_longitudinal


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


# FirstDifferencesFormat is a near replica of AlphaDiversityFormat in q2_types
FirstDifferences = SemanticType(
    'FirstDifferences', variant_of=SampleData.field['type'])


class FirstDifferencesFormat(model.TextFileFormat):
    def _check_n_records(self, n_lines=None):
        with self.open() as fh:
            line = fh.readline()
            for line, _ in zip(fh, range(n_lines)):
                cells = line.strip().split('\t')
                if len(cells) != 2:
                    raise ValidationError('File is not FirstDifferencesFormat')

    def _validate_(self, level):
        record_count_map = {'min': 5, 'max': 1000}
        self._check_n_records(record_count_map[level])


FirstDifferencesDirectoryFormat = model.SingleFileDirectoryFormat(
    'FirstDifferencesDirectoryFormat', 'FirstDifferences.tsv',
    FirstDifferencesFormat)


# borrowed from q2_types
def _read_dataframe(fh):
    # Using `dtype=object` and `set_index` to avoid type casting/inference
    # of any columns or the index.
    df = pd.read_csv(fh, sep='\t', header=0, dtype=object)
    df.set_index(df.columns[0], drop=True, append=False, inplace=True)
    df.index.name = None
    return df


@plugin.register_transformer
def _4(data: pd.Series) -> (FirstDifferencesFormat):
    ff = FirstDifferencesFormat()
    with ff.open() as fh:
        data.to_csv(fh, sep='\t', header=True)
    return ff


@plugin.register_transformer
def _5(ff: FirstDifferencesFormat) -> (pd.Series):
    with ff.open() as fh:
        df = _read_dataframe(fh)
        return df.iloc[:, 0]


@plugin.register_transformer
def _6(ff: FirstDifferencesFormat) -> (qiime2.Metadata):
    with ff.open() as fh:
        return qiime2.Metadata(_read_dataframe(fh))


plugin.register_formats(
    FirstDifferencesFormat, FirstDifferencesDirectoryFormat)

plugin.register_semantic_types(FirstDifferences)

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
        'inferno', 'magma', 'terrain', 'rainbow']),
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
    'state_column': ('Metadata column containing state (e.g., Time) '
                     'across which samples are paired.'),
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
        'state_column': miscellaneous_parameter_descriptions['state_column'],
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
            miscellaneous_parameter_descriptions['replicate_handling']),
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


plugin.methods.register_function(
    function=first_differences,
    inputs={'table': FeatureTable[RelativeFrequency]},
    parameters={**miscellaneous_parameters,
                **shared_parameters,
                'metric': Str},
    outputs=[('first_differences', SampleData[FirstDifferences])],
    input_descriptions={
        'table': ('Feature table to optionally use for computing first '
                  'differences.')},
    parameter_descriptions={
        **miscellaneous_parameter_descriptions,
        **shared_parameter_descriptions,
        'metric': 'Numerical metadata or artifact column to test.',
    },
    output_descriptions={'first_differences': 'Series of first differences.'},
    name='First difference computation between sequential states',
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
        'groups of subjects.')
)


plugin.methods.register_function(
    function=first_distances,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={**miscellaneous_parameters,
                **shared_parameters},
    outputs=[('first_distances', SampleData[FirstDifferences])],
    input_descriptions={
        'distance_matrix': 'Matrix of distances between pairs of samples.'},
    parameter_descriptions={
        **miscellaneous_parameter_descriptions,
        **shared_parameter_descriptions,
    },
    output_descriptions={'first_distances': 'Series of first distances.'},
    name='First distance computation between sequential states',
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
        'groups of subjects.')
)

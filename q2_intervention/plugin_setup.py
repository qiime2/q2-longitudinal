# ----------------------------------------------------------------------------
# Copyright (c) 2017-, q2-intervention development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from qiime2.plugin import Str, Bool, Plugin, Metadata, Choices
from q2_types.feature_table import FeatureTable, Frequency
from ._intervention import paired_differences
import q2_intervention


plugin = Plugin(
    name='intervention',
    version=q2_intervention.__version__,
    website="https://github.com/nbokulich/q2-intervention",
    package='q2_intervention'
)


plugin.visualizers.register_function(
    function=paired_differences,
    inputs={'table': FeatureTable[Frequency]},
    parameters={
        'metadata': Metadata,
        'metric': Str,
        'group_category': Str,
        'state_category': Str,
        'state_pre': Str,
        'state_post': Str,
        'individual_id_category': Str,
        'parametric': Bool,
        'pallette': Str % Choices([
            'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'tab10', 'tab20', 'tab20b', 'tab20c', 'viridis', 'plasma',
            'inferno', 'magma', 'terrain', 'rainbow']),
        'drop_duplicates': Bool,
    },
    input_descriptions={'table': (
        'Feature table to optionally use for paired comparisons.')},
    parameter_descriptions={
        'metadata': (
            'Sample metadata containing group_category,  state_category, '
            'individual_id_category, and optionally metric values.'),
        'metric': 'Numerical metadata or artifact category to test.',
        'group_category': (
            'Metadata category on which to separate groups for comparison'),
        'state_category': ('Metadata category containing state (e.g., Time) '
                           'across which samples are paired.'),
        'state_pre': 'Baseline state category value.',
        'state_post': 'State category value to pair with baseline.',
        'individual_id_category': (
            'Metadata category containing subject IDs  to use for pairing '
            'samples. WARNING: if duplicates exist for an individual ID at '
            'either state_pre or state_post, that subject will be dropped and '
            'reported in standard output by default. Set duplicates="ignore" '
            'to instead randomly select one member, and use --verbose to list '
            'conflicts.'),
        'parametric': ('Perform parametric (ANOVA and t-tests) or non-'
                       'parametric (Kruskal-Wallis, Wilcoxon, and Mann-'
                       'Whitney U tests) statistical tests.'),
        'pallette': 'Color palette to use for generating boxplots.',
        'drop_duplicates': (
            'If True, will discard all subject IDs with duplicate samples '
            'at either state_pre or state_post. If False, will instead '
            'choose one representative at random from among duplicates.')
    },
    name='Paired difference testing and boxplots',
    description=(
        'Performs paired difference testing between samples from each subject'
        '. Sample pairs may represent a typical intervention study, e.g., '
        'samples collected pre- and post-treatment; paired samples from two '
        'different timepoints (e.g., in a longitudinal study design), or '
        'identical samples receiving different two different treatments. '
        'This action tests whether the change in a continuous metadata value '
        '"metric" differs from zero and differs between groups (e.g., groups '
        'of subjects receiving different treatments), and produces boxplots '
        'of paired difference distributions for each group. A feature table '
        'artifact is required input, though whether "metric" is derived from '
        'the feature table or metadata is optional.')
)

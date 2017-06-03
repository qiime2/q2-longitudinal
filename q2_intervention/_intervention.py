# ----------------------------------------------------------------------------
# Copyright (c) 2017-, q2-intervention development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import qiime2
import pandas as pd
from skbio.stats.distance import MissingIDError
from skbio import DistanceMatrix

from ._utilities import (_get_group_pairs, _extract_distance_distribution,
                         _between_subject_distance_distribution,
                         _get_paired_differences, test_paired_differences,
                         _multiple_group_difference, per_method_pairwise_tests,
                         boxplot_from_dict, _visualize)




def paired_differences(output_dir: str, table: pd.DataFrame,
                       metadata: qiime2.Metadata, metric: str,
                       group_category: str, state_category: str='Time',
                       state_pre: str='Pre', state_post: str='post',
                       individual_id_category: str='SubjectID',
                       parametric: bool=True, pallette: str='Set1',
                       drop_duplicates: bool=True):

    # find metric in metadata or derive from table and merge into metadata
    metadata = metadata.to_dataframe()
    if metric not in metadata.columns:
        if metric in table.columns:
            metadata = pd.concat(
                [metadata, pd.DataFrame(table[metric])], axis=1, join='inner')
        else:
            raise ValueError(
                'metric must be a valid metadata or feature table category.')

    # calculate paired difference distributions
    pairs = {}
    group_names = metadata[group_category].unique()
    for group in group_names:
        group_pairs = _get_group_pairs(
            metadata, group_value=group,
            individual_id_category=individual_id_category,
            group_category=group_category, state_category=state_category,
            state_values=[state_pre, state_post],
            drop_duplicates=drop_duplicates)
        pairs[group] = _get_paired_differences(metadata, group_pairs, metric)

    # kruskal test or ANOVA between groups
    multiple_group_test = _multiple_group_difference(
        pairs.values(), parametric=parametric)

    # pairwise testing
    pairwise_tests = per_method_pairwise_tests(
        pairs, paired=False, parametric=parametric)

    # paired difference tests
    paired_difference_tests = test_paired_differences(
        pairs, parametric=parametric)

    # boxplots
    ax = boxplot_from_dict(
        pairs, pallette=pallette, y_label=metric, x_label=group_category)

    summary = pd.Series(
        [metric, group_category, state_category, state_pre, state_post,
         individual_id_category, parametric, drop_duplicates],
        index=['Metric', 'Group category', 'State category', 'State pre',
               'State post', 'Individual ID category', 'Parametric',
               'Drop duplicates'],
        name='Paired difference tests')

    _visualize(output_dir, multiple_group_test, pairwise_tests,
               paired_difference_tests, ax.get_figure(), summary=summary)

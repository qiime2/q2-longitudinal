# ----------------------------------------------------------------------------
# Copyright (c) 2017-, q2-intervention development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import qiime2
import pandas as pd
from skbio import DistanceMatrix

from ._utilities import (_get_group_pairs, _extract_distance_distribution,
                         _between_subject_distance_distribution,
                         _get_paired_differences, _stats_and_visuals)


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

    # Calculate test statistics and generate boxplots
    _stats_and_visuals(
        output_dir, pairs, metric, group_category, state_category, state_pre,
        state_post, individual_id_category, parametric, pallette,
        drop_duplicates, multiple_group_test=True, pairwise_tests=True,
        paired_difference_tests=True, boxplot=True)


def pairwise_distance(output_dir: str, distance_matrix: DistanceMatrix,
                      metadata: qiime2.Metadata,
                      group_category: str, state_category: str='Time',
                      state_pre: str='Pre', state_post: str='post',
                      individual_id_category: str='SubjectID',
                      parametric: bool=True, pallette: str='Set1',
                      drop_duplicates: bool=True, between_group_distance=True):

    metadata = metadata.to_dataframe()

    # calculate pairwise distance distributions
    pairs = {}
    group_names = metadata[group_category].unique()
    for group in group_names:
        group_pairs = _get_group_pairs(
            metadata, group_value=group,
            individual_id_category=individual_id_category,
            group_category=group_category, state_category=state_category,
            state_values=[state_pre, state_post],
            drop_duplicates=drop_duplicates)
        pairs[group] = _extract_distance_distribution(
            distance_matrix, group_pairs)
        if between_group_distance:
            between = group + '_between_subject'
            pairs[between] = _between_subject_distance_distribution(
                distance_matrix, group_pairs, metadata, group_category, group)

    # Calculate test statistics and generate boxplots
    _stats_and_visuals(
        output_dir, pairs, 'distance', group_category,
        state_category, state_pre, state_post, individual_id_category,
        parametric, pallette, drop_duplicates, multiple_group_test=True,
        pairwise_tests=True, paired_difference_tests=False, boxplot=True)

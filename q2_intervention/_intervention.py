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
                         _between_subject_distance_distribution, _visualize,
                         _get_paired_differences, _stats_and_visuals,
                         _add_metric_to_metadata, _linear_effects,
                         _regplot_subplots_from_dataframe, _load_metadata)


def paired_differences(output_dir: str, metadata: qiime2.Metadata,
                       group_category: str, metric: str,
                       table: pd.DataFrame=None,
                       state_category: str='Time',
                       state_pre: str='Pre', state_post: str='post',
                       individual_id_category: str='SubjectID',
                       parametric: bool=True, palette: str='Set1',
                       drop_duplicates: bool=True) -> None:

    # find metric in metadata or derive from table and merge into metadata
    metadata = _add_metric_to_metadata(table, metadata, metric)

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
        state_post, individual_id_category, parametric, palette,
        drop_duplicates, multiple_group_test=True, pairwise_tests=True,
        paired_difference_tests=True, boxplot=True)


def pairwise_distance(output_dir: str, distance_matrix: DistanceMatrix,
                      metadata: qiime2.Metadata,
                      group_category: str, state_category: str='Time',
                      state_pre: str='Pre', state_post: str='post',
                      individual_id_category: str='SubjectID',
                      parametric: bool=True, palette: str='Set1',
                      drop_duplicates: bool=True,
                      between_group_distance=True) -> None:

    metadata = _load_metadata(metadata)

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
        parametric, palette, drop_duplicates, multiple_group_test=True,
        pairwise_tests=True, paired_difference_tests=False, boxplot=True)


def linear_mixed_effects(output_dir: str, metadata: qiime2.Metadata,
                         group_categories: str, metric: str,
                         table: pd.DataFrame=None, state_category: str='Time',
                         individual_id_category: str='SubjectID',
                         palette: str='Set1', lowess=False, ci=95) -> None:
    # split group_categories into list of categories
    group_categories = group_categories.split(",")

    # find metric in metadata or derive from table and merge into metadata
    metadata = _add_metric_to_metadata(table, metadata, metric)

    # Generate LME model summary
    model_summary, model_results = _linear_effects(
        metadata, metric, state_category, group_categories,
        individual_id_category)

    # Plot dependent variable as function of independent variables
    g = _regplot_subplots_from_dataframe(
        state_category, metric, metadata, group_categories, lowess=lowess,
        ci=ci, palette=palette)

    # summarize parameters and visualize
    summary = pd.Series(
        [metric, ", ".join(group_categories), state_category,
         individual_id_category],
        index=['Metric', 'Group category', 'State category',
               'Individual ID category'],
        name='Linear mixed effects parameters')

    _visualize(output_dir, model_summary=model_summary,
               model_results=model_results, plot=g, summary=summary)

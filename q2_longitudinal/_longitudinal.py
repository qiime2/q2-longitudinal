# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path

import pandas as pd
import skbio
import qiime2

from ._utilities import (_get_group_pairs, _extract_distance_distribution,
                         _visualize, _validate_metadata_is_superset,
                         _get_pairwise_differences, _stats_and_visuals,
                         _add_metric_to_metadata, _linear_effects,
                         _regplot_subplots_from_dataframe, _load_metadata,
                         _validate_input_values, _validate_input_columns,
                         _control_chart_subplots, _nmit,
                         _validate_is_numeric_column, _tabulate_matrix_ids,
                         _first_differences)


def pairwise_differences(output_dir: str, metadata: qiime2.Metadata,
                         group_column: str, metric: str, state_column: str,
                         state_1: str, state_2: str, individual_id_column: str,
                         parametric: bool=False, palette: str='Set1',
                         replicate_handling: str='error',
                         table: pd.DataFrame=None) -> None:

    # find metric in metadata or derive from table and merge into metadata
    metadata = _add_metric_to_metadata(table, metadata, metric)

    _validate_input_values(metadata, metric, individual_id_column,
                           group_column, state_column, state_1, state_2)

    # calculate paired difference distributions
    pairs = {}
    pairs_summaries = {}
    errors = []
    pairs_summary = pd.DataFrame()
    group_names = metadata[group_column].unique()
    for group in group_names:
        group_pairs, error = _get_group_pairs(
            metadata, group_value=group,
            individual_id_column=individual_id_column,
            group_column=group_column, state_column=state_column,
            state_values=[state_1, state_2],
            replicate_handling=replicate_handling)
        pairs[group], pairs_summaries[group] = _get_pairwise_differences(
            metadata, group_pairs, metric, individual_id_column, group_column)
        pairs_summary = pd.concat([pairs_summary, pairs_summaries[group]])
        errors.extend(error)
    pairs_summary.to_csv(os.path.join(output_dir, 'pairs.tsv'), sep='\t')

    # Calculate test statistics and generate boxplots
    y_label = 'Difference in {0} ({1} {2} - {1} {3})'.format(
        metric, state_column, state_2, state_1)

    _stats_and_visuals(
        output_dir, pairs, y_label, group_column, state_column, state_1,
        state_2, individual_id_column, errors, parametric, palette,
        replicate_handling, multiple_group_test=True, pairwise_tests=True,
        paired_difference_tests=True, boxplot=True)


def pairwise_distances(output_dir: str, distance_matrix: skbio.DistanceMatrix,
                       metadata: qiime2.Metadata, group_column: str,
                       state_column: str, state_1: str, state_2: str,
                       individual_id_column: str, parametric: bool=False,
                       palette: str='Set1', replicate_handling: str='error',
                       ) -> None:

    metadata = _load_metadata(metadata)

    _validate_input_values(metadata, None, individual_id_column, group_column,
                           state_column, state_1, state_2)

    # calculate pairwise distance distributions
    pairs = {}
    pairs_summaries = {}
    errors = []
    pairs_summary = pd.DataFrame()
    group_names = metadata[group_column].unique()
    for group in group_names:
        group_pairs, error = _get_group_pairs(
            metadata, group_value=group,
            individual_id_column=individual_id_column,
            group_column=group_column, state_column=state_column,
            state_values=[state_1, state_2],
            replicate_handling=replicate_handling)
        pairs[group], pairs_summaries[group] = _extract_distance_distribution(
            distance_matrix, group_pairs, metadata, individual_id_column,
            group_column)
        pairs_summary = pd.concat([pairs_summary, pairs_summaries[group]])
        errors.extend(error)
    pairs_summary.to_csv(os.path.join(output_dir, 'pairs.tsv'), sep='\t')

    # Calculate test statistics and generate boxplots
    _stats_and_visuals(
        output_dir, pairs, 'distance', group_column,
        state_column, state_1, state_2, individual_id_column, errors,
        parametric, palette, replicate_handling, multiple_group_test=True,
        pairwise_tests=True, paired_difference_tests=False, boxplot=True,
        plot_name='Pairwise distance boxplot')


def linear_mixed_effects(output_dir: str, metadata: qiime2.Metadata,
                         group_categories: str, metric: str, state_column: str,
                         individual_id_column: str, table: pd.DataFrame=None,
                         palette: str='Set1', lowess: bool=False, ci: int=95
                         ) -> None:

    # split group_categories into list of categories
    group_categories = group_categories.split(",")

    # find metric in metadata or derive from table and merge into metadata
    metadata = _add_metric_to_metadata(table, metadata, metric)

    _validate_input_columns(metadata, individual_id_column, group_categories,
                            state_column, metric)

    # let's force states to be numeric
    _validate_is_numeric_column(metadata, state_column)

    # Generate LME model summary
    model_summary, model_results = _linear_effects(
        metadata, metric, state_column, group_categories,
        individual_id_column)

    # Plot dependent variable as function of independent variables
    g = _regplot_subplots_from_dataframe(
        state_column, metric, metadata, group_categories, lowess=lowess,
        ci=ci, palette=palette)

    # summarize parameters and visualize
    summary = pd.Series(
        [metric, ", ".join(group_categories), state_column,
         individual_id_column],
        index=['Metric', 'Group column', 'State column',
               'Individual ID column'],
        name='Linear mixed effects parameters')

    raw_data = metadata[[
        metric, state_column, individual_id_column, *group_categories]]

    _visualize(output_dir, model_summary=model_summary,
               model_results=model_results, plot=g, summary=summary,
               raw_data=raw_data,
               plot_name='Regression scatterplots')


def volatility(output_dir: str, metadata: qiime2.Metadata, group_column: str,
               metric: str, state_column: str, individual_id_column: str,
               table: pd.DataFrame=None, palette: str='Set1', ci: int=95,
               plot_control_limits: bool=True, xtick_interval: int=None,
               yscale: str='linear',  spaghetti: str='no') -> None:

    # find metric in metadata or derive from table and merge into metadata
    metadata = _add_metric_to_metadata(table, metadata, metric)

    _validate_input_columns(metadata, individual_id_column, group_column,
                            state_column, metric)

    # let's force states to be numeric
    _validate_is_numeric_column(metadata, state_column)

    # plot control charts
    chart, global_mean, global_std = _control_chart_subplots(
        state_column, metric, metadata, group_column, individual_id_column,
        ci=ci, palette=palette, plot_control_limits=plot_control_limits,
        xtick_interval=xtick_interval, yscale=yscale, spaghetti=spaghetti)

    # summarize parameters and visualize
    summary = pd.Series(
        [metric, group_column, state_column, individual_id_column, global_mean,
         global_std],
        index=['Metric', 'Group column', 'State column',
               'Individual ID column', 'Global mean',
               'Global standard deviation'],
        name='Volatility test parameters')

    raw_data = metadata[[
        metric, state_column, individual_id_column, group_column]]

    _visualize(output_dir, plot=chart, summary=summary, raw_data=raw_data,
               plot_name='Control charts')


def nmit(table: pd.DataFrame, metadata: qiime2.Metadata,
         individual_id_column: str, corr_method: str="kendall",
         dist_method: str="fro") -> skbio.DistanceMatrix:

    # load and prep metadata
    metadata = _load_metadata(metadata)
    _validate_metadata_is_superset(metadata, table)
    metadata = metadata[metadata.index.isin(table.index)]

    # validate id column
    _validate_input_columns(metadata, individual_id_column, None, None, None)

    # run NMIT
    _dist = _nmit(
        table, metadata, individual_id_column=individual_id_column,
        corr_method=corr_method, dist_method=dist_method)

    return _dist


def first_differences(metadata: qiime2.Metadata, state_column: str,
                      individual_id_column: str, metric: str,
                      replicate_handling: str='error',
                      table: pd.DataFrame=None) -> pd.Series:

    # find metric in metadata or derive from table and merge into metadata
    if table is not None:
        _validate_metadata_is_superset(metadata.to_dataframe(), table)
        metadata = _add_metric_to_metadata(table, metadata, metric)
    else:
        metadata = _load_metadata(metadata)
        _validate_is_numeric_column(metadata, metric)

    # validate columns
    _validate_input_columns(
        metadata, individual_id_column, None, state_column, metric)

    return _first_differences(
        metadata, state_column, individual_id_column, metric,
        replicate_handling, distance_matrix=None)


def first_distances(distance_matrix: skbio.DistanceMatrix,
                    metadata: qiime2.Metadata, state_column: str,
                    individual_id_column: str,
                    replicate_handling: str='error') -> pd.Series:

    # load and validate metadata
    metadata = _load_metadata(metadata)
    _validate_metadata_is_superset(
        metadata, _tabulate_matrix_ids(distance_matrix))

    # validate columns
    # "Distance" is not actually a metadata value, so don't validate metric!
    _validate_input_columns(
        metadata, individual_id_column, None, state_column, None)

    return _first_differences(
        metadata, state_column, individual_id_column, metric=None,
        replicate_handling=replicate_handling, distance_matrix=distance_matrix)

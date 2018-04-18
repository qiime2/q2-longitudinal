# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import pkg_resources
from distutils.dir_util import copy_tree

import pandas as pd
import skbio
import qiime2
import q2templates

from ._utilities import (_get_group_pairs, _extract_distance_distribution,
                         _visualize, _validate_metadata_is_superset,
                         _get_pairwise_differences, _stats_and_visuals,
                         _add_metric_to_metadata, _linear_effects,
                         _regplot_subplots_from_dataframe, _load_metadata,
                         _validate_input_values, _validate_input_columns,
                         _nmit, _validate_is_numeric_column,
                         _tabulate_matrix_ids, _first_differences)
from ._vega import _render_volatility_spec


TEMPLATES = pkg_resources.resource_filename('q2_longitudinal', 'assets')


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
                         metric: str, state_column: str,
                         individual_id_column: str, group_columns: str=None,
                         random_effects: str=None, table: pd.DataFrame=None,
                         palette: str='Set1', lowess: bool=False, ci: int=95
                         ) -> None:

    raw_data_columns = [metric, state_column, individual_id_column]

    # split group_columns into list of columns
    if group_columns is not None:
        group_columns = group_columns.split(",")
        raw_data_columns.extend(group_columns)
    if random_effects is not None:
        random_effects = random_effects.split(",")
        raw_data_columns.extend(random_effects)

    # find metric in metadata or derive from table and merge into metadata
    metadata = _add_metric_to_metadata(table, metadata, metric)

    _validate_input_columns(metadata, individual_id_column, group_columns,
                            state_column, metric)
    # separately validate random_effects, since these can recycle state_column
    # and individual_id_column and group_column values, but not metric
    _validate_input_columns(metadata, None, random_effects, None, metric)

    # let's force states to be numeric
    _validate_is_numeric_column(metadata, state_column)

    # Generate LME model summary
    model_summary, model_results = _linear_effects(
        metadata, metric, state_column, group_columns,
        individual_id_column, random_effects=random_effects)

    # Plot dependent variable as function of independent variables
    g = _regplot_subplots_from_dataframe(
        state_column, metric, metadata, group_columns, lowess=lowess,
        ci=ci, palette=palette)

    # summarize parameters and visualize
    summary = pd.Series(
        [metric, group_columns, state_column,
         individual_id_column, random_effects],
        index=['Metric', 'Group column', 'State column',
               'Individual ID column', 'Random effects'],
        name='Linear mixed effects parameters')

    raw_data = metadata[list(set(raw_data_columns))]

    _visualize(output_dir, model_summary=model_summary,
               model_results=model_results, plot=g, summary=summary,
               raw_data=raw_data,
               plot_name='Regression scatterplots')


def volatility(output_dir: str, metadata: qiime2.Metadata,
               default_group_column: str, default_metric: str,
               state_column: str, individual_id_column: str,
               table: pd.DataFrame=None, yscale: str='linear') -> None:
    data = _add_metric_to_metadata(table, metadata, default_metric)

    _validate_input_columns(data, individual_id_column, default_group_column,
                            state_column, default_metric)
    _validate_is_numeric_column(data, state_column)

    # Convert back to metadata so that we can save these data in a
    # consistent way.
    data.index.name = 'id'
    metadata = qiime2.Metadata(data)
    metadata.save(os.path.join(output_dir, 'data.tsv'))

    categorical = metadata.filter_columns(column_type='categorical')
    group_columns = list(categorical.columns.keys())
    numeric = metadata.filter_columns(column_type='numeric')
    metric_columns = list(numeric.columns.keys())

    vega_spec = _render_volatility_spec(data, individual_id_column,
                                        state_column, default_group_column,
                                        group_columns, default_metric,
                                        metric_columns, yscale)

    # Order matters here - need to render the template *after* copying the
    # directory tree, otherwise we will overwrite the index.html
    copy_tree(os.path.join(TEMPLATES, 'volatility'), output_dir)
    index = os.path.join(TEMPLATES, 'volatility', 'index.html')
    q2templates.render(index, output_dir, context={'vega_spec': vega_spec})


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
                      replicate_handling: str='error', baseline: float=None,
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
        replicate_handling, baseline=baseline, distance_matrix=None)


def first_distances(distance_matrix: skbio.DistanceMatrix,
                    metadata: qiime2.Metadata, state_column: str,
                    individual_id_column: str, baseline: float=None,
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
        replicate_handling=replicate_handling, baseline=baseline,
        distance_matrix=distance_matrix)

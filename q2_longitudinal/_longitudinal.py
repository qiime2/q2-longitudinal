# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
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
import warnings
import biom
import statsmodels.api as sm
from statsmodels.formula.api import ols

from ._utilities import (_get_group_pairs, _extract_distance_distribution,
                         _visualize, _validate_metadata_is_superset,
                         _get_pairwise_differences, _stats_and_visuals,
                         _add_metric_to_metadata, _linear_effects,
                         _regplot_subplots_from_dataframe, _load_metadata,
                         _validate_input_values, _validate_input_columns,
                         _nmit, _validate_is_numeric_column, _maz_score,
                         _first_differences, _importance_filtering,
                         _summarize_feature_stats, _convert_nan_to_none,
                         _parse_formula, _visualize_anova)
from ._vega_specs import render_spec_volatility


TEMPLATES = pkg_resources.resource_filename('q2_longitudinal', 'assets')


def pairwise_differences(output_dir: str, metadata: qiime2.Metadata,
                         metric: str, state_column: str, state_1: str,
                         state_2: str, individual_id_column: str,
                         group_column: str = None, parametric: bool = False,
                         palette: str = 'Set1',
                         replicate_handling: str = 'error',
                         table: pd.DataFrame = None) -> None:

    # find metric in metadata or derive from table and merge into metadata
    metadata = _add_metric_to_metadata(table, _load_metadata(metadata), metric)

    _validate_input_values(metadata, metric, individual_id_column,
                           group_column, state_column, state_1, state_2)

    # calculate paired difference distributions
    pairs = {}
    pairs_summaries = {}
    errors = []
    pairs_summary = pd.DataFrame()
    if group_column is not None:
        group_names = metadata[group_column].unique()
    else:
        # if we do not stratify by groups, find all pairs.
        group_names = ['all subjects']
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

    if group_column is not None:
        multiple_group_test = pairwise_tests = True
    else:
        multiple_group_test = pairwise_tests = False

    _stats_and_visuals(
        output_dir, pairs, y_label, group_column, state_column, state_1,
        state_2, individual_id_column, errors, parametric, palette,
        replicate_handling, multiple_group_test=multiple_group_test,
        pairwise_tests=pairwise_tests, paired_difference_tests=True,
        boxplot=True)


def pairwise_distances(output_dir: str, distance_matrix: skbio.DistanceMatrix,
                       metadata: qiime2.Metadata, group_column: str,
                       state_column: str, state_1: str, state_2: str,
                       individual_id_column: str, parametric: bool = False,
                       palette: str = 'Set1',
                       replicate_handling: str = 'error',) -> None:

    metadata = _load_metadata(metadata)

    _validate_input_values(metadata, None, individual_id_column, group_column,
                           state_column, state_1, state_2)

    # subset metadata to match distance_matrix ids
    metadata = metadata.filter(distance_matrix.ids, axis=0)

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
                         state_column: str, individual_id_column: str,
                         metric: str = None, group_columns: str = None,
                         random_effects: str = None,
                         table: pd.DataFrame = None,
                         palette: str = 'Set1', lowess: bool = False,
                         ci: int = 95, formula: str = None) -> None:

    metadata = _load_metadata(metadata)

    # Must use formula and/or metric.
    if metric is None and formula is None:
        raise ValueError('Must specify either a metric or a formula that '
                         'contains a valid metric to use as a dependent '
                         'variable.')

    # Formulae must contain state_column and designate a metric.
    # (state_column is used separately for plotting but is ambiguous in the
    # formula, so it does make sense to require even when formula is passed.)
    if formula is not None:
        if '~' not in formula:
            raise ValueError('Formula must be in format "metric ~ independent '
                             'variables".')
        if state_column not in formula:
            raise ValueError(
                '"formula" must contain the "state_column" value: '
                '{0} does not contain {1}'.format(formula, state_column))

    # optionally parse R-style formula for validation. Note that this will
    # override "metric" and "group_columns" parameters if input separately.
    # "formula" is meant to be a "secret" feature for power users familiar with
    # R-style formulae, so I think it is okay to just clarify this in the docs
    # instead of putting in too many safety features (e.g., to prevent this
    # override behavior).
    if formula is not None:
        metric, group_columns = _parse_formula(formula)
        group_columns.remove(state_column)
        group_columns = ','.join(list(group_columns))

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
    model_summary, model_results, model_fit, formula = _linear_effects(
        metadata, metric, state_column, group_columns,
        individual_id_column, random_effects=random_effects, formula=formula)

    # Plot dependent variable as function of independent variables
    g = _regplot_subplots_from_dataframe(
        state_column, metric, metadata, group_columns, lowess=lowess,
        ci=ci, palette=palette)

    # Plot fit vs. residuals
    metadata['residual'] = model_fit.resid
    predicted = 'predicted {0}'.format(metric)
    metadata[predicted] = model_fit.predict()
    res = _regplot_subplots_from_dataframe(
        predicted, 'residual', metadata, group_columns, lowess=lowess,
        ci=ci, palette=palette)

    # add predicted/residual values to "raw" data just for fun
    raw_data_columns.extend([predicted, 'residual'])

    # if the name of predicted/residual is already in metadata, warn users that
    # predicted column is overwritten in the viz raw data download
    for term in [predicted, 'residual']:
        if term in metadata.columns:
            _warn_column_name_exists(predicted)

    # summarize parameters and visualize
    summary = pd.Series(
        [formula, metric, group_columns, state_column,
         individual_id_column, random_effects],
        index=['Fixed Effects formula', 'Metric', 'Group column',
               'State column', 'Individual ID column', 'Random effects'],
        name='Linear mixed effects parameters')

    raw_data = metadata[list(set(raw_data_columns))]

    _visualize(output_dir, model_summary=model_summary,
               model_results=model_results, plot=g, summary=summary,
               raw_data=raw_data,
               plot_name='Regression scatterplots', residuals=res)


def anova(output_dir: str,
          metadata: qiime2.Metadata,
          formula: str,
          sstype: str = 'II') -> None:

    # Grab metric and covariate names from formula
    metric, group_columns = _parse_formula(formula)
    columns = [metric] + list(group_columns)

    # Validate formula (columns are in metadata, etc)
    for col in columns:
        metadata.get_column(col)
    # store categorical column names for later use
    cats = metadata.filter_columns(column_type='categorical').columns.keys()
    metadata = metadata.to_dataframe()[columns].dropna()

    # Run anova
    lm = ols(formula, metadata).fit()
    results = pd.DataFrame(sm.stats.anova_lm(lm, typ=sstype)).fillna('')
    results.to_csv(os.path.join(output_dir, 'anova.tsv'), sep='\t')

    # Run pairwise t-tests with multiple test correction
    pairwise_tests = pd.DataFrame()
    for group in group_columns:
        # only run on categorical columns — numeric columns raise error
        if group in cats:
            ttests = lm.t_test_pairwise(group, method='fdr_bh').result_frame
            pairwise_tests = pd.concat([pairwise_tests, pd.DataFrame(ttests)])
    if pairwise_tests.empty:
        pairwise_tests = False

    # Plot fit vs. residuals
    metadata['residual'] = lm.resid
    metadata['fitted_values'] = lm.fittedvalues
    res = _regplot_subplots_from_dataframe(
        'fitted_values', 'residual', metadata, group_columns, lowess=False,
        ci=95, palette='Set1', fit_reg=False)

    # Visualize results
    _visualize_anova(output_dir, pairwise_tests=pairwise_tests,
                     model_results=results, residuals=res,
                     pairwise_test_name='Pairwise t-tests')


def _warn_column_name_exists(column_name):
    warning = (
        'This is only a warning, and the results of this action are still '
        'valid. The column name "{0}" already exists in your metadata file. '
        'Any "raw" metadata that can be downloaded from the resulting '
        'visualization will contain overwritten values for this metadata '
        'column, not the original values.'.format(column_name))
    warnings.warn(warning, UserWarning)


def _volatility(output_dir, metadata, state_column, individual_id_column,
                default_group_column, default_metric, table, yscale,
                importances):
    if individual_id_column == state_column:
        raise ValueError('individual_id_column & state_column must be set to '
                         'unique values.')

    # verify that individual_id_column exists in metadata
    # other metadata columns are validated later (to ensure correct types)
    if individual_id_column is not None:
        individual_ids = metadata.get_column(
            individual_id_column).to_dataframe()

    is_feat_vol_plot = importances is not None

    if is_feat_vol_plot:
        # We don't want to include any MD columns in the metric select in the
        # feature volatility variant of the viz, except for the state vo
        state_md_col = metadata.get_column(state_column).to_dataframe()
        metadata = metadata.filter_columns(column_type='categorical')
        metadata = metadata.merge(qiime2.Metadata(state_md_col))

        # Compile first differences and other stats on feature data
        stats_chart_data = _summarize_feature_stats(table, state_md_col)
        stats_chart_data = importances.join(stats_chart_data, how='inner')
        qiime2.Metadata(stats_chart_data).save(
            os.path.join(output_dir, 'feature_metadata.tsv'))
        stats_chart_data = stats_chart_data.reset_index(drop=False)
        # convert np.nan to None (nans and vega don't mix)
        stats_chart_data = _convert_nan_to_none(stats_chart_data)

    # Convert table to metadata and merge, if present.
    if table is not None:
        table.index.name = 'id'
        table_md = qiime2.Metadata(table)
        metadata = metadata.merge(table_md)

    # Partition the metadata into constituent types and assign defaults.
    categorical = metadata.filter_columns(column_type='categorical')
    numeric = metadata.filter_columns(column_type='numeric')
    if default_group_column is None:
        default_group_column = list(categorical.columns.keys())[0]
    if default_metric is None:
        default_metric = list(numeric.columns.keys())[0]

    # Ensure the default_* columns are members of their respective groups.
    # This will raise a uniform framework error on our behalf if necessary.
    categorical.get_column(default_group_column)
    numeric.get_column(default_metric)

    # We don't need to do any additional validation on the
    # individual_id_column after this point, since it doesn't matter if it is
    # categorical, numeric, only one value, etc.

    # Verify states column is numeric
    states = metadata.get_column(state_column)
    if not isinstance(states, qiime2.NumericMetadataColumn):
        raise TypeError('state_column must be numeric.')

    # Verify that the state column has more than one value present
    uniq_states = states.to_series().unique()
    if len(uniq_states) < 2:
        raise ValueError('state_column must contain at least two unique '
                         'values.')

    group_columns = list(categorical.columns.keys())
    if individual_id_column and individual_id_column not in group_columns:
        group_columns += [individual_id_column]
        if individual_id_column not in metadata.columns.keys():
            metadata = metadata.merge(qiime2.Metadata(individual_ids))
    metric_columns = list(numeric.columns.keys())
    control_chart_data = metadata.to_dataframe()
    # convert np.nan to None (nans and vega don't mix)
    control_chart_data = _convert_nan_to_none(control_chart_data)

    if is_feat_vol_plot:
        metric_columns.remove(state_column)

    # If we made it this far that means we can let Vega do it's thing!
    vega_spec = render_spec_volatility(control_chart_data,
                                       (stats_chart_data if is_feat_vol_plot
                                        else None),
                                       individual_id_column,
                                       state_column, default_group_column,
                                       group_columns, default_metric,
                                       metric_columns, yscale,
                                       is_feat_vol_plot)

    # Order matters here - need to render the template *after* copying the
    # directory tree, otherwise we will overwrite the index.html
    metadata.save(os.path.join(output_dir, 'data.tsv'))
    copy_tree(os.path.join(TEMPLATES, 'volatility'), output_dir)
    index = os.path.join(TEMPLATES, 'volatility', 'index.html')
    q2templates.render(index, output_dir,
                       context={'vega_spec': vega_spec,
                                'is_feat_vol_plot': is_feat_vol_plot})


def volatility(output_dir: str,
               metadata: qiime2.Metadata,
               state_column: str,
               individual_id_column: str = None,
               default_group_column: str = None,
               default_metric: str = None,
               table: pd.DataFrame = None,
               yscale: str = 'linear') -> None:
    _volatility(output_dir, metadata, state_column, individual_id_column,
                default_group_column, default_metric, table, yscale,
                importances=None)


def plot_feature_volatility(output_dir: str,
                            table: pd.DataFrame,
                            importances: pd.DataFrame,
                            metadata: qiime2.Metadata,
                            state_column: str,
                            individual_id_column: str = None,
                            default_group_column: str = None,
                            yscale: str = 'linear',
                            importance_threshold: float = None,
                            feature_count: int = 100,
                            missing_samples: str = 'error') -> None:
    # validate importances index is superset of table columns
    if missing_samples == 'error':
        _validate_metadata_is_superset(importances, table.T)

    # filter table, importances based on importance threshold / feature count
    table, importances = _importance_filtering(
        table, importances, importance_threshold, feature_count)

    # default_metric should be whatever the most important feature is
    default_metric = importances.index[0]

    _volatility(output_dir, metadata, state_column, individual_id_column,
                default_group_column, default_metric, table, yscale,
                importances)


def nmit(table: pd.DataFrame, metadata: qiime2.Metadata,
         individual_id_column: str, corr_method: str = "kendall",
         dist_method: str = "fro") -> skbio.DistanceMatrix:

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


def maturity_index(ctx,
                   table,
                   metadata,
                   state_column,
                   group_by,
                   control,
                   individual_id_column=None,
                   estimator='RandomForestRegressor',
                   n_estimators=100,
                   test_size=0.5,
                   step=0.05,
                   cv=5,
                   random_state=None,
                   n_jobs=1,
                   parameter_tuning=False,
                   optimize_feature_selection=False,
                   stratify=False,
                   missing_samples='error',
                   feature_count=50):

    filter_samples = ctx.get_action('feature_table', 'filter_samples')
    filter_features = ctx.get_action('feature_table', 'filter_features')
    group_table = ctx.get_action('feature_table', 'group')
    heatmap = ctx.get_action('feature_table', 'heatmap')
    split = ctx.get_action('sample_classifier', 'split_table')
    fit = ctx.get_action('sample_classifier', 'fit_regressor')
    predict_test = ctx.get_action('sample_classifier', 'predict_regression')
    summarize_estimator = ctx.get_action('sample_classifier', 'summarize')
    scatter = ctx.get_action('sample_classifier', 'scatterplot')
    volatility = ctx.get_action('longitudinal', 'volatility')

    # we must perform metadata superset validation here before we start
    # slicing and dicing.
    md_as_frame = metadata.to_dataframe()
    if missing_samples == 'error':
        _validate_metadata_is_superset(md_as_frame, table.view(biom.Table))

    # Let's also validate metadata columns before we get busy
    _validate_input_columns(
        md_as_frame, individual_id_column, group_by, state_column, None)

    # train regressor on subset of control samples
    control_table, = filter_samples(
        table, metadata=metadata, where="{0}='{1}'".format(group_by, control))

    md_column = metadata.get_column(state_column)
    X_train, X_test = split(control_table, md_column, test_size, random_state,
                            stratify, missing_samples='ignore')

    sample_estimator, importance = fit(
        X_train, md_column, step, cv, random_state, n_jobs, n_estimators,
        estimator, optimize_feature_selection, parameter_tuning,
        missing_samples='ignore')

    # drop training samples from rest of dataset; we will predict all others
    control_ids = pd.DataFrame(index=X_train.view(biom.Table).ids())
    control_ids.index.name = 'id'
    control_ids = qiime2.Metadata(control_ids)
    test_table, = filter_samples(table, metadata=control_ids, exclude_ids=True)

    # predict test samples
    predictions, = predict_test(test_table, sample_estimator, n_jobs)

    # summarize estimator params
    summary, = summarize_estimator(sample_estimator)

    # only report accuracy on control test samples
    test_ids = X_test.view(biom.Table).ids()
    accuracy_md = metadata.filter_ids(test_ids).get_column(state_column)
    accuracy_results, = scatter(predictions, accuracy_md, 'ignore')

    # calculate MAZ score
    # merge is inner join by default, so training samples are dropped (good!)
    pred_md = metadata.merge(predictions.view(qiime2.Metadata)).to_dataframe()
    pred_md['prediction'] = pd.to_numeric(pred_md['prediction'])
    pred_md = _maz_score(
        pred_md, 'prediction', state_column, group_by, control)
    maz = '{0} MAZ score'.format(state_column)
    maz_scores = ctx.make_artifact('SampleData[RegressorPredictions]',
                                   pred_md[maz])

    # make heatmap
    # load importance data and sum rows (to average importances if there are
    # multiple scores).
    imp = importance.view(pd.DataFrame).sum(1)

    # filter importances by user criteria
    imp = imp.sort_values(ascending=False)
    if feature_count > 0:
        imp = imp[:feature_count]
    imp.name = 'importance'
    imp = qiime2.Metadata(imp.to_frame())

    # filter table to important features for viewing as heatmap
    table, = filter_features(table, metadata=imp)
    # make sure IDs match between table and metadata
    cluster_table, = filter_samples(table, metadata=metadata)
    # need to group table by two columns together, so do this ugly hack
    cluster_by = group_by + '-' + state_column
    md_as_frame[cluster_by] = (md_as_frame[group_by].astype(str) + '-' +
                               md_as_frame[state_column].astype(str))
    cluster_md = qiime2.CategoricalMetadataColumn(md_as_frame[cluster_by])
    cluster_table, = group_table(cluster_table, axis='sample',
                                 metadata=cluster_md, mode='median-ceiling')
    # group metadata to match grouped sample IDs and sort by group/column
    clust_md = md_as_frame.groupby(cluster_by).first()
    clust_md = clust_md.sort_values([group_by, state_column])
    # sort table using clustered/sorted metadata as guide
    sorted_table = cluster_table.view(biom.Table).sort_order(clust_md.index)
    sorted_table = ctx.make_artifact('FeatureTable[Frequency]', sorted_table)
    clustermap, = heatmap(sorted_table, cluster='features')

    # visualize MAZ vs. time (column)
    lineplots, = volatility(
        qiime2.Metadata(pred_md), state_column=state_column,
        individual_id_column=individual_id_column,
        default_group_column=group_by, default_metric=maz, yscale='linear')

    return (
        sample_estimator, importance, predictions, summary, accuracy_results,
        maz_scores, clustermap, lineplots)


def first_differences(metadata: qiime2.Metadata, state_column: str,
                      individual_id_column: str, metric: str,
                      replicate_handling: str = 'error',
                      baseline: float = None,
                      table: pd.DataFrame = None) -> pd.Series:

    # find metric in metadata or derive from table and merge into metadata
    metadata = _load_metadata(metadata)
    if table is not None:
        _validate_metadata_is_superset(metadata, table)
        metadata = _add_metric_to_metadata(table, metadata, metric)
    else:
        _validate_is_numeric_column(metadata, metric)

    # validate columns
    _validate_input_columns(
        metadata, individual_id_column, None, state_column, metric)

    return _first_differences(
        metadata, state_column, individual_id_column, metric,
        replicate_handling, baseline=baseline, distance_matrix=None)


def first_distances(distance_matrix: skbio.DistanceMatrix,
                    metadata: qiime2.Metadata, state_column: str,
                    individual_id_column: str, baseline: float = None,
                    replicate_handling: str = 'error') -> pd.Series:

    # load and validate metadata
    metadata = _load_metadata(metadata)
    _validate_metadata_is_superset(metadata, distance_matrix)

    # validate columns
    # "Distance" is not actually a metadata value, so don't validate metric!
    _validate_input_columns(
        metadata, individual_id_column, None, state_column, None)

    return _first_differences(
        metadata, state_column, individual_id_column, metric=None,
        replicate_handling=replicate_handling, baseline=baseline,
        distance_matrix=distance_matrix)


def feature_volatility(ctx,
                       table,
                       metadata,
                       state_column,
                       individual_id_column=None,
                       cv=5,
                       random_state=None,
                       n_jobs=1,
                       n_estimators=100,
                       estimator='RandomForestRegressor',
                       parameter_tuning=False,
                       missing_samples='error',
                       importance_threshold=None,
                       feature_count=100):
    regress = ctx.get_action('sample_classifier', 'regress_samples')
    # TODO: Add this back once filter_features can operate on a
    # FeatureTable[RelativeFrequency] artifact (see notes below)
    # filter_tab = ctx.get_action('feature_table', 'filter_features')
    relative = ctx.get_action('feature_table', 'relative_frequency')
    volatility = ctx.get_action('longitudinal', 'plot_feature_volatility')

    # this validation must be tested here ahead of supervised regression
    states = metadata.get_column(state_column)
    if not isinstance(states, qiime2.NumericMetadataColumn):
        raise TypeError('state_column must be numeric.')

    estimator, importances, predictions, summary, accuracy = regress(
        table, metadata=states, cv=cv, random_state=random_state,
        n_jobs=n_jobs, n_estimators=n_estimators, estimator=estimator,
        parameter_tuning=parameter_tuning, optimize_feature_selection=True,
        missing_samples=missing_samples)

    # filter table to important features and convert to relative frequency
    feature_md = importances.view(qiime2.Metadata)
    filtered_table, = relative(table=table)
    # TODO: use feature_table.relative_frequency to convert to transform
    # once feature_table.filter_features can accept a relative frequency table.
    # We must transform and then filter, because otherwise relative frequencies
    # are based only on filtered features, which can seriously distort
    # frequency if a large number of features is filtered out. At which point
    # we can just filter with this code:
    # filtered_table, = filter_tab(table=filtered_table, metadata=feature_md)
    # filter_features also seems problematic, since it drops SAMPLES that have
    # none of the features to keep... distorting averages in the control chart.
    # For now let's use biom for this:
    filtered_table = filtered_table.view(biom.Table).filter(
        ids_to_keep=feature_md.get_ids(), axis='observation', inplace=False)
    filtered_table = ctx.make_artifact('FeatureTable[RelativeFrequency]',
                                       filtered_table)

    volatility_plot, = volatility(metadata=metadata,
                                  table=filtered_table,
                                  importances=importances,
                                  state_column=state_column,
                                  individual_id_column=individual_id_column,
                                  default_group_column=None,
                                  yscale='linear',
                                  importance_threshold=importance_threshold,
                                  feature_count=feature_count,
                                  missing_samples='ignore')

    return filtered_table, importances, volatility_plot, accuracy, estimator

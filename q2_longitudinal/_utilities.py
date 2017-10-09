# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from statsmodels.sandbox.stats.multicomp import multipletests
from itertools import combinations, cycle
import pkg_resources
import q2templates
from random import choice
from statsmodels.formula.api import mixedlm
from patsy import PatsyError
from math import ceil

import numpy as np
from scipy import linalg
from numpy.linalg.linalg import LinAlgError
from skbio.stats.distance import MissingIDError
from skbio import DistanceMatrix
from scipy.stats import (kruskal, mannwhitneyu, wilcoxon, ttest_ind, ttest_rel,
                         ttest_1samp, f_oneway)


TEMPLATES = pkg_resources.resource_filename('q2_longitudinal', 'assets')


def _validate_input_values(df, individual_id_column, group_column,
                           state_column, state_1, state_2):
    # confirm that different state values are input
    if state_1 is not None and state_1 == state_2:
        raise ValueError((
            'You have chosen the same value for state_1 and state_2. These '
            'parameters must be given different values.'))
    # confirm that individual, group, and state columns are in metadata
    _validate_input_columns(
        df, individual_id_column, group_column, state_column)
    # confirm that state_1 and state_2 exist in metadata and in each group in
    # group_column. Both checks are performed to give specific error messages.
    if state_1 is not None:
        for state in [state_1, state_2]:
            state = df[state_column].dtype.type(state)
            if state not in df[state_column].values:
                raise ValueError((
                    'State {0} not present in column {1} of metadata'.format(
                        state, state_column)))
            for group in df[group_column].unique():
                df1 = df[df[group_column] == group]
                if state not in df1[state_column].values:
                    raise ValueError((
                        'State {0} is not represented by any members of group '
                        '{1} in metadata. Consider using a different '
                        'group_column or state value.'.format(state, group)))


def _validate_input_columns(df, individual_id_column, group_column,
                            state_column):
    # confirm that individual, group, and state columns are in metadata
    if isinstance(group_column, list):
        cols = [individual_id_column, state_column] + group_column
    else:
        cols = [individual_id_column, group_column, state_column]
    for column in cols:
        if column is not None and column not in df.columns:
            raise ValueError('{0} is not a column in your metadata'.format(
                             column))


def _get_group_pairs(df, group_value, individual_id_column='SubjectID',
                     group_column='Group', state_column='time_point',
                     state_values=['1', '2'], replicate_handling='error'):
    results = []
    errors = []
    group_members = df[group_column] == group_value
    group_md = df[group_members]
    for individual_id in set(group_md[individual_id_column]):
        result = []
        for state_value in state_values:
            state_value = df[state_column].dtype.type(state_value)
            individual_id = df[individual_id_column].dtype.type(individual_id)
            _state = df[state_column] == state_value
            _ind = df[individual_id_column] == individual_id
            individual_at_state_idx = group_md[_state & _ind].index
            if len(individual_at_state_idx) > 1:
                errors.append(
                    "Multiple values for {0} {1} at {2} {3} ({4})".format(
                        individual_id_column, individual_id, state_column,
                        state_value,
                        ' '.join(map(str, individual_at_state_idx))))
                if replicate_handling == 'error':
                    raise ValueError((
                        'Replicate values for individual {0} at state {1}. '
                        'Remove replicate values from input files or set '
                        'replicate_handling parameter to select how '
                        'replicates are handled.'))
                elif replicate_handling == 'random':
                    result.append(choice(individual_at_state_idx))
                elif replicate_handling == 'drop':
                    pass
            elif len(individual_at_state_idx) == 0:
                errors.append("No values for {0} {1} at {2} {3}".format(
                    individual_id_column, individual_id, state_column,
                    state_value))
                pass
            else:
                result.append(individual_at_state_idx[0])
        if len(result) == len(state_values):
            results.append(tuple(result))
    return results, errors


def _extract_distance_distribution(distance_matrix: DistanceMatrix, pairs,
                                   df, individual_id_column, group_column):
    result = []
    pairs_summary = []
    for p in pairs:
        try:
            dist = distance_matrix[p]
            result.append(dist)
            individual_id = df[individual_id_column][p[0]]
            group = df[group_column][p[0]]
            pairs_summary.append((individual_id, dist, group))
        except MissingIDError:
            pass
    pairs_summary = pd.DataFrame(
        pairs_summary, columns=['SubjectID', 'Distance', 'Group'])
    pairs_summary.set_index('SubjectID', inplace=True)
    return result, pairs_summary


def _between_subject_distance_distribution(
        distance_matrix: DistanceMatrix, pairs, metadata=None,
        group_column=None, group_values=None):
    '''Extract distance distribution between all sample pairs in distance
    matrix, if samples are not listed in pairs; optionally, also if both
    samples are members of group_column's group_values in metadata
    pd.DataFrame.
    pairs: list of tuples
        paired samples to exclude from "between" distribution
    metadata: pd.DataFrame
        sample metadata
    group_column: str
        dataframe column to search
    group_values: list
        groups to include
    '''
    results = dict()
    for i in distance_matrix.ids:
        # ignore sample if belongs to group
        if group_column and group_values:
            if metadata.loc[i][group_column] not in group_values:
                continue
        for j in distance_matrix.ids:
            # ignore sample if belongs to group
            if group_column and group_values:
                if metadata.loc[j][group_column] not in group_values:
                    continue
            # ignore sample if i = j (self match)
            if i == j:
                continue
            if ((i, j) not in pairs and
                    (j, i) not in pairs and
                    (j, i) not in results):
                try:
                    results[i, j] = distance_matrix[i, j]
                except MissingIDError:
                    pass
    return list(results.values())


def _get_pairwise_differences(df, pairs, category, individual_id_column,
                              group_column):
    result = []
    pairs_summary = []
    for pre_idx, post_idx in pairs:
        individual_id = df[individual_id_column][pre_idx]
        group = df[group_column][pre_idx]
        pre_value = float(df[category][pre_idx])
        post_value = float(df[category][post_idx])
        paired_difference = post_value - pre_value
        if not np.isnan(paired_difference):
            result.append(paired_difference)
            pairs_summary.append((individual_id, paired_difference, group))
    pairs_summary = pd.DataFrame(
        pairs_summary, columns=['SubjectID', 'Difference', 'Group'])
    pairs_summary.set_index('SubjectID', inplace=True)
    return result, pairs_summary


def _compare_pairwise_differences(groups, parametric=False):
    pvals = []
    stat = 'W (wilcoxon signed-rank test)'
    for name, values in groups.items():
        try:
            if parametric:
                t, p = ttest_1samp(values, 0.0)
                stat = 't (one-sample t-test)'
            else:
                t, p = wilcoxon(values)
        except ValueError:
            # if test fails (e.g., because of zero variance), just skip
            pass
        pvals.append((name, t, p))
    result = pd.DataFrame(pvals, columns=["Group", stat, "P-value"])
    result.set_index(["Group"], inplace=True)
    result = _multiple_tests_correction(result, method='fdr_bh')
    return result


def _multiple_group_difference(groups, parametric=False):
    '''groups: list of lists of values.'''
    if parametric:
        stat, p_val = f_oneway(*groups)
        stat_name = 'F'
        test_name = 'One-way ANOVA'
    else:
        stat, p_val = kruskal(*groups, nan_policy='omit')
        stat_name = 'H'
        test_name = 'Kruskal Wallis test'
    multiple_group_test = pd.Series(
        [stat, p_val], index=[stat_name, 'P value'], name=test_name)
    return multiple_group_test


def _per_method_pairwise_stats(groups, paired=False, parametric=True):
    '''Perform mann whitney U tests between group distance distributions,
    followed by FDR correction. Returns pandas dataframe of p-values.
    groups: dict
        {group_names: [distribution of values to compare]}
    paired: bool
        Perform Wilcoxon signed rank test instead of Mann Whitney U. groups
        values must be ordered such that paired samples appear in the same
        order in each group.
    '''
    pvals = []
    combos = [a for a in combinations(groups.keys(), 2)]
    for a in combos:
        try:
            if not paired and not parametric:
                stat_name = 'Mann-Whitney U'
                u, p = mannwhitneyu(
                    groups[a[0]], groups[a[1]], alternative='two-sided')

            elif not paired and parametric:
                stat_name = 't (two-sample t-test)'
                u, p = ttest_ind(
                    groups[a[0]], groups[a[1]], nan_policy='raise')

            elif paired and not parametric:
                stat_name = 'W (wilcoxon signed-rank test)'
                u, p = wilcoxon(groups[a[0]], groups[a[1]])

            else:
                stat_name = 't (paired t-test)'
                u, p = ttest_rel(
                    groups[a[0]], groups[a[1]], nan_policy='raise')

            pvals.append((a[0], a[1], u, p))
        except ValueError:
            # if test fails (e.g., because of zero variance), just skip
            pass

    result = pd.DataFrame(
        pvals, columns=["Group A", "Group B", stat_name, "P-value"])
    result.set_index(['Group A', 'Group B'], inplace=True)
    result = _multiple_tests_correction(result, method='fdr_bh')

    return result


def _linear_effects(metadata, metric, state_column, group_categories,
                    individual_id_column):
    # format formula
    formula = "{0} ~ {1} * {2}".format(
        metric, state_column, " * ".join(group_categories))

    # generate model
    try:
        mlm = mixedlm(
            formula, metadata, groups=metadata[individual_id_column])

    # semicolon-delimited taxonomies cause an error; rename and run
    except (SyntaxError, PatsyError):
        new_metric = metric.split(';')[-1]
        metadata[new_metric] = metadata[metric]
        formula = "{0} ~ {1} * {2}".format(
            new_metric, state_column, " * ".join(group_categories))
        mlm = mixedlm(
            formula, metadata, groups=metadata[individual_id_column])

    # numpy.linalg.linalg.LinAlgError appears to raise
    # See https://github.com/qiime2/q2-longitudinal/issues/39
    try:
        model_fit = mlm.fit()
    except LinAlgError:
        raise ValueError(
            'Linear model will not compute due to singular matrix error. '
            'This may occur if input variables correlate closely or exhibit '
            'zero variance. Please check your input variables. Removing '
            'potential covariates may resolve this issue.')

    # summarize model and format prettily
    model_summary, model_results = model_fit.summary().tables
    model_summary = pd.Series(
        data=list(model_summary[1].values) + list(model_summary[3].values),
        index=list(model_summary[0].values) + list(model_summary[2].values),
        name='model summary').to_frame()

    return model_summary, model_results


def _boxplot_from_dict(groups, hue=None, y_label=None, x_label=None,
                       y_min=None, y_max=None, palette=None,
                       label_rotation=45):
    """Generate boxplot of metric (y) by groups (x), input as dict of lists of
    values.

    hue, color variables all pass directly to equivalently named
        variables in seaborn.boxplot().
    """
    x_tick_labels = [v for k, v in sorted(
        _add_sample_size_to_xtick_labels(groups).items())]
    vals = [v for k, v in sorted(groups.items())]

    ax = sns.boxplot(data=vals, hue=hue, palette=palette)
    ax.set_ylim(bottom=y_min, top=y_max)
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    ax.set_xticklabels(x_tick_labels, rotation=label_rotation)
    return ax


def _add_sample_size_to_xtick_labels(groups):
    '''groups is a dict of lists of values.'''
    x_tick_labels = {k: '{0} (n={1})'.format(k, len(v))
                     for k, v in groups.items()}
    return x_tick_labels


def _lmplot_from_dataframe(state_column, metric, metadata, group_by,
                           lowess=False, ci=95, palette='Set1'):
    g = sns.lmplot(state_column, metric, data=metadata, hue=group_by,
                   fit_reg=True, scatter_kws={"marker": ".", "s": 100},
                   legend=False, lowess=lowess, ci=ci, palette=palette)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    return g


def _regplot_subplots_from_dataframe(state_column, metric, metadata,
                                     group_by, lowess=False, ci=95,
                                     palette='Set1'):
    '''plot a single regplot for each group in group_by.'''
    num_groups = len(group_by)
    height = 6 * num_groups
    f, axes = plt.subplots(num_groups, figsize=(6, height))
    for num in range(num_groups):
        if num_groups > 1:
            ax = axes[num]
        else:
            ax = axes
        sns.set_palette(palette)
        for group in metadata[group_by[num]].unique():
            subset = metadata[metadata[group_by[num]] == group]
            sns.regplot(state_column, metric, data=subset, fit_reg=True,
                        scatter_kws={"marker": ".", "s": 100}, label=group,
                        ax=ax, lowess=lowess, ci=ci)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    return f


def _control_chart_subplots(state_column, metric, metadata, group_column,
                            individual_id_column, ci=95, palette='Set1',
                            plot_control_limits=True, xtick_interval=None,
                            yscale='linear', spaghetti=False):

    groups = metadata[group_column].unique()
    states = metadata[state_column].unique()
    fig_count = len(groups) + 1
    fig_height = fig_count * 6
    chart, axes = plt.subplots(fig_count, figsize=(6, fig_height))

    # determine x tick interval: autoscale so that ≤ 30 labels appear
    xtick_interval = _set_xtick_interval(xtick_interval, states)

    # plot individual groups' control charts
    colors = cycle(sns.color_palette(palette, n_colors=len(groups)))
    num = 1
    for group, group_md in metadata.groupby(group_column):
        color=next(colors)
        c, gm, gs = _control_chart(
            state_column, metric, group_md, None, ci=ci, legend=False,
            color=color, plot_control_limits=True, ax=axes[num],
            palette=None, xtick_interval=xtick_interval)
        c.set_title('{0}: {1}'.format(group_column, group))
        if spaghetti:
            # plot group's sphaghetti on main plot and current subplot
            for ax in [0, num]:
                c = _make_spaghetti(
                    group_md, state_column, metric, individual_id_column,
                    ax=axes[ax], color=color, alpha=0.3)
        c = _set_xticks(c, group_md, state_column, states, xtick_interval)
        num += 1

    # plot all groups together, compare variances
    c, global_mean, global_std = _control_chart(
        state_column, metric, metadata, group_column, ci=ci, palette=palette,
        plot_control_limits=plot_control_limits, ax=axes[0],
        xtick_interval=xtick_interval)
    c.set_title('Group volatility comparison plot')
    c = _set_xticks(c, metadata, state_column, states, xtick_interval)

    plt.yscale(yscale)

    plt.tight_layout()

    return chart, global_mean, global_std


def _make_spaghetti(metadata, state_column, metric, individual_id_column,
                    ax, color=None, alpha=1.0):
    for ind, ind_data in metadata.groupby(individual_id_column):
        ax.plot(ind_data[state_column], ind_data[metric], alpha=alpha, c=color,
                label='_nolegend_')
    return ax


def _set_xtick_interval(xtick_interval, states):
    if xtick_interval is None:
        if len(states) > 30:
            xtick_interval = ceil(len(states) / 30)
        else:
            xtick_interval = 1
    return xtick_interval


def _set_xtick_labels(metadata, state_column, states, xtick_interval):
    # pull x-axis labels from array of unique states, slice xtick_interval
    x_tick_labels = {state: metadata[metadata[state_column] == state]
                     for state in states[::xtick_interval]}
    # sort labels by key (state), add sample size to labels.
    # x_tick_labels[0] must be empty because this label does not appear in plot
    x_tick_labels = [''] + [v for k, v in sorted(
        _add_sample_size_to_xtick_labels(x_tick_labels).items())]
    return x_tick_labels


def _set_xticks(ax, metadata, state_column, states, xtick_interval):
    # generate x-axis tick labels
    x_tick_labels = _set_xtick_labels(
        metadata, state_column, states, xtick_interval)
    # Find tick locations at set at xtick_interval
    ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_interval))
    # add new labels and rotate
    ax.set_xticklabels(x_tick_labels, rotation=90)
    return ax


def _pointplot_from_dataframe(state_column, metric, metadata, group_by,
                              ci=95, palette='Set1', ax=None, legend=True,
                              color=None, xtick_interval=None):

    g = sns.pointplot(data=metadata, x=state_column, y=metric, hue=group_by,
                      ci=ci, palette=palette, color=color, ax=ax, markers=".")

    # place legend to right side of plot
    if legend:
        if ax is not None:
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        else:
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    return g


def _calculate_variability(metadata, metric):
    global_mean = metadata[metric].mean()
    global_std = metadata[metric].std()
    upper_limit = global_mean + global_std * 3
    lower_limit = global_mean - global_std * 3
    upper_warning = global_mean + global_std * 2
    lower_warning = global_mean - global_std * 2
    return (global_mean, global_std, upper_limit, lower_limit, upper_warning,
            lower_warning)


def _multiple_tests_correction(df, method='fdr_bh', sort=True):
    try:
        df['FDR P-value'] = multipletests(df['P-value'], method=method)[1]
    # zero division error will occur if the P-value series is empty. Ignore.
    except ZeroDivisionError:
        pass
    if sort:
        df.sort_index(inplace=True)
    return df


def _control_chart(state_column, metric, metadata, group_by, ci=95,
                   palette='Set1', plot_control_limits=True, ax=None,
                   color=None, legend=True, xtick_interval=None):
    g = _pointplot_from_dataframe(state_column, metric, metadata, group_by,
                                  ci=95, palette=palette, ax=ax, legend=legend,
                                  color=color, xtick_interval=xtick_interval)

    m, stdev, ul, ll, uw, lw = _calculate_variability(metadata, metric)
    if plot_control_limits:
        for lm, ls in [(m, '-'), (ul, '--'), (ll, '--'), (uw, ':'), (lw, ':')]:
            g.axes.plot(g.get_xlim(), [lm, lm], ls=ls, c='grey')
    return g, m, stdev


def _load_metadata(metadata):
    metadata = metadata.to_dataframe()
    metadata = metadata.apply(lambda x: pd.to_numeric(x, errors='ignore'))
    return metadata


def _add_metric_to_metadata(table, metadata, metric):
    '''find metric in metadata or derive from table and merge into metadata.'''
    metadata = _load_metadata(metadata)
    if metric not in metadata.columns:
        if table is not None and metric in table.columns:
            metadata = pd.concat(
                [metadata, pd.DataFrame(table[metric])], axis=1, join='inner')
        else:
            raise ValueError(
                'metric must be a valid metadata or feature table column.')
    return metadata


def _visualize(output_dir, multiple_group_test=False, pairwise_tests=False,
               paired_difference_tests=False, plot=False, summary=False,
               errors=False, model_summary=False, model_results=False,
               raw_data=False, plot_name='Pairwise difference boxplot',
               pairwise_test_name='Pairwise group comparison tests'):

    pd.set_option('display.max_colwidth', -1)

    if summary is not False:
        summary = q2templates.df_to_html(summary.to_frame())

    if multiple_group_test is not False:
        multiple_group_test = multiple_group_test.to_frame().transpose()
        multiple_group_test = q2templates.df_to_html(multiple_group_test)

    if pairwise_tests is not False:
        pairwise_tests.to_csv(join(output_dir, 'pairwise_tests.tsv'), sep='\t')
        pairwise_tests = q2templates.df_to_html(pairwise_tests)

    if raw_data is not False:
        raw_data.to_csv(join(output_dir, 'raw-data.tsv'), sep='\t')
        raw_data = True

    if paired_difference_tests is not False:
        paired_difference_tests.to_csv(join(
            output_dir, 'paired_difference_tests.tsv'), sep='\t')
        paired_difference_tests = q2templates.df_to_html(
            paired_difference_tests)

    if model_summary is not False:
        model_summary.to_csv(join(output_dir, 'model_summary.tsv'), sep='\t')
        model_summary = q2templates.df_to_html(model_summary)

    if model_results is not False:
        model_results.to_csv(join(output_dir, 'model_results.tsv'), sep='\t')
        model_results = q2templates.df_to_html(model_results)

    if plot is not False:
        plot.savefig(join(output_dir, 'plot.png'), bbox_inches='tight')
        plot.savefig(join(output_dir, 'plot.pdf'), bbox_inches='tight')
        plt.close('all')

    index = join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context={
        'errors': errors,
        'summary': summary,
        'model_summary': model_summary,
        'model_results': model_results,
        'multiple_group_test': multiple_group_test,
        'pairwise_tests': pairwise_tests,
        'paired_difference_tests': paired_difference_tests,
        'plot': plot,
        'plot_name': plot_name,
        'raw_data': raw_data,
        'pairwise_test_name': pairwise_test_name,
    })


def _stats_and_visuals(output_dir, pairs, metric, group_column,
                       state_column, state_1, state_2,
                       individual_id_column, errors, parametric, palette,
                       replicate_handling,
                       multiple_group_test=True, pairwise_tests=True,
                       paired_difference_tests=True, boxplot=True,
                       plot_name='Pairwise difference boxplot',
                       pairwise_test_name='Pairwise group comparison tests'):
    # kruskal test or ANOVA between groups
    if multiple_group_test:
        multiple_group_test = _multiple_group_difference(
            pairs.values(), parametric=parametric)

    # pairwise testing
    if pairwise_tests:
        pairwise_tests = _per_method_pairwise_stats(
            pairs, paired=False, parametric=parametric)

    # paired difference tests
    if paired_difference_tests:
        paired_difference_tests = _compare_pairwise_differences(
            pairs, parametric=parametric)

    # boxplots
    if boxplot:
        boxplot = _boxplot_from_dict(
            pairs, palette=palette, y_label=metric, x_label=group_column)
    boxplot = boxplot.get_figure()

    summary = pd.Series(
        [metric, group_column, state_column, state_1, state_2,
         individual_id_column, parametric, replicate_handling],
        index=['Metric', 'Group column', 'State column', 'State 1',
               'State 2', 'Individual ID column', 'Parametric',
               'Replicates handling'],
        name='Paired difference tests')

    _visualize(output_dir, multiple_group_test, pairwise_tests,
               paired_difference_tests, boxplot, summary=summary,
               errors=errors, plot_name=plot_name,
               pairwise_test_name=pairwise_test_name)


def _temporal_corr(table, individual_id, corr_method="kendall"):
    '''Create Temporal correlation from a feature table
    containing repeated measures samples.
    table: pd.DataFrame
        Feature table, rows are samples, columns are features
    individual_id: pd.Series
        sample ids, with the same length as table
    corr_method: str
        temporal correlation method, "kendall", "pearson", "spearman"
    '''

    table["individual_id"] = individual_id
    results = table.groupby(["individual_id"]).corr(method=corr_method)
    results = results.fillna(0)

    return results


def _temporal_distance(corr, id_set, dist_method="fro"):
    '''Calculate Distance Matrix from temporal correlation data.
    corr: pd.DataFrame
        table grouped by individual ids, this is the output from _temporal_corr
    id_set: pd.Series
        unique subject ids from individual_id with index attached
    dist_method: str
        method supported by scipy.linalg.norm parameter ord
    '''
    id_n = len(id_set)

    dist = np.zeros((id_n, id_n))
    for i, id_i in enumerate(id_set):
        for j, id_j in enumerate(id_set[:i]):
            dist[i, j] = dist[j, i] = linalg.norm(
                corr.loc[id_i] - corr.loc[id_j], ord=dist_method)
    return DistanceMatrix(dist, ids=id_set.index)


def _nmit(table, sample_md, individual_id_column, corr_method="kendall",
          dist_method="fro"):
    '''Function to perform nonparametric microbial interdependence test (nmit)
    test.
    table: pd.DataFrame
        Feature table, rows are samples, columns are features
    sample_md: pd.DataFrame
        Sample metadata
    individual_id_column: str
        Metadata column containing IDs for individual subjects
    corr_method: str
        temporal correlation method
    dist_method: str
        temporal distance method from numpy.linalg.norm, default is "fro"
    '''
    # full series of individual ids
    individual_id = sample_md[individual_id_column]

    # compile series of unique individuals with index retained
    # the goal here is to have ordered lists of ids that we eventually append
    # to distance matrix (for metadata extraction during permanova or other
    # follow-up tests/plotting), hence we do not use pd.unique()
    id_set = individual_id.drop_duplicates()

    # calculate correlation in individuals
    _corr = _temporal_corr(table, individual_id,  corr_method)

    # calculate distance between individuals
    _dist = _temporal_distance(_corr, id_set, dist_method)

    return _dist


def _validate_metadata_is_superset(metadata, table):
    metadata_ids = set(metadata.index.tolist())
    table_ids = set(table.index.tolist())
    if not table_ids.issubset(metadata_ids):
        raise ValueError('Missing samples in metadata: %r' %
                         table_ids.difference(metadata_ids))

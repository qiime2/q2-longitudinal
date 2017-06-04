# ----------------------------------------------------------------------------
# Copyright (c) 2017-, q2-intervention development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from os.path import join
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.sandbox.stats.multicomp import multipletests
from itertools import combinations
import pkg_resources
import q2templates
from random import choice

import numpy as np
from skbio.stats.distance import MissingIDError
from skbio import DistanceMatrix
from scipy.stats import (kruskal, mannwhitneyu, wilcoxon, ttest_ind, ttest_rel,
                         ttest_1samp, f_oneway)


TEMPLATES = pkg_resources.resource_filename('q2_intervention', 'assets')


def _get_group_pairs(df, group_value, individual_id_category='SubjectID',
                     group_category='Group', state_category='time_point',
                     state_values=['pre', 'post'], drop_duplicates=True):
    results = []
    group_members = df[group_category] == group_value
    group_md = df[group_members]
    for individual_id in set(group_md[individual_id_category]):
        result = []
        for state_value in state_values:
            _state = df[state_category] == state_value
            _ind = df[individual_id_category] == individual_id
            individual_at_state_idx = group_md[_state & _ind].index
            if len(individual_at_state_idx) > 1:
                print("Multiple values for {0} {1} at {2} {3} ({4})".format(
                    individual_id_category, individual_id, state_category,
                    state_value, ' '.join(individual_at_state_idx)))
                if drop_duplicates:
                    break
                else:
                    individual_at_state_idx = [choice(individual_at_state_idx)]
            elif len(individual_at_state_idx) == 0:
                print("No values for {0} {1} at {2} {3}".format(
                    individual_id_category, individual_id, state_category,
                    state_value))
                break
            result.append(individual_at_state_idx[0])
        if len(result) == len(state_values):
            results.append(tuple(result))
    return results


def _extract_distance_distribution(distance_matrix: DistanceMatrix, pairs):
    result = []
    for p in pairs:
        try:
            result.append(distance_matrix[p])
        except MissingIDError:
            pass
    return result


def _between_subject_distance_distribution(
        distance_matrix: DistanceMatrix, pairs, metadata=None,
        group_category=None, group_values=None):
    '''Extract distance distribution between all sample pairs in distance
    matrix, if samples are not listed in pairs; optionally, also if neither
    sample is a not member of group_category's group_values in metadata
    pd.DataFrame.
    pairs: list of tuples
        paired samples to exclude from "between" distribution
    metadata: pd.DataFrame
        sample metadata
    group_category: str
        dataframe category to search
    group_values: list
        groups to exclude
    '''
    results = dict()
    for i in distance_matrix.ids:
        # ignore sample if belongs to group
        if group_category and group_values:
            if metadata.loc[i][group_category] in group_values:
                continue
        for j in distance_matrix.ids:
            # ignore sample if belongs to group
            if group_category and group_values:
                if metadata.loc[j][group_category] in group_values:
                    continue
            if ((i, j) not in pairs and
                    (j, i) not in pairs and
                    (j, i) not in results):
                try:
                    results[i, j] = distance_matrix[i, j]
                except MissingIDError:
                    pass
    return list(results.values())


def _get_paired_differences(df, pairs, category):
    result = []
    for pre_idx, post_idx in pairs:
        pre_value = float(df[category][pre_idx])
        post_value = float(df[category][post_idx])
        paired_difference = post_value - pre_value
        if not np.isnan(paired_difference):
            result.append(paired_difference)
    return result


def test_paired_differences(groups, parametric=True):
    pvals = []
    for name, values in groups.items():
        try:
            if parametric:
                t, p = ttest_1samp(values, 0.0)
            else:
                t, p = wilcoxon(values)
        except ValueError:
            # if test fails (e.g., because of zero variance), just skip
            pass
        pvals.append((name, t, p))
    result = pd.DataFrame(pvals, columns=["Group", "stat", "P"])
    result.set_index(["Group"], inplace=True)
    result['FDR P'] = multipletests(result['P'], method='fdr_bh')[1]
    result.sort_index(inplace=True)
    return result


def _multiple_group_difference(groups, parametric=True):
    '''groups: list of lists of values.'''
    if parametric:
        stat, p_val = f_oneway(*groups)
    else:
        stat, p_val = kruskal(*groups, nan_policy='omit')
    multiple_group_test = pd.Series(
        [stat, p_val], index=['test statistic', 'P value'],
        name='Multiple group test')
    return multiple_group_test


def per_method_pairwise_tests(groups, paired=False, parametric=True):
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
                u, p = mannwhitneyu(
                    groups[a[0]], groups[a[1]], alternative='two-sided')

            elif not paired and parametric:
                u, p = ttest_ind(
                    groups[a[0]], groups[a[1]], nan_policy='raise')

            elif paired and not parametric:
                u, p = wilcoxon(groups[a[0]], groups[a[1]])

            else:
                u, p = ttest_rel(
                    groups[a[0]], groups[a[1]], nan_policy='raise')

            pvals.append((a[0], a[1], u, p))
        except ValueError:
            # if test fails (e.g., because of zero variance), just skip
            pass

    result = pd.DataFrame(pvals, columns=["Group A", "Group B", "stat", "P"])
    result.set_index(['Group A', 'Group B'], inplace=True)
    try:
        result['FDR P'] = multipletests(result['P'], method='fdr_bh')[1]
    except ZeroDivisionError:
        pass
    result.sort_index(inplace=True)

    return result


def boxplot_from_dict(groups, hue=None, y_label=None, x_label=None, y_min=None,
                      y_max=None, pallette=None, label_rotation=45):
    """Generate boxplot of metric (y) by groups (x), input as dict of lists of
    values.

    hue, color variables all pass directly to equivalently named
        variables in seaborn.boxplot().
    """
    x_tick_labels = [k for k, v in sorted(groups.items())]
    vals = [v for k, v in sorted(groups.items())]

    ax = sns.boxplot(data=vals, hue=hue, palette=pallette)
    ax.set_ylim(bottom=y_min, top=y_max)
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    ax.set_xticklabels(x_tick_labels, rotation=label_rotation)
    return ax


def _visualize(output_dir, multiple_group_test=None, pairwise_tests=None,
               paired_difference_tests=None, plot=None, summary=None):

    pd.set_option('display.max_colwidth', -1)

    if summary is not None:
        summary = summary.to_frame().to_html(classes=(
            "table table-striped table-hover")).replace(
                'border="1"', 'border="0"')

    if multiple_group_test is not None:
        multiple_group_test = multiple_group_test.to_frame().to_html(classes=(
            "table table-striped table-hover")).replace(
                'border="1"', 'border="0"')

    if pairwise_tests is not None:
        pairwise_tests.to_csv(join(output_dir, 'pairwise_tests.tsv'), sep='\t')
        pairwise_tests = pairwise_tests.to_html(classes=(
            "table table-striped table-hover")).replace(
                'border="1"', 'border="0"')

    if paired_difference_tests is not None:
        paired_difference_tests.to_csv(join(
            output_dir, 'paired_difference_tests.tsv'), sep='\t')
        paired_difference_tests = paired_difference_tests.to_html(classes=(
            "table table-striped table-hover")).replace(
                'border="1"', 'border="0"')

    if plot is not None:
        plot.savefig(join(output_dir, 'plot.png'), bbox_inches='tight')
        plot.savefig(join(output_dir, 'plot.pdf'), bbox_inches='tight')
        plt.close('all')

    index = join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context={
        'summary': summary,
        'multiple_group_test': multiple_group_test,
        'pairwise_tests': pairwise_tests,
        'paired_difference_tests': paired_difference_tests,
        'plot': plot,
    })

# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .const import (
    LEG_CTRL_LABEL, LEG_CTRL_SYMBOL, SIG_CTRL_CHART_HEIGHT, SIG_COLOR_SCHEME,
    SIG_CTRL_CHART_WIDTH, SIG_CTRL_MEAN_LINE_THICKNESS,
    SIG_CTRL_MEAN_LINE_OPACITY, SIG_CTRL_MEAN_SYMBOL_SIZE,
    SIG_CTRL_MEAN_SYMBOL_OPACITY, SIG_CTRL_SPG_LINE_THICKNESS,
    SIG_CTRL_SPG_LINE_OPACITY, SIG_CTRL_SPG_SYMBOL_SIZE,
    SIG_CTRL_SPG_SYMBOL_OPACITY, SIG_WIDTH, SIG_SHOW_ERROR_BARS, SIG_METRIC,
    SIG_GROUP, SIG_SHOW_GLOBAL_MEAN, SIG_SHOW_GLOBAL_CTRL_LIMS,
    SIG_STATS_CHART_WIDTH, SIG_STATS_CHART_HEIGHT, DAT_STATS, SIG_STATS_LEFT,
    SIG_STATS_RIGHT, SIG_STATS_SORT, SIG_STATS_SORT_DIR, MRK_STATS,
    MRK_STATS_CIRCLES)


def render_signals_ctrl(default_group, group_columns, default_metric,
                        metric_columns, is_feat_vol_plot):
    _signals_ctrl = [
        # LAYOUT/DIMENSIONS
        {'name': SIG_WIDTH, 'value': '', 'bind': {'input': 'text'},
         'on': [{'events': {'source': 'window', 'type': 'resize'},
                 'update': 'containerSize()[0]'}]},
        {'name': SIG_CTRL_CHART_HEIGHT, 'value': 400},
        {'name': SIG_CTRL_CHART_WIDTH, 'update': '[0, %s]' % SIG_WIDTH},

        # UI WIDGETS
        {'name': SIG_SHOW_ERROR_BARS, 'value': False,
         'bind': {'input': 'checkbox', 'element': '#toggle-error-bars'}},
        {'name': SIG_GROUP,
         'value': default_group,
         'bind': {'input': 'select', 'element': '#group-column',
                  'options': group_columns}},
        {'name': SIG_METRIC, 'value': default_metric,
         'bind': {'input': 'select', 'element': '#metric-column',
                  'options': metric_columns}},
        {'name': SIG_SHOW_GLOBAL_MEAN, 'value': False,
         'bind': {'input': 'checkbox', 'element': '#toggle-global-mean'}},
        {'name': SIG_SHOW_GLOBAL_CTRL_LIMS, 'value': False,
         'bind': {'input': 'checkbox',
                  'element': '#toggle-global-control-limits'}},
        {'name': SIG_CTRL_MEAN_LINE_THICKNESS, 'value': 3,
         'bind': {'input': 'range', 'min': 0.1, 'max': 10, 'step': 0.1,
                  'element': '#mean-line-thickness'}},
        {'name': SIG_CTRL_MEAN_LINE_OPACITY, 'value': 1.0,
         'bind': {'input': 'range', 'min': 0.0, 'max': 1.0, 'step': 0.01,
                  'element': '#mean-line-opacity'}},
        {'name': SIG_CTRL_MEAN_SYMBOL_SIZE, 'value': 50.0,
         'bind': {'input': 'range', 'min': 0.0, 'max': 500.0, 'step': 1.0,
                  'element': '#mean-symbol-size'}},
        {'name': SIG_CTRL_MEAN_SYMBOL_OPACITY, 'value': 0.0,
         'bind': {'input': 'range', 'min': 0.0, 'max': 1.0, 'step': 0.01,
                  'element': '#mean-symbol-opacity'}},
        {'name': SIG_COLOR_SCHEME, 'value': 'category10',
         'bind': {'input': 'select', 'element': '#color-scheme',
                  'options': ['accent', 'category10', 'category20',
                              'category20b', 'category20c', 'dark2', 'paired',
                              'pastel1', 'pastel2', 'set1', 'set2', 'set3',
                              'tableau10', 'tableau20']}},

        # LEGEND EVENTS
        {'name': 'clear', 'value': True,
         'on': [{'events': 'mouseup[!event.item]', 'update': 'true',
                 'force': True}]},
        {'name': 'shift', 'value': False,
         'on': [{'events': '@%s:click, @%s:click' %
                           (LEG_CTRL_SYMBOL, LEG_CTRL_LABEL),
                 'update': 'event.shiftKey', 'force': True}]},
        {'name': 'clicked', 'value': None,
         'on': [{'events': '@%s:click, @%s:click' %
                           (LEG_CTRL_SYMBOL, LEG_CTRL_LABEL),
                 'update': '{value: datum.value}', 'force': True}]},
    ]
    # only render feature metadata stats if feat vol
    if is_feat_vol_plot:
        _signals_ctrl[5]['on'] = [
            {'events': '@%s:click' % MRK_STATS,
             'update': 'datum.id', 'force': True},
            {'events': '@%s:click' % MRK_STATS_CIRCLES,
             'update': 'datum.id', 'force': True}]
    return _signals_ctrl


def render_signals_ctrl_individual():
    return [
        {'name': SIG_CTRL_SPG_LINE_THICKNESS, 'value': 0.5,
         'bind': {'input': 'range', 'min': 0.1, 'max': 10, 'step': 0.1,
                  'element': '#spaghetti-line-thickness'}},
        {'name': SIG_CTRL_SPG_LINE_OPACITY, 'value': 0.5,
         'bind': {'input': 'range', 'min': 0.0, 'max': 1.0, 'step': 0.01,
                  'element': '#spaghetti-line-opacity'}},
        {'name': SIG_CTRL_SPG_SYMBOL_SIZE, 'value': 50.0,
         'bind': {'input': 'range', 'min': 0.0, 'max': 500.0, 'step': 1.0,
                  'element': '#spaghetti-symbol-size'}},
        {'name': SIG_CTRL_SPG_SYMBOL_OPACITY, 'value': 0.0,
         'bind': {'input': 'range', 'min': 0.0, 'max': 1.0, 'step': 0.01,
                  'element': '#spaghetti-symbol-opacity'}}]


def render_signals_stats(stat_opts, sort_opts):
    return [
        # LAYOUT/DIMENSIONS
        {'name': SIG_STATS_CHART_WIDTH,
         'update': '(%s / 2) - 25' % SIG_WIDTH},
        {'name': SIG_STATS_CHART_HEIGHT,
         'update': '10 * length(data("%s"))' % DAT_STATS},

        # UI WIDGETS
        {'name': SIG_STATS_LEFT,
         'value': 'importance',
         'bind': {
             'input': 'select', 'element': '#metric-stats-left',
             'options': stat_opts}},
        {'name': SIG_STATS_RIGHT,
         'value': 'Global Mean',
         'bind': {
             'input': 'select', 'element': '#metric-stats-right',
             'options': stat_opts}},
        {'name': SIG_STATS_SORT,
         'value': 'importance',
         'bind': {
             'input': 'select', 'element': '#sort-stats',
             'options': sort_opts}},
        {'name': SIG_STATS_SORT_DIR,
         'value': 'descending',
         'bind': {
             'input': 'select', 'element': '#sort-stats-dir',
             'options': ['ascending', 'descending']}}]

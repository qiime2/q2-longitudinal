# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .const import (
    DAT_AGG_BY, FLD_GROUP_BY, SCL_CTRL_COLOR, DAT_GLOBAL_VALS, FLD_MIN_X,
    DAT_INDIVIDUAL, FLD_CTRL_CL1, FLD_CTRL_CL2, FLD_CTRL_CL3, SCL_CTRL_X,
    DAT_SERIES, SIG_SHOW_GLOBAL_MEAN, STY_STROKE_2,
    FLD_MAX_X, FLD_CTRL_MEAN, FLD_CTRL_CL0, FLD_CTRL_CI1, DAT_SPAGHETTIS,
    FLD_STATS_AVG_CHANGE, MRK_STATS, MRK_STATS_CIRCLES,
    SCL_CTRL_Y, FLD_METRIC, SIG_CTRL_MEAN_LINE_OPACITY, SIG_CTRL_CHART_HEIGHT,
    SCL_STATS_Y, FLD_STATS_AVG_INC, FLD_STATS_AVG_DEC, SIG_METRIC,
    SIG_CTRL_MEAN_SYMBOL_OPACITY, FLD_CTRL_CI0, SIG_CTRL_MEAN_LINE_THICKNESS,
    SIG_CTRL_MEAN_SYMBOL_SIZE, TST_GROUP, SIG_CTRL_SPG_LINE_OPACITY,
    SIG_CTRL_SPG_SYMBOL_OPACITY, FLD_CTRL_COUNT, SIG_SHOW_ERROR_BARS,
    SIG_CTRL_SPG_SYMBOL_SIZE, SIG_CTRL_SPG_LINE_THICKNESS, SIG_GROUP,
    SIG_SHOW_GLOBAL_CTRL_LIMS, SIG_CTRL_CHART_WIDTH, STY_DASH_A, STY_DASH_B,
    SIG_STATS_CHART_WIDTH, SIG_STATS_CHART_HEIGHT, FLD_STATS_ID, DAT_STATS)


def render_marks_ctrl(yscale):
    return \
        {'description': 'Control Chart',
         'name': 'controlChart',
         'type': 'group',
         'encode': {
             'enter': {
                 'x': {'value': 0},
                 'y': {'value': 0},
                 'width': {'signal': SIG_CTRL_CHART_WIDTH},
                 'height': {'signal': SIG_CTRL_CHART_HEIGHT},
                }},
         'marks': [],
         'scales': [],
         'axes': [],
         'legends': []}


def render_marks_stats():
    return [
        {'description': 'Descriptive Stats (Left)',
         'name': 'statsChartLeft',
         'type': 'group',
         'encode': {
             'enter': {
                 'x': {'value': 0},
                 'y': {'signal': SIG_CTRL_CHART_HEIGHT, 'offset': 75},
                 'width': {'signal': SIG_STATS_CHART_WIDTH},
                 'height': {'signal': SIG_STATS_CHART_HEIGHT}}},
         'marks': [],
         'scales': [],
         'axes': [],
         'legends': []},
        {'description': 'Descriptive Stats (Right)',
         'name': 'statsChartRight',
         'type': 'group',
         'encode': {
             'enter': {
                 'x': {'signal': SIG_STATS_CHART_WIDTH, 'offset': 50},
                 'y': {'signal': SIG_CTRL_CHART_HEIGHT, 'offset': 75},
                 'width': {'signal': SIG_STATS_CHART_WIDTH},
                 'height': {'signal': SIG_STATS_CHART_HEIGHT}}},
         'marks': [],
         'scales': [],
         'axes': [],
         'legends': []}]


def render_marks_ctrl_global():
    return [
        # Global Mean
        {'type': 'rule',
         'from': {'data': DAT_GLOBAL_VALS},
         'encode': {
             'update': {
                 'strokeWidth': {'value': STY_STROKE_2},
                 'x': {'scale': SCL_CTRL_X, 'field': FLD_MIN_X},
                 'x2': {'scale': SCL_CTRL_X, 'field': FLD_MAX_X},
                 'y': {'scale': SCL_CTRL_Y, 'field': FLD_CTRL_MEAN},
                 'strokeOpacity': [
                     {'test': SIG_SHOW_GLOBAL_MEAN, 'value': 1.0},
                     {'value': 0.0}]}}},
        # Global confidence limit, -3x std dev
        {'type': 'rule',
         'from': {'data': DAT_GLOBAL_VALS},
         'encode': {
             'update': {
                 'strokeWidth': {'value': STY_STROKE_2},
                 'strokeDash': {'value': STY_DASH_A},
                 'x': {'scale': SCL_CTRL_X, 'field': FLD_MIN_X},
                 'x2': {'scale': SCL_CTRL_X, 'field': FLD_MAX_X},
                 'y': {'scale': SCL_CTRL_Y, 'field': FLD_CTRL_CL0},
                 'strokeOpacity': [
                     {'test': SIG_SHOW_GLOBAL_CTRL_LIMS,
                      'value': 1.0},
                     {'value': 0.0},
                 ]}}},
        # Global confidence limit, -2x std dev
        {'type': 'rule',
         'from': {'data': DAT_GLOBAL_VALS},
         'encode': {
             'update': {
                 'strokeWidth': {'value': STY_STROKE_2},
                 'strokeDash': {'value': STY_DASH_B},
                 'x': {'scale': SCL_CTRL_X, 'field': FLD_MIN_X},
                 'x2': {'scale': SCL_CTRL_X, 'field': FLD_MAX_X},
                 'y': {'scale': SCL_CTRL_Y, 'field': FLD_CTRL_CL1},
                 'strokeOpacity': [
                     {'test': SIG_SHOW_GLOBAL_CTRL_LIMS,
                      'value': 1.0},
                     {'value': 0.0},
                 ]}}},
        # Global confidence limit, +2x std dev
        {'type': 'rule',
         'from': {'data': DAT_GLOBAL_VALS},
         'encode': {
             'update': {
                 'strokeWidth': {'value': STY_STROKE_2},
                 'strokeDash': {'value': STY_DASH_A},
                 'x': {'scale': SCL_CTRL_X, 'field': FLD_MIN_X},
                 'x2': {'scale': SCL_CTRL_X, 'field': FLD_MAX_X},
                 'y': {'scale': SCL_CTRL_Y, 'field': FLD_CTRL_CL2},
                 'strokeOpacity': [
                     {'test': SIG_SHOW_GLOBAL_CTRL_LIMS,
                      'value': 1.0},
                     {'value': 0.0},
                 ]}}},
        # Global confidence limit, +3x std dev
        {'type': 'rule',
         'from': {'data': DAT_GLOBAL_VALS},
         'encode': {
             'update': {
                 'strokeWidth': {'value': STY_STROKE_2},
                 'strokeDash': {'value': STY_DASH_B},
                 'x': {'scale': SCL_CTRL_X, 'field': FLD_MIN_X},
                 'x2': {'scale': SCL_CTRL_X, 'field': FLD_MAX_X},
                 'y': {'scale': SCL_CTRL_Y, 'field': FLD_CTRL_CL3},
                 'strokeOpacity': [
                     {'test': SIG_SHOW_GLOBAL_CTRL_LIMS,
                      'value': 1.0},
                     {'value': 0.0},
                 ]}}}]


def render_marks_ctrl_grouped(state):
    datum_state = "datum['%s']" % state
    mean_signal = ('{"title": "group mean", "group": datum.%s,'
                   ' "state": datum["%s"], "count": datum.%s,'
                   ' "mean": datum.%s, "ci0": datum.%s, "ci1": datum.%s}'
                   % (FLD_GROUP_BY, state, FLD_CTRL_COUNT, FLD_CTRL_MEAN,
                      FLD_CTRL_CI0, FLD_CTRL_CI1))
    return [
        {'type': 'group',
         'from': {
             # Regroup by "group" column
             'facet': {
                 'name': DAT_SERIES,
                 'data': DAT_AGG_BY,
                 'groupby': FLD_GROUP_BY}},
         'marks': [
             # Per-group mean lines
             {'type': 'line',
              'from': {'data': DAT_SERIES},
              'sort': {'field': datum_state, 'order': 'ascending'},
              'encode': {
                  'update': {
                      'x': {'scale': SCL_CTRL_X, 'field': state},
                      'y': {'scale': SCL_CTRL_Y, 'field': FLD_CTRL_MEAN},
                      'stroke': {'scale': SCL_CTRL_COLOR,
                                 'field': FLD_GROUP_BY},
                      'strokeWidth': {'signal':
                                      SIG_CTRL_MEAN_LINE_THICKNESS},
                      'opacity': [
                          {'test': TST_GROUP,
                           'signal': SIG_CTRL_MEAN_LINE_OPACITY},
                          {'value': 0.0},
                      ]}}},
             # Per-group symbols
             {'type': 'symbol',
              'from': {'data': DAT_SERIES},
              'encode': {
                  'update': {
                      'tooltip': {'signal': mean_signal},
                      'x': {'scale': SCL_CTRL_X, 'field': state},
                      'y': {'scale': SCL_CTRL_Y, 'field': FLD_CTRL_MEAN},
                      'stroke': {'scale': SCL_CTRL_COLOR,
                                 'field': FLD_GROUP_BY},
                      'fill': {'scale': SCL_CTRL_COLOR,
                               'field': FLD_GROUP_BY},
                      'size': {'signal': SIG_CTRL_MEAN_SYMBOL_SIZE},
                      'opacity': [
                          {'test': TST_GROUP,
                           'signal': SIG_CTRL_MEAN_SYMBOL_OPACITY},
                          {'value': 0.0}]}}},
             # Per-group error bars
             {'type': 'rect',
              'from': {'data': DAT_SERIES},
              'encode': {
                  'update': {
                      'width': {'value': 2.0},
                      'x': {'scale': SCL_CTRL_X, 'field': state,
                            'band': 0.5},
                      'y': {'scale': SCL_CTRL_Y, 'field': FLD_CTRL_CI0},
                      'y2': {'scale': SCL_CTRL_Y, 'field': FLD_CTRL_CI1},
                      'fill': {'scale': SCL_CTRL_COLOR,
                               'field': FLD_GROUP_BY},
                      'opacity': [
                          {'test': '%s && (%s)' % (SIG_SHOW_ERROR_BARS,
                                                   TST_GROUP), 'value': 1.0},
                          {'value': 0.0}]}}}]}]


def render_marks_ctrl_individual(individual_id, state):
    datum_state = 'datum["%s"]' % state
    tooltip_expr = ('{"title": "spaghetti", "individual_id": '
                    ' datum["%s"], "group": datum.%s, "state": '
                    ' datum["%s"], "metric": datum.%s}' %
                    (individual_id, FLD_GROUP_BY, state, FLD_METRIC))
    return \
        {'type': 'group',
         'from': {
             # Regroup by "individual_id" column
             'facet': {
                 'name': DAT_SPAGHETTIS,
                 'data': DAT_INDIVIDUAL,
                 'groupby': individual_id}},
         'marks': [
             {'type': 'line', 'from': {'data': DAT_SPAGHETTIS},
              'sort': {'field': datum_state, 'order': 'ascending'},
              'encode': {
                  'update': {
                      'strokeWidth': {'signal':
                                      SIG_CTRL_SPG_LINE_THICKNESS},
                      'x': {'scale': SCL_CTRL_X, 'field': state},
                      'y': {'scale': SCL_CTRL_Y,
                            'field': FLD_METRIC},
                      'stroke': {'scale': SCL_CTRL_COLOR,
                                 'field': {'signal': SIG_GROUP}},
                      'opacity': [
                          {'test': TST_GROUP,
                           'signal': SIG_CTRL_SPG_LINE_OPACITY},
                          {'value': 0.0}]}}},
             # Need to add symbols into plot for mouseover
             # https://github.com/vega/vega-tooltip/issues/120
             {'type': 'symbol',
              'from': {'data': DAT_SPAGHETTIS},
              'encode': {
                  'update': {
                      'tooltip': {'signal': tooltip_expr},
                      'size': {'signal': SIG_CTRL_SPG_SYMBOL_SIZE},
                      'x': {'scale': SCL_CTRL_X, 'field': state},
                      'y': {'scale': SCL_CTRL_Y,
                            'field': FLD_METRIC},
                      'stroke': {'scale': SCL_CTRL_COLOR,
                                 'field': {'signal': SIG_GROUP}},
                      'fill': {'scale': SCL_CTRL_COLOR,
                               'field': {'signal': SIG_GROUP}},
                      'opacity': [
                          {'test': TST_GROUP,
                           'signal': SIG_CTRL_SPG_SYMBOL_OPACITY},
                          {'value': 0.0}]}}}]}


def render_marks_stats_bars(x_scale, selected_stat):
    test = '%s === "%s"' % (selected_stat, FLD_STATS_AVG_CHANGE)
    return [
        {'name': MRK_STATS,
         'type': 'rect',
         'from': {'data': DAT_STATS},
         'encode': {
             'enter': {'height': {'scale': SCL_STATS_Y, 'band': 1}},
             'hover': {'fill': {'value': '#f7f591'}},
             'update': {
                 'tooltip': {'signal': 'datum'},
                 'x': [{'test': test, 'scale': x_scale,
                        'field': FLD_STATS_AVG_DEC},
                       {'scale': x_scale, 'value': 0}],
                 'x2': [{'test': test, 'scale': x_scale,
                         'field': FLD_STATS_AVG_INC},
                        {'scale': x_scale,
                         'field': {'signal': selected_stat}}],
                 'y': {'scale': SCL_STATS_Y, 'field': FLD_STATS_ID},
                 'fill': [{'test': '%s === datum.%s' %
                           (SIG_METRIC, FLD_STATS_ID), 'value': '#59bbe5'},
                          {'value': '#d3d3d3'}]}}},
        {'name': MRK_STATS_CIRCLES,
         'type': 'symbol',
         'from': {'data': DAT_STATS},
         'encode': {
             'enter': {'size': 50},
             'hover': {'fill': {'value': '#f7f591'}},
             'update': {
                 'tooltip': {'signal': 'datum'},
                 'x': {'value': -10},
                 'y': {'scale': SCL_STATS_Y, 'field': FLD_STATS_ID,
                       'offset': 5},
                 'fill': [{'test': '%s === datum.%s' %
                           (SIG_METRIC, FLD_STATS_ID), 'value': '#59bbe5'},
                          {'value': '#ededed'}]}}},
        {'type': 'rule',
         'encode': {
             'update': {
                 'y': {'value': 0, 'scale': SCL_STATS_Y},
                 'y2': {'signal': SIG_STATS_CHART_HEIGHT},
                 'x': {'value': 0, 'scale': x_scale, 'offset': 0.5},
                 'x2': {'value': 0, 'scale': x_scale, 'offset': 0.5},
                 'strokeWidth': {'value': 0.25},
                }}}]

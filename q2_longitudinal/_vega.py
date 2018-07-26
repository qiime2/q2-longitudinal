# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json

import pandas as pd


def _render_volatility_spec(data: pd.DataFrame, individual_id: str, state: str,
                            default_group: str, group_columns: list,
                            default_metric: str, metric_columns: list,
                            yscale: str) -> str:
    # Double-quotes for many of the strings below, just so that we don't have
    # to escape the single quotes - it is already hard enough to read when
    # wrapped.
    opacity_test = ("!length(data('selected')) || "
                    "indata('selected', 'value', datum.value)")
    group_test = ("!length(data('selected')) || "
                  "indata('selected', 'value', datum.groupByVal)")
    error_bar_test = 'showErrorBars && (%s)' % group_test
    metric_signal = {'signal': 'metric'}
    group_signal = {'signal': 'grouper'}
    # This looks grosser than it is (you can't do variable assignment in a
    # vega expr, so no temp helper vars) - basically find the min and max
    # extents of the metric in question for the y-axis rendering, including
    # the 3x stdev (depending on the spread this could be beyond the metric's
    # own limits.
    domain_expr = ("[min(data('globalVals')[0].cl0,"
                   "data('globalVals')[0].minY),"
                   "max(data('globalVals')[0].cl3,"
                   "data('globalVals')[0].maxY)]")

    # These templates customize the tooltips
    mean_signal = ('{"title": "group mean", "group": datum.groupByVal,'
                   ' "state": datum["%s"], "count": datum.count,'
                   ' "mean": datum.mean, "ci0": datum.ci0, "ci1": datum.ci1}'
                   % state)

    marks = [
        {
            'type': 'rule',
            'from': {
                'data': 'globalVals',
            },
            'encode': {
                'update': {
                    'strokeWidth': {
                        'value': 2,
                    },
                    'x': {
                        'scale': 'x',
                        'field': 'minX',
                    },
                    'x2': {
                        'scale': 'x',
                        'field': 'maxX',
                    },
                    'y': {
                        'scale': 'y',
                        'field': 'mean',
                    },
                    'strokeOpacity': [
                        {
                            'test': 'showGlobalMean',
                            'value': 1.0,
                        },
                        {
                            'value': 0.0,
                        },
                    ],
                },
            },
        },
        {
            'type': 'rule',
            'from': {
                'data': 'globalVals',
            },
            'encode': {
                'update': {
                    'strokeWidth': {
                        'value': 2,
                    },
                    'strokeDash': {
                        'value': [8, 8],
                    },
                    'x': {
                        'scale': 'x',
                        'field': 'minX',
                    },
                    'x2': {
                        'scale': 'x',
                        'field': 'maxX',
                    },
                    'y': {
                        'scale': 'y',
                        'field': 'cl0',
                    },
                    'strokeOpacity': [
                        {
                            'test': 'showGlobalControlLimits',
                            'value': 1.0,
                        },
                        {
                            'value': 0.0,
                        },
                    ],
                },
            },
        },
        {
            'type': 'rule',
            'from': {
                'data': 'globalVals',
            },
            'encode': {
                'update': {
                    'strokeWidth': {
                        'value': 2,
                    },
                    'strokeDash': {
                        'value': [6, 2],
                    },
                    'x': {
                        'scale': 'x',
                        'field': 'minX',
                    },
                    'x2': {
                        'scale': 'x',
                        'field': 'maxX',
                    },
                    'y': {
                        'scale': 'y',
                        'field': 'cl1',
                    },
                    'strokeOpacity': [
                        {
                            'test': 'showGlobalControlLimits',
                            'value': 1.0,
                        },
                        {
                            'value': 0.0,
                        },
                    ],
                },
            },
        },
        {
            'type': 'rule',
            'from': {
                'data': 'globalVals',
            },
            'encode': {
                'update': {
                    'strokeWidth': {
                        'value': 2,
                    },
                    'strokeDash': {
                        'value': [6, 2],
                    },
                    'x': {
                        'scale': 'x',
                        'field': 'minX',
                    },
                    'x2': {
                        'scale': 'x',
                        'field': 'maxX',
                    },
                    'y': {
                        'scale': 'y',
                        'field': 'cl2',
                    },
                    'strokeOpacity': [
                        {
                            'test': 'showGlobalControlLimits',
                            'value': 1.0,
                        },
                        {
                            'value': 0.0,
                        },
                    ],
                },
            },
        },
        {
            'type': 'rule',
            'from': {
                'data': 'globalVals',
            },
            'encode': {
                'update': {
                    'strokeWidth': {
                        'value': 2,
                    },
                    'strokeDash': {
                        'value': [8, 8],
                    },
                    'x': {
                        'scale': 'x',
                        'field': 'minX',
                    },
                    'x2': {
                        'scale': 'x',
                        'field': 'maxX',
                    },
                    'y': {
                        'scale': 'y',
                        'field': 'cl3',
                    },
                    'strokeOpacity': [
                        {
                            'test': 'showGlobalControlLimits',
                            'value': 1.0,
                        },
                        {
                            'value': 0.0,
                        },
                    ],
                },
            },
        },
        {
            'type': 'group',
            'from': {
                'facet': {
                    'name': 'series',
                    'data': 'aggBy',
                    'groupby': 'groupByVal',
                },
            },
            'marks': [
                {
                    'type': 'line',
                    'from': {
                        'data': 'series',
                    },
                    'sort': {
                        'field': 'datum.%s' % state,
                        'order': 'ascending',
                    },
                    'encode': {
                        'update': {
                            'x': {
                                'scale': 'x',
                                'field': state,
                            },
                            'y': {
                                'scale': 'y',
                                'field': 'mean',
                            },
                            'stroke': {
                                'scale': 'color',
                                'field': 'groupByVal',
                            },
                            'strokeWidth': {
                                'signal': 'meanLineThickness',
                            },
                            'opacity': [
                                {
                                    'test': group_test,
                                    'signal': 'meanLineOpacity',
                                },
                                {
                                    'value': 0.0,
                                },
                            ],
                        },
                    },
                },
                # Need to add symbols into plot for mouseover
                # https://github.com/vega/vega-tooltip/issues/120
                {
                    'type': 'symbol',
                    'from': {
                        'data': 'series',
                    },
                    'encode': {
                        'update': {
                            'tooltip': {
                                'signal': mean_signal,
                            },
                            'x': {
                                'scale': 'x',
                                'field': state,
                            },
                            'y': {
                                'scale': 'y',
                                'field': 'mean',
                            },
                            'stroke': {
                                'scale': 'color',
                                'field': 'groupByVal',
                            },
                            'fill': {
                                'scale': 'color',
                                'field': 'groupByVal',
                            },
                            'size': {
                                'signal': 'meanSymbolSize',
                            },
                            'opacity': [
                                {
                                    'test': group_test,
                                    'signal': 'meanSymbolOpacity',
                                },
                                {
                                    'value': 0.0,
                                },
                            ],
                        },
                    },
                },
                {
                    'type': 'rect',
                    'from': {
                        'data': 'series',
                    },
                    'encode': {
                        'update': {
                            'width': {
                                'value': 2.0,
                            },
                            'x': {
                                'scale': 'x',
                                'field': state,
                                'band': 0.5,
                            },
                            'y': {
                                'scale': 'y',
                                'field': 'ci0',
                            },
                            'y2': {
                                'scale': 'y',
                                'field': 'ci1',
                            },
                            'fill': {
                                'scale': 'color',
                                'field': 'groupByVal',
                            },
                            'opacity': [
                                {
                                    'test': error_bar_test,
                                    'value': 1.0,
                                },
                                {
                                    'value': 0.0,
                                },
                            ],
                        },
                    },
                },
            ],
        },
    ]

    signals = [
        {
            'name': 'grouper',
            'value': default_group,
            'bind': {
                'input': 'select',
                'element': '#group-column',
                'options': group_columns,
            }
        },
        {
            'name': 'metric',
            'value': default_metric,
            'bind': {
                'input': 'select',
                'element': '#metric-column',
                'options': metric_columns,
            },
        },
        {
            'name': 'width',
            'value': '',
            'bind': {
                'input': 'text',
            },
            'on': [
                {
                    'events': {
                        'source': 'window',
                        'type': 'resize',
                    },
                    'update': 'containerSize()[0]',
                },
            ],
        },
        {
            'name': 'showErrorBars',
            'value': False,
            'bind': {
                'input': 'checkbox',
                'element': '#toggle-error-bars',
            },
        },
        {
            'name': 'showGlobalMean',
            'value': False,
            'bind': {
                'input': 'checkbox',
                'element': '#toggle-global-mean',
            },
        },
        {
            'name': 'showGlobalControlLimits',
            'value': False,
            'bind': {
                'input': 'checkbox',
                'element': '#toggle-global-control-limits',
            },
        },
        {
            'name': 'meanLineThickness',
            'value': 3,
            'bind': {
                'input': 'range',
                'min': 0.1,
                'max': 10,
                'step': 0.1,
                'element': '#mean-line-thickness',
            },
        },
        {
            'name': 'meanLineOpacity',
            'value': 1.0,
            'bind': {
                'input': 'range',
                'min': 0.0,
                'max': 1.0,
                'step': 0.01,
                'element': '#mean-line-opacity',
            },
        },
        {
            'name': 'meanSymbolSize',
            'value': 50.0,
            'bind': {
                'input': 'range',
                'min': 0.0,
                'max': 500.0,
                'step': 1.0,
                'element': '#mean-symbol-size',
            },
        },
        {
            'name': 'meanSymbolOpacity',
            'value': 0.0,
            'bind': {
                'input': 'range',
                'min': 0.0,
                'max': 1.0,
                'step': 0.01,
                'element': '#mean-symbol-opacity',
            },
        },
        {
            'name': 'colorScheme',
            'value': 'category10',
            'bind': {
                'input': 'select',
                'element': '#color-scheme',
                'options': [
                    'accent',
                    'category10',
                    'category20',
                    'category20b',
                    'category20c',
                    'dark2',
                    'paired',
                    'pastel1',
                    'pastel2',
                    'set1',
                    'set2',
                    'set3',
                    'tableau10',
                    'tableau20',
                ],
            }
        },
        {
            'name': 'clear',
            'value': True,
            'on': [
                {
                    'events': 'mouseup[!event.item]',
                    'update': 'true',
                    'force': True,
                },
            ],
        },
        {
            'name': 'shift',
            'value': False,
            'on': [
                {
                    'events': '@legendSymbol:click, @legendLabel:click',
                    'update': 'event.shiftKey',
                    'force': True,
                },
            ],
        },
        {
            'name': 'clicked',
            'value': None,
            'on': [
                {
                    'events': '@legendSymbol:click, @legendLabel:click',
                    'update': '{value: datum.value}',
                    'force': True,
                },
            ],
        }]

    if individual_id:
        spaghetti_signal = ('{"title": "spaghetti", "individual_id": '
                            'datum["%s"], "group": datum.groupByVal, "state": '
                            'datum["%s"], "metric": datum.metricVal}' %
                            (individual_id, state))
        marks.append({
                'type': 'group',
                'from': {
                    'facet': {
                        'name': 'spaghettis',
                        'data': 'individual',
                        'groupby': individual_id,
                    },
                },
                'marks': [
                    {
                        'type': 'line',
                        'from': {
                            'data': 'spaghettis',
                        },
                        'sort': {
                            'field': 'datum.%s' % state,
                            'order': 'ascending',
                        },
                        'encode': {
                            'update': {
                                'strokeWidth': {
                                    'signal': 'spaghettiLineThickness',
                                },
                                'x': {
                                    'scale': 'x',
                                    'field': state,
                                },
                                'y': {
                                    'scale': 'y',
                                    'field': metric_signal,
                                },
                                'stroke': {
                                    'scale': 'color',
                                    'field': group_signal,
                                },
                                'opacity': [
                                    {
                                        'test': group_test,
                                        'signal': 'spaghettiLineOpacity',
                                    },
                                    {
                                        'value': 0.0,
                                    },
                                ],
                            },
                        },
                    },
                    # Need to add symbols into plot for mouseover
                    # https://github.com/vega/vega-tooltip/issues/120
                    {
                        'type': 'symbol',
                        'from': {
                            'data': 'spaghettis',
                        },
                        'encode': {
                            'update': {
                                'tooltip': {
                                    'signal': spaghetti_signal,
                                },
                                'size': {
                                    'signal': 'spaghettiSymbolSize',
                                },
                                'x': {
                                    'scale': 'x',
                                    'field': state,
                                },
                                'y': {
                                    'scale': 'y',
                                    'field': metric_signal,
                                },
                                'stroke': {
                                    'scale': 'color',
                                    'field': group_signal,
                                },
                                'fill': {
                                    'scale': 'color',
                                    'field': group_signal,
                                },
                                'opacity': [
                                    {
                                        'test': group_test,
                                        'signal': 'spaghettiSymbolOpacity',
                                    },
                                    {
                                        'value': 0.0,
                                    },
                                ],
                            },
                        },
                    },
                ],
            })
        signals.extend([
            {
                'name': 'spaghettiLineThickness',
                'value': 0.5,
                'bind': {
                    'input': 'range',
                    'min': 0.1,
                    'max': 10,
                    'step': 0.1,
                    'element': '#spaghetti-line-thickness',
                },
            },
            {
                'name': 'spaghettiLineOpacity',
                'value': 0.5,
                'bind': {
                    'input': 'range',
                    'min': 0.0,
                    'max': 1.0,
                    'step': 0.01,
                    'element': '#spaghetti-line-opacity',
                },
            },
            {
                'name': 'spaghettiSymbolSize',
                'value': 50.0,
                'bind': {
                    'input': 'range',
                    'min': 0.0,
                    'max': 500.0,
                    'step': 1.0,
                    'element': '#spaghetti-symbol-size',
                },
            },
            {
                'name': 'spaghettiSymbolOpacity',
                'value': 0.0,
                'bind': {
                    'input': 'range',
                    'min': 0.0,
                    'max': 1.0,
                    'step': 0.01,
                    'element': '#spaghetti-symbol-opacity',
                },
            }])

    # This looks weird, but the idea is to render a unique string into the JSON
    # spec, that way later on when we render the datatable to JSON, we can
    # just inject the datatable JSON wholesale.
    DATA_PLACEHOLDER = '26b16f18-fd66-4531-bbcb-080beba01086'

    # Just a quick note, order doesn't matter here (JSON documents are not
    # ordered) - this will render out stochastically, which is fine - vega
    # knows what to do.
    spec = {
        # This `$schema` is only fetched as part of the interactive vega
        # editor, which opens up outside of the visualization - this doesn't
        # appear to create any kind of XHR side-effect when loading the
        # visualization in an offline context.
        '$schema': 'https://vega.github.io/schema/vega/v3.0.json',
        'autosize': {
            'type': 'fit-x',
            'contains': 'padding',
            'resize': True,
        },
        # These dimensions are here for when the viz is opened in the
        # Vega Editor.
        'width': 800,
        'height': 400,
        'signals': signals,
        'scales': [
            {
                'name': 'x',
                'type': 'linear',
                'range': 'width',
                'nice': True,
                'domain': {
                    'data': 'individual',
                    'field': state,
                    'sort': True,
                },
            },
            {
                'name': 'y',
                # Signal registration on this param is currently blocked by
                # https://github.com/vega/vega/issues/525, which is why this
                # setting is still a QIIME 2 param to this viz.
                'type': yscale,
                'range': 'height',
                'nice': True,
                'domain': {
                    'signal': domain_expr,
                    'sort': True,
                },
            },
            {
                'name': 'color',
                'type': 'ordinal',
                'range': {
                    'scheme': {
                        'signal': 'colorScheme',
                    },
                },
                'domain': {
                    'data': 'individual',
                    'field': group_signal,
                },
            },
        ],
        'axes': [
            {
                'orient': 'bottom',
                'scale': 'x',
                'title': state,
            },
            {
                'orient': 'left',
                'scale': 'y',
                'title': metric_signal,
            },
        ],
        'legends': [
            {
                'stroke': 'color',
                'title': group_signal,
                'encode': {
                    'symbols': {
                        'name': 'legendSymbol',
                        'interactive': True,
                        'update': {
                            'fill': {
                                'value': 'transparent',
                            },
                            'strokeWidth': {
                                'value': 2,
                            },
                            'opacity': [
                                {
                                    'test': opacity_test,
                                    'value': 1.0,
                                },
                                {
                                    'value': 0.15,
                                },
                            ],
                            'size': {
                                'value': 100,
                            },
                        },
                    },
                    'labels': {
                        'name': 'legendLabel',
                        'interactive': True,
                        'update': {
                            'opacity': [
                                {
                                    'test': opacity_test,
                                    'value': 1,
                                },
                                {
                                    'value': 0.25,
                                },
                            ],
                        },
                    },
                },
            },
        ],
        'marks': marks,
        'data': [
            {
                'name': 'individual',
                'values': DATA_PLACEHOLDER,
                'transform': [
                    {
                        'type': 'formula',
                        'as': 'groupByVal',
                        'expr': 'datum[grouper]',
                    },
                    {
                        'type': 'formula',
                        'as': 'metricVal',
                        'expr': 'datum[metric]',
                    },
                ],
            },
            {
                'name': 'globalVals',
                'source': 'individual',
                'transform': [
                    {
                        'type': 'aggregate',
                        'ops': [
                            'mean',
                            'min',
                            'max',
                            'stdev',
                            'min',
                            'max',
                        ],
                        'fields': [
                            metric_signal,
                            state,
                            state,
                            metric_signal,
                            metric_signal,
                            metric_signal,
                        ],
                        'as': [
                            'mean',
                            'minX',
                            'maxX',
                            'stdev',
                            'minY',
                            'maxY',
                        ]
                    },
                    {
                        'type': 'formula',
                        'as': 'cl0',
                        'expr': 'datum.mean - (3 * datum.stdev)'
                    },
                    {
                        'type': 'formula',
                        'as': 'cl1',
                        'expr': 'datum.mean - (2 * datum.stdev)'
                    },
                    {
                        'type': 'formula',
                        'as': 'cl2',
                        'expr': 'datum.mean + (2 * datum.stdev)'
                    },
                    {
                        'type': 'formula',
                        'as': 'cl3',
                        'expr': 'datum.mean + (3 * datum.stdev)'
                    },
                    {
                        'type': 'formula',
                        'as': 'ext',
                        'expr': '[datum.cl0, datum.cl3]',
                    },
                ],
            },
            {
                'name': 'aggBy',
                'source': 'individual',
                'transform': [
                    {
                        'type': 'aggregate',
                        'groupby': [
                            'groupByVal',
                            state,
                        ],
                        'ops': [
                            'mean',
                            # TODO: parameterize these intervals
                            # I don't see an easy way at the moment to define
                            # your own confidence interval in vega.
                            'ci0',
                            'ci1',
                            'count',
                        ],
                        'fields': [
                            metric_signal,
                            metric_signal,
                            metric_signal,
                            metric_signal,
                        ],
                        'as': [
                            'mean',
                            'ci0',
                            'ci1',
                            'count',
                        ],
                    },
                ],
            },
            {
                'name': 'selected',
                'on': [
                    {
                        'trigger': 'clear',
                        'remove': True
                    },
                    {
                        'trigger': '!shift',
                        'remove': True
                    },
                    {
                        'trigger': '!shift && clicked',
                        'insert': 'clicked'
                    },
                    {
                        'trigger': 'shift && clicked',
                        'toggle': 'clicked'
                    },
                ],
            },
        ],
    }
    rendered_spec = json.dumps(spec)
    rendered_spec = rendered_spec.replace("'", r"\'")
    rendered_data = data.to_json(orient='records')
    # convert metadata single quotes to unicode
    rendered_data = rendered_data.replace("'", r'\u0027')
    return rendered_spec.replace('"%s"' % DATA_PLACEHOLDER, rendered_data)

# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json

import pandas as pd

from .const import (
    FLD_STATS_AVG_CHANGE, FLD_STATS_AVG_DEC, FLD_STATS_AVG_INC, FLD_STATS_ID,
    SIG_METRIC, SCL_STATS_X_LEFT, SIG_STATS_LEFT, DAT_STATS_GLOB_EXT_LEFT,
    SCL_STATS_X_RIGHT, SIG_STATS_RIGHT, DAT_STATS_GLOB_EXT_RIGHT)
from .axis import render_axes_ctrl, render_axes_stats
from .legend import render_legends_ctrl
from .mark import (
    render_marks_ctrl, render_marks_stats, render_marks_ctrl_global,
    render_marks_ctrl_grouped, render_marks_ctrl_individual,
    render_marks_stats_bars)
from .signal import (
    render_signals_ctrl, render_signals_ctrl_individual, render_signals_stats)
from .scale import render_scales_ctrl, render_scales_stats
from .data import render_data_ctrl, render_data_stats


def render_subplot_ctrl(yscale, state):
    control_chart = render_marks_ctrl(yscale)
    control_chart['marks'] = render_marks_ctrl_global() + \
        render_marks_ctrl_grouped(state)
    control_chart['scales'] = render_scales_ctrl(state, yscale)
    control_chart['axes'] = render_axes_ctrl(state)
    control_chart['legends'] = render_legends_ctrl()

    return control_chart


def render_subplot_stats():
    # add in scales, axes, and marks for left plot
    stats_chart = render_marks_stats()
    stats_chart[0]['scales'] = render_scales_stats(
        SCL_STATS_X_LEFT, SIG_STATS_LEFT, DAT_STATS_GLOB_EXT_LEFT)
    stats_chart[0]['axes'] = render_axes_stats(SCL_STATS_X_LEFT,
                                               SIG_STATS_LEFT)
    stats_chart[0]['marks'] = render_marks_stats_bars(SCL_STATS_X_LEFT,
                                                      SIG_STATS_LEFT)
    # add in scales, axes, and marks for right plot
    stats_chart[1]['scales'] = render_scales_stats(
        SCL_STATS_X_RIGHT, SIG_STATS_RIGHT, DAT_STATS_GLOB_EXT_RIGHT)
    stats_chart[1]['axes'] = render_axes_stats(SCL_STATS_X_RIGHT,
                                               SIG_STATS_RIGHT)
    stats_chart[1]['marks'] = render_marks_stats_bars(SCL_STATS_X_RIGHT,
                                                      SIG_STATS_RIGHT)
    return stats_chart


def render_spec_volatility(control_chart_data: pd.DataFrame,
                           stats_chart_data: pd.DataFrame,
                           individual_id: str, state: str,
                           default_group: str, group_columns: list,
                           default_metric: str, metric_columns: list,
                           yscale: str,
                           is_feat_vol_plot: bool) -> str:
    spec = {
        # This `$schema` is only fetched as part of the interactive vega
        # editor, which opens up outside of the visualization - this doesn't
        # appear to create any kind of XHR side-effect when loading the
        # visualization in an offline context.
        '$schema': 'https://vega.github.io/schema/vega/v4.0.json',
        'autosize': {'type': 'fit-x', 'contains': 'padding', 'resize': True},
        # This dimension is here for when the viz is opened in the online
        # Vega Editor.
        'width': 800,
        'title': {'text': {'signal': SIG_METRIC}},
        'background': '#FFFFFF',
        # Not registering signals on a subplot level since they aren't super
        # helpful (e.g. can't `bind` in a subplot signal).
        'signals': render_signals_ctrl(default_group, group_columns,
                                       default_metric, metric_columns,
                                       is_feat_vol_plot),
        'scales': [],
        'marks': [],
        # Add data at root of plot, it is easier to use the built in view
        # accessor to debug values
        'data': render_data_ctrl(control_chart_data, state),
    }
    control_chart = render_subplot_ctrl(yscale, state)
    if individual_id:
        control_chart['marks'].append(
            render_marks_ctrl_individual(individual_id, state))
        spec['signals'].extend(render_signals_ctrl_individual())
    spec['marks'].append(control_chart)

    if stats_chart_data is not None:
        stat_opts, sort_opts = make_stats_opts(stats_chart_data)
        spec['signals'].extend(render_signals_stats(stat_opts, sort_opts))
        spec['marks'].extend(render_subplot_stats())
        spec['data'].extend(render_data_stats(stats_chart_data))

    return json.dumps(spec, indent=2, sort_keys=True)


def make_stats_opts(data):
    sort_opts = data.columns.values.tolist()
    sort_opts.remove(FLD_STATS_ID)

    stat_opts = data.columns.values.tolist()
    stat_opts.append(FLD_STATS_AVG_CHANGE)
    stat_opts.remove(FLD_STATS_AVG_INC)
    stat_opts.remove(FLD_STATS_AVG_DEC)
    stat_opts.remove(FLD_STATS_ID)
    return (stat_opts, sort_opts)

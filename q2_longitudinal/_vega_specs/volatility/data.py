# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .const import (
    DAT_INDIVIDUAL, FLD_GROUP_BY, FLD_METRIC, DAT_GLOBAL_VALS, FLD_CTRL_CL3,
    FLD_CTRL_EXT, DAT_AGG_BY, FLD_CTRL_CI0, FLD_CTRL_CI1, FLD_CTRL_COUNT,
    DAT_SELECTED, FLD_CTRL_MEAN, FLD_CTRL_STDEV, SIG_METRIC, SIG_GROUP,
    FLD_MIN_X, FLD_MAX_X, FLD_MIN_Y, FLD_MAX_Y, FLD_CTRL_CL0, FLD_CTRL_CL1,
    FLD_CTRL_CL2)


def render_data_ctrl(control_chart_data, state):
    return [
        {'name': DAT_INDIVIDUAL,
         'values': control_chart_data.to_dict('record'),
         'transform': [
             {'type': 'formula', 'as': FLD_GROUP_BY,
              'expr': 'datum[%s]' % SIG_GROUP},
             {'type': 'formula', 'as': FLD_METRIC,
              'expr': 'datum[%s]' % SIG_METRIC}]},
        {'name': DAT_GLOBAL_VALS,
         'source': DAT_INDIVIDUAL,
         'transform': [
             {'type': 'aggregate',
              'ops': ['mean', 'min', 'max', 'stdev', 'min', 'max'],
              'fields': [FLD_METRIC, state, state, FLD_METRIC, FLD_METRIC,
                         FLD_METRIC],
              'as': [FLD_CTRL_MEAN, FLD_MIN_X, FLD_MAX_X, FLD_CTRL_STDEV,
                     FLD_MIN_Y, FLD_MAX_Y]},
             {'type': 'formula', 'as': FLD_CTRL_CL0,
              'expr': 'datum.%s - (3 * datum.%s)' % (FLD_CTRL_MEAN,
                                                     FLD_CTRL_STDEV)},
             {'type': 'formula', 'as': FLD_CTRL_CL1,
              'expr': 'datum.%s - (2 * datum.%s)' % (FLD_CTRL_MEAN,
                                                     FLD_CTRL_STDEV)},
             {'type': 'formula', 'as': FLD_CTRL_CL2,
              'expr': 'datum.%s + (2 * datum.%s)' % (FLD_CTRL_MEAN,
                                                     FLD_CTRL_STDEV)},
             {'type': 'formula', 'as': FLD_CTRL_CL3,
              'expr': 'datum.%s + (3 * datum.%s)' % (FLD_CTRL_MEAN,
                                                     FLD_CTRL_STDEV)},
             {'type': 'formula', 'as': FLD_CTRL_EXT,
              'expr': '[datum.%s, datum.%s]' % (FLD_CTRL_CL0, FLD_CTRL_CL3)}]},
        {'name': DAT_AGG_BY,
         'source': DAT_INDIVIDUAL,
         'transform': [
             {'type': 'aggregate',
              'groupby': [FLD_GROUP_BY, state],
              # TODO: parameterize these intervals
              # I don't see an easy way at the moment to define
              # your own confidence interval in vega.
              'ops': ['mean', 'ci0', 'ci1', 'count'],
              'fields': [FLD_METRIC, FLD_METRIC, FLD_METRIC, FLD_METRIC],
              'as': [FLD_CTRL_MEAN, FLD_CTRL_CI0, FLD_CTRL_CI1,
                     FLD_CTRL_COUNT]}]},
        # These are just UI state vars to keep track of what has been clicked
        # in the legend.
        {'name': DAT_SELECTED,
         'on': [
             {'trigger': 'clear', 'remove': True},
             {'trigger': '!shift', 'remove': True},
             {'trigger': '!shift && clicked', 'insert': 'clicked'},
             {'trigger': 'shift && clicked', 'toggle': 'clicked'}]}]

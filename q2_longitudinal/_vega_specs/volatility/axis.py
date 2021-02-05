# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .const import SCL_CTRL_X, SCL_CTRL_Y


def render_axes_ctrl(state):
    return [
        {'orient': 'bottom', 'scale': SCL_CTRL_X,
         'title': state},
        {'orient': 'left', 'scale': SCL_CTRL_Y,
         'title': 'Metric Column'}]


def render_axes_stats(scale, selected_stat):
    return [
        {'orient': 'top',
         'scale': scale,
         'title': {'signal': selected_stat},
         'zindex': 1}]

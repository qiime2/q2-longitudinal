# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .const import (
    LEG_CTRL_SYMBOL, LEG_CTRL_LABEL, SIG_GROUP, SCL_CTRL_COLOR, TST_OPACITY,
    STY_STROKE_2)


def render_legends_ctrl():
    return [
        {'stroke': SCL_CTRL_COLOR,
         'title': {'signal': SIG_GROUP},
         'encode': {
             'symbols': {
                 'name': LEG_CTRL_SYMBOL,
                 'interactive': True,
                 'update': {
                     'fill': {'value': 'transparent'},
                     'strokeWidth': {'value': STY_STROKE_2},
                     'opacity': [
                         {'test': TST_OPACITY, 'value': 1.0},
                         {'value': 0.15},
                     ],
                     'size': {'value': 100}}},
             'labels': {
                 'name': LEG_CTRL_LABEL,
                 'interactive': True,
                 'update': {
                     'opacity': [
                         {'test': TST_OPACITY, 'value': 1.0},
                         {'value': 0.25}]}}}}]

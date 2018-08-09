# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# 1. These constants should be used for identifiers only (don't worry about
#    values or things like mark types). The exception is when we are
#    coordinating styles across multiple components.
# 2. Names should be prefixed with a 3 character abbreviation:
#    SIG - signal
#    SCL - scale
#    DAT - data
#    FLD - field
#    TST - test
#    EXP - expression
#    STY - style
#    LEG - legend
#    VAR - any other variable that doesn't fit into the above categories

# SIGNALS
###############################################################################
SIG_METRIC = 'metric'
SIG_GROUP = 'grouper'
SIG_SHOW_GLOBAL_MEAN = 'showGlobalMean'
SIG_SHOW_GLOBAL_CTRL_LIMS = 'showGlobalControlLimits'
SIG_SHOW_ERROR_BARS = 'showErrorBars'
SIG_COLOR_SCHEME = 'colorScheme'
SIG_CTRL_MEAN_LINE_THICKNESS = 'meanLineThickness'
SIG_CTRL_MEAN_LINE_OPACITY = 'meanLineOpacity'
SIG_CTRL_MEAN_SYMBOL_SIZE = 'meanSymbolSize'
SIG_CTRL_MEAN_SYMBOL_OPACITY = 'meanSymbolOpacity'
SIG_CTRL_SPG_LINE_THICKNESS = 'spaghettiLineThickness'
SIG_CTRL_SPG_LINE_OPACITY = 'spaghettiLineOpacity'
SIG_CTRL_SPG_SYMBOL_SIZE = 'spaghettiSymbolSize'
SIG_CTRL_SPG_SYMBOL_OPACITY = 'spaghettiSymbolOpacity'
SIG_WIDTH = 'width'  # this is a special signal in vega
SIG_HEIGHT = 'height'  # this is a special signal in vega
SIG_CTRL_CHART_WIDTH = 'controlChartWidth'
SIG_CTRL_CHART_HEIGHT = 'controlChartHeight'

# STYLES
###############################################################################
STY_STROKE_2 = 2
STY_DASH_A = [8, 8]
STY_DASH_B = [6, 2]

# DATA
###############################################################################
DAT_GLOBAL_VALS = 'globalVals'
DAT_INDIVIDUAL = 'individual'
DAT_SERIES = 'series'
DAT_SPAGHETTIS = 'spaghettis'
DAT_SELECTED = 'selected'
DAT_AGG_BY = 'aggBy'

# FIELDS
###############################################################################
FLD_MIN_X = 'minX'
FLD_MAX_X = 'maxX'
FLD_MIN_Y = 'minY'
FLD_MAX_Y = 'maxY'
FLD_CTRL_MEAN = 'mean'
FLD_CTRL_MIN = 'min'
FLD_CTRL_STDEV = 'stdev'
FLD_CTRL_COUNT = 'count'
FLD_CTRL_CL0 = 'cl0'
FLD_CTRL_CL1 = 'cl1'
FLD_CTRL_CL2 = 'cl2'
FLD_CTRL_CL3 = 'cl3'
FLD_CTRL_CI0 = 'ci0'
FLD_CTRL_CI1 = 'ci1'
FLD_CTRL_EXT = 'ext'
FLD_GROUP_BY = 'groupByVal'
FLD_METRIC = 'metricVal'

# SCALES
###############################################################################
SCL_CTRL_X = 'x'
SCL_CTRL_Y = 'y'
SCL_CTRL_COLOR = 'color'

# TESTS
###############################################################################
TST_OPACITY = "!length(data('{0}')) || indata('{0}', 'value', datum.value)".format(DAT_SELECTED)  # noqa: E501
TST_GROUP = "!length(data('{0}')) || indata('{0}', 'value', datum.{1})".format(DAT_SELECTED, FLD_GROUP_BY)  # noqa: E501

# LEGEND
###############################################################################
LEG_CTRL_SYMBOL = 'legendSymbol'
LEG_CTRL_LABEL = 'legendLabel'

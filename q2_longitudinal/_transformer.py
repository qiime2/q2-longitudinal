# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import qiime2

from .plugin_setup import plugin
from ._format import FirstDifferencesFormat


# modified from q2_types/sample_data/_transformer.py
def _read_first_differences(fh):
    return qiime2.Metadata.load(fh)


@plugin.register_transformer
def _1(data: pd.Series) -> FirstDifferencesFormat:
    ff = FirstDifferencesFormat()
    with ff.open() as fh:
        data.to_csv(fh, sep='\t', header=True)
    return ff


@plugin.register_transformer
def _2(ff: FirstDifferencesFormat) -> pd.Series:
    return _read_first_differences(str(ff)).to_dataframe().iloc[:, 0]


@plugin.register_transformer
def _3(ff: FirstDifferencesFormat) -> qiime2.Metadata:
    return _read_first_differences(str(ff))

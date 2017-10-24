# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
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
    # Using `dtype=object` and `set_index` to avoid type casting/inference
    # of any columns or the index.
    df = pd.read_csv(fh, sep='\t', header=0, dtype=object)
    df.set_index(df.columns[0], drop=True, append=False, inplace=True)
    # The only column in the dataframe should be numeric, so try to cast the
    # entire dataframe.
    df = df.astype(float)
    return df


@plugin.register_transformer
def _1(data: pd.Series) -> FirstDifferencesFormat:
    ff = FirstDifferencesFormat()
    with ff.open() as fh:
        data.to_csv(fh, sep='\t', header=True)
    return ff


@plugin.register_transformer
def _2(ff: FirstDifferencesFormat) -> pd.Series:
    with ff.open() as fh:
        df = _read_first_differences(fh)
        return df.iloc[:, 0]


@plugin.register_transformer
def _3(ff: FirstDifferencesFormat) -> qiime2.Metadata:
    with ff.open() as fh:
        return qiime2.Metadata(_read_first_differences(fh))

# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
import pandas.util.testing as pdt
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_longitudinal._format import FirstDifferencesFormat


class TestFirstDifferencesTransformers(TestPluginBase):
    package = "q2_longitudinal.tests"

    def test_pd_series_to_first_differences_format(self):
        transformer = self.get_transformer(pd.Series, FirstDifferencesFormat)

        index = pd.Index(['a', 'b', 'c', 'd'], name='#SampleID', dtype=object)
        series = pd.Series([1, 2, 3, 4], name='Difference', index=index)

        result = transformer(series)

        with open(str(result), 'r') as fh:
            obs = fh.read()

        self.assertEqual(
            obs, '#SampleID\tDifference\na\t1\nb\t2\nc\t3\nd\t4\n')

    def test_first_differences_format_to_pd_series(self):
        _, obs = self.transform_format(
            FirstDifferencesFormat, pd.Series, 'first-differences-simple.tsv')

        exp = pd.Series([1.0, 2.0, 3.0, 4.0], name='Difference',
                        index=pd.Index(['a', 'b', 'c', 'd'], name='#SampleID'))

        pdt.assert_series_equal(obs, exp)

    def test_first_differences_format_to_metadata(self):
        _, obs = self.transform_format(
            FirstDifferencesFormat, qiime2.Metadata,
            'first-differences-simple.tsv')

        exp_index = pd.Index(['a', 'b', 'c', 'd'], name='#SampleID',
                             dtype=str)
        exp_df = pd.DataFrame({'Difference': [1.0, 2.0, 3.0, 4.0]},
                              index=exp_index)
        exp = qiime2.Metadata(exp_df)

        self.assertEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()

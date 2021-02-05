# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase

from q2_longitudinal._format import FirstDifferencesFormat


class TestFirstDifferencesFormat(TestPluginBase):
    package = "q2_longitudinal.tests"

    def test_valid_simple(self):
        filepath = self.get_data_path('first-differences-simple.tsv')
        format = FirstDifferencesFormat(filepath, mode='r')

        format.validate('min')
        format.validate('max')

    def test_valid_real_data(self):
        filepath = self.get_data_path('ecam-taxa-first-differences.tsv')
        format = FirstDifferencesFormat(filepath, mode='r')

        format.validate('min')
        format.validate('max')

    def test_valid_alternate_header(self):
        # `Distance` instead of `Difference` for second column name
        filepath = self.get_data_path('ecam-first-distances.tsv')
        format = FirstDifferencesFormat(filepath, mode='r')

        format.validate('min')
        format.validate('max')

    def test_invalid_header(self):
        filepath = self.get_data_path('first-differences-invalid-header.tsv')
        format = FirstDifferencesFormat(filepath, mode='r')

        for level in 'min', 'max':
            with self.assertRaisesRegex(ValidationError,
                                        'Header.*\n\n.*#SampleID.*Peanut'):
                format.validate(level)

    def test_invalid_not_two_fields(self):
        filepath = self.get_data_path('first-differences-not-two-fields.tsv')
        format = FirstDifferencesFormat(filepath, mode='r')

        for level in 'min', 'max':
            with self.assertRaisesRegex(ValidationError,
                                        'two fields.*3.*line 2.*\n\n.*a'):
                format.validate(level)

    def test_invalid_non_numeric_column(self):
        filepath = self.get_data_path('first-differences-non-numeric.tsv')
        format = FirstDifferencesFormat(filepath, mode='r')

        for level in 'min', 'max':
            with self.assertRaisesRegex(ValidationError,
                                        "non-numeric.*'e'.*line 4"):
                format.validate(level)

    def test_invalid_header_only(self):
        filepath = self.get_data_path('first-differences-header-only.tsv')
        format = FirstDifferencesFormat(filepath, mode='r')

        for level in 'min', 'max':
            with self.assertRaisesRegex(ValidationError,
                                        "at least one data record"):
                format.validate(level)

    def test_invalid_empty_file(self):
        filepath = self.get_data_path('empty')
        format = FirstDifferencesFormat(filepath, mode='r')

        for level in 'min', 'max':
            with self.assertRaisesRegex(ValidationError,
                                        'Header line must be TSV'):
                format.validate(level)

    def test_valid_min_invalid_max(self):
        filepath = self.get_data_path(
                'first-differences-valid-min-invalid-max.tsv')
        format = FirstDifferencesFormat(filepath, mode='r')

        # Error not raised with min validation.
        format.validate('min')

        # Error raised with max validation.
        with self.assertRaisesRegex(ValidationError,
                                    "non-numeric.*'milo'.*line 7"):
            format.validate('max')


if __name__ == '__main__':
    unittest.main()

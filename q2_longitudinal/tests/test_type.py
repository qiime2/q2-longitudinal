# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.plugin.testing import TestPluginBase
from q2_types.sample_data import SampleData

from q2_longitudinal._type import FirstDifferences
from q2_longitudinal._format import FirstDifferencesDirectoryFormat


class TestTypes(TestPluginBase):
    package = "q2_longitudinal.tests"

    def test_first_differences_semantic_type_registration(self):
        self.assertRegisteredSemanticType(FirstDifferences)

    def test_first_differences_semantic_type_to_format_registration(self):
        self.assertSemanticTypeRegisteredToFormat(
            SampleData[FirstDifferences],
            FirstDifferencesDirectoryFormat
        )


if __name__ == '__main__':
    unittest.main()

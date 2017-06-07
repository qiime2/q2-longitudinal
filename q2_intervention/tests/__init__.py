# ----------------------------------------------------------------------------
# Copyright (c) 2017--, q2-intervention development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import pkg_resources

from qiime2.plugin.testing import TestPluginBase


class InterventionTestPluginBase(TestPluginBase):
    package = 'q2_intervention'

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory(
            prefix='q2-intervention-test-temp-')

    def tearDown(self):
        self.temp_dir.cleanup()

    def get_data_path(self, filename):
        return pkg_resources.resource_filename(self.package,
                                               'test_data/%s' % filename)

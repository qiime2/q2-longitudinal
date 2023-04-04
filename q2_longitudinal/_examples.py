# ----------------------------------------------------------------------------
# Copyright (c) 2017-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources

import qiime2


def _get_data_from_tests(path):
    return pkg_resources.resource_filename('q2_longitudinal.tests',
                                           os.path.join('data', path))


def sample_md_blanks_factory():
    return qiime2.Metadata.load(
        _get_data_from_tests('ecam-sample-md-blanks.tsv')
    )


def volatility_missing_md_in_state_column(use):
    md = use.init_metadata('md', sample_md_blanks_factory)

    volatility_viz, = use.action(
        use.UsageAction('longitudinal', 'volatility'),
        use.UsageInputs(
            metadata=md,
            state_column='month'
        ),
        use.UsageOutputNames(
            visualization='volatility_plot'
        )
    )
    # This is testing the expected behavior for _convert_nan_to_none()
    # Assuring us that NaNs are successfully converted to nulls in the viz
    volatility_viz.assert_has_line_matching('index.html',
                                            r'^\s*"month": null,$')

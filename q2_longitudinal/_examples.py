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


def ecam_missing_numerical_md_factory():
    return qiime2.Metadata.load(
        _get_data_from_tests('ecam_sample_md_blanks.tsv'))


def longitudinal_volatility_missing_numerical_md(use):
    metadata = use.init_metadata('metadata', ecam_missing_numerical_md_factory)

    volatility, = use.action(
        use.UsageAction('longitudinal', 'volatility'),
        use.UsageInputs(
            metadata=metadata,
            state_column='month'
        ),
        use.UsageOutputNames(
            visualization='volatility_plot',

        )
    )

# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup
import versioneer


setup(
    name='q2-longitudinal',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    packages=find_packages(),
    author="Nicholas Bokulich",
    author_email="nbokulich@gmail.com",
    description=(
        "QIIME2 plugin for longitudinal studies and paired comparisons."),
    url="https://github.com/qiime2/q2-longitudinal",
    entry_points={
        'qiime2.plugins':
        ['q2-longitudinal=q2_longitudinal.plugin_setup:plugin']
    },
    package_data={
        'q2_longitudinal.tests': ['data/*'],
        'q2_longitudinal': [
            'citations.bib',
            'assets/index.html',
            'assets/volatility/index.html',
            'assets/volatility/js/*',
            'assets/volatility/css/*',
            'assets/volatility/licenses/*',
        ],
    },
    zip_safe=False,
)

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, q2-longitudinal development team.
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
    install_requires=['biom-format', 'pandas', 'scipy',
                      'scikit-bio', 'seaborn', 'statsmodels'],
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
        'q2_longitudinal.tests': ['tests/data/*'],
        'q2_longitudinal': ['assets/index.html'],
    },
    zip_safe=False,
)

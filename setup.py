# ----------------------------------------------------------------------------
# Copyright (c) 2017--, q2-intervention development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup
import versioneer


setup(
    name='q2-intervention',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    packages=find_packages(),
    install_requires=['biom-format', 'pandas', 'scipy',
                      'scikit-bio', 'seaborn', 'statsmodels'],
    author="Nicholas Bokulich",
    author_email="nbokulich@gmail.com",
    description=(
        "QIIME2 plugin for intervention studies and paired comparisons."),
    url="https://github.com/nbokulich/q2-intervention",
    entry_points={
        'qiime2.plugins':
        ['q2-intervention=q2_intervention.plugin_setup:plugin']
    },
    package_data={
        'q2_intervention.tests': ['test_data/*'],
        'q2_intervention': ['assets/index.html'],
    },
    zip_safe=False,
)

# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model
from qiime2.plugin import ValidationError


class FirstDifferencesFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        with self.open() as fh:
            line = fh.readline().rstrip()
            if line not in ['#SampleID\tDifference', '#SampleID\tDistance']:
                raise ValidationError(
                    "Header line must be TSV with column names '#SampleID' "
                    "and either 'Difference' or 'Distance'. Found the "
                    "following header:\n\n{0!r}".format(line))

            has_data = False
            for line_number, line in enumerate(fh, start=2):
                cells = line.strip().split('\t')
                if len(cells) != 2:
                    # TODO indicate tab separated
                    raise ValidationError(
                        "Expected data record to be TSV with two fields, "
                        "detected {0} fields at line {1}:\n\n{2!r}"
                        .format(len(cells), line_number, cells))
                try:
                    float(cells[1])
                except ValueError:
                    raise ValidationError(
                        "Second column must contain only numeric values. "
                        "A non-numeric value ({0!r}) was detected at line "
                        "{1}.".format(cells[1], line_number))

                has_data = True
                if n_records is not None and (line_number - 1) >= n_records:
                    break

            if not has_data:
                raise ValidationError(
                    "There must be at least one data record present in the "
                    "file in addition to the header line.")

    def _validate_(self, level):
        record_count_map = {'min': 5, 'max': None}
        self._validate(record_count_map[level])


FirstDifferencesDirectoryFormat = model.SingleFileDirectoryFormat(
    'FirstDifferencesDirectoryFormat', 'FirstDifferences.tsv',
    FirstDifferencesFormat)

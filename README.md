# q2-intervention
QIIME2 plugin for paired sample comparisons

## Examples

### Paired difference testing

Paired differences in metadata
```
cd ~/Desktop/projects/q2-intervention/q2_intervention/test_data

qiime intervention paired-differences --i-table ecam-table-taxa.qza --m-metadata-file ecam_map_maturity.txt --p-metric 'k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__' --p-group-category delivery --p-state-category month --p-state-pre 1 --p-state-post 12 --p-individual-id-category studyid --o-visualization ecam-delivery --p-no-drop-duplicates
qiime intervention paired-differences --i-table ecam-table-taxa.qza --m-metadata-file ecam_map_maturity.txt --p-metric 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__;s__' --p-group-category diet_3 --p-state-category month --p-state-pre 0 --p-state-post 12 --p-individual-id-category studyid --o-visualization ecam-diet --p-no-drop-duplicates --verbose
```

Paired differences in feature table
```
```
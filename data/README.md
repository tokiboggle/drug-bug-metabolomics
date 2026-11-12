# Data files

All three files live in `raw/` and come from the Dataset 2 handout (`group2_drugs_*.csv`). Put anything you derive later in `processed/`.

The sample design is summarised in the main README. Facts below were read from the files themselves, and from `docs/data_description.pdf`.

## `group2_drugs_metadata.csv`

69 rows (one per injection) and 8 columns.

| Column | What it is |
| --- | --- |
| `filename` | The `.mzML` file name. This is what joins the metadata to the `* Peak area` columns in the feature table. |
| `ATRRIBUTE_SampleType` | Time point (`T0`–`T8`) or `Blank` / `QC` / `Pooled_QC`. The name is misspelled in the file: three R's (`ATRRIBUTE`). |
| `ATTRIBUTE_Replicate` | `1`, `2` or `3` for the experimental injections, or the control label (`Blank`, `QC`, `Pooled_QC`). |
| `ATTRIBUTE_Bacteria` | `Yes` (community present), `No` (medium only), or the control label. |
| `mzmine_sample_type` | `sample`, `blank` or `qc`. |
| `ATTRIBUTE_time_run` | A time stamp for the run. A few QC rows are empty here. |
| `ATTRIBUTE_Drug` | `Omeprazole`, `DMSO`, `Medium`, `Blank`, `QC` or `Pooled_QC`. |
| `run_date` | When the injection was run on the instrument. |

## `group2_drugs_iimn_gnps_quant.csv`

The feature intensity table (MZmine / IIMN / GNPS export). **4,525 rows** (one per feature) and 83 columns. The data-description PDF says the MZmine table has 6,993 features; this file is smaller, so do not assume the two numbers refer to the same table.

| Column | What it is |
| --- | --- |
| `row ID` | Feature identifier. This is what `id` in the annotation file refers to. |
| `row m/z` | Mass-to-charge ratio of the feature. |
| `row retention time` | When the feature came off the chromatography column, in minutes. |
| `row ion mobility`, `row ion mobility unit`, `row CCS` | Ion-mobility fields. Often empty in this export. |
| `correlation group ID`, `annotation network number`, `best ion`, `auto MS2 verify`, `identified by n=`, `partners`, `neutral M mass` | Fields from the IIMN / GNPS networking step. |
| `… Peak area` | One column per injection (69 in total). The column name is the metadata `filename` plus ` Peak area`. |
| (trailing empty column) | An unnamed column left by the export. Ignore it. |

## `group2_drugs_annotations.csv`

2,713 rows and 20 columns, covering 335 distinct `id` values. Every `id` is present as a `row ID` in the feature table. A feature can appear more than once if it has several suggested names.

| Column | What it is |
| --- | --- |
| `id` | Feature identifier; matches `row ID` in the feature table. |
| `compound_name` | Suggested name. |
| `adduct` | Ion form, for example `[M+H]+`. |
| `score` | Match score from the annotation method. |
| `precursor_mz` | Precursor mass-to-charge ratio used in the match. |
| `ion_mobility`, `ccs`, `rt` | Ion mobility, collision cross-section, and retention time. Many rows leave these empty. |
| `mol_formula`, `smiles`, `isomeric_smiles`, `inchi`, `inchi_key` | Chemical identifiers for the suggested compound, when present. |
| `method` | How the annotation was made. In this file: mostly `spectral_library_annotation`, a small number `matchedlipid`. |
| `usi`, `entry_id`, `iupac_name`, `cas`, `internal_id`, `query_spectrum_usi` | Library entry and spectrum identifiers. Several of these are empty or mostly empty. |

## Things to watch for when you join the tables

These are facts about the files, not analysis results.

- The metadata column is spelled `ATRRIBUTE_SampleType` (three R's).
- Pooled QC file names use `Polled_Qcs` (not "Pooled").
- The T7 injections are named with a lowercase `t7_…mzML` in both the metadata and the feature table. The `ATRRIBUTE_SampleType` value for those rows is still `T7`.
- Medium controls appear only at T0 (`T0_Medium_noBac_rep1/2/3.mzML`).
- Join peak-area columns to metadata by stripping the ` Peak area` suffix and matching `filename`. All 69 names match.

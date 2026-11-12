# Omeprazole and a bacterial community: an untargeted metabolomics study

We look at how a bacterial culture or synthetic community (a set of bacterial strains grown together on purpose) changes its metabolism over time when treated with the drug omeprazole, compared with a DMSO control. We'll present what we find.

Data source: Pakkir Shah et al., *A Functional Metabolomics Framework to Track Microbiome Drug Metabolism*, bioRxiv. The data description is in `docs/data_description.pdf`.

## Question

> How do omeprazole-treated cultures change over time compared with DMSO controls?

## The data

Three tables are in `data/raw/`. They are small enough to keep in git. Column details are in `data/README.md`.

| File | What it is |
| --- | --- |
| `group2_drugs_metadata.csv` | One row per injection: treatment, time point, replicate, and whether it is a real sample, a blank, or a QC (quality-control) injection. 69 rows. |
| `group2_drugs_iimn_gnps_quant.csv` | The feature table: one row per metabolite feature, with a peak area (how much of that feature was seen) in each sample. **4,525 features** in this file. The data description says the MZmine table has 6,993 features. The counts don't match, so check this before you start. |
| `group2_drugs_annotations.csv` | Suggested compound names for some features, mostly from spectral-library matching. 2,713 annotation rows covering 335 of the 4,525 features. One feature can have more than one suggested name. |

Sample design, from `docs/data_description.pdf` and checked against the metadata file:

| Kind | Count | Detail |
| --- | --- | --- |
| Bacterial samples | 54 | 2 treatments (omeprazole and DMSO) × 9 time points (T0–T8) × 3 replicates |
| Medium controls | 3 | Medium with no bacteria, at T0 only |
| Blanks | 4 | Background / solvent injections |
| QC injections | 8 | 4 pooled QCs + 4 QC-mix |

A few file quirks (a typo in a column name, lowercase `t7` in some file names, and `Polled_Qcs` instead of "Pooled") are listed in `data/README.md`.

## How to start

```bash
git clone <this-repo-url>
cd <repo-folder>
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
jupyter notebook notebooks/01_load_and_inspect.ipynb
```

The first notebook loads the three files and prints their size and a count of each sample type. The later steps are empty headings to fill in.

## Folder map

| Folder | What's in it |
| --- | --- |
| `data/raw/` | The three CSV files above |
| `data/processed/` | Cleaned tables made later (large outputs are kept out of git) |
| `docs/` | Project PDFs, including the data description |
| `notebooks/` | Analysis notebooks, starting with `01_load_and_inspect.ipynb` |
| `results/figures/` and `results/tables/` | Plots and tables for the presentation |
| `notes/` | Decision log and AI-use log |
| `slides/` | Presentation slides |

## Analysis plan

Tick items off as you go. This is a plan, not a finished analysis, and you don't have to do every step.

- [ ] **1. Know the design.** Groups, replicates, and anything that could confound the results (for example when each sample was run).
- [ ] **2. Clean and transform.** Filter features, including blank and QC filtering; log-transform; normalise; impute missing values if needed.
- [ ] **3. Test.** Find features whose abundance changes with treatment, with time, or with the treatment × time interaction (a change over time that differs between omeprazole and DMSO). Use a multiple-testing correction such as FDR (a way to limit false hits when testing thousands of features).
- [ ] **4. Show.** PCA (a plot that summarises how samples group), volcano plots, and a heatmap.
- [ ] **5. Interpret the list.** Use the annotation table, then enrichment (do the changing features cluster in a pathway or compound class?). Rank unknown features too, not only annotated ones.
- [ ] **6. Connect.** Interaction networks, hubs and modules.
- [ ] **7. Zoom in on hits.** Take one or two interesting features and look them up more closely (identity, structure, what is already known).

## Presenting the findings

- [ ] A clear line from the question, through the analysis, to the conclusions.
- [ ] The limits of the method: what LC–MS/MS metabolomics can and cannot see, and how that limits the conclusions.
- [ ] What AI tools we used, what for, and how we checked them (see `notes/ai-use.md`).
- [ ] Last slide: each member's name, student number, and who did what.

## Who does what

| Member | Student number | Part of the project |
| --- | --- | --- |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |

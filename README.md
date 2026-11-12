# 27200 group project — ME2: omeprazole and a bacterial community

**Keep this repository private.** It holds course handouts.

**Group __** (ME2 was allocated to groups 9, 12 and 22; fill in your number.)

## What this project is

This is untargeted metabolomics of a bacterial culture or synthetic community treated with omeprazole, compared with a DMSO control, and followed over time. The starting question, from the group-projects introduction, is:

> How do omeprazole-treated cultures change over time compared with DMSO controls?

There is no written report and nothing to hand in before the exam. The only deliverable is the oral exam on Thursday 10 December 2026. The project is 70% of the course grade (guidelines PDF).

## The data

Three tables are in `data/raw/`, from the Dataset 2 handout. They are small enough to keep in git. Column details are in `data/README.md`.

| File | What it is |
| --- | --- |
| `group2_drugs_metadata.csv` | One row per injection: treatment, time point, replicate, and whether it is a real sample, a blank, or a QC (quality-control) injection. 69 rows. |
| `group2_drugs_iimn_gnps_quant.csv` | The feature table: one row per metabolite feature, with a peak area (how much of that feature was seen) in each sample. **4,525 features** in this file. The data-description PDF says the MZmine table has 6,993 features — the counts do not match, so check this before you start. |
| `group2_drugs_annotations.csv` | Suggested compound names for some features, mostly from spectral-library matching. 2,713 annotation rows covering 335 of the 4,525 features. One feature can have more than one suggested name. |

Sample design, from `data_description.pdf` and checked against the metadata file:

| Kind | Count | Detail |
| --- | --- | --- |
| Bacterial samples | 54 | 2 treatments (omeprazole and DMSO) × 9 time points (T0–T8) × 3 replicates |
| Medium controls | 3 | Medium with no bacteria, at T0 only |
| Blanks | 4 | Background / solvent injections |
| QC injections | 8 | 4 pooled QCs + 4 QC-mix |

Source of the dataset, from the introduction and the data description: Pakkir Shah et al., *A Functional Metabolomics Framework to Track Microbiome Drug Metabolism*, bioRxiv.

A few file quirks (a typo in a column name, lowercase `t7` in some file names, and `Polled_Qcs` instead of "Pooled") are listed in `data/README.md`.

## How to start

```bash
git clone <this-repo-url>
cd 27200-ME2-project
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
jupyter notebook notebooks/01_load_and_inspect.ipynb
```

If you received this folder as a zip and it already contains a `.git` directory, skip the clone. Create an empty **private** repository on GitHub, then from this folder:

```bash
git remote add origin <your-private-repo-url>
git push -u origin main
```

The first notebook only loads the three files and prints their size and a count of each sample type. Later steps are empty headings for the group to fill in.

## Folder map

| Folder | What's in it |
| --- | --- |
| `data/raw/` | The three CSV files above |
| `data/processed/` | Cleaned tables you make later (large outputs are kept out of git) |
| `docs/` | Course PDFs: guidelines, allocation, introduction, and data description |
| `notebooks/` | Analysis notebooks, starting with `01_load_and_inspect.ipynb` |
| `results/figures/` and `results/tables/` | Plots and tables for the talk |
| `notes/` | Decision log and AI-use log |
| `slides/` | Presentation slides, when you make them |

## Analysis plan

This follows the usual route in `docs/27200_project_guidelines.pdf`. Tick items off as you go. This is a checklist, not a finished analysis, and it does not tell you what you will find.

- [ ] **1. Know your design.** Write down the groups, replicates, and anything that could confound the results (for example when each sample was run). Practised in Week 4.
- [ ] **2. Clean and transform.** Filter features, including blank and QC filtering; log-transform; normalise; impute missing values if needed. Weeks 4, 6 and 7.
- [ ] **3. Test.** Find features whose abundance changes with treatment, with time, or with the treatment × time interaction (a change over time that differs between omeprazole and DMSO). Use a multiple-testing correction such as FDR (a way to limit false hits when you test thousands of features). Weeks 4, 6 and 7.
- [ ] **4. Show.** PCA (a plot that summarises how samples group), volcano plots, and a heatmap. Weeks 6, 7 and 9.
- [ ] **5. Interpret the list.** Use the annotation table, then enrichment (do the changing features cluster in a pathway or compound class?). Weeks 6 and 7.
- [ ] **6. Connect.** Interaction networks, hubs and modules. Week 3.
- [ ] **7. Zoom in on hits.** Take one or two interesting features and look them up more closely (identity, structure, what is already known). Week 2.

The data-description PDF also suggests starting with a treatment × time model, and ranking both annotated and unknown features.

You need **at least three** tools or representations from the course. The list above is a menu, not a requirement to do every step.

## Exam requirements

From the project guidelines:

- [ ] At least three course tools or representations (for example statistics, enrichment, networks, or sequence and structure).
- [ ] Discuss measurement constraints: what LC–MS/MS metabolomics can and cannot see, and how that limits what you can conclude.
- [ ] Reflect on AI use: what you used, what for, and how you checked it. Keep a running log in `notes/ai-use.md`.
- [ ] Last slide of the talk: each member's name, student number, and who did what.

The exam is taken as a group, but **grades are individual**. Each person is graded on the part they presented and on their answers in the discussion.

What the examiners look for (short labels from the guidelines): *coherence* (a clear line from the question through the analysis to the conclusion), *systems reasoning* (going beyond a list of hits, to pathways, networks and mechanisms), *model awareness* (knowing the assumptions and limits of your methods and of what the technology can measure), and *interpretation quality* (conclusions that follow from the data, with caveats named).

Format: 15–20 minute presentation (every member presents a part) plus 15–20 minute discussion. Two instructors from metabolomics examine you. One other group sits in and may ask questions; those questions do not affect your grade. You do not hand in slides or a report beforehand. The time slot will be published on DTU Learn.

## Key dates

| When | What |
| --- | --- |
| Thu 12 Nov, Thu 19 Nov, Thu 26 Nov, Thu 3 Dec 2026 | Sessions 10–13: class time to work on the project |
| Thu 26 Nov 2026 (Week 12) | Peer feedback between groups |
| Thu 10 Dec 2026 | Oral exam (group exam). Exact time on DTU Learn |

## Who does what

Fill this in, and copy it onto the last slide of the talk.

| Member | Student number | Part of the project |
| --- | --- | --- |
|  |  |  |
|  |  |  |
|  |  |  |
|  |  |  |

# Omeprazole and a bacterial community: an untargeted metabolomics study

We investigate how the Com20 synthetic human gut bacterial community (20 commensal strains grown anaerobically in mGAM medium) changes its metabolism over an 8-hour time course (`T0–T8`, hourly) when treated with the proton-pump inhibitor **omeprazole**, compared with a **DMSO** vehicle control.

Data source: Pakkir Shah et al. (2026), *A Functional Metabolomics Framework to Track Microbiome Drug Metabolism*, bioRxiv (`10.64898/2026.01.30.702925v1` / PMC12889426). The data description is in `docs/data_description.pdf`.

## Question

> How do omeprazole-treated cultures change over time compared with DMSO controls?

## Key Findings at a Glance

1. **QC & Confounder Audit (Weeks 4 & 8):**
   - **4,525 vs 6,993 feature count:** `row ID` spans `16–11,816` because this table is a subset of the 12-drug Com20 study; consequently, **Pooled QCs** contain features from other drugs (e.g., `erythromycin`, `ID 4135`) and cannot be used as a strict presence filter (`629/1,060` omeprazole-specific features are absent from Pooled QCs).
   - **4 failed blank-like injections removed** (`T4_DMSO_Bac_rep3`, `T4_Omepr_Bac_rep3`, `T5_DMSO_Bac_rep3`, `T5_Omepr_Bac_rep2`): TIC `4.4–4.7 × 10⁷` (below all 4 Blanks), `>81%` zeros, `r = 0.83–0.89` with Blanks vs `r = 0.21–0.40` with true replicates (`n = 50` valid bacterial samples retained).
   - **Run-order confounding & autosampler instability:** Sampling timepoint correlates with 5.3-day LC-MS acquisition order (`r = −0.47` in DMSO, `−0.43` in Omeprazole). Acid-labile parent **omeprazole (`ID 2828`)** decays in the autosampler queue (`r = −0.62, p = 9.4 × 10⁻⁴`, dropping to `0` in T1/T2 run on Days 5–6), whereas microbial biotransformation products track biological time cleanly (`Ufiprazole ID 2962`, `r_time = +0.976`).
2. **Normalization & Measurement Constraints (Weeks 4 & 8):**
   - **TIC vs PQN:** Ufiprazole (`ID 2962`) + its 9 co-eluting IIMN fragments/isotopes (`cg=220`) surge 9-fold from T0→T8 (~47% of total TIC), doubling Omeprazole TIC (`1.05× → 1.99×`) and causing blind TIC normalization to fabricate **`1,780/2,146` (83%)** false down-regulated endogenous features. **PQN** eliminates this compositional artifact.
   - **ESI co-elution ion suppression (`RT 1.80–2.25 min`):** Even after PQN, **`114/172` (66.3%)** of significant negative `Trt × Time` shared interactions fall inside Ufiprazole's elution/tailing window (`9.0%` of features; `OR = 46.5`, Fisher `p = 2.1 × 10⁻⁹¹`), flagging features like `N-Phenethylacetamide (ID 2972)` as ion-suppression artifacts.
3. **Biological Conclusion (Weeks 2, 3, 4, 6–7, 9, 13):**
   - **Global PERMANOVA (`3,414` filtered features):** `Treatment` explains **50.8%**, `Time` **28.4%**, `Trt × Time` **11.7%** (`all p = 0.001`).
   - **Clean Endogenous PERMANOVA (`1,915` shared features outside `RT 1.80–2.25 min`):** `Time` dominates at **59.1%** (`Treatment` `16.1%`, `Trt × Time` `7.1%`) — Com20 preserves its core 0–8 h metabolic growth program (`Bile Acids` 5.5× and `N-Acyl Amino Acids` 2.9× enriched in produced metabolites, `q = 0.012`; `Cyclic Dipeptides` 3.7× enriched in stable mGAM background, `q = 1.1 × 10⁻¹⁸`), with selective modest perturbations under omeprazole (e.g., ~2× attenuated `Indole-3-propionic acid ID 4553`, delayed `PE 13:0 ID 7477`).
   - **Dominant drug transformation cascade:** Anaerobic bacterial sulfoxide reduction converts omeprazole to **Ufiprazole (`ID 2962`, `m/z 330.1270`, +9×)**, branching into **O-demethylation (`ID 2465`)**, **C–S cleavage (`5-MeO-2-mercaptobenzimidazole ID 2480`, `³⁴S` confirmed)**, **regioselective methylation (`ID 3605` benzimidazole-Me vs `ID 2500` pyridine-Me, localized via co-eluting fragments)**, and sulfur-verified dark-metabolome conjugates (notably **`ID 2368`, `m/z 274.0565`, `C₁₁H₁₅NO₃S₂`, `−0.46 ppm`, 2S confirmed by `8.2% ³⁴S` + intact pyridine `151.0992`, +95×**).

## The data

Three raw tables live in `data/raw/` (column details in `data/README.md`):

| File | What it is |
| --- | --- |
| `group2_drugs_metadata.csv` | One row per injection: treatment, time point, replicate, run date, and sample/QC/blank labels (`69` rows). |
| `group2_drugs_iimn_gnps_quant.csv` | MZmine/IIMN/GNPS feature table (`4,525` features × `69` injections; subset of the 12-drug study table). |
| `group2_drugs_annotations.csv` | Spectral-library & lipid matches (`2,713` rows covering `335` distinct features; `259` retained after filtering). |

## How to reproduce the analysis

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# 1. Run end-to-end analysis (writes data/processed/, results/tables/, results/networks/)
python src/pipeline.py

# 2. Generate all 9 publication-quality figures (writes results/figures/)
python src/make_figures.py

# 3. Explore executed walkthrough notebooks
jupyter notebook notebooks/
```

## Folder map

| Folder | What's in it |
| --- | --- |
| `data/raw/` | Original CSV files (`group2_drugs_*.csv`) |
| `data/processed/` | Cleaned PQN & `log₂(PQN)` matrices, sample QC table, curated RDKit annotations, and `analysis_summary.json` |
| `docs/` | Course project PDFs (`27200_project_guidelines.pdf`, `data_description.pdf`, etc.) |
| `src/` | Reproducible analysis pipeline (`pipeline.py`) and figure generator (`make_figures.py`) |
| `notebooks/` | Starter inspection notebook (`01`) + 3 fully executed analysis notebooks (`02_qc_normalization_and_statistics.ipynb`, `03_enrichment_and_endogenous_metabolism.ipynb`, `04_networks_and_hit_zoomin.ipynb`) |
| `results/figures/` | 9 high-DPI (300 DPI) slide-ready figures (`fig01`–`fig09`) |
| `results/tables/` | 10 structured CSV result tables (`table01`–`table09` + ML CV metrics) |
| `results/networks/` | Cytoscape-ready `.graphml` files (`omeprazole_biotransformation_network.graphml`, `endogenous_coabundance_network.graphml`) |
| `notes/` | Full analytical `decisions.md` and `ai-use.md` logs |
| `slides/` | Presentation slides |

## Analysis plan

- [x] **1. Know your design.** Audited groups, replicates, 4 failed injections, Pooled QC/QC-Mix/Blank/Medium roles, and 5.3-day LC-MS run-order confounding (`fig01`, `table01`).
- [x] **2. Clean and transform.** Filtered to `3,414` features (`≥2` reps in `≥1` group, `>3×` blank); proved PQN prevents TIC's 83% false-down artifact; left-censored `0.2×min` imputation + `log₂`; diagnosed `RT 1.80–2.25 min` ESI ion suppression (`fig02`).
- [x] **3. Test.** Run-order-adjusted linear models, two-way categorical ANOVA (`+774` non-linear transient hits), within-treatment temporal tests, BH-FDR (`q < 0.05`), and 999-permutation PERMANOVA (`fig03`, `fig04`, `table02`, `table09`).
- [x] **4. Show.** Global & clean-endogenous PCA, volcano plots, IIMN-collapsed temporal heatmap, 4 k-means kinetic modules (`n = 606`), and 5-fold CV Random Forest / Logistic Regression prioritization of known & unknown hits (`fig03`–`fig06`, `table08`).
- [x] **5. Interpret the list.** Audited all `259` annotations (flagging carryover drugs & false-positive `Sertindole`), ran RDKit chemical-class Fisher's exact ORA (`q < 0.05`), and mapped endogenous Com20 trajectories (`fig07`, `table04`, `table05`).
- [x] **6. Connect.** Built IIMN-collapsed exact-mass biotransformation network (`57` nodes, `82` edges) and endogenous Com20 co-abundance network (`91` nodes, `834` edges, `8` modules, hubs/betweenness) exported to Cytoscape `.graphml` (`fig08`, `table06`, `table07`).
- [x] **7. Zoom in on hits.** RDKit 2D structures, Morgan Tanimoto similarity, `³⁴S` M+2 isotope verification, regioisomeric methylation localization (`ID 3605` vs `ID 2500`), and dark-metabolome candidate `ID 2368` (`C₁₁H₁₅NO₃S₂`) (`fig09`, `table03`).

## Presenting the findings

- [x] A clear line from the question, through the analysis, to the conclusions (`notebooks/02–04`, `fig01–fig09`).
- [x] The limits of the method: positive-mode ESI/EtOAc bias, ionization efficiency vs concentration, ESI co-elution suppression, autosampler instability, lack of abiotic omeprazole controls, dark metabolome (`92.4%` unannotated), and `n = 2–3` power (`notebooks/04`, `notes/decisions.md`).
- [x] What AI tools we used, what for, and how we checked them (`notes/ai-use.md`).
- [ ] Last slide: each member's name, student number, and who did what.

## Who does what

| Member | Student number | Part of the project |
| --- | --- | --- |
| Tóki |  | Design, QC, normalization & statistical modeling (`notebooks/02`, `fig01`–`fig04`) |
| Magda |  | PCA, kinetic clustering, ML prioritization & statistics (`notebooks/02`, `fig03`–`fig06`) |
| Caroline |  | Annotation audit, chemical-class enrichment & endogenous biology (`notebooks/03`, `fig07`) |
| Irene |  | Multi-layer networks, hit zoom-in & measurement constraints (`notebooks/04`, `fig08`–`fig09`) |

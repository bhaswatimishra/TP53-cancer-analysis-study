# Exploratory Analysis of TP53 Mutations Across TCGA Studies

This repository contains a beginner-level exploratory project analyzing **TP53** mutation frequencies, codon hotspots, and 3D structural positioning using publicly available data from cBioPortal.

> **Project Intent:** This project was created primarily for hands-on learning and practice with Python-based data processing, visualization, and basic molecular modeling in ChimeraX. It is a simple exploration rather than a formal academic study.

---

## Overview of Workflow & Work Included

1. **Data Collection:** Public TCGA Pan-Cancer Atlas dataset retrieved from cBioPortal (10,967 samples across 32 study types).
2. **Data Processing & Plotting:** Self-written Python scripts (`Pandas`, `NumPy`, `Matplotlib`) used to compute sample counts, mutation percentages, hotspot distributions, and heatmap matrices.
3. **3D Structural Mapping:** Manual labeling and visualization of identified hotspot residues (positions 175, 248, 273) on the p53 core domain–DNA complex (PDB: 1TUP) using **UCSF ChimeraX**.
4. **Report Generation:** Written analysis text and formatting in the PDF report were drafted with assistance from LLM/AI tooling, then reviewed for alignment with the underlying data plots.

---

## Repository Structure

```text
TP53-cancer-analysis-study/
│
├── 3D_structure/         # ChimeraX visualization scripts and rendered structure images
├── data/                 # Raw/processed mutation data CSVs
├── graphical_data/       # Python-generated plots and heatmap images
├── README.md             # Project overview
└── TP53 analysis report.pdf  # Summary PDF report

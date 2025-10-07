# 📘 Thesis Project: "Analysing trans-chromosomal contacts in mouse embryonic stem cells using machine learning"

This project was completed as a part of my MSc Bioinformatics degree at the University of Birmingham. It was supervised by Dr Csilla Varnai from the Institute of Cancer and Genomic Sciences at the University of Birmingham.

---

## Abstract:

The 3D organisation of the genome within the nucleus plays a critical role in gene regulation, with chromatin loops and Topologically Associating Domains (TADs) being essential for bringing distant genomic regions into close proximity. This study investigates the significance of pluripotency factors (Oct4, Sox2, Nanog) and repetitive elements in the 3D genome architecture of mouse embryonic stem cells (ESCs) using machine learning models, specifically Linear Regression and Random Forest Regression. We constructed feature matrices from single-cell Hi-C data, incorporating chromosome lengths, repeat element counts, and their densities, to predict inter-chromosomal interactions.

Our analysis reveals that OSN12, a feature representing regions with high pluripotency factor binding, is consistently important in determining ESC chromosomal interactions. Additionally, SINE repeats emerged as significant, likely due to their abundance. The Random Forest Regression model incorporating topospace interactions of repeats provided the highest predictive accuracy, suggesting that spatial proximity within the 3D genome is crucial for chromosomal interactions. These findings highlight the roles of both pluripotency factors and repetitive elements in the 3D genome organisation of ESCs. This work provides a framework for further exploration of genomic features influencing 3D chromosomal architecture.

---

## 📓 Notebooks Overview:

This repository contains several Jupyter notebooks:

### 1️⃣`first_feature_matrix.ipynb`

🔍 **Goal**: Flatten individual count matrices as a vector and pass them through regression models with ESC data as the dependent variable and pluripotency factors and chromosome lengths data as the independent variables.

**Data Used**: 
- Embryonic Stem Cell CSV data
- Pluripotency factors (OSN0, OSN1, OSN4, OSN12, SE105) CSV data
- Chromosome lengths CSV data

# 🧬 Predicting Neuroblastoma Clinical Outcomes from Gene Expression Data

This project predicts key clinical outcomes in neuroblastoma patients using gene 
expression data from two platforms: **RNA-Seq (log2 FPKM)** and **Microarray probe 
intensities**. Four machine learning models were developed, evaluated, and compared 
across four clinical endpoints to identify the best performing approach for each 
outcome and platform.

## Project Aims

- **Explore and preprocess** gene expression data from both platforms using Scanpy and AnnData, 
  including quality control filtering, highly variable gene selection, PCA, Nearest neighbor graph construction and visualisation and 
  train/test splitting to prevent data leakage.

- **Predict four clinical endpoints** using gene expression features:

  - **Death from Disease** — occurrence of death from disease (yes = 1, no = 0)
  - **High Risk** — clinically considered high-risk neuroblastoma (yes = 1, no = 0)
  - **INSS Stage** — disease stage per International Neuroblastoma Staging System (1, 2, 3, 4, 4S)
  - **Progression** — occurrence of tumour progression event (yes = 1, no = 0)

- **Train and compare four classification models** across both platforms:
  - Logistic Regression
  - Random Forest
  - XGBoost
  - Neural Network

- **Evaluate model performance** using accuracy, precision, recall, F1-score and 
  ROC-AUC, with SMOTE applied to address class imbalance.

- **Identify the best performing model** per outcome and platform, and interpret 
  key predictive genes driving model decisions.

- **Generate predictions** for missing clinical outcomes using the best performing 
  trained models.

## Project Structure
```
neuroblastoma-gene-expression/
│
├── exploratory_data_analysis.ipynb   # Data loading, QC, scanpy preprocessing
├── logistic_regression.ipynb         # LR model training and evaluation
├── random_forest.ipynb               # RF model training and evaluation
├── xgboost.ipynb                     # XGBoost model training and evaluation
├── neural_network.ipynb              # Neural network training and evaluation
├── plot_results.ipynb                # Cross-model comparison and visualisation
├── predictive_modelling.ipynb        # Initial practice notebook (not directly part of analysis)
│
├── utils.py                          # Data loading, QC filtering, HVG selection
│
└── README.md
```

## Notebook Overview

## 📊 Model Performance Summary:
### RNA-Seq

| Outcome | Logistic Regression | Random Forest | XGBoost | Neural Network | **Best** |
|---------|--------------------:|-------------:|--------:|---------------:|----------|
| Death From Disease | **0.877** | 0.851 | 0.849 | 0.869 | **Logistic Regression** |
| High Risk | **0.994** | 0.985 | 0.973 | 0.992 | **Logistic Regression** |
| Progression | 0.797 | 0.832 | 0.737 | **0.833** | **Neural Network** |
| Inss Stage | N/A | **0.797** | 0.769 | 0.783 | **Random Forest** |

### Microarray

| Outcome | Logistic Regression | Random Forest | XGBoost | Neural Network | **Best** |
|---------|--------------------:|-------------:|--------:|---------------:|----------|
| Death From Disease | 0.883 | 0.895 | **0.927** | 0.790 | **XGBoost** |
| High Risk | 0.989 | 0.993 | **0.998** | 0.974 | **XGBoost** |
| Progression | 0.731 | 0.759 | **0.786** | 0.740 | **XGBoost** |
| Inss Stage | N/A | **0.790** | 0.767 | 0.744 | **Random Forest** |

## Key Findings

- **High Risk is the most predictable outcome** across all models and platforms 
  (AUC: 0.973–0.998), suggesting neuroblastoma risk classification has a strong 
  and consistent transcriptional signature detectable by all model types.

- **Platform influences model choice** — Logistic Regression dominates on RNA-Seq 
  while XGBoost consistently outperforms all models on microarray data, indicating 
  the two platforms have fundamentally different data characteristics.

- **Simpler models generalise better on RNA-Seq** — Logistic Regression achieves 
  the best or near-best performance for death from disease, high risk and progression 
  on RNA-Seq, suggesting gene expression relationships with clinical outcomes are 
  largely linear in this platform.

- **Progression and INSS Stage are the hardest outcomes to predict** (AUC: 0.73–0.83), 
  reflecting the biological complexity of disease progression and staging which may 
  require additional clinical or genomic features beyond gene expression alone.

- **Neural networks offer no consistent advantage** over traditional models on this 
  dataset, likely due to the small sample size (~200 patients), where simpler models 
  are less prone to overfitting and generalise more reliably.

## 🧠 Future Work
- **External validation** — These models were trained and evaluated on a single dataset, and 
  validation on independent neuroblastoma cohorts is essential before any clinical 
  application.

- **Feature interpretation** — identifying and validating the biological relevance 
  of top predictive genes (e.g. via pathway enrichment analysis) could provide 
  mechanistic insights into neuroblastoma outcomes.

- **Multi-modal integration** — combining RNA-Seq and microarray features, or 
  incorporating additional data types such as DNA methylation, copy number variation 
  or clinical variables, may improve prediction of harder outcomes like progression 
  and INSS stage.

- **Survival analysis** — extending binary outcome prediction to time-to-event 
  modelling (e.g. using Cox regression) would provide more clinically 
  actionable prognostic information than binary classification alone.
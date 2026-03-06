# Imports
import os
from google.colab import userdata

import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

file_path = userdata.get('DATA_PATH')

META_COLS = [
    'FactorValue..Sex.',
    'FactorValue..age.at.diagnosis.',
    'FactorValue..death.from.disease.',
    'FactorValue..high.risk.',
    'FactorValue..inss.stage.',
    'FactorValue..progression.'
]

OUTCOME_COLS = {
    'death_from_disease': 'FactorValue..death.from.disease.',
    'high_risk': 'FactorValue..high.risk.',
    'inss_stage': 'FactorValue..inss.stage.',
    'progression': 'FactorValue..progression.'
}

# Functions
def _qc_filter_rnaseq(df_genes, df_meta):
    """
    Apply scanpy QC filtering for RNA-seq data.
    Only sample/gene filtering without log transform and HVG selection.
    These steps must happen after train/test split.
    """
    adata = sc.AnnData(
        X=df_genes.values,
        obs=pd.DataFrame(index=df_genes.index),
        var=pd.DataFrame(index=df_genes.columns)
    )
    adata.obs_names = df_genes.index
    adata.var_names = df_genes.columns

    sc.pp.filter_cells(adata, min_counts=44000)
    sc.pp.filter_genes(adata, min_cells=10)

    print(f"  After QC filtering: {adata.n_obs} samples, {adata.n_vars} genes")

    df_filtered = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    df_meta_filtered = df_meta.loc[adata.obs_names]

    return df_filtered, df_meta_filtered


def _qc_filter_microarray(df_genes, df_meta):
    """
    Apply scanpy QC filtering for microarray data.
    Only sample filtering without log transform and HVG selection.
    These steps must happen after train/test split.
    """
    # Convert column names to strings to avoid AnnData warning
    df_genes.columns = df_genes.columns.astype(str)

    adata = sc.AnnData(
        X=df_genes.values,
        obs=pd.DataFrame(index=df_genes.index),
        var=pd.DataFrame(index=df_genes.columns)
    )
    adata.obs_names = df_genes.index
    adata.var_names = df_genes.columns

    sc.pp.filter_cells(adata, min_genes=43285)

    print(f"  After QC filtering: {adata.n_obs} samples, {adata.n_vars} genes")

    df_filtered = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    df_meta_filtered = df_meta.loc[adata.obs_names]

    return df_filtered, df_meta_filtered


def _log_transform_and_hvg(df_train_genes, df_test_genes, n_top_genes=2000, flavor="seurat"):
    """
    Apply log transformation and HVG selection on train set only.
    Then apply the same genes to test set.
    """
    # --- Train ---
    adata_train = sc.AnnData(
        X=df_train_genes.values,
        var=pd.DataFrame(index=df_train_genes.columns)
    )
    adata_train.obs_names = df_train_genes.index
    adata_train.var_names = df_train_genes.columns

    # sc.pp.log1p(adata_train)

    sc.pp.highly_variable_genes(adata_train, flavor=flavor, n_top_genes=n_top_genes)
    hvg_genes = adata_train.var_names[adata_train.var.highly_variable].tolist()
    adata_train = adata_train[:, adata_train.var.highly_variable]

    print(f"  HVG selected from train: {len(hvg_genes)} genes")

    df_train_filtered = pd.DataFrame(adata_train.X, index=adata_train.obs_names, columns=adata_train.var_names)

    # --- Test: apply same HVG genes ---
    adata_test = sc.AnnData(
        X=df_test_genes.values,
        var=pd.DataFrame(index=df_test_genes.columns)
    )
    adata_test.obs_names = df_test_genes.index
    adata_test.var_names = df_test_genes.columns

    # sc.pp.log1p(adata_test)
    adata_test = adata_test[:, hvg_genes]

    df_test_filtered = pd.DataFrame(adata_test.X, index=adata_test.obs_names, columns=adata_test.var_names)

    return df_train_filtered, df_test_filtered


def get_data(platform, outcome, test_size=0.2, random_state=42, n_top_genes=2000):
    """
    Get X_train, X_test, y_train, y_test for modeling.
    Features are gene expression only (age and sex excluded).
    For multiclass outcomes (inss_stage), labels are integer encoded.

    Pipeline:
      1. Load data
      2. QC filtering (sample/gene level) on full data 
      3. Train/test split (stratified)
      4. HVG selection on train only, apply same genes to test

    Parameters:
    - platform: 'rnaseq' or 'microarray'
    - outcome: 'death_from_disease', 'high_risk', 'inss_stage', 'progression'
    - test_size: fraction of data for test set (default 0.2)
    - random_state: random seed (default 42)
    - n_top_genes: number of highly variable genes to select (default 2000)

    Returns:
    - X_train, X_test, y_train, y_test
    - For inss_stage: y values are integer encoded (0, 1, 2, ...)
    """

    outcome_col = OUTCOME_COLS[outcome]

    # Load patient metadata
    df_patient_info = pd.read_csv(f"{file_path}/patientInfo.tsv", sep="\t").set_index('ID')
    df_patient_info.columns.name = 'FactorValues'
    df_patient_info = df_patient_info.sort_index(ascending=True)

    if platform == 'rnaseq':
        df_rnaseq = pd.read_csv(f"{file_path}/log2FPKM.tsv", sep="\t").rename({'00gene_id': 'gene_id'}, axis=1)
        df_rnaseq = df_rnaseq.set_index(['gene_id'])
        df_rnaseq.columns.name = 'ID'

        patient_ids = df_patient_info.index.values
        df_rnaseq.columns = patient_ids
        df_rnaseq = df_rnaseq.T

        df_merged = pd.concat([df_rnaseq, df_patient_info], axis=1)
        df_merged = df_merged[df_merged['FactorValue..death.from.disease.'].notna()]

        df_genes = df_merged.drop(columns=META_COLS)
        df_meta = df_merged[META_COLS]

        df_genes, df_meta = _qc_filter_rnaseq(df_genes, df_meta)

    elif platform == 'microarray':
        df_microarray = pd.read_csv(f"{file_path}/allProbIntensities.tsv", sep="\t").set_index(['Reporter.Identifier'])
        df_microarray.columns.name = 'ID'
        df_microarray = df_microarray.dropna(how='all')
        df_microarray = df_microarray.drop(columns='GeneSymbols')

        patient_ids = df_patient_info.index.values
        df_microarray.columns = patient_ids
        df_microarray = df_microarray.T

        df_merged = pd.concat([df_microarray, df_patient_info], axis=1)
        df_merged = df_merged[df_merged['FactorValue..death.from.disease.'].notna()]

        df_genes = df_merged.drop(columns=META_COLS)
        df_meta = df_merged[META_COLS]

        df_genes, df_meta = _qc_filter_microarray(df_genes, df_meta)

    else:
        raise ValueError("platform must be 'rnaseq' or 'microarray'")

    # Keep only the target outcome from metadata
    y = df_meta[outcome_col]

    # Remove samples with missing target outcome
    valid_mask = y.notna()
    df_genes = df_genes[valid_mask]
    y = y[valid_mask]

    # Encode labels for multiclass outcomes (e.g. inss_stage has '1','2','3','4','4S')
    # Also encode binary outcomes to ensure consistent 0/1 integer format
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    y_encoded = pd.Series(y_encoded, index=y.index)

    # Print encoding map for multiclass so you know what each integer means
    n_classes = len(le.classes_)
    if n_classes > 2:
        encoding_map = {cls: i for i, cls in enumerate(le.classes_)}
        print(f"  Label encoding: {encoding_map}")

    # Train/test split BEFORE HVG selection
    df_genes_train, df_genes_test, y_train, y_test = train_test_split(
        df_genes, y_encoded,
        test_size=test_size,
        stratify=y_encoded,
        random_state=random_state
    )

    # HVG on train only, apply same genes to test
    df_genes_train, df_genes_test = _log_transform_and_hvg(
        df_genes_train, df_genes_test, n_top_genes=n_top_genes
    )

    X_train = df_genes_train.copy()
    X_test = df_genes_test.copy()

    X_train.columns = X_train.columns.astype(str)
    X_test.columns = X_test.columns.astype(str)

    print(f"\n{platform.upper()} - {outcome.replace('_', ' ').title()}")
    print(f"  Train samples: {X_train.shape[0]} | Test samples: {X_test.shape[0]}")
    print(f"  Total features: {X_train.shape[1]}")
    print(f"  Train outcome distribution: {y_train.value_counts().to_dict()}")
    print(f"  Test outcome distribution:  {y_test.value_counts().to_dict()}")

    return X_train, X_test, y_train, y_test
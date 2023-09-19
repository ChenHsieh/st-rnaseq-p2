import streamlit as st
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2 import DeseqDataSet, DeseqStats

# Separate functions for better organization

def upload_data():
    counts_df = st.file_uploader("Upload Read Counts (CSV)", type=['csv'])
    sample_info = st.file_uploader("Upload Sample Info (CSV)", type=['csv'])
    
    if counts_df and sample_info:
        counts_df = pd.read_csv(counts_df, index_col=0)
        sample_info = pd.read_csv(sample_info, index_col=0)
        return counts_df, sample_info
    return None, None

def remove_outliers(counts_df, sample_info):
    # Placeholder for outlier detection. You can expand this with specific methods for outlier detection.
    # Here, I'll provide a simplistic way to remove samples.
    outliers = st.multiselect("Select Outliers to Remove (if any)", sample_info.index.tolist())
    if outliers:
        counts_df = counts_df.drop(outliers, axis=0)
        sample_info = sample_info.drop(outliers, axis=0)
    return counts_df, sample_info

def run_DE_analysis(counts_df, sample_info, contrast):
    # Your provided DE analysis code, slightly modified for Streamlit
    counts_df = counts_df[counts_df.sum(axis=1) >= 10].transpose().round(0)
    
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=sample_info,
        design_factors=contrast[0],
        refit_cooks=True
    )
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=contrast)
    stat_res.summary()

    return stat_res.results_df

def plot_MA(results):
    """Plot MA plot."""
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=results['log2FoldChange'], y=results['baseMean'], alpha=0.6)
    plt.axvline(0, color='red', linestyle='--')
    plt.title("MA Plot")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("Base Mean")
    st.pyplot(plt.clf())

def plot_volcano(results, alpha=0.05):
    """Plot Volcano plot."""
    plt.figure(figsize=(10, 6))
    # Assuming results have a column 'padj' for adjusted p-values
    results['-log10(p-value)'] = -np.log10(results['padj'])
    sns.scatterplot(x=results['log2FoldChange'], y=results['-log10(p-value)'], alpha=0.6)
    plt.axvline(0, color='red', linestyle='--')
    plt.axhline(-np.log10(alpha), color='blue', linestyle='--')
    plt.title("Volcano Plot")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10(p-value)")
    st.pyplot(plt.clf())

def filter_and_download_results(results):
    # Filter options
    padj_threshold = st.slider("Adjusted p-value threshold", 0.0, 1.0, 0.05)
    logfc_threshold = st.slider("Log2 Fold Change threshold", 0.0, 5.0, 1.0)
    
    filtered_results = results[(results['padj'] < padj_threshold) & (abs(results['log2FoldChange']) > logfc_threshold)]
    st.write(filtered_results)
    
    # Download option
    if st.button('Download Filtered Results'):
        csv = filtered_results.to_csv(index=True)
        b64 = base64.b64encode(csv.encode()).decode()
        href = f'<a href="data:file/csv;base64,{b64}" download="filtered_results.csv">Download CSV File</a>'
        st.markdown(href, unsafe_allow_html=True)

# Streamlit App Execution

st.title('Differential Expression Analysis with PyDESeq2')

counts_df, sample_info = upload_data()

if counts_df is not None and sample_info is not None:
    counts_df, sample_info = remove_outliers(counts_df, sample_info)
    
    # Run DE Analysis
    if st.button('Run DE Analysis'):
        contrast = ["genotype4testing", "4cl1", "control"]  # Example contrast, modify if needed
        results = run_DE_analysis(counts_df, sample_info, contrast)
        
        st.subheader("MA Plot")
        plot_MA(results)
        
        st.subheader("Volcano Plot")
        plot_volcano(results)
        
        st.subheader("Filter and Download Results")
        filter_and_download_results(results)
